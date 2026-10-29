/**
 * Stern/Dumer information-set decoder for deqode, over the GF(4) (Pauli) weight.
 *
 * - One information set is built per code and reused for a batch of syndromes, since the
 *   code is fixed for a block and every nonce gives a new syndrome.
 * - Pattern lists are bucketed by their window signature (a counting sort over 2^l
 *   buckets), so a collision lookup is an array index.
 * - Errors are stored interleaved (X word, Z word) and a candidate is built one qubit
 *   word at a time, stopping as soon as its partial weight exceeds w.
 * - A candidate is P(L) ^ P(R) ^ P(s), where P(entry) = err ^ represent(syn) is
 *   precomputed for every list entry and syndrome (represent() is linear).
 * - For a batch larger than the left list the join is done the other way round: the
 *   batch is bucketed by signature and every (left, right) pair is looked up in it. The
 *   candidates are the same, only cheaper.
 *
 * Every witness found is reported, not only the first.
 */

#include <stdlib.h>
#include <string.h>

#include "malloc-huge.h"
#include "deqode-gate.h"

// Stern's p is split between two lists, so at most this many symbols per side.
#define DQ_MAX_ORDER 4

// Upper bound on the collision window.
#define DQ_MAX_WINDOW 24

// Upper bound on the automatic batch of 2^(l+1) syndromes. The batch's P's take 16 bytes
// per qubit word per syndrome.
#define DQ_MAX_AUTO_BATCH ( (size_t)1 << 20 )

//////////////////////////////////////////////////////////////////////////////////////////
// State
//
// Errors are stored interleaved: word 2j is the X mask's word j and word 2j+1 the Z
// mask's. Witnesses passed to callers are de-interleaved.
//////////////////////////////////////////////////////////////////////////////////////////

/** One free single-qubit Pauli: a column outside the information set. */
typedef struct
{
   int qubit;
   int type;              // 1 = X, 2 = Z, 3 = Y
   uint32_t signature;    // window bits of represent(this symbol's column syndrome)
   const uint64_t *syn;   // its column syndrome, rw words, borrowed from the code
   const uint64_t *p;     // err ^ represent(syn), interleaved, pw words
} dq_symbol;

struct dq_solver
{
   dq_tuning tuning;

   int n, r, w;
   size_t nw;        // DQ_WORDS(n), words per Pauli mask
   size_t pw;        // 2 * nw, words per interleaved error
   size_t rw;        // DQ_WORDS(r), words per syndrome
   int chunks;       // ceil(r/8), representer table chunks
   int window_bits;

   uint64_t rng;

   const dq_code *code;   // borrowed for the lifetime of one information set

   // Elimination, keyed by a vector's highest set bit.
   uint8_t  *basis_used;      // r
   uint64_t *basis_vec;       // r * rw
   uint64_t *basis_wit;       // r * 2nw   (x mask then z mask)
   uint8_t  *independent;     // 2n
   int      *order;           // 2n, the shuffled column order

   // Representer tables, word-major and interleaved: 2 words per (qubit-word, chunk, byte).
   uint64_t *tab;             // nw * chunks * 256 * 2
   size_t    tab_words;
   uint32_t *sig_tab;         // chunks * 256, the window signature of each chunk value

   int *window_columns;       // window_bits

   // Free symbols, grouped one group per qubit, groups shuffled then split in half.
   dq_symbol *symbols;        // up to 3 * n
   uint64_t  *symbol_p;       // 3n * pw, each symbol's P
   uint64_t  *y_syn;          // n * rw, the Y symbols' column syndromes
   int *group_start;          // group_count + 1, indexed by the group's original number
   int *group_order;          // group_count, the shuffled visiting order
   int group_count;
   int split;                 // order positions [0,split) are left, [split,count) right

   // Pattern lists: one signature and one P per entry, the left list bucketed by
   // signature. All four arrays live in `pool`, P's first, signatures after them.
   uint64_t *left_p, *right_p;
   uint32_t *left_sig, *right_sig;
   size_t left_count, right_count, max_entries;
   uint32_t *bucket;          // (1 << window_bits) + 1 offsets into the left list
   size_t bucket_slots;
   void *pool;
   size_t pool_bytes;
   bool pool_huge;

   // The syndrome batch of the swapped join, grown on demand.
   uint32_t *batch_sig;       // batch_cap
   uint32_t *batch_order;     // batch_cap, sorted position -> index in the caller's batch
   uint64_t *batch_p;         // batch_cap * pw, P(s) in sorted order
   uint32_t *batch_bucket;    // batch_slots, offsets into the sorted batch
   size_t batch_cap, batch_slots;

   // Per-pass scratch.
   uint64_t *scratch_syn;     // rw
   uint64_t *target_p;        // pw, one syndrome's P
   uint64_t *cand_x;          // nw
   uint64_t *cand_z;          // nw
   uint64_t *unit;            // 2nw, one represented unit syndrome

   bool have_set;
};

//////////////////////////////////////////////////////////////////////////////////////////
// Small helpers
//////////////////////////////////////////////////////////////////////////////////////////

static inline uint64_t dq_rand( struct dq_solver *s )
{
   // splitmix64
   uint64_t z = ( s->rng += 0x9E3779B97F4A7C15ULL );
   z = ( z ^ ( z >> 30 ) ) * 0xBF58476D1CE4E5B9ULL;
   z = ( z ^ ( z >> 27 ) ) * 0x94D049BB133111EBULL;
   return z ^ ( z >> 31 );
}

static inline void xor_words( uint64_t *dst, const uint64_t *src, size_t words )
{
   for ( size_t i = 0; i < words; i++ ) dst[i] ^= src[i];
}

static inline int highest_bit( const uint64_t *bits, size_t words )
{
   for ( size_t i = words; i-- > 0; )
      if ( bits[i] ) return (int)( 64 * i ) + 63 - __builtin_clzll( bits[i] );
   return -1;
}

/**
 * Allocate an array that is read at random, asking for transparent huge pages to avoid
 * TLB misses. Freed with free().
 */
static void *batch_alloc( size_t bytes )
{
   const size_t huge = (size_t)2 << 20;
   void *p = NULL;
   if ( bytes >= huge )
   {
      if ( posix_memalign( &p, huge, ( bytes + huge - 1 ) & ~( huge - 1 ) ) ) p = NULL;
#if defined(MADV_HUGEPAGE)
      if ( p ) madvise( p, ( bytes + huge - 1 ) & ~( huge - 1 ), MADV_HUGEPAGE );
#endif
   }
   return p ? p : malloc( bytes );
}

//////////////////////////////////////////////////////////////////////////////////////////
// Tuning
//////////////////////////////////////////////////////////////////////////////////////////

void dq_solver_autotune( dq_tuning *tuning, const dq_params *params )
{
   // p = 2 + 2 measured best from n = 112 to 160.
   if ( tuning->left_order <= 0 && tuning->right_order <= 0 )
   {
      tuning->left_order  = 2;
      tuning->right_order = 2;
   }
   if ( tuning->left_order  <= 0 ) tuning->left_order  = tuning->right_order;
   if ( tuning->right_order <= 0 ) tuning->right_order = tuning->left_order;

   // Each pattern symbol costs one unit of weight, so keep p at most w / 2.
   int budget = params->w / 2;
   if ( budget < 1 ) budget = 1;
   while ( tuning->left_order + tuning->right_order > budget )
   {
      if ( tuning->left_order >= tuning->right_order ) tuning->left_order--;
      else                                            tuning->right_order--;
      if ( tuning->left_order + tuning->right_order <= 1 ) break;
   }
   if ( tuning->left_order  < 0 ) tuning->left_order  = 0;
   if ( tuning->right_order < 0 ) tuning->right_order = 0;
   if ( tuning->left_order  > DQ_MAX_ORDER ) tuning->left_order  = DQ_MAX_ORDER;
   if ( tuning->right_order > DQ_MAX_ORDER ) tuning->right_order = DQ_MAX_ORDER;

   // window_bits 0 means choose it per information set (choose_window_bits).
   if ( tuning->window_bits > DQ_MAX_WINDOW ) tuning->window_bits = DQ_MAX_WINDOW;
   if ( tuning->window_bits > params->r )     tuning->window_bits = params->r;
   if ( tuning->window_bits < 0 )             tuning->window_bits = 0;

   // Per thread pattern list budget; the default orders need a few MiB.
   if ( !tuning->memory_bytes ) tuning->memory_bytes = (size_t)64 << 20;

   // Nonces per information set; 0 means automatic (dq_solver_batch_size).
   if ( tuning->batch < 0 ) tuning->batch = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
// Construction
//////////////////////////////////////////////////////////////////////////////////////////

void dq_solver_free( struct dq_solver *s )
{
   if ( !s ) return;
   free( s->basis_used );  free( s->basis_vec );   free( s->basis_wit );
   free( s->independent ); free( s->order );       free( s->tab );
   free( s->sig_tab );
   free( s->window_columns ); free( s->symbols );  free( s->symbol_p );
   free( s->group_start ); free( s->group_order ); free( s->y_syn );
   free( s->bucket );
   free( s->batch_sig );   free( s->batch_order ); free( s->batch_p );
   free( s->batch_bucket );
   free( s->scratch_syn ); free( s->target_p );
   free( s->cand_x );      free( s->cand_z );      free( s->unit );
   if ( s->pool )
   {
      // malloc-huge.c has no matching free; unmap with the same 2 MiB rounding.
#if defined(MAP_HUGETLB) && defined(MAP_ANON)
      if ( s->pool_huge )
      {
         const size_t mask = (size_t)( 2 * 1024 * 1024 ) - 1;
         munmap( s->pool, ( s->pool_bytes + mask ) & ~mask );
      }
      else
#endif
         free( s->pool );
   }
   free( s );
}

struct dq_solver *dq_solver_new( const dq_params *params, const dq_tuning *tuning,
                                 uint64_t seed )
{
   if ( !dq_params_are_sane( params ) ) return NULL;

   struct dq_solver *s = (struct dq_solver*) calloc( 1, sizeof(*s) );
   if ( !s ) return NULL;

   s->tuning = *tuning;
   dq_solver_autotune( &s->tuning, params );

   s->n  = params->n;
   s->r  = params->r;
   s->w  = params->w;
   s->nw = DQ_WORDS( params->n );
   s->pw = 2 * s->nw;
   s->rw = DQ_WORDS( params->r );
   s->chunks = ( params->r + 7 ) / 8;
   s->window_bits = 0;   // chosen per information set, see choose_window_bits
   s->rng = seed ? seed : 0x243F6A8885A308D3ULL;

   const size_t n2 = 2 * (size_t)s->n;
   s->basis_used  = (uint8_t*)  calloc( (size_t)s->r, 1 );
   s->basis_vec   = (uint64_t*) calloc( (size_t)s->r * s->rw, 8 );
   s->basis_wit   = (uint64_t*) calloc( (size_t)s->r * s->pw, 8 );
   s->independent = (uint8_t*)  calloc( n2, 1 );
   s->order       = (int*)      calloc( n2, sizeof(int) );
   s->window_columns = (int*)   calloc( (size_t)DQ_MAX_WINDOW, sizeof(int) );
   s->symbols     = (dq_symbol*) calloc( 3 * (size_t)s->n, sizeof(dq_symbol) );
   s->symbol_p    = (uint64_t*) calloc( 3 * (size_t)s->n * s->pw, 8 );
   s->group_start = (int*)      calloc( (size_t)s->n + 1, sizeof(int) );
   s->group_order = (int*)      calloc( (size_t)s->n, sizeof(int) );
   s->y_syn       = (uint64_t*) calloc( (size_t)s->n * s->rw, 8 );
   s->scratch_syn = (uint64_t*) calloc( s->rw, 8 );
   s->target_p    = (uint64_t*) calloc( s->pw, 8 );
   s->cand_x      = (uint64_t*) calloc( s->nw, 8 );
   s->cand_z      = (uint64_t*) calloc( s->nw, 8 );
   s->unit        = (uint64_t*) calloc( s->pw, 8 );
   s->sig_tab     = (uint32_t*) calloc( (size_t)s->chunks * 256, sizeof(uint32_t) );

   s->tab_words = (size_t)s->nw * s->chunks * 256 * 2;
   s->tab = (uint64_t*) calloc( s->tab_words, 8 );

   if ( !s->basis_used || !s->basis_vec || !s->basis_wit || !s->independent || !s->order
     || !s->window_columns || !s->symbols || !s->symbol_p || !s->group_start
     || !s->group_order || !s->y_syn || !s->scratch_syn || !s->target_p
     || !s->cand_x || !s->cand_z || !s->unit || !s->sig_tab || !s->tab )
   {
      dq_solver_free( s );
      return NULL;
   }

   // One entry is a P and a signature.
   const size_t entry_bytes = s->pw * 8 + sizeof(uint32_t);

   // The pool holds both lists' P's and then both lists' signatures. The bucket array is
   // separate because its length depends on the window chosen per information set.
   s->max_entries = s->tuning.memory_bytes / entry_bytes;
   if ( s->max_entries < 64 ) { dq_solver_free( s ); return NULL; }

   s->pool_bytes = s->max_entries * entry_bytes;
   s->pool = malloc_hugepages( s->pool_bytes );
   s->pool_huge = s->pool != NULL;
   if ( !s->pool ) s->pool = malloc( s->pool_bytes );
   if ( !s->pool ) { dq_solver_free( s ); return NULL; }
   s->left_p = (uint64_t*) s->pool;
   s->left_sig = (uint32_t*)( s->left_p + s->max_entries * s->pw );
   return s;
}

size_t dq_solver_batch_size( const struct dq_solver *s )
{
   if ( s->tuning.batch > 0 ) return (size_t)s->tuning.batch;
   const size_t batch = (size_t)2 << s->window_bits;
   return batch < DQ_MAX_AUTO_BATCH ? batch : DQ_MAX_AUTO_BATCH;
}

//////////////////////////////////////////////////////////////////////////////////////////
// Information set: everything that does not depend on the syndrome
//////////////////////////////////////////////////////////////////////////////////////////

/**
 * Add a column to the elimination, keyed by highest set bit, keeping the Pauli that
 * produced each pivot. Returns true if the column was independent.
 */
static bool basis_add( struct dq_solver *s, int column )
{
   const size_t rw = s->rw, pw = s->pw;
   uint64_t *value   = s->scratch_syn;
   uint64_t *witness = s->unit;

   memcpy( value, dq_column_syndrome( s->code, column ), rw * 8 );
   memset( witness, 0, pw * 8 );
   // Column q is an X on qubit q, column n+q a Z on qubit q.
   if ( column < s->n ) witness[ (size_t)column / 64 ] |= (uint64_t)1 << ( column % 64 );
   else
   {
      const int q = column - s->n;
      witness[ s->nw + (size_t)q / 64 ] |= (uint64_t)1 << ( q % 64 );
   }

   for ( ;; )
   {
      const int pivot = highest_bit( value, rw );
      if ( pivot < 0 ) return false;              // reduced to zero: linearly dependent
      if ( !s->basis_used[ pivot ] )
      {
         s->basis_used[ pivot ] = 1;
         memcpy( s->basis_vec + (size_t)pivot * rw, value,   rw * 8 );
         memcpy( s->basis_wit + (size_t)pivot * pw, witness, pw * 8 );
         return true;
      }
      xor_words( value,   s->basis_vec + (size_t)pivot * rw, rw );
      xor_words( witness, s->basis_wit + (size_t)pivot * pw, pw );
   }
}

/** The witness for a syndrome, if the basis spans it. `out` holds 2nw words. */
static bool basis_represent( const struct dq_solver *s, const uint64_t *target,
                             uint64_t *scratch, uint64_t *out )
{
   const size_t rw = s->rw, pw = s->pw;
   memcpy( scratch, target, rw * 8 );
   memset( out, 0, pw * 8 );
   for ( ;; )
   {
      const int pivot = highest_bit( scratch, rw );
      if ( pivot < 0 ) return true;
      if ( !s->basis_used[ pivot ] ) return false;
      xor_words( scratch, s->basis_vec + (size_t)pivot * rw, rw );
      xor_words( out,     s->basis_wit + (size_t)pivot * pw, pw );
   }
}

/** The window bits of a represented witness: its coefficients on the window's columns. */
static uint64_t basis_signature( const struct dq_solver *s, const uint64_t *witness )
{
   uint64_t signature = 0;
   for ( int i = 0; i < s->window_bits; i++ )
   {
      const int column = s->window_columns[i];
      const int q      = column < s->n ? column : column - s->n;
      const size_t off = column < s->n ? 0 : s->nw;
      if ( ( witness[ off + (size_t)q / 64 ] >> ( q % 64 ) ) & 1 )
         signature |= (uint64_t)1 << i;
   }
   return signature;
}

/** Build the word-major representer tables from the r represented unit syndromes. */
static bool build_tables( struct dq_solver *s )
{
   const size_t rw = s->rw, pw = s->pw, nw = s->nw;
   uint64_t *units = (uint64_t*) calloc( (size_t)s->r * pw, 8 );
   uint64_t *unit_syn = (uint64_t*) calloc( rw, 8 );
   uint64_t *scratch  = (uint64_t*) calloc( rw, 8 );
   if ( !units || !unit_syn || !scratch )
   {
      free( units ); free( unit_syn ); free( scratch );
      return false;
   }

   for ( int bit = 0; bit < s->r; bit++ )
   {
      memset( unit_syn, 0, rw * 8 );
      unit_syn[ (size_t)bit / 64 ] = (uint64_t)1 << ( bit % 64 );
      if ( !basis_represent( s, unit_syn, scratch, units + (size_t)bit * pw ) )
      {
         free( units ); free( unit_syn ); free( scratch );
         return false;   // not full rank; the caller draws a different set or gives up
      }
   }

   memset( s->tab, 0, s->tab_words * 8 );
   // tab[(j * chunks + c) * 256 + byte] holds the pair (x_j, z_j). Each entry is the entry
   // without its lowest set bit XOR that bit's row.
   for ( int c = 0; c < s->chunks; c++ )
   {
      const int start = 8 * c;
      const int width = ( s->r - start ) < 8 ? ( s->r - start ) : 8;
      const int values = 1 << width;
      for ( int value = 1; value < values; value++ )
      {
         const int low    = value & -value;
         const int offset = __builtin_ctz( (unsigned)low );
         const uint64_t *u = units + (size_t)( start + offset ) * pw;
         for ( size_t j = 0; j < nw; j++ )
         {
            uint64_t *dst  = s->tab + ( ( j * s->chunks + c ) * 256 + (size_t)value ) * 2;
            const uint64_t *prev = s->tab
                                 + ( ( j * s->chunks + c ) * 256 + (size_t)( value ^ low ) ) * 2;
            dst[0] = prev[0] ^ u[j];
            dst[1] = prev[1] ^ u[ nw + j ];
         }
      }
   }

   free( units ); free( unit_syn ); free( scratch );
   return true;
}

/**
 * The window signature of every chunk value. A syndrome's signature is the XOR of its
 * chunks' signatures, since represent() is linear.
 */
static void build_signature_table( struct dq_solver *s )
{
   const int chunks = s->chunks;
   for ( int c = 0; c < chunks; c++ )
      for ( int value = 0; value < 256; value++ )
      {
         uint32_t signature = 0;
         for ( int i = 0; i < s->window_bits; i++ )
         {
            const int column = s->window_columns[i];
            const int q      = column < s->n ? column : column - s->n;
            const uint64_t *e = s->tab
                              + ( ( (size_t)q / 64 * chunks + c ) * 256 + (size_t)value ) * 2;
            if ( ( e[ column < s->n ? 0 : 1 ] >> ( q % 64 ) ) & 1 )
               signature |= (uint32_t)1 << i;
         }
         s->sig_tab[ (size_t)c * 256 + value ] = signature;
      }
}

/** P(s) = represent(s), interleaved, from the chunk tables. */
static inline void represent_syndrome( const struct dq_solver *s, const uint64_t *syndrome,
                                       uint64_t *p )
{
   const int chunks = s->chunks;
   for ( size_t j = 0; j < s->nw; j++ )
   {
      const uint64_t *tj = s->tab + j * (size_t)chunks * 512;
      uint64_t xw = 0, zw = 0;
      for ( int c = 0; c < chunks; c++ )
      {
         const unsigned byte = (unsigned)( syndrome[ c / 8 ] >> ( 8 * ( c % 8 ) ) ) & 0xFF;
         const uint64_t *e = tj + ( (size_t)c * 256 + byte ) * 2;
         xw ^= e[0];
         zw ^= e[1];
      }
      p[ 2 * j ]     = xw;
      p[ 2 * j + 1 ] = zw;
   }
}

static inline uint32_t syndrome_signature( const struct dq_solver *s,
                                           const uint64_t *syndrome )
{
   uint32_t signature = 0;
   for ( int c = 0; c < s->chunks; c++ )
      signature ^= s->sig_tab[ (size_t)c * 256
                             + ( ( syndrome[ c / 8 ] >> ( 8 * ( c % 8 ) ) ) & 0xFF ) ];
   return signature;
}

/** Pick the l window columns: paired X/Z on one qubit first, then whatever is left. */
static bool build_window( struct dq_solver *s )
{
   const int n = s->n;
   int *paired  = (int*) malloc( (size_t)n * sizeof(int) );
   int *singles = (int*) malloc( (size_t)2 * n * sizeof(int) );
   if ( !paired || !singles ) { free( paired ); free( singles ); return false; }

   int np = 0, ns = 0;
   for ( int q = 0; q < n; q++ )
   {
      const bool xin = s->independent[ q ];
      const bool zin = s->independent[ n + q ];
      if      ( xin && zin ) paired[ np++ ]  = q;
      else if ( xin )        singles[ ns++ ] = q;
      else if ( zin )        singles[ ns++ ] = n + q;
   }
   for ( int i = np; i > 1; i-- )
   {
      const int j = (int)( dq_rand( s ) % (uint64_t)i );
      const int t = paired[i-1]; paired[i-1] = paired[j]; paired[j] = t;
   }
   for ( int i = ns; i > 1; i-- )
   {
      const int j = (int)( dq_rand( s ) % (uint64_t)i );
      const int t = singles[i-1]; singles[i-1] = singles[j]; singles[j] = t;
   }

   int filled = 0;
   for ( int i = 0; i < np && filled < s->window_bits; i++ )
   {
      s->window_columns[ filled++ ] = paired[i];
      if ( filled < s->window_bits ) s->window_columns[ filled++ ] = n + paired[i];
   }
   for ( int i = 0; i < ns && filled < s->window_bits; i++ )
      s->window_columns[ filled++ ] = singles[i];

   free( paired ); free( singles );
   return filled == s->window_bits;
}

/** Group the free columns by qubit. */
static bool build_symbols( struct dq_solver *s )
{
   const int n = s->n;
   const size_t rw = s->rw;

   int count = 0;
   s->group_count = 0;
   s->group_start[0] = 0;
   for ( int q = 0; q < n; q++ )
   {
      const bool xfree = !s->independent[ q ];
      const bool zfree = !s->independent[ n + q ];
      if ( !xfree && !zfree ) continue;
      // X, Z and Y on one qubit each count as weight 1.
      if ( xfree )
         s->symbols[ count++ ] = (dq_symbol){ q, 1, 0, dq_column_syndrome( s->code, q ), NULL };
      if ( zfree )
         s->symbols[ count++ ] = (dq_symbol){ q, 2, 0, dq_column_syndrome( s->code, n + q ), NULL };
      if ( xfree && zfree )
      {
         // A Y's column syndrome is the XOR of the X and Z columns.
         uint64_t *ys = s->y_syn + (size_t)q * rw;
         memcpy( ys, dq_column_syndrome( s->code, q ), rw * 8 );
         xor_words( ys, dq_column_syndrome( s->code, n + q ), rw );
         s->symbols[ count++ ] = (dq_symbol){ q, 3, 0, ys, NULL };
      }
      s->group_start[ ++s->group_count ] = count;
   }
   if ( s->group_count < 2 ) return false;

   // Shuffle the groups and split them in half, so the two lists use disjoint qubits and
   // the partition changes with every information set.
   for ( int i = 0; i < s->group_count; i++ ) s->group_order[i] = i;
   for ( int i = s->group_count; i > 1; i-- )
   {
      const int j = (int)( dq_rand( s ) % (uint64_t)i );
      const int t = s->group_order[i-1];
      s->group_order[i-1] = s->group_order[j];
      s->group_order[j] = t;
   }

   s->split = s->group_count / 2;
   return true;
}

/**
 * Give every free symbol its window signature and its P. A list entry's P is the XOR of
 * its symbols' P's.
 */
static bool build_signatures( struct dq_solver *s )
{
   const size_t nw = s->nw, pw = s->pw;
   uint64_t *witness = (uint64_t*) malloc( pw * 8 );
   if ( !witness ) return false;
   const int count = s->group_start[ s->group_count ];
   for ( int i = 0; i < count; i++ )
   {
      dq_symbol *sym = &s->symbols[i];
      if ( !basis_represent( s, sym->syn, s->scratch_syn, witness ) )
      {
         free( witness );
         return false;
      }
      sym->signature = (uint32_t)basis_signature( s, witness );

      uint64_t *p = s->symbol_p + (size_t)i * pw;
      for ( size_t j = 0; j < nw; j++ )
      {
         p[ 2 * j ]     = witness[ j ];
         p[ 2 * j + 1 ] = witness[ nw + j ];
      }
      const size_t word = (size_t)sym->qubit / 64;
      const uint64_t bit = (uint64_t)1 << ( sym->qubit % 64 );
      if ( sym->type & 1 ) p[ 2 * word ]     ^= bit;
      if ( sym->type & 2 ) p[ 2 * word + 1 ] ^= bit;
      sym->p = p;
   }
   free( witness );
   return true;
}

/**
 * The number of patterns a side would enumerate: the elementary symmetric polynomial of
 * degree `order` in the group sizes, saturating at `cap`.
 */
static size_t count_side( const struct dq_solver *s, int group_lo, int group_hi,
                          int order, size_t cap )
{
   if ( order < 0 ) return 0;
   if ( order > group_hi - group_lo ) order = group_hi - group_lo;
   size_t dp[ DQ_MAX_ORDER + 1 ] = { 0 };
   dp[0] = 1;
   for ( int g = group_lo; g < group_hi; g++ )
   {
      const int group = s->group_order[g];
      const size_t size = (size_t)( s->group_start[group+1] - s->group_start[group] );
      for ( int k = order; k >= 1; k-- )
      {
         const size_t add = dp[k-1] > cap / ( size ? size : 1 ) ? cap : dp[k-1] * size;
         dp[k] = dp[k] > cap - add ? cap : dp[k] + add;
      }
   }
   return dp[ order ];
}

/**
 * The window width l. At least the balanced Stern choice (2^l near the list size). With
 * the swapped join a wider window makes candidates more likely to be light at the same
 * cost, so l grows while 2^l <= 2T for a batch of T syndromes and a nonce still yields
 * DQ_MIN_CANDIDATES_PER_NONCE candidates.
 */
#define DQ_MIN_CANDIDATES_PER_NONCE 128

static int choose_window_bits( const struct dq_solver *s, size_t left, size_t right )
{
   size_t smaller = left < right ? left : right;
   if ( smaller < 2 ) return 1;
   int bits = 63 - __builtin_clzll( (uint64_t)smaller );

   const size_t batch = s->tuning.batch > 0 ? (size_t)s->tuning.batch : DQ_MAX_AUTO_BATCH;
   const int batch_bits = 64 - __builtin_clzll( (uint64_t)batch );     // 2^l <= 2T
   const double per_nonce = (double)left * (double)right / DQ_MIN_CANDIDATES_PER_NONCE;
   int wide = per_nonce >= 2 ? 63 - __builtin_clzll( (uint64_t)per_nonce ) : 1;
   if ( wide > batch_bits ) wide = batch_bits;
   if ( wide > bits ) bits = wide;

   if ( bits > DQ_MAX_WINDOW ) bits = DQ_MAX_WINDOW;
   if ( bits > s->r )          bits = s->r;
   if ( bits < 1 )             bits = 1;
   return bits;
}

/** Recursive enumeration of one side's pattern list. */
typedef struct
{
   struct dq_solver *s;
   uint32_t *sig;      // out, one per entry
   uint64_t *p;        // out, pw words per entry
   size_t count;
   size_t cap;
   int group_hi;
   uint32_t signature;
   uint64_t *run;      // pw, the running P
} dq_enum;

static void enumerate( dq_enum *e, int start, int remaining )
{
   struct dq_solver *s = e->s;
   const size_t pw = s->pw;

   if ( e->count >= e->cap ) return;
   if ( remaining == 0 )
   {
      e->sig[ e->count ] = e->signature;
      memcpy( e->p + e->count * pw, e->run, pw * 8 );
      e->count++;
      return;
   }

   for ( int g = start; g + remaining <= e->group_hi; g++ )
   {
      const int group = s->group_order[g];
      for ( int k = s->group_start[group]; k < s->group_start[group+1]; k++ )
      {
         const dq_symbol *sym = &s->symbols[k];

         xor_words( e->run, sym->p, pw );
         e->signature ^= sym->signature;

         enumerate( e, g + 1, remaining - 1 );

         e->signature ^= sym->signature;
         xor_words( e->run, sym->p, pw );

         if ( e->count >= e->cap ) return;
      }
   }
}

static size_t enumerate_side( struct dq_solver *s, uint32_t *sig, uint64_t *p, size_t cap,
                              int group_lo, int group_hi, int order )
{
   uint64_t *run = (uint64_t*) calloc( s->pw, 8 );
   if ( !run ) return 0;

   // A side with fewer groups than its order would have no patterns.
   if ( order > group_hi - group_lo ) order = group_hi - group_lo;

   dq_enum e = { s, sig, p, 0, cap, group_hi, 0, run };
   if ( order >= 0 ) enumerate( &e, group_lo, order );

   free( run );
   return e.count;
}

bool dq_solver_new_information_set( struct dq_solver *s, const dq_code *code )
{
   if ( code->n != s->n || code->r != s->r ) return false;
   s->code = code;
   s->have_set = false;

   const size_t n2 = 2 * (size_t)s->n;
   for ( size_t i = 0; i < n2; i++ ) s->order[i] = (int)i;
   for ( size_t i = n2; i > 1; i-- )
   {
      const size_t j = (size_t)( dq_rand( s ) % (uint64_t)i );
      const int t = s->order[i-1]; s->order[i-1] = s->order[j]; s->order[j] = t;
   }

   memset( s->basis_used, 0, (size_t)s->r );
   memset( s->independent, 0, n2 );
   int rank = 0;
   // Offer every column; the dependent ones are the free columns for the pattern lists.
   for ( size_t i = 0; i < n2; i++ )
   {
      const int column = s->order[i];
      if ( basis_add( s, column ) ) { s->independent[ column ] = 1; rank++; }
   }
   if ( rank != s->r ) return false;   // the code's syndrome map is not onto

   if ( !build_tables( s ) )  return false;
   if ( !build_symbols( s ) ) return false;

   // Size the lists before building them, so the window can be chosen to match.
   const size_t half = s->max_entries / 2;
   size_t want_left  = count_side( s, 0, s->split, s->tuning.left_order, half );
   size_t want_right = count_side( s, s->split, s->group_count, s->tuning.right_order,
                                   s->max_entries - half );
   if ( !want_left || !want_right ) return false;

   s->window_bits = s->tuning.window_bits > 0 ? s->tuning.window_bits
                                              : choose_window_bits( s, want_left, want_right );

   const size_t slots = ( (size_t)1 << s->window_bits ) + 1;
   if ( slots != s->bucket_slots )
   {
      free( s->bucket );
      s->bucket = (uint32_t*) malloc( slots * sizeof(uint32_t) );
      if ( !s->bucket ) { s->bucket_slots = 0; return false; }
      s->bucket_slots = slots;
   }

   if ( !build_window( s ) )     return false;
   build_signature_table( s );
   if ( !build_signatures( s ) ) return false;

   s->left_count = enumerate_side( s, s->left_sig, s->left_p, want_left, 0, s->split,
                                   s->tuning.left_order );
   s->right_p   = s->left_p + s->left_count * s->pw;
   s->right_sig = s->left_sig + s->left_count;
   const size_t right_cap = s->max_entries - s->left_count;
   s->right_count = enumerate_side( s, s->right_sig, s->right_p,
                                    want_right < right_cap ? want_right : right_cap,
                                    s->split, s->group_count, s->tuning.right_order );
   if ( !s->left_count || !s->right_count ) return false;

   // Bucket the left list by signature (counting sort over 2^l buckets).
   memset( s->bucket, 0, s->bucket_slots * sizeof(uint32_t) );
   for ( size_t i = 0; i < s->left_count; i++ ) s->bucket[ s->left_sig[i] + 1 ]++;
   for ( size_t i = 1; i < s->bucket_slots; i++ ) s->bucket[i] += s->bucket[i-1];

   {
      const size_t pw = s->pw;
      uint64_t *tmp_p   = (uint64_t*) malloc( s->left_count * pw * 8 );
      uint32_t *tmp_sig = (uint32_t*) malloc( s->left_count * sizeof(uint32_t) );
      uint32_t *at      = (uint32_t*) malloc( s->bucket_slots * sizeof(uint32_t) );
      if ( !tmp_p || !tmp_sig || !at ) { free( tmp_p ); free( tmp_sig ); free( at ); return false; }
      memcpy( at, s->bucket, s->bucket_slots * sizeof(uint32_t) );
      for ( size_t i = 0; i < s->left_count; i++ )
      {
         const uint32_t k = at[ s->left_sig[i] ]++;
         tmp_sig[k] = s->left_sig[i];
         memcpy( tmp_p + (size_t)k * pw, s->left_p + i * pw, pw * 8 );
      }
      memcpy( s->left_p, tmp_p, s->left_count * pw * 8 );
      memcpy( s->left_sig, tmp_sig, s->left_count * sizeof(uint32_t) );
      free( at ); free( tmp_sig ); free( tmp_p );
   }

   s->have_set = true;
   return true;
}

//////////////////////////////////////////////////////////////////////////////////////////
// Join: syndromes against the current information set
//////////////////////////////////////////////////////////////////////////////////////////

/**
 * One candidate, P(L) ^ P(R) ^ P(s), built one qubit word at a time and abandoned as soon
 * as its weight exceeds w. On success the error is in cand_x / cand_z.
 */
static inline bool candidate( struct dq_solver *s, const uint64_t *a, const uint64_t *b,
                              const uint64_t *c )
{
   const size_t nw = s->nw;
   const int w = s->w;
   int weight = 0;
   for ( size_t j = 0; j < nw; j++ )
   {
      const uint64_t xw = a[ 2 * j ]     ^ b[ 2 * j ]     ^ c[ 2 * j ];
      const uint64_t zw = a[ 2 * j + 1 ] ^ b[ 2 * j + 1 ] ^ c[ 2 * j + 1 ];
      weight += __builtin_popcountll( xw | zw );
      if ( weight > w ) return false;
      s->cand_x[j] = xw;
      s->cand_z[j] = zw;
   }
   return true;
}

/** Grow the swapped join's batch arrays to `count` syndromes at the current window. */
static bool reserve_batch( struct dq_solver *s, size_t count )
{
   if ( count > s->batch_cap )
   {
      free( s->batch_sig ); free( s->batch_order ); free( s->batch_p );
      s->batch_sig   = (uint32_t*) batch_alloc( count * sizeof(uint32_t) );
      s->batch_order = (uint32_t*) batch_alloc( count * sizeof(uint32_t) );
      s->batch_p     = (uint64_t*) batch_alloc( count * s->pw * 8 );
      if ( !s->batch_sig || !s->batch_order || !s->batch_p )
      {
         free( s->batch_sig ); free( s->batch_order ); free( s->batch_p );
         s->batch_sig = s->batch_order = NULL;
         s->batch_p = NULL;
         s->batch_cap = 0;
         return false;
      }
      s->batch_cap = count;
   }
   if ( s->bucket_slots > s->batch_slots )
   {
      free( s->batch_bucket );
      s->batch_bucket = (uint32_t*) batch_alloc( s->bucket_slots * sizeof(uint32_t) );
      if ( !s->batch_bucket ) { s->batch_slots = 0; return false; }
      s->batch_slots = s->bucket_slots;
   }
   return true;
}

bool dq_solver_join_batch( struct dq_solver *s, const uint64_t *syndromes, size_t count,
                           dq_batch_witness_fn on_witness, void *context,
                           const volatile uint8_t *abort_flag, double *progress )
{
   if ( progress ) *progress = 0;
   if ( !s->have_set || !count ) return false;

   const size_t pw = s->pw, rw = s->rw;
   const size_t L = s->left_count, R = s->right_count;

   if ( count <= L )
   {
      // Batch smaller than the left list: each syndrome probes every right entry against
      // the left buckets.
      uint64_t *ps = s->target_p;
      for ( size_t t = 0; t < count; t++ )
      {
         const uint64_t *syndrome = syndromes + t * rw;
         represent_syndrome( s, syndrome, ps );
         const uint32_t target_sig = syndrome_signature( s, syndrome );

         for ( size_t ri = 0; ri < R; ri++ )
         {
            if ( abort_flag && ( ri & 0x3FF ) == 0 && *abort_flag )
            {
               if ( progress ) *progress = (double)t / (double)count;
               return false;
            }
            const uint32_t key = s->right_sig[ri] ^ target_sig;
            const uint32_t lo = s->bucket[ key ], hi = s->bucket[ key + 1 ];
            if ( lo == hi ) continue;
            const uint64_t *pr = s->right_p + ri * pw;
            for ( uint32_t li = lo; li < hi; li++ )
            {
               if ( candidate( s, s->left_p + (size_t)li * pw, pr, ps )
                 && on_witness( context, t, s->cand_x, s->cand_z ) )
               {
                  if ( progress ) *progress = (double)( t + 1 ) / (double)count;
                  return true;
               }
            }
         }
      }
      if ( progress ) *progress = 1;
      return false;
   }

   // Otherwise bucket the batch by signature and look up every (left, right) pair in it.
   if ( !reserve_batch( s, count ) ) return false;
   uint32_t *sig = s->batch_sig, *order = s->batch_order, *bucket = s->batch_bucket;
   uint64_t *bp = s->batch_p;
   const size_t slots = s->bucket_slots;

   memset( bucket, 0, slots * sizeof(uint32_t) );
   for ( size_t t = 0; t < count; t++ )
   {
      sig[t] = syndrome_signature( s, syndromes + t * rw );
      bucket[ sig[t] + 1 ]++;
   }
   for ( size_t i = 1; i < slots; i++ ) bucket[i] += bucket[i-1];
   // After placing through bucket[sig]++, shift the starts back.
   for ( size_t t = 0; t < count; t++ )
   {
      const uint32_t k = bucket[ sig[t] ]++;
      order[k] = (uint32_t)t;
      represent_syndrome( s, syndromes + t * rw, bp + (size_t)k * pw );
   }
   memmove( bucket + 1, bucket, ( slots - 1 ) * sizeof(uint32_t) );
   bucket[0] = 0;

   for ( size_t li = 0; li < L; li++ )
   {
      if ( abort_flag && *abort_flag )
      {
         if ( progress ) *progress = (double)li / (double)L;
         return false;
      }
      const uint32_t lsig = s->left_sig[li];
      const uint64_t *pl = s->left_p + li * pw;
      for ( size_t ri = 0; ri < R; ri++ )
      {
         const uint32_t key = lsig ^ s->right_sig[ri];
         const uint32_t lo = bucket[ key ], hi = bucket[ key + 1 ];
         if ( lo == hi ) continue;
         const uint64_t *pr = s->right_p + ri * pw;
         for ( uint32_t k = lo; k < hi; k++ )
         {
            if ( candidate( s, pl, pr, bp + (size_t)k * pw )
              && on_witness( context, order[k], s->cand_x, s->cand_z ) )
            {
               if ( progress ) *progress = (double)( li + 1 ) / (double)L;
               return true;
            }
         }
      }
   }
   if ( progress ) *progress = 1;
   return false;
}
