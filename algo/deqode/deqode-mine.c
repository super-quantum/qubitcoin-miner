/**
 * The deqode mining loop.
 *
 * The code depends only on hashPrevBlock, so it is derived once per work item. Each
 * nonce gives a new syndrome. One information set is built per batch of nonces and all
 * of the batch's syndromes are decoded against it, and every witness found is tried
 * against the acceptance target.
 */

#include <stdlib.h>
#include <string.h>

#include "algo/sha/sha256-hash.h"
#include "deqode-gate.h"

// Per-thread state; nothing is shared between threads.
typedef struct
{
   bool ready;
   dq_params params;
   dq_tuning tuning;
   dq_code code;
   struct dq_solver *solver;
   bool have_code;
   unsigned char code_seed[32];   // the hashPrevBlock the current code was derived from

   // Aligned for v128_bswap32_80.
   unsigned char prefix[ DQ_PREFIX_SIZE ] __attribute__((aligned(64)));
   uint64_t *syndromes;           // one batch, DQ_WORDS(r) words per nonce
   size_t syndrome_cap;
   unsigned char *witness;
   size_t witness_bytes;

   // End of this thread's nonce range. A batch may run past max_nonce but not past this.
   uint32_t end_nonce;

   // SHA-256 context with the prefix's first 64 bytes absorbed.
   sha256_context midstate;
} dq_thread;

static __thread dq_thread *tls = NULL;

/** SHA-256 of the 80 byte prefix and an optional witness, starting from the midstate. */
static void hash_from_midstate( const dq_thread *t, const unsigned char *tail,
                                size_t tail_len, void *out )
{
   sha256_context ctx = t->midstate;
   sha256_update( &ctx, t->prefix + 64, 16 );
   if ( tail_len ) sha256_update( &ctx, tail, tail_len );
   sha256_final( &ctx, out );
}

typedef struct
{
   dq_thread *t;
   struct work *work;
   struct thr_info *mythr;
   uint32_t first_nonce;          // nonce of the batch's first syndrome
   const uint32_t *target;
   bool solved;
} dq_attempt;

/** Try a witness against the acceptance target; submit it if it wins. */
static bool on_witness( void *context, size_t index, const uint64_t *x, const uint64_t *z )
{
   dq_attempt *a = (dq_attempt*) context;
   dq_thread *t = a->t;
   uint32_t hash[8] __attribute__((aligned(32)));
   const uint32_t nonce = a->first_nonce + (uint32_t)index;

   t->prefix[ DQ_NONCE_OFFSET     ] = (unsigned char)( nonce       );
   t->prefix[ DQ_NONCE_OFFSET + 1 ] = (unsigned char)( nonce >>  8 );
   t->prefix[ DQ_NONCE_OFFSET + 2 ] = (unsigned char)( nonce >> 16 );
   t->prefix[ DQ_NONCE_OFFSET + 3 ] = (unsigned char)( nonce >> 24 );
   dq_encode_witness( x, z, t->params.n, t->witness );
   // SHA256(prefix80 || witness), without the length byte the serialization adds.
   hash_from_midstate( t, t->witness, t->witness_bytes, hash );

   if ( !valid_hash( hash, a->target ) ) return false;
   if ( opt_benchmark ) return false;

   a->work->data[19] = bswap_32( nonce );
   a->work->taskdata_len = (uint16_t) t->witness_bytes;
   memcpy( a->work->taskdata, t->witness, t->witness_bytes );
   submit_solution( a->work, hash, a->mythr );
   a->solved = true;
   return true;
}

bool deqode_thread_init( int thr_id )
{
   tls = (dq_thread*) aligned_alloc( 64, ( sizeof(dq_thread) + 63 ) & ~(size_t)63 );
   if ( tls ) memset( tls, 0, ( sizeof(dq_thread) + 63 ) & ~(size_t)63 );
   if ( !tls )
   {
      applog( LOG_ERR, "deqode: thread %d out of memory", thr_id );
      return false;
   }
   return true;
}

void deqode_get_new_work( struct work *work, struct work *g_work, int thr_id,
                          uint32_t *end_nonce_ptr )
{
   std_get_new_work( work, g_work, thr_id, end_nonce_ptr );
   if ( tls ) tls->end_nonce = *end_nonce_ptr;
}

/** Back off for a few seconds, or until new work arrives. */
static void wait_for_work( int thr_id )
{
   for ( int i = 0; i < 5 && !work_restart[ thr_id ].restart; i++ ) sleep( 1 );
}

/** Room for `count` syndromes. */
static bool reserve_syndromes( dq_thread *t, size_t count )
{
   if ( count <= t->syndrome_cap ) return true;
   free( t->syndromes );
   t->syndromes = (uint64_t*) malloc( count * DQ_WORDS( t->params.r ) * 8 );
   t->syndrome_cap = t->syndromes ? count : 0;
   return t->syndromes != NULL;
}

/** Rebuild whatever the new sizing invalidated. Returns false if the sizing is unusable. */
static bool ensure_sizing( dq_thread *t, const dq_params *params, int thr_id )
{
   if ( t->ready && !memcmp( &t->params, params, sizeof(*params) ) ) return true;

   if ( !dq_params_are_sane( params ) )
   {
      applog( LOG_ERR, "deqode: unusable sizing n=%d r=%d w=%d gate_multiplier=%d",
              params->n, params->r, params->w, params->gate_multiplier );
      return false;
   }

   dq_solver_free( t->solver );  t->solver = NULL;
   dq_code_free( &t->code );
   free( t->syndromes ); t->syndromes = NULL; t->syndrome_cap = 0;
   free( t->witness );   t->witness = NULL;
   t->ready = t->have_code = false;

   t->params  = *params;
   t->tuning  = deqode_opt_tuning;
   dq_solver_autotune( &t->tuning, params );

   if ( !dq_code_alloc( &t->code, params->n, params->r ) ) return false;
   t->witness_bytes = dq_witness_bytes( params->n );
   t->witness = (unsigned char*) calloc( t->witness_bytes, 1 );
   // A different seed per thread.
   t->solver = dq_solver_new( params, &t->tuning,
                              0x9E3779B97F4A7C15ULL * (uint64_t)( thr_id + 1 ) + 1 );
   if ( !t->witness || !t->solver )
   {
      applog( LOG_ERR, "deqode: thread %d could not allocate solver state", thr_id );
      return false;
   }

   if ( thr_id == 0 )
   {
      char batch[24];
      if ( t->tuning.batch > 0 ) snprintf( batch, sizeof batch, "%d", t->tuning.batch );
      else                       snprintf( batch, sizeof batch, "auto" );
      applog( LOG_INFO, "deqode n=%d r=%d w=%d gates=%d/qubit, witness %zu bytes, "
                        "Stern p=%d+%d batch=%s",
              params->n, params->r, params->w, params->gate_multiplier,
              t->witness_bytes, t->tuning.left_order, t->tuning.right_order, batch );
   }
   t->ready = true;
   return true;
}

int scanhash_deqode( struct work *work, uint32_t max_nonce, uint64_t *hashes_done,
                     struct thr_info *mythr )
{
   dq_thread *t = tls;
   const int thr_id = mythr->id;
   uint32_t *pdata = work->data;
   uint32_t n = pdata[19];

   *hashes_done = 0;
   if ( !t ) return 0;

   dq_params params = { (int)work->task_params[0], (int)work->task_params[1],
                        (int)work->task_params[2], (int)work->task_params[3] };
   if ( opt_benchmark ) params = deqode_opt_params;
   if ( !ensure_sizing( t, &params, thr_id ) )
   {
      wait_for_work( thr_id );
      return 0;
   }

   // The header in serialization order.
   v128_bswap32_80( t->prefix, pdata );

   // The code is seeded by hashPrevBlock.
   if ( !t->have_code
     || memcmp( t->code_seed, t->prefix + DQ_PREVHASH_OFFSET, 32 ) )
   {
      memcpy( t->code_seed, t->prefix + DQ_PREVHASH_OFFSET, 32 );
      dq_derive_columns( &t->code, t->code_seed, 32, t->params.gate_multiplier );
      t->have_code = true;
   }

   // Only the nonce changes within a work item, so the first 64 bytes are hashed once.
   sha256_ctx_init( &t->midstate );
   sha256_update( &t->midstate, t->prefix, 64 );

   // Witnesses are checked against the template's acceptance target, not work->target.
   dq_attempt attempt = { t, work, mythr, 0, work->task_target, false };
   const uint32_t last_nonce = max_nonce - 1;
   double searched = 0;
   const uint32_t range_end = t->end_nonce > n ? t->end_nonce : max_nonce;

   while ( n < last_nonce && !work_restart[ thr_id ].restart )
   {
      // One information set per batch of nonces; each nonce is tried once.
      if ( !dq_solver_new_information_set( t->solver, &t->code ) )
      {
         applog( LOG_WARNING, "deqode: no information set at height %d (degenerate code "
                              "or out of memory), waiting for work", work->height );
         wait_for_work( thr_id );
         break;
      }

      // Not cut at max_nonce: a scanhash call runs at least one whole batch.
      size_t batch = dq_solver_batch_size( t->solver );
      if ( batch > (size_t)( range_end - n ) ) batch = (size_t)( range_end - n );
      if ( !batch ) break;
      if ( !reserve_syndromes( t, batch ) )
      {
         applog( LOG_ERR, "deqode: thread %d cannot hold a batch of %zu syndromes",
                 thr_id, batch );
         wait_for_work( thr_id );
         break;
      }

      const size_t rw = DQ_WORDS( t->params.r );
      size_t derived = 0;
      for ( ; derived < batch; derived++ )
      {
         if ( ( derived & 0xFFF ) == 0 && work_restart[ thr_id ].restart ) break;
         const uint32_t nonce = n + (uint32_t)derived;
         unsigned char seed[32];
         t->prefix[ DQ_NONCE_OFFSET     ] = (unsigned char)( nonce       );
         t->prefix[ DQ_NONCE_OFFSET + 1 ] = (unsigned char)( nonce >>  8 );
         t->prefix[ DQ_NONCE_OFFSET + 2 ] = (unsigned char)( nonce >> 16 );
         t->prefix[ DQ_NONCE_OFFSET + 3 ] = (unsigned char)( nonce >> 24 );
         hash_from_midstate( t, NULL, 0, seed );
         dq_derive_syndrome( t->params.r, seed, 32, t->syndromes + derived * rw );
      }
      if ( derived < batch ) break;

      attempt.first_nonce = n;
      double progress = 0;
      if ( !dq_solver_join_batch( t->solver, t->syndromes, batch, on_witness, &attempt,
                                  &work_restart[ thr_id ].restart, &progress )
        && progress == 0 && !work_restart[ thr_id ].restart )
      {
         applog( LOG_ERR, "deqode: thread %d cannot hold a batch of %zu syndromes",
                 thr_id, batch );
         wait_for_work( thr_id );
         break;
      }
      searched += progress * (double)batch;
      n += (uint32_t)batch;
      if ( attempt.solved ) goto done;
   }

done:
   // Counted in nonces (decoding attempts), not hashes.
   *hashes_done = (uint64_t)( searched + 0.5 );
   pdata[19] = n;
   return 0;
}
