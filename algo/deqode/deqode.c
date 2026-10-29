/**
 * deqode instance derivation and witness encoding, as in the node's
 * src/pow/tasks/deqode.cpp. See doc/deqode.md in the node repository.
 */

#include <stdlib.h>
#include <string.h>

#include "algo/sha/sha256-hash.h"
#include "deqode-gate.h"

//////////////////////////////////////////////////////////////////////////////////////////
// Bit vector helpers
//////////////////////////////////////////////////////////////////////////////////////////

static inline bool test_bit( const uint64_t *bits, int index )
{
   return ( bits[ (size_t)index / 64 ] >> ( index % 64 ) ) & 1;
}

static inline void set_bit( uint64_t *bits, int index )
{
   bits[ (size_t)index / 64 ] |= (uint64_t)1 << ( index % 64 );
}

static inline void xor_into( uint64_t *target, const uint64_t *source, size_t words )
{
   for ( size_t i = 0; i < words; i++ ) target[i] ^= source[i];
}

static const char DQ_CODE_TAG[]     = "QubitcoinDeqode/code";
static const char DQ_SYNDROME_TAG[] = "QubitcoinDeqode/syndrome";

#define DQ_TAG_PREFIX_LEN 64

// SHA-256 in counter mode, read as big endian 64 bit words. Each block is
// SHA256(SHA256(tag) || SHA256(tag) || seed || counter) with a big endian 64 bit counter
// starting at 1.
typedef struct
{
   size_t seed_len;
   uint64_t counter;
   int used;
   unsigned char block[32];
   // tag_prefix || seed || counter, hashed in one call
   unsigned char buf[ DQ_TAG_PREFIX_LEN + DQ_MAX_QUBITS + 8 ];
} dq_hash_stream;

static void hash_stream_init( dq_hash_stream *stream, const char *tag,
                              const void *seed, size_t seed_len )
{
   unsigned char tag_hash[32];
   sha256_full( tag_hash, (const unsigned char*)tag, strlen( tag ) );
   memcpy( stream->buf,      tag_hash, 32 );
   memcpy( stream->buf + 32, tag_hash, 32 );

   stream->seed_len = seed_len;
   stream->counter = 0;
   stream->used = 4;
   memcpy( stream->buf + DQ_TAG_PREFIX_LEN, seed, seed_len );
}

static uint64_t hash_stream_next( dq_hash_stream *stream )
{
   if ( stream->used == 4 )
   {
      stream->counter++;
      stream->used = 0;
   }
   if ( stream->used == 0 )
   {
      unsigned char *counter = stream->buf + DQ_TAG_PREFIX_LEN + stream->seed_len;
      const uint64_t c = stream->counter;
      counter[0] = (unsigned char)( c >> 56 );  counter[1] = (unsigned char)( c >> 48 );
      counter[2] = (unsigned char)( c >> 40 );  counter[3] = (unsigned char)( c >> 32 );
      counter[4] = (unsigned char)( c >> 24 );  counter[5] = (unsigned char)( c >> 16 );
      counter[6] = (unsigned char)( c >>  8 );  counter[7] = (unsigned char)( c       );
      sha256_full( stream->block, stream->buf,
                   DQ_TAG_PREFIX_LEN + stream->seed_len + 8 );
   }
   const unsigned char *p = stream->block + 8 * stream->used;
   stream->used++;
   return   ( (uint64_t)p[0] << 56 ) | ( (uint64_t)p[1] << 48 )
          | ( (uint64_t)p[2] << 40 ) | ( (uint64_t)p[3] << 32 )
          | ( (uint64_t)p[4] << 24 ) | ( (uint64_t)p[5] << 16 )
          | ( (uint64_t)p[6] <<  8 ) | ( (uint64_t)p[7]       );
}

// Bits off a hash stream, low end first. A draw wider than what is left discards the
// remainder and starts a fresh word.
typedef struct
{
   dq_hash_stream stream;
   uint64_t reservoir;
   int available;
} dq_bit_reservoir;

static void reservoir_init( dq_bit_reservoir *res, const char *tag,
                            const void *seed, size_t seed_len )
{
   hash_stream_init( &res->stream, tag, seed, seed_len );
   res->reservoir = 0;
   res->available = 0;
}

static uint64_t reservoir_take( dq_bit_reservoir *res, int width )
{
   if ( width <= 0 ) return 0;
   if ( res->available < width )
   {
      res->reservoir = hash_stream_next( &res->stream );
      res->available = 64;
   }
   const uint64_t value = res->reservoir & ( ( (uint64_t)1 << width ) - 1 );
   res->reservoir >>= width;
   res->available -= width;
   return value;
}

// Uniform on [0, bound), by rejection sampling.
static uint64_t reservoir_below( dq_bit_reservoir *res, uint64_t bound )
{
   int width = 0;
   while ( ( (uint64_t)1 << width ) < bound ) width++;
   for ( ;; )
   {
      const uint64_t value = reservoir_take( res, width );
      if ( value < bound ) return value;
   }
}

// One gate: 0 is H and 1 is S on `first`, 2 is CNOT from `first` to `second`.
static void reservoir_next_gate( dq_bit_reservoir *res, int n, dq_gate_draw *draw )
{
   draw->gate   = (int)reservoir_below( res, 3 );
   draw->first  = (int)reservoir_below( res, (uint64_t)n );
   draw->second = 0;
   if ( draw->gate == 2 )
   {
      draw->second = (int)reservoir_below( res, (uint64_t)n - 1 );
      if ( draw->second >= draw->first ) draw->second++;
   }
}

//////////////////////////////////////////////////////////////////////////////////////////
// The derivation
//////////////////////////////////////////////////////////////////////////////////////////

bool dq_params_are_sane( const dq_params *params )
{
   if ( params->n <= 1 || params->n > DQ_MAX_QUBITS ) return false;
   if ( params->r <= 0 || params->r > params->n )     return false;
   if ( params->w <  0 || params->w > params->n )     return false;
   if ( params->gate_multiplier <= 0 )                return false;
   return true;
}

bool dq_code_alloc( dq_code *code, int n, int r )
{
   memset( code, 0, sizeof(*code) );
   if ( n <= 1 || n > DQ_MAX_QUBITS || r <= 0 || r > n ) return false;
   code->n  = n;
   code->r  = r;
   code->rw = DQ_WORDS( r );
   code->col_x = (uint64_t*) calloc( (size_t)n * code->rw, sizeof(uint64_t) );
   code->col_z = (uint64_t*) calloc( (size_t)n * code->rw, sizeof(uint64_t) );
   if ( !code->col_x || !code->col_z ) { dq_code_free( code ); return false; }
   return true;
}

void dq_code_free( dq_code *code )
{
   free( code->col_x );
   free( code->col_z );
   code->col_x = code->col_z = NULL;
}

void dq_derive_columns( dq_code *code, const void *seed, size_t seed_len,
                        int gate_multiplier )
{
   const size_t rw = code->rw;
   const int n = code->n;

   memset( code->col_x, 0, (size_t)n * rw * sizeof(uint64_t) );
   memset( code->col_z, 0, (size_t)n * rw * sizeof(uint64_t) );

   // Row i starts as Z on qubit i.
   for ( int row = 0; row < code->r; row++ )
      set_bit( code->col_z + (size_t)row * rw, row );

   dq_bit_reservoir res;
   reservoir_init( &res, DQ_CODE_TAG, seed, seed_len );

   const int64_t gates = (int64_t)gate_multiplier * n;
   for ( int64_t i = 0; i < gates; i++ )
   {
      dq_gate_draw draw;
      reservoir_next_gate( &res, n, &draw );
      uint64_t *fx = code->col_x + (size_t)draw.first * rw;
      uint64_t *fz = code->col_z + (size_t)draw.first * rw;
      switch ( draw.gate )
      {
         case 0:  // H
            for ( size_t k = 0; k < rw; k++ )
            {
               const uint64_t t = fx[k];  fx[k] = fz[k];  fz[k] = t;
            }
            break;
         case 1:  // S
            xor_into( fz, fx, rw );
            break;
         default: // CNOT
         {
            uint64_t *sx = code->col_x + (size_t)draw.second * rw;
            uint64_t *sz = code->col_z + (size_t)draw.second * rw;
            xor_into( sx, fx, rw );
            xor_into( fz, sz, rw );
            break;
         }
      }
   }
}

void dq_derive_syndrome( int r, const void *seed, size_t seed_len, uint64_t *out )
{
   dq_hash_stream stream;
   hash_stream_init( &stream, DQ_SYNDROME_TAG, seed, seed_len );
   const size_t words = DQ_WORDS( r );
   for ( size_t i = 0; i < words; i++ ) out[i] = hash_stream_next( &stream );
   if ( r % 64 ) out[ words - 1 ] &= ( (uint64_t)1 << ( r % 64 ) ) - 1;
}

//////////////////////////////////////////////////////////////////////////////////////////
// Acceptance primitives
//////////////////////////////////////////////////////////////////////////////////////////

void dq_encode_witness( const uint64_t *x, const uint64_t *z, int n, unsigned char *out )
{
   const size_t mask_bytes = (size_t)( ( n + 7 ) / 8 );
   memset( out, 0, 2 * mask_bytes );
   for ( int q = 0; q < n; q++ )
   {
      if ( test_bit( x, q ) ) out[ (size_t)q / 8 ]              |= 1 << ( q % 8 );
      if ( test_bit( z, q ) ) out[ mask_bytes + (size_t)q / 8 ] |= 1 << ( q % 8 );
   }
}
