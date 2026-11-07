/**
 * deqode algo gate: command line options, reading the task fields from a block template
 * and building the submitblock payload with the witness in the header.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "deqode-gate.h"

// Command line options. Zero means not given, or choose automatically.
dq_params deqode_opt_params = { 0, 0, 0, 0 };
dq_tuning deqode_opt_tuning = { 0, 0, 0, 0, 0 };

//////////////////////////////////////////////////////////////////////////////////////////
// Command line
//////////////////////////////////////////////////////////////////////////////////////////

bool deqode_parse_sizing( const char *arg )
{
   dq_params p = { 0, 0, 0, 24 };
   int end = 0;
   if ( sscanf( arg, "%d:%d:%d%n:%d%n", &p.n, &p.r, &p.w, &end, &p.gate_multiplier, &end ) < 3
     || arg[end] ) return false;
   if ( !dq_params_are_sane( &p ) ) return false;
   deqode_opt_params = p;
   return true;
}

bool deqode_parse_orders( const char *arg )
{
   int left = 0, right = 0, end = 0;
   if ( sscanf( arg, "%d:%d%n", &left, &right, &end ) != 2 || arg[end] ) return false;
   if ( left < 0 || right < 0 ) return false;
   deqode_opt_tuning.left_order  = left;
   deqode_opt_tuning.right_order = right;
   return true;
}

//////////////////////////////////////////////////////////////////////////////////////////
// getblocktemplate
//////////////////////////////////////////////////////////////////////////////////////////

static bool json_int_field( const json_t *obj, const char *name, int *out )
{
   const json_t *v = json_object_get( obj, name );
   if ( !v || !json_is_integer( v ) ) return false;
   const json_int_t value = json_integer_value( v );
   if ( value < INT_MIN || value > INT_MAX ) return false;
   *out = (int) value;
   return true;
}

/** Read the task 1 fields of a getblocktemplate response. */
static bool deqode_gbt_decode_task( const json_t *val, struct work *work )
{
   const json_t *task = json_object_get( val, "task" );
   if ( !task || !json_is_object( task ) )
   {
      applog( LOG_ERR, "deqode: the template has no task object. This node predates "
                       "SuperDense Consensus, or does not know task %u.",
              (unsigned) DEQODE_TASK_ID );
      return false;
   }

   const json_t *unsupported = json_object_get( task, "unsupported" );
   if ( unsupported && json_is_true( unsupported ) )
   {
      applog( LOG_ERR, "deqode: the node reports task %u as unsupported, so it cannot be "
                       "mined from outside it", (unsigned) DEQODE_TASK_ID );
      return false;
   }

   int id = -1;
   if ( !json_int_field( task, "id", &id ) || id != DEQODE_TASK_ID )
   {
      applog( LOG_ERR, "deqode: the template is for task %d, not task %u",
              id, (unsigned) DEQODE_TASK_ID );
      return false;
   }

   // A version carries a task id in bits 20-27 only when its top three bits are 010.
   const uint32_t version = bswap_32( work->data[0] );
   const uint8_t header_task = ( version & 0xE0000000 ) == 0x40000000
                             ? (uint8_t)( ( version >> 20 ) & 0xFF ) : 0;
   if ( header_task != DEQODE_TASK_ID )
   {
      applog( LOG_ERR, "deqode: the template's version 0x%08x carries task %u, not %u",
              version, (unsigned) header_task, (unsigned) DEQODE_TASK_ID );
      return false;
   }

   const json_t *sizing = json_object_get( task, "deqode" );
   if ( !sizing || !json_is_object( sizing ) )
   {
      applog( LOG_ERR, "deqode: the task object carries no deqode sizing" );
      return false;
   }

   dq_params params = { 0, 0, 0, 0 };
   if ( !json_int_field( sizing, "n", &params.n )
     || !json_int_field( sizing, "r", &params.r )
     || !json_int_field( sizing, "gate_multiplier", &params.gate_multiplier ) )
   {
      applog( LOG_ERR, "deqode: the sizing is missing n, r or gate_multiplier" );
      return false;
   }

   // w and the acceptance target depend on the template's nBits.
   int expected_witnesses = 0;
   if ( !json_int_field( sizing, "w", &params.w )
     || !json_int_field( sizing, "expected_witnesses", &expected_witnesses ) )
   {
      applog( LOG_ERR, "deqode: the template has no weight bound" );
      return false;
   }

   if ( !dq_params_are_sane( &params ) )
   {
      applog( LOG_ERR, "deqode: the node published an unusable sizing n=%d r=%d w=%d gm=%d",
              params.n, params.r, params.w, params.gate_multiplier );
      return false;
   }

   unsigned char raw[32];
   const char *target_hex =
      json_string_value( json_object_get( sizing, "acceptance_target" ) );
   if ( !target_hex || strlen( target_hex ) != 2 * sizeof(raw)
     || !hex2bin( raw, target_hex, sizeof(raw) ) )
   {
      applog( LOG_ERR, "deqode: the template carries no valid acceptance target" );
      return false;
   }
   // Byte order as for "target" in cpu-miner.c.
   for ( int i = 0; i < 32; i++ )
      ( (unsigned char*) work->task_target )[i] = raw[ 31 - i ];

   int data_bytes = 0;
   if ( json_object_get( task, "data_bytes" )
     && ( !json_int_field( task, "data_bytes", &data_bytes )
       || (size_t)data_bytes != dq_witness_bytes( params.n ) ) )
   {
      applog( LOG_ERR, "deqode: the template has invalid data_bytes for n=%d; expected %zu",
              params.n, dq_witness_bytes( params.n ) );
      return false;
   }

   work->task_params[0] = (uint32_t) params.n;
   work->task_params[1] = (uint32_t) params.r;
   work->task_params[2] = (uint32_t) params.w;
   work->task_params[3] = (uint32_t) params.gate_multiplier;

   if ( opt_debug )
      applog( LOG_DEBUG, "deqode: template wants w=%d with %d expected witnesses",
              params.w, expected_witnesses );
   return true;
}

//////////////////////////////////////////////////////////////////////////////////////////
// submitblock
//////////////////////////////////////////////////////////////////////////////////////////

/**
 * submitblock with the witness in the header: the 80 byte header, the witness with a
 * CompactSize length, then the transactions.
 */
static char *deqode_malloc_txs_request( struct work *work )
{
   char *req, *params = NULL;
   unsigned char header[ DQ_PREFIX_SIZE + 3 + MAX_TASK_DATA_SIZE ];
   size_t at = 0;

   // Into a copy: a failed submission is retried with the same work.
   for ( int i = 0; i < DQ_PREFIX_SIZE / 4; i++ )
      be32enc( header + 4 * i, work->data[i] );
   at = DQ_PREFIX_SIZE;

   const size_t len = work->taskdata_len;
   if ( len < 253 )
      header[ at++ ] = (unsigned char) len;
   else
   {
      header[ at++ ] = 0xFD;
      header[ at++ ] = (unsigned char)( len       );
      header[ at++ ] = (unsigned char)( len >>  8 );
   }
   memcpy( header + at, work->taskdata, len );
   at += len;

   char *hex = (char*) malloc( 2 * at + 1 );
   if ( !hex ) return NULL;
   bin2hex( hex, header, at );

   if ( work->workid )
   {
      json_t *v = json_object();
      json_object_set_new( v, "workid", json_string( work->workid ) );
      params = json_dumps( v, 0 );
      json_decref( v );
      req = (char*) malloc( 128 + strlen( hex ) + strlen( work->txs )
                                + ( params ? strlen( params ) : 0 ) );
      if ( req )
         sprintf( req,
            "{\"method\": \"submitblock\", \"params\": [\"%s%s\", %s], \"id\":4}\r\n",
            hex, work->txs, params );
      free( params );
   }
   else
   {
      req = (char*) malloc( 128 + strlen( hex ) + strlen( work->txs ) );
      if ( req )
         sprintf( req,
            "{\"method\": \"submitblock\", \"params\": [\"%s%s\"], \"id\":4}\r\n",
            hex, work->txs );
   }
   free( hex );
   return req;
}

//////////////////////////////////////////////////////////////////////////////////////////
// Registration
//////////////////////////////////////////////////////////////////////////////////////////

bool register_deqode_algo( algo_gate_t *gate )
{
   gate->task_id            = DEQODE_TASK_ID;
   gate->scanhash           = &scanhash_deqode;
   gate->miner_thread_init  = &deqode_thread_init;
   gate->get_new_work       = &deqode_get_new_work;
   gate->gbt_decode_task    = &deqode_gbt_decode_task;
   gate->malloc_txs_request = &deqode_malloc_txs_request;

   if ( opt_benchmark && !deqode_opt_params.n )
   {
      applog( LOG_ERR, "deqode: --benchmark needs --deqode-sizing=n:r:w[:gate_multiplier], "
                       "for example --deqode-sizing=128:82:15:24" );
      return false;
   }
   if ( !opt_benchmark && deqode_opt_params.n )
   {
      applog( LOG_ERR, "deqode: --deqode-sizing is only for --benchmark; when mining, the "
                       "block template supplies the sizing" );
      return false;
   }
   return true;
}
