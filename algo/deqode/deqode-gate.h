#ifndef DEQODE_GATE_H__
#define DEQODE_GATE_H__

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "algo-gate-api.h"

/**
 * deqode, proof-of-work task 1. A solution is a Pauli error of weight at most w whose
 * syndrome under the block's random stabilizer code is the block's syndrome, carried in
 * the header after the nonce. See doc/deqode.md in the node repository.
 */

#define DEQODE_TASK_ID 1

// The node's MAX_DEQODE_QUBITS.
#define DQ_MAX_QUBITS 1024

// Number of 64 bit words needed to hold `bits` bits.
#define DQ_WORDS(bits) ( ( (size_t)(bits) + 63 ) / 64 )

// Byte offsets into the 80 byte header prefix.
#define DQ_PREFIX_SIZE      80
#define DQ_PREVHASH_OFFSET   4
#define DQ_NONCE_OFFSET     76

//////////////////////////////////////////////////////////////////////////////////////////
// Instance derivation: deqode.c
//////////////////////////////////////////////////////////////////////////////////////////

/** An instance's sizing, as the template publishes it. */
typedef struct
{
   int n;                 // qubits
   int r;                 // stabilizer generators (syndrome bits)
   int w;                 // largest Pauli weight a solution may have
   int gate_multiplier;   // Clifford gates per qubit in the code circuit
} dq_params;

/**
 * A stabilizer code in column-major form: col_x[q] and col_z[q] are r-bit vectors holding
 * each stabilizer's X and Z bit on qubit q.
 */
typedef struct
{
   int n;
   int r;
   size_t rw;         // DQ_WORDS(r), the stride of both column arrays
   uint64_t *col_x;   // n * rw words
   uint64_t *col_z;   // n * rw words
} dq_code;

/** One gate of the code circuit. */
typedef struct
{
   int gate;     // 0 = H on `first`, 1 = S on `first`, 2 = CNOT control `first` target `second`
   int first;
   int second;
} dq_gate_draw;

/** True if the sizing is one this implementation can derive an instance for. */
bool dq_params_are_sane( const dq_params *params );

/** Allocate a code's column arrays. Returns false on a bad sizing or out of memory. */
bool dq_code_alloc( dq_code *code, int n, int r );
void dq_code_free( dq_code *code );

/** Run the Clifford circuit for a seed, yielding the code's columns. */
void dq_derive_columns( dq_code *code, const void *seed, size_t seed_len,
                        int gate_multiplier );

/** The r-bit syndrome a seed determines. `out` holds DQ_WORDS(r) words. */
void dq_derive_syndrome( int r, const void *seed, size_t seed_len, uint64_t *out );

/**
 * The syndrome of a single-qubit error. Column q is an X on qubit q, column n+q a Z on
 * qubit q; an X on q flips the rows in col_z[q] and a Z those in col_x[q].
 */
static inline const uint64_t *dq_column_syndrome( const dq_code *code, int column )
{
   return column < code->n ? code->col_z + (size_t)column * code->rw
                           : code->col_x + (size_t)( column - code->n ) * code->rw;
}

/** The canonical 2 * ceil(n/8) byte encoding: the X mask, then the Z mask. */
void dq_encode_witness( const uint64_t *x, const uint64_t *z, int n, unsigned char *out );

/** 2 * ceil(n/8): the exact witness length this sizing expects. */
static inline size_t dq_witness_bytes( int n ) { return 2 * (size_t)( ( n + 7 ) / 8 ); }

//////////////////////////////////////////////////////////////////////////////////////////
// Solver: deqode-stern.c, a Stern decoder. It knows nothing about headers or the chain.
//////////////////////////////////////////////////////////////////////////////////////////

/** Solver tuning, from the command line or dq_solver_autotune. */
typedef struct
{
   int left_order;        // Stern's p, split between the two lists
   int right_order;
   int window_bits;       // Stern's l, the collision window
   size_t memory_bytes;   // per thread budget for the pattern lists
   int batch;             // syndromes run against one information set; 0 = per set
} dq_tuning;

/** Fill in whatever the caller left at zero, from the sizing. */
void dq_solver_autotune( dq_tuning *tuning, const dq_params *params );

struct dq_solver;

/**
 * Called with each witness of weight at most w whose syndrome is a target, as packed masks
 * of DQ_WORDS(n) words; `index` is the position of that syndrome in the batch. Returns
 * true to stop.
 */
typedef bool ( *dq_batch_witness_fn )( void *context, size_t index, const uint64_t *x,
                                       const uint64_t *z );

/** Allocate solver state for one thread. Returns NULL on a bad sizing or out of memory. */
struct dq_solver *dq_solver_new( const dq_params *params, const dq_tuning *tuning,
                                 uint64_t seed );
void dq_solver_free( struct dq_solver *solver );

/**
 * Draw a fresh information set for `code` and build everything that does not depend on
 * the syndrome (the elimination, the two pattern lists and the hash table over the left
 * one), so it can be reused for many syndromes. Returns false if the code is degenerate
 * or the solver's tables cannot be allocated.
 */
bool dq_solver_new_information_set( struct dq_solver *solver, const dq_code *code );

/**
 * Run `count` syndromes (DQ_WORDS(r) words each) against the current information set,
 * reporting every witness found. Returns true if a callback asked to stop. `abort_flag`,
 * if not NULL, is polled during the pass. `progress`, if not NULL, receives the fraction
 * of the batch done when the call returns; it is 0 if the batch could not be allocated.
 */
bool dq_solver_join_batch( struct dq_solver *solver, const uint64_t *syndromes,
                           size_t count, dq_batch_witness_fn on_witness, void *context,
                           const volatile uint8_t *abort_flag, double *progress );

/** Syndromes per batch: --deqode-batch, or 2^(l+1) for window l, at most 2^20. */
size_t dq_solver_batch_size( const struct dq_solver *solver );

//////////////////////////////////////////////////////////////////////////////////////////
// Algo gate: deqode-gate.c and deqode-mine.c
//////////////////////////////////////////////////////////////////////////////////////////

bool register_deqode_algo( algo_gate_t *gate );

bool deqode_thread_init( int thr_id );

/** std_get_new_work, also recording where this thread's nonce range ends. */
void deqode_get_new_work( struct work *work, struct work *g_work, int thr_id,
                          uint32_t *end_nonce_ptr );

int scanhash_deqode( struct work *work, uint32_t max_nonce, uint64_t *hashes_done,
                     struct thr_info *mythr );

/** Command line options. The sizing is only for --benchmark, which requires it. */
extern dq_params  deqode_opt_params;    // --deqode-sizing=n:r:w:g
extern dq_tuning  deqode_opt_tuning;    // --deqode-p, --deqode-l, --deqode-mem, --deqode-batch

/** Parse one of the --deqode-* option values. Returns false on a malformed value. */
bool deqode_parse_sizing( const char *arg );
bool deqode_parse_orders( const char *arg );

#endif
