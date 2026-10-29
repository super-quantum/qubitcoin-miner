#ifndef QHASH_GATE_H__
#define QHASH_GATE_H__

#include "algo-gate-api.h"

#define NUM_QUBITS 16
#define NUM_LAYERS 2
#define SHA256_BLOCK_SIZE 32
#define INPUT_SIZE 80

// Offset of nTime in the block header.
#define NTIME_OFFSET 68

bool register_qhash_algo(algo_gate_t *gate);
bool register_qhash_cpu_algo(algo_gate_t *gate);

typedef void (*qhash_simulation_fn)(const unsigned char data[2 * SHA256_BLOCK_SIZE],
                                    int angle_offset, double expectations[NUM_QUBITS]);

/** Hashes an 80-byte header, running the circuit with the given simulator. */
int qhash_digest(void *output, const void *input, qhash_simulation_fn simulate);

// Simulator-specific functions

/**
 * Initializes all necessary simulation variables for each mining thread.
 * @param thr_id id of the thread being initialized
 * @return       true if the initialization was successful and false otherwise
 */
bool qhash_thread_init(int thr_id);

/**
 * Runs the quantum simulation and computes qubit expectations.
 * @param data         input hash for the circuit (each byte represents a single nibble for angle parametrization)
 * @param angle_offset 0 or 1; the rotation angle for a data byte d is -(2 * d + angle_offset) * pi / 32
 * @param expectations output parameter for Z-basis expectations for each circuit qubit
 */
void run_simulation(const unsigned char data[2 * SHA256_BLOCK_SIZE], int angle_offset,
                    double expectations[NUM_QUBITS]);

/** CPU counterpart of qhash_thread_init. */
bool qhash_cpu_thread_init(int thr_id);

/** CPU counterpart of run_simulation, in double precision like the node. */
void run_simulation_cpu(const unsigned char data[2 * SHA256_BLOCK_SIZE], int angle_offset,
                        double expectations[NUM_QUBITS]);

#endif
