#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "qhash-gate.h"

#define STATE_SIZE ((size_t)1 << NUM_QUBITS)

// Amplitude i is the basis state whose bit q is the value of qubit q, as in cuStateVec.
static __thread double complex *state = NULL;

bool qhash_cpu_thread_init(int thr_id)
{
    state = malloc(STATE_SIZE * sizeof(*state));
    return state != NULL;
}

// a * exp(i * theta), given cos(theta) and sin(theta).
static inline double complex rotate_phase(double complex a, double c, double s)
{
    return CMPLX(c * creal(a) - s * cimag(a), s * creal(a) + c * cimag(a));
}

// The rotations apply exp(i * theta * P), like custatevecApplyPauliRotation.
static void apply_rotation_y(int qubit, double theta)
{
    const double c = cos(theta), s = sin(theta);
    const size_t bit = (size_t)1 << qubit;
    for (size_t hi = 0; hi < STATE_SIZE; hi += 2 * bit)
    {
        for (size_t i = hi; i < hi + bit; ++i)
        {
            const double complex a0 = state[i], a1 = state[i + bit];
            state[i] = c * a0 + s * a1;
            state[i + bit] = c * a1 - s * a0;
        }
    }
}

static void apply_rotation_z(int qubit, double theta)
{
    const double c = cos(theta), s = sin(theta);
    const size_t bit = (size_t)1 << qubit;
    for (size_t hi = 0; hi < STATE_SIZE; hi += 2 * bit)
    {
        for (size_t i = hi; i < hi + bit; ++i)
        {
            state[i] = rotate_phase(state[i], c, s);
            state[i + bit] = rotate_phase(state[i + bit], c, -s);
        }
    }
}

// Requires control < target, which holds for the ladder.
static void apply_cnot(int control, int target)
{
    const size_t control_bit = (size_t)1 << control, target_bit = (size_t)1 << target;
    for (size_t hi = 0; hi < STATE_SIZE; hi += 2 * target_bit)
    {
        for (size_t mid = hi + control_bit; mid < hi + target_bit; mid += 2 * control_bit)
        {
            for (size_t i = mid; i < mid + control_bit; ++i)
            {
                const double complex a = state[i];
                state[i] = state[i + target_bit];
                state[i + target_bit] = a;
            }
        }
    }
}

static void get_expectations(double expectations[NUM_QUBITS])
{
    for (size_t q = 0; q < NUM_QUBITS; ++q)
        expectations[q] = 0.0;
    for (size_t i = 0; i < STATE_SIZE; ++i)
    {
        const double p = creal(state[i]) * creal(state[i]) + cimag(state[i]) * cimag(state[i]);
        for (size_t q = 0; q < NUM_QUBITS; ++q)
            expectations[q] += ((i >> q) & 1) ? -p : p;
    }
}

void run_simulation_cpu(const unsigned char data[2 * SHA256_BLOCK_SIZE], int angle_offset,
                        double expectations[NUM_QUBITS])
{
    memset(state, 0, STATE_SIZE * sizeof(*state));
    state[0] = 1.0;

    for (size_t l = 0; l < NUM_LAYERS; ++l)
    {
        for (size_t i = 0; i < NUM_QUBITS; ++i)
        {
            const unsigned char y = data[(2 * l * NUM_QUBITS + i) % (2 * SHA256_BLOCK_SIZE)];
            const unsigned char z = data[((2 * l + 1) * NUM_QUBITS + i) % (2 * SHA256_BLOCK_SIZE)];
            apply_rotation_y(i, -(2 * y + angle_offset) * M_PI / 32);
            apply_rotation_z(i, -(2 * z + angle_offset) * M_PI / 32);
        }
        for (size_t i = 0; i < NUM_QUBITS - 1; ++i)
            apply_cnot(i, i + 1);
    }

    get_expectations(expectations);
}
