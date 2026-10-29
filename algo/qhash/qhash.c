#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "algo/sha/sha256-hash.h"

#include "qhash-gate.h"

#define FIXED_FRACTION int16_t
#define FIXED_INTERMEDIATE int32_t
#define FRACTION_BITS 15

// Number of fixed point expectation bytes appended to the initial hash.
#define EXPECTATION_BYTES (NUM_QUBITS * sizeof(FIXED_FRACTION))

// Activation times of qhash rule changes, as in the node's src/crypto/qhash.cpp.
#define ANGLE_OFFSET_TIME   1758762000  // rotation angles move to odd multiples of pi/32
#define DEGENERATE_ALL_TIME 1753105444  // all expectation bytes zero
#define DEGENERATE_3_4_TIME 1753305380  // three quarters of them zero
#define DEGENERATE_1_4_TIME 1754220531  // one quarter of them zero

static FIXED_FRACTION toFixed(double x)
{
    static_assert(FRACTION_BITS <= (sizeof(FIXED_FRACTION) * 8 - 1));
    const FIXED_INTERMEDIATE fractionMult = 1 << FRACTION_BITS;
    return (x >= 0.0) ? (x * fractionMult + 0.5) : (x * fractionMult - 0.5);
}

#define NIBBLE_MASK (unsigned char)0xF
#define NIBBLE_SIZE 4

static void split_nibbles(const unsigned char input[SHA256_BLOCK_SIZE],
                          unsigned char output[2 * SHA256_BLOCK_SIZE])
{
    for (size_t i = 0; i < SHA256_BLOCK_SIZE; ++i)
    {
        output[2 * i] = (input[i] >> NIBBLE_SIZE) & NIBBLE_MASK;
        output[2 * i + 1] = input[i] & NIBBLE_MASK;
    }
}

int qhash_hash(void *output, const void *input, int thr_id)
{
    const uint32_t ntime = le32dec((const unsigned char *)input + NTIME_OFFSET);

    unsigned char buf[SHA256_BLOCK_SIZE + EXPECTATION_BYTES];
    sha256_full(buf, input, INPUT_SIZE);

    unsigned char nibbles[2 * SHA256_BLOCK_SIZE];
    split_nibbles(buf, nibbles);

    double expectations[NUM_QUBITS];
    run_simulation(nibbles, ntime >= ANGLE_OFFSET_TIME, expectations);

    size_t zeroes = 0;
    for (size_t i = 0; i < NUM_QUBITS; ++i)
    {
        const FIXED_FRACTION fixed = toFixed(expectations[i]);
        const size_t j = SHA256_BLOCK_SIZE + i * sizeof(FIXED_FRACTION);
        for (size_t k = 0; k < sizeof(FIXED_FRACTION); ++k)
        {
            // Little endian representation
            buf[j + k] = fixed >> (k * 8);
            if (buf[j + k] == 0)
                ++zeroes;
        }
    }

    // Too many zero expectation bytes: the hash is all ones, which no target accepts.
    if ((zeroes == EXPECTATION_BYTES && ntime >= DEGENERATE_ALL_TIME) ||
        (zeroes >= EXPECTATION_BYTES * 3 / 4 && ntime >= DEGENERATE_3_4_TIME) ||
        (zeroes >= EXPECTATION_BYTES * 1 / 4 && ntime >= DEGENERATE_1_4_TIME))
    {
        memset(output, 0xFF, SHA256_BLOCK_SIZE);
        return 1;
    }

    sha256_full(output, buf, sizeof buf);
    return 1;
}
