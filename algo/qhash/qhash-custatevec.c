#include <complex.h>
#include <cuComplex.h>
#include <cuda_runtime_api.h>
#include <custatevec.h>
#include <math_constants.h>
#include <stdbool.h>

#include "qhash-gate.h"

#define HANDLE_CUSTATEVEC_ERROR(x)                                                              \
    {                                                                                           \
        const custatevecStatus_t err = x;                                                       \
        if unlikely (err != CUSTATEVEC_STATUS_SUCCESS)                                          \
        {                                                                                       \
            fprintf(stdout, "Error: %s in line %d\n", custatevecGetErrorString(err), __LINE__); \
        }                                                                                       \
    };

#define HANDLE_CUDA_ERROR(x)                                                              \
    {                                                                                     \
        const cudaError_t err = x;                                                        \
        if unlikely (err != cudaSuccess)                                                  \
        {                                                                                 \
            fprintf(stdout, "Error: %s in line %d\n", cudaGetErrorString(err), __LINE__); \
        }                                                                                 \
    };


// Calculation constants
static const complex float matrixX[] = {0.0f, 1.0f, 1.0f, 0.0f};
static const custatevecPauli_t pauliY[] = {CUSTATEVEC_PAULI_Y};
static const custatevecPauli_t pauliZ[] = {CUSTATEVEC_PAULI_Z};
static const custatevecPauli_t *const pauliExpectations[NUM_QUBITS] = {
    pauliZ, pauliZ, pauliZ, pauliZ, pauliZ, pauliZ, pauliZ, pauliZ,
    pauliZ, pauliZ, pauliZ, pauliZ, pauliZ, pauliZ, pauliZ, pauliZ};
static const int32_t basisBits[NUM_QUBITS] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
static const int32_t *const basisBitsArr[NUM_QUBITS] = {
    basisBits, basisBits + 1, basisBits + 2, basisBits + 3, basisBits + 4, basisBits + 5,
    basisBits + 6, basisBits + 7, basisBits + 8, basisBits + 9, basisBits + 10, basisBits + 11,
    basisBits + 12, basisBits + 13, basisBits + 14, basisBits + 15};
static const uint32_t nBasisBits[NUM_QUBITS] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

// GPU variables
static __thread custatevecHandle_t handle;
static __thread cuComplex *dStateVec = NULL;
static __thread size_t extraSize = 0;
static __thread void *extra = NULL;
bool qhash_thread_init(int)
{
    custatevecStatus_t custatevecErr;
    cudaError_t cudaErr;
    custatevecErr = custatevecCreate(&handle);
    if (custatevecErr != CUSTATEVEC_STATUS_SUCCESS)
        return false;
    const size_t stateVecSizeBytes = (1 << NUM_QUBITS) * sizeof(cuComplex);
    cudaErr = cudaMalloc((void **)&dStateVec, stateVecSizeBytes);
    if (cudaErr != cudaSuccess)
        return false;
    custatevecErr = custatevecApplyMatrixGetWorkspaceSize(handle, CUDA_C_32F, NUM_QUBITS,
                                                          matrixX, CUDA_C_32F,
                                                          CUSTATEVEC_MATRIX_LAYOUT_ROW, 0,
                                                          1, 1, CUSTATEVEC_COMPUTE_DEFAULT,
                                                          &extraSize);
    if (custatevecErr != CUSTATEVEC_STATUS_SUCCESS)
        return false;
    if (extraSize)
    {
        cudaErr = cudaMalloc(&extra, extraSize);
        if (cudaErr != cudaSuccess)
            return false;
    }
    return true;
}

static void get_expectations(double expectations[NUM_QUBITS])
{
    HANDLE_CUSTATEVEC_ERROR(custatevecComputeExpectationsOnPauliBasis(
        handle, dStateVec, CUDA_C_32F, NUM_QUBITS, expectations,
        (const custatevecPauli_t **)pauliExpectations, NUM_QUBITS, (const int32_t **)basisBitsArr,
        nBasisBits));
}

static void main_circuit(const unsigned char data[2 * SHA256_BLOCK_SIZE], int angle_offset)
{
    for (size_t l = 0; l < NUM_LAYERS; ++l)
    {
        for (size_t i = 0; i < NUM_QUBITS; ++i)
        {
            const int32_t target = i;
            // RY gates
            HANDLE_CUSTATEVEC_ERROR(custatevecApplyPauliRotation(
                handle, dStateVec, CUDA_C_32F, NUM_QUBITS,
                -(2 * data[(2 * l * NUM_QUBITS + i) % (2 * SHA256_BLOCK_SIZE)] + angle_offset) *
                    CUDART_PI / 32,
                pauliY, &target, 1, NULL, NULL, 0));
            // RZ gates
            HANDLE_CUSTATEVEC_ERROR(custatevecApplyPauliRotation(
                handle, dStateVec, CUDA_C_32F, NUM_QUBITS,
                -(2 * data[((2 * l + 1) * NUM_QUBITS + i) % (2 * SHA256_BLOCK_SIZE)] + angle_offset) *
                    CUDART_PI / 32,
                pauliZ, &target, 1, NULL, NULL, 0));
        }
        for (size_t i = 0; i < NUM_QUBITS - 1; ++i)
        {
            const int32_t control = i;
            const int32_t target = control + 1;

            HANDLE_CUSTATEVEC_ERROR(custatevecApplyMatrix(handle, dStateVec, CUDA_C_32F, NUM_QUBITS,
                                                          matrixX, CUDA_C_32F,
                                                          CUSTATEVEC_MATRIX_LAYOUT_ROW, 0, &target,
                                                          1, &control, NULL, 1,
                                                          CUSTATEVEC_COMPUTE_DEFAULT, extra,
                                                          extraSize));
        }
    }
}

void run_simulation(const unsigned char data[2 * SHA256_BLOCK_SIZE], int angle_offset,
                    double expectations[NUM_QUBITS])
{
    HANDLE_CUSTATEVEC_ERROR(custatevecInitializeStateVector(handle, dStateVec, CUDA_C_32F,
                                                            NUM_QUBITS,
                                                            CUSTATEVEC_STATE_VECTOR_TYPE_ZERO));
    main_circuit(data, angle_offset);
    get_expectations(expectations);
}