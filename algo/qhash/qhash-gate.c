#include "qhash-gate.h"

#if defined(TASK_QHASH_ENABLED)
static int qhash_hash(void *output, const void *input, int thr_id)
{
    return qhash_digest(output, input, run_simulation);
}

bool register_qhash_algo(algo_gate_t *gate)
{
    gate->hash = &qhash_hash;
    gate->miner_thread_init = &qhash_thread_init;
    // gate->optimizations = ()
    return true;
}
#endif

#if defined(TASK_QHASH_CPU_ENABLED)
static int qhash_cpu_hash(void *output, const void *input, int thr_id)
{
    return qhash_digest(output, input, run_simulation_cpu);
}

bool register_qhash_cpu_algo(algo_gate_t *gate)
{
    gate->hash = &qhash_cpu_hash;
    gate->miner_thread_init = &qhash_cpu_thread_init;
    return true;
}
#endif
