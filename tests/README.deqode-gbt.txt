Offline deqode getblocktemplate integer-decoding regression tests

Purpose
-------
Exercise the production static json_int_field and deqode_gbt_decode_task
functions against valid, missing, wrongly typed and out-of-int-range values.
Large JSON integers must be rejected before conversion to int. An optional
data_bytes value may be absent, but a present invalid value must be rejected.

Run on Linux with GCC (or a compatible compiler and ELF linker), Autotools,
the normal project header dependencies, pkg-config and development Jansson:

  sh autogen.sh
  sh configure --enable-tasks=deqode --disable-assembly
  sh tests/run-deqode-gbt.sh

The final command builds and runs a small offline test executable in a temporary
directory. It does not build or start the miner, connect to a node, or make RPC
requests. CUDA and cuQuantum are not required for these parser tests.

Test design and scope
---------------------
The test translation unit includes deqode-gate.c directly to exercise its real
static functions. deqode.c supplies the real dq_params_are_sane implementation.
Unused registration, mining and derivation functions are discarded at link time.
The test provides a no-op applog and an explicitly limited hex2bin test double
for the fixed all-zero acceptance target. It does not validate target decoding,
proof-of-work solving, block submission, or live mining behavior.

The regression fixtures require a 32-bit int and a wider json_int_t. The runner
uses -ffunction-sections/-fdata-sections and GNU-compatible --gc-sections; this
is not a portable replacement for a complete platform build or full miner test.

Expected result: 35 checks pass. The unpatched source at revision
62b3d51395ce2e2bda0d4ca56368842f4c11b6a9 fails 14 of these checks on Linux x86_64
with GCC and Jansson 2.14. In particular, 4294967297 must not become task id 1,
4294967424 must not become n=128, and 4294967328 must not become data_bytes=32.
