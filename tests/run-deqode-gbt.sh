#!/bin/sh
# Offline parser regression tests. Does not start a miner or make RPC requests.
set -eu
cd "$(dirname "$0")/.."
if [ ! -f cpuminer-config.h ]; then
   echo "Generate project headers first: sh autogen.sh && sh configure --enable-tasks=deqode --disable-assembly" >&2
   exit 2
fi
: "${CC:=cc}"
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/deqode-gbt-test.XXXXXX")
cleanup() {
   rm -f "$test_dir/test-deqode-gbt"
   rmdir "$test_dir"
}
trap cleanup 0
trap 'exit 130' INT
trap 'exit 143' TERM
# Section garbage collection excludes unrelated mining/registration functions.
# pkg-config output is intentionally split into compiler arguments.
"$CC" -std=gnu11 -O2 -I. -ffunction-sections -fdata-sections \
   tests/test-deqode-gbt.c algo/deqode/deqode.c \
   -Wl,--gc-sections $(pkg-config --cflags --libs jansson) \
   -o "$test_dir/test-deqode-gbt"
"$test_dir/test-deqode-gbt"
