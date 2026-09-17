#!/bin/bash
# Quick integration test for libsais-packed.
#
# Runs the program against example_data/*.txt with a set of flag/sparseness
# combinations and checks that the output hash matches the value recorded in
# expected_outputs.sha256. If the hash still matches, the program still
# behaves the same as when the test was last recorded.
#
# Usage:
#   tests/run_tests.sh            check output against expected_outputs.sha256
#   tests/run_tests.sh --update   rebuild and overwrite expected_outputs.sha256
#                                 with the current output. Only do this after
#                                 confirming a change in output is intentional.

set -u

# Find the project root relative to this script, so it works no matter which
# directory you run it from.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

BIN="$ROOT_DIR/build/libsais-packed"
OUT_DIR="$SCRIPT_DIR/output"
EXPECTED_FILE="$SCRIPT_DIR/expected_outputs.sha256"
HG="$ROOT_DIR/example_data/human_genome.1000.txt"
UP="$ROOT_DIR/example_data/uniprot_entries.1000.txt"

UPDATE=false
[ "${1:-}" = "--update" ] && UPDATE=true

echo "Building..."
mkdir -p "$ROOT_DIR/build"
( cd "$ROOT_DIR/build" && cmake .. && make ) > "$SCRIPT_DIR/build.log" 2>&1
if [ ! -x "$BIN" ]; then
    echo "Build failed, see tests/build.log"
    exit 1
fi

mkdir -p "$OUT_DIR"
[ "$UPDATE" = true ] && > "$EXPECTED_FILE"   # start the file fresh when recording

PASS=0
FAIL=0

# Runs the program once, then either checks its output hash against the
# recorded one, or (in --update mode) records it.
#   $1 = test label (used as the output filename and the key in the hash file)
#   $2 = flags to pass to the program, e.g. "-s 4" or "-u -s 3"
#   $3 = input file
check() {
    label=$1
    flags=$2
    input=$3
    output="$OUT_DIR/$label.ssa"

    if ! $BIN $flags "$input" "$output" > "$OUT_DIR/$label.log" 2>&1; then
        echo "FAIL $label (program exited with an error, see tests/output/$label.log)"
        FAIL=$((FAIL + 1))
        return
    fi

    hash=$(shasum -a 256 "$output" | cut -d ' ' -f 1)

    if [ "$UPDATE" = true ]; then
        echo "$label $hash" >> "$EXPECTED_FILE"
        echo "recorded $label"
        return
    fi

    expected=$(grep "^$label " "$EXPECTED_FILE" | cut -d ' ' -f 2)
    if [ "$hash" = "$expected" ]; then
        echo "PASS $label"
        PASS=$((PASS + 1))
    else
        echo "FAIL $label (output changed)"
        FAIL=$((FAIL + 1))
    fi
}

# One line per test case. Together they cover every code path in the
# program: sparseness 1 (raw, no packing), the 8-bit, 16-bit and 32-bit
# packed paths, -u (unoptimized), -d (dynamic sparseness) and -c (compressed
# output).
check hg_raw_s1         "-s 1"    "$HG"
check hg_8bit_s2        "-s 2"    "$HG"
check hg_16bit_s4       "-s 4"    "$HG"
check hg_32bit_s8       "-s 8"    "$HG"
check hg_unoptimized_s3 "-u -s 3" "$HG"
check hg_dynamic_s8     "-d -s 8" "$HG"
check hg_compressed_s2  "-c -s 2" "$HG"
check up_raw_s1         "-s 1"    "$UP"
check up_16bit_s2       "-s 2"    "$UP"
check up_32bit_s5       "-s 5"    "$UP"
check up_unoptimized_s2 "-u -s 2" "$UP"
check up_dynamic_s5     "-d -s 5" "$UP"
check up_compressed_s3  "-c -s 3" "$UP"

if [ "$UPDATE" = true ]; then
    echo ""
    echo "Updated $EXPECTED_FILE"
    exit 0
fi

echo ""
echo "$PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ]
