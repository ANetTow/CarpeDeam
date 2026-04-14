#!/bin/bash
set -e

# Base directory
BASE_DIR=$(dirname "$0")/..
BUILD_DIR="$BASE_DIR/build"
EXE="$BUILD_DIR/bin/carpedeam"

if [ ! -f "$EXE" ]; then
    echo "Error: carpedeam binary not found at $EXE. Please build it first."
    exit 1
fi

echo "Running test with Pen1-A10_S10_L001_CarpeDeam profiles..."

rm -rf "$BASE_DIR/test_tmp" "$BASE_DIR/test_output.fasta"

"$EXE" ancient_assemble \
    "$BASE_DIR/example/test_data.fq.gz" \
    "$BASE_DIR/test_output.fasta" \
    "$BASE_DIR/test_tmp" \
    --ancient-damage "$BASE_DIR/example/Pen1-A10_S10_L001_CarpeDeam_"

if [ -f "$BASE_DIR/test_output.fasta" ] && [ -s "$BASE_DIR/test_output.fasta" ]; then
    echo "Success: Assembly completed and output generated."
    # Clean up
    rm -rf "$BASE_DIR/test_tmp" "$BASE_DIR/test_output.fasta"
    exit 0
else
    echo "Failure: Assembly did not produce output.fasta or file is empty."
    exit 1
fi
