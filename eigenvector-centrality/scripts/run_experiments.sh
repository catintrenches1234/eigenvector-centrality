#!/usr/bin/env bash

set -e

BUILD_DIR=build
INPUT_DIR=data/input
OUTPUT_DIR=data/output

mkdir -p "$OUTPUT_DIR"

EXEC="$BUILD_DIR/main"

echo "Running Power Iteration and Direct Eigenvalue Decomposition..."

for input_file in "$INPUT_DIR"/*.txt; do
    filename=$(basename "$input_file")

    power_output="$OUTPUT_DIR/${filename%.txt}_power.out"
    direct_output="$OUTPUT_DIR/${filename%.txt}_direct.out"

    echo "Processing: $filename"

    echo "  -> Power Iteration"
    "$EXEC" "$input_file" "$power_output" power

    echo "  -> Direct Eigenvalue Decomposition"
    "$EXEC" "$input_file" "$direct_output" direct
done

echo "All experiments completed."
