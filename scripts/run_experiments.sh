#!/usr/bin/env bash

set -e

BUILD_DIR=build
INPUT_DIR=data/input
OUTPUT_DIR=data/output

mkdir -p $OUTPUT_DIR

cmake -S . -B $BUILD_DIR
cmake --build $BUILD_DIR

EXEC="$BUILD_DIR/main"

for input_file in $INPUT_DIR/*.txt; do
    filename=$(basename "$input_file")
    output_file="$OUTPUT_DIR/$filename.out"

    echo "running $filename"
    $EXEC "$input_file" "$output_file"
done

echo "done"
