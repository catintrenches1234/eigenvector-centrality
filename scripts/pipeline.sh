#!/usr/bin/env bash

set -e

echo "[1/4] Generating graphs"
python3 scripts/generate_graphs.py

echo "[2/4] Building project"
cmake -S . -B build
cmake --build build

echo "[3/4] Running experiments"
bash scripts/run_experiments.sh

echo "[4/4] Plotting results"
python3 scripts/plot.py

echo "Pipeline complete"
