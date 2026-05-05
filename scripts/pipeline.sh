#!/usr/bin/env bash

set -e

echo "[1/5] Generating graphs"
python3 scripts/generate_graphs.py

echo "[2/5] Building project"
cmake -S . -B build
cmake --build build

echo "[3/5] Running experiments"
bash scripts/run_experiments.sh

echo "[4/5] Generating tables"
python3 scripts/generate_tables.py

echo "[5/5] Generating plots"
python3 scripts/plot.py

echo "Pipeline complete."
