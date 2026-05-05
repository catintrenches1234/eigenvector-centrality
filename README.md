# Eigenvector Centrality

Implementation of eigenvector centrality for **undirected, unweighted graphs** using both **C++17** (general-purpose) and **Python 3** (molecular graph analysis) implementations.

## Introduction to Algorithm Engineering Project

- Ayush Subham Moharana   2024101074
- Ishaan Shiv Kumar       2024101098

---

## Project Structure

```
.
├── eigenvector-centrality/            # C++ Subproject
│   ├── CMakeLists.txt                 # C++ build configuration
│   ├── include/                       # C++ headers
│   │   ├── graph.hpp
│   │   └── eigenvector_centrality.hpp
│   ├── src/                           # C++ source
│   │   ├── graph.cpp
│   │   ├── eigenvector_centrality.cpp
│   │   └── main.cpp
│   ├── scripts/                       # Utility scripts
│   │   ├── generate_graphs.py         # Generate random graphs
│   │   ├── generate_tables.py         # Generate result tables
│   │   ├── run_experiments.sh         # Run C++ experiments
│   │   ├── plot.py                    # Plot C++ results
│   │   └── pipeline.sh                # Full pipeline
│   ├── data/                          # C++ experiment data
│   │   ├── input/                     # Generated graphs
│   │   └── output/                    # Centrality results
│   └── results/                       # C++ result output
│       ├── figures/                   # Generated plots
│       └── tables/                    # Generated tables
├── molecular-graphs/                  # Python Subproject
│   ├── src/                           # Python source
│   │   ├── main.py                    # Molecular analysis pipeline
│   │   ├── algorithms.py              # Direct & Power Iteration methods
│   │   ├── analysis.py                # Topological indices
│   │   └── molecular_graphs.py        # Alkane molecules (C4-C8)
│   ├── plots/                         # Generated visualizations
│   └── data/                          # Analysis output
│       ├── analysis_results.json      # Molecular analysis
│       └── benchmark_results.csv      # Algorithm benchmarks
├── .gitignore
├── requirements.txt                   # Python dependencies
└── README.md                          # This file
```

---

## Two Implementations

### C++ Implementation (General Purpose)

**Location:** `eigenvector-centrality/` subdirectory

- High-performance eigenvector centrality for arbitrary graphs
- Direct eigenvalue decomposition using Eigen library
- Generates random graphs and computes centrality
- Benchmarks on graphs with 100 to 2000 vertices

**Build & Run:**

```bash
cd eigenvector-centrality
cmake -S . -B build
cmake --build build
./build/main <input_file> <output_file>
```

### Python Implementation (Molecular Analysis)

**Location:** `molecular-graphs/` subdirectory

- Analysis of **11 alkane isomers** (C4 to C8) using spectral graph theory
- Compares two algorithms: Direct eigendecomposition vs. Power Iteration
- Computes 7 topological indices from eigenvector centrality
- Validates against NetworkX and published chemical data

**Run:**

```bash
cd molecular-graphs
python src/main.py
```

---

## Build (C++ only)

```bash
cd eigenvector-centrality
cmake -S . -B build
cmake --build build
```

Binary: `eigenvector-centrality/build/main`

---

## Usage

### C++ Pipeline (Single Graph)

```bash
cd eigenvector-centrality
./build/main <input_file> <output_file>
```

**Input Format:**
```
n m
u1 v1
u2 v2
...
```

Where:
- `n`: number of vertices (0-indexed)
- `m`: number of edges

**Output Format:**
```
vertex_id centrality_score
```

### Full C++ Pipeline

```bash
# One-time setup
pip install -r requirements.txt
cd eigenvector-centrality
chmod +x scripts/*.sh

# Run everything
./scripts/pipeline.sh
```

Or manually:

```bash
cd eigenvector-centrality
python3 scripts/generate_graphs.py
cmake -S . -B build
cmake --build build
bash scripts/run_experiments.sh
python3 scripts/plot.py
python3 scripts/generate_tables.py
```

### Python Pipeline (Molecular Analysis)

```bash
cd molecular-graphs
python src/main.py
```

Produces:
- 5 publication-ready plots in `plots/`
- JSON analysis results in `data/analysis_results.json`
- CSV benchmark comparison in `data/benchmark_results.csv`

---

## Dependencies

### C++ Build
- C++17 compiler (gcc, clang, MSVC)
- CMake ≥ 3.10
- Eigen3 library

### Python Runtime
- Python 3.10+
- numpy ≥ 1.26.0
- scipy ≥ 1.11.0
- matplotlib ≥ 3.8.0
- networkx ≥ 3.2

**Install all Python dependencies:**

```bash
pip install -r requirements.txt
```

---

## Results & Findings

### C++ Results
- Tested on random sparse graphs (100-2000 vertices, 0.2%-0.5% edge density)
- Output saved to `eigenvector-centrality/report/figures/` and `eigenvector-centrality/report/tables/`

### Python Results
- Eigenvector centrality computed for 11 alkane structures
- Algorithm benchmarks: Direct vs. Power Iteration
  - Direct eigendecomposition faster for n < 500
  - Power Iteration competitive at n ≈ 200+
  - Both methods converge correctly (L₂ error < 10⁻⁴)
- 7 topological indices computed and validated
- Boiling point correlation analysis

---

## Reproducibility

- C++ graphs generated deterministically (seed = hardcoded)
- Python analyses use fixed seeds (42) for reproducibility
- All results logged to output files

---

## References

1. **Eigenvector Centrality**: Bonacich (1987), "Power and Centrality"
2. **Power Iteration**: Kaczynski et al. (2014)
3. **Molecular Analysis**: Rani & Balamurugan (2025), "Molecular descriptors from graph topology"
4. **Alkane Structures**: IUPAC nomenclature and experimental boiling points
* Pipeline is reproducible via `pipeline.sh`
* Generated data is not stored in the repository
