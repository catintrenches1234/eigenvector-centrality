# Eigenvector Centrality

Implementation of eigenvector centrality for **undirected, unweighted graphs** using a C++17 codebase and a simple experiment pipeline.

---

## Project Structure

```
.
├── CMakeLists.txt
├── data
│   ├── input        # generated graphs (ignored)
│   └── output       # centrality results (ignored)
├── include
│   ├── graph.hpp
│   └── eigenvector_centrality.hpp
├── src
│   ├── graph.cpp
│   ├── eigenvector_centrality.cpp
│   └── main.cpp
├── scripts
│   ├── generate_graphs.py
│   ├── run_experiments.sh
│   ├── plot.py
│   └── pipeline.sh
├── report
│   ├── figures
│   └── report.pdf
├── requirements.txt
└── README.md
```

---

## Build

```bash
cmake -S . -B build
cmake --build build
```

Binary:

```
build/main
```

---

## Usage

### Run on a single graph

```bash
./build/main <input_file> <output_file>
```

Example:

```bash
./build/main data/input/graph_n100_p0.05.txt data/output/out.txt
```

---

## Input Format

```
n m
u1 v1
u2 v2
...
```

* `n`: number of vertices (0-indexed)
* `m`: number of edges

---

## Full Pipeline

### One-time setup

```bash
pip install -r requirements.txt
chmod +x scripts/*.sh
```

### Run everything

```bash
./scripts/pipeline.sh
```

---

## Manual Steps (optional)

```bash
python3 scripts/generate_graphs.py
cmake -S . -B build
cmake --build build
bash scripts/run_experiments.sh
python3 scripts/plot.py
```

---

## Output

Each output file contains:

```
vertex_id centrality_score
```

---

## Visualization

Plots are saved to:

```
report/figures/
```

---

## Dependencies

* C++17 compiler
* CMake ≥ 3.10
* Python 3
* Python packages:

  ```
  matplotlib
  ```

---

## Reproducibility

* Graphs are generated deterministically
* Pipeline is reproducible via `pipeline.sh`
* Generated data is not stored in the repository
