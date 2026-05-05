# Eigenvector Centrality for Molecular Graphs

**Project 20** — Intro to Algorithms Engineering Course Project

**Authors:** Ayush Subham Moharana (2024101074) | Ishaan Shiv Kumar (2024101098)

## Overview

This project implements and analyzes **eigenvector centrality** for carbon skeleton graphs of alkane molecules, comparing two computational methods and deriving topological indices based on spectral graph theory.

### Key Features

- **Two Algorithms**: Direct eigenvalue decomposition (exact) vs. Power Iteration method (scalable)
- **Molecular Analysis**: 11 alkane isomers (Butane, pentane, hexane, octane) with diverse structures
- **Topological Indices**: Seven eigenvector-based indices from Rani & Balamurugan (2025)
- **Benchmarking**: Comprehensive performance comparison across graphs (n = 10 to 500)
- **Validation**: Results verified against NetworkX reference implementation and published values
- **Visualization**: Five high-quality publication-ready plots
- **Theory**: Complete mathematical foundation with spectral graph theory and Perron–Frobenius theorem

## Project Structure

```
eigenvector-centrality/
├── README.md                          # This file
├── src/
│   ├── main.py                        # Main execution script
│   ├── molecular_graphs.py            # Molecule definitions (C4-C8)
│   ├── algorithms.py                  # Direct & Power Iteration methods
│   └── analysis.py                    # Topological indices & utilities
├── plots/
│   ├── 01_centrality_profiles.png     # Eigenvector centrality for isomers
│   ├── 02_topological_indices.png     # Seven indices comparison
│   ├── 03_algorithm_benchmark.png     # Performance scaling
│   ├── 04_convergence_curves.png      # Power Iteration convergence
│   └── 05_boiling_point_correlation.png # Structure-property relationships
├── data/
│   ├── analysis_results.json          # Complete analysis output (JSON)
│   └── benchmark_results.csv          # Algorithm benchmarking (CSV)
├── report/
│   ├── project_report.tex             # Full LaTeX document
│   └── project_report.pdf             # Compiled PDF (if available)
└── requirements.txt                   # Python dependencies

```


## Quick Start

### Prerequisites

- Python 3.8+
- Virtual environment (recommended)

### Installation & Execution

```bash
# Navigate to project
cd eigenvector-centrality

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run analysis
cd src
python main.py
```

**Expected Runtime:** ~30–60 seconds | **Expected Output:** 5 plots + 2 data files + console report

---

## Project Description

### Problem Statement

Given an undirected, connected graph $G = (V, E)$ with adjacency matrix $A \in \mathbb{R}^{n \times n}$, compute the centrality vector $\mathbf{x} \in \mathbb{R}^n$ satisfying:

$$A\mathbf{x} = \lambda_1 \mathbf{x}$$

where $\lambda_1$ is the dominant eigenvalue. Each component $x_i$ represents the eigenvector centrality score of node $i$.

### Algorithms

#### 1. Direct Eigenvalue Decomposition

**Method**: Full spectral decomposition using `numpy.linalg.eigh()`

- Computes all eigenvalues and eigenvectors
- Extracts dominant eigenvector
- Exact solution (numerical precision only)

**Advantages**:
- ✓ Exact baseline
- ✓ Guaranteed convergence

**Disadvantages**:
- ✗ $O(n^3)$ time complexity
- ✗ $O(n^2)$ space complexity
- ✗ Prohibitive for large graphs

**Use cases**: Small graphs ($n < 5000$), ground truth validation

#### 2. Power Iteration Method

**Algorithm**:
```
Input: Adjacency matrix A, tolerance ε, max iterations K
Initialize: v₀ ← random unit vector
for k = 1, 2, ..., K:
    v ← A·vₖ₋₁
    λ ← v·vₖ₋₁  (Rayleigh quotient)
    vₖ ← v / ||v||
    if ||vₖ - vₖ₋₁|| < ε: return (vₖ, λ)
return (vₖ, λ)
```

**Properties**:
- Iterative convergence: $||x_k - x^*|| \propto \lambda_2^k / \lambda_1^k$
- Cost per iteration: $O(m)$ where $m = $ number of edges
- Memory: $O(n)$ for sparse graphs

**Advantages**:
- ✓ $O(m \cdot k)$ total cost (sparse graphs)
- ✓ $O(n)$ memory footprint
- ✓ Scalable to large networks

**Disadvantages**:
- ✗ Iterative (not exact)
- ✗ Convergence dependent on $\lambda_1 / \lambda_2$ ratio

**Use cases**: Large, sparse graphs; networks with thousands to millions of nodes

### Molecular Graph Representation

Molecules are represented as undirected graphs where:
- **Nodes**: Carbon atoms (indices 0 to n-1)
- **Edges**: C-C bonds
- **Hydrogen atoms**: Omitted (skeletal graph)

**Examples**:

| Molecule | Structure | λ₁ | Most Central | Notes |
|----------|-----------|----|-----------|----|
| **n-butane** | Linear C-C-C-C | 1.618 | Central atoms | Even distribution |
| **isobutane** | Star C(C)(C)C | 1.732 | Central C | High centralization |
| **n-pentane** | Linear chain (5C) | 1.809 | Middle atoms | Gradual decrease |
| **neopentane** | Star with 4 branches | 2.0 | Center atom | Maximum star structure |

### Normalization

Eigenvectors are normalized to ensure comparability:

1. **L₂ (Euclidean) Normalization**: $||x||_2 = 1$
   - Enables node-centrality status (Ruhnau 2000)
   - Scaled: $nc_i = \sqrt{2} \cdot x_i$ ensures star center = 1

2. **L∞ (Maximum) Normalization**: $||x||_∞ = 1$ (relative centrality)

3. **L₁ (Sum) Normalization**: $\sum x_i = 1$ (proportional centrality)

### Topological Indices

Seven eigenvector-based indices (Rani & Balamurugan, 2025):

| Index | Formula | Meaning |
|-------|---------|---------|
| **M1x** | $\sum_{(i,j) \in E} (x_i + x_j)$ | First Zagreb-like index |
| **M2x** | $\sum_{(i,j) \in E} x_i x_j$ | Second Zagreb-like index |
| **Fx** | $\sum_{(i,j) \in E} (x_i^2 + x_j^2)$ | Forgotten index |
| **Hx** | $\sum_{(i,j) \in E} \frac{2}{x_i + x_j}$ | Harmonic index |
| **ISx** | $\sum_{(i,j) \in E} \frac{x_i x_j}{x_i + x_j}$ | Inverse sum index |
| **RRx** | $\sum_{(i,j) \in E} \sqrt{x_i x_j}$ | Reciprocal Randić index |
| **ABCx** | $\sum_{(i,j) \in E} \sqrt{\left\|\frac{x_i + x_j - 2}{x_i x_j}\right\|}$ | Atom-bond connectivity |

## Key Results

### Molecular Analyses

**Example: Butane Isomers**

| Property | n-butane | isobutane |
|----------|----------|-----------|
| λ₁ | 1.618 | 1.732 |
| Most Central | C2, C3 | C2 |
| M1x | 2.0 | 2.185 |
| Boiling Point | -0.5°C | -11.7°C |
| Structure Effect | Linear → even centrality | Branched → concentrated |

**Key Observation**: Higher centralization in branched isomers correlates with lower boiling points (fewer degrees of freedom).

### Algorithm Comparison

**Benchmark on sparse graphs (chain + random edges)**:

| n | m | Direct (ms) | Power Iter (ms) | Winner | Speedup |
|---|---|-------------|-----------------|--------|---------|
| 10 | 11 | 0.099 | 2.400 | Direct | 24.4× |
| 20 | 25 | 0.093 | 0.577 | Direct | 6.2× |
| 50 | 64 | 0.558 | 1.097 | Direct | 1.96× |
| 100 | 131 | 0.944 | 2.173 | Direct | 2.30× |
| 200 | 264 | 3.894 | 2.554 | **Power Iteration** | 1.5× faster |
| 500 | 665 | 21.946 | 29.815 | Direct | 1.36× |

**Key Findings**:

1. **Why speedup is non-monotonic**: Power Iteration reaches competitive speed at n=200 (graph large enough to offset overhead), but gets slower at n=500 due to more edges increasing convergence iterations.

2. **Error metric**: All solutions (n≥20) converge correctly. The previous L2 differences of 2.0 were due to sign ambiguity in eigenvectors (Power Iteration converged to `-v` instead of `v`). Fixed by using sign-aware metric.

3. **Practical recommendation**: Use Direct for most applications (n < 1000). Power Iteration advantage emerges at very large scales (n > 5000) with sparse matrices where dense eigendecomposition becomes impractical.

### Verification

All results verified against NetworkX 2.6 reference implementation:

```
n-octane Topological Indices (Rani & Balamurugan 2025):
Index    Published    Computed    Match
M1x      5.0240       5.0240      ✓
M2x      0.9395       0.9395      ✓
Fx       1.9476       1.9476      ✓
...
```

## Code Modules

### `molecular_graphs.py`

Defines molecule structures:
- **MOLECULES**: Dictionary of 11 alkane isomers (C4-C8)
- **BOILING_POINTS**: Known experimental values

### `algorithms.py`

Core computational methods:
- **EigenvectorCentralityDirect**: Direct eigendecomposition
- **EigenvectorCentralityPowerIteration**: Power Iteration algorithm
- **compare_methods()**: Head-to-head comparison

### `analysis.py`

Analysis utilities:
- **build_adjacency_matrix()**: Create sparse adjacency matrix
- **topological_indices()**: Compute seven indices
- **MoleculeAnalysis**: Unified analysis wrapper
- **format_result()**: Pretty-print results

### `main.py`

Orchestration and plotting:
- **run_molecular_analysis()**: Analyze all molecules
- **benchmark_algorithms()**: Compare methods on synthetic graphs
- **convergence_analysis()**: Study Power Iteration convergence
- **generate_all_plots()**: Create 5 visualization plots
- **save_results_json()** / **save_benchmark_csv()**: Export data

## Theory & References

### Perron-Frobenius Theorem

For connected graphs with non-negative adjacency matrices:
- $\lambda_1$ is real, positive, and strictly dominant
- Corresponding eigenvector has all positive entries
- Suitable as centrality measure

### Node-Centrality Requirements (Ruhnau 2000)

Valid centrality measure must satisfy:
1. Values in $[0, 1]$
2. Maximum value 1 achieved only by star graph center

Only **Euclidean normalization** with scaling $nc_i = \sqrt{2} \cdot x_i$ meets both criteria.

### Power Iteration Convergence

Residual norm decreases as:
$$||v_k - v^*||_2 \leq C \left|\frac{\lambda_2}{\lambda_1}\right|^k$$

Convergence rate depends on spectral gap $(\lambda_1 - \lambda_2) / \lambda_1$.

## References

- **Ruhnau, B.** (2000). Eigenvector-centrality — a node-centrality?. *Social Networks*, 22(4), 357–365.
- **Rani, A.S.J. & Balamurugan, B.J.** (2025). Novel eigenvector centrality indices for octane isomers to explore their physicochemical properties. *Scientific Reports*, 15, 34730.
- **NetworkX Documentation**: https://networkx.org/documentation/
- **Cvetković, D., Doob, M., & Sachs, H.** (1995). *Spectra of Graphs: Theory and Application*. Academic Press.

## Troubleshooting

### Issue: ImportError for numpy or matplotlib

**Solution**:
```bash
pip install --upgrade numpy matplotlib
```

### Issue: Plots not appearing

**Solution**:
- Ensure `plots/` directory exists
- Check write permissions: `touch plots/test.txt`
- Verify matplotlib backend

### Issue: "ConnectionError" when benchmarking large graphs

**Solution**:
- Reduce `num_nodes_range` in benchmark call
- Increase `tolerance` in Power Iteration (trades accuracy for speed)

## Performance Tips

### For Large Graphs

1. Use **Power Iteration** (not Direct method)
2. Increase `tolerance` parameter (default 1e-6)
3. Limit `max_iterations` based on required accuracy
4. Use sparse adjacency matrix (scipy.sparse) for modifications

### For Production Use

```python
from algorithms import EigenvectorCentralityPowerIteration
import scipy.sparse as sp

# Load/build sparse matrix
A_sparse = sp.csr_matrix(...)

# Create solver
solver = EigenvectorCentralityPowerIteration(
    A_sparse, 
    max_iterations=500,
    tolerance=1e-8
)

# Compute
x = solver.compute()
print(f"Converged in {solver.iterations} iterations")
```

## License & Authors

**Course**: Intro to Algorithms Engineering  
**Authors**: 
- Ayush Subham Moharana (2024101074)
- Ishaan Shiv Kumar (2024101098)

**Project Year**: 2024-2025

## Contact & Support

For issues, questions, or suggestions:
- Review the project report in `report/project_report.pdf`
- Check inline code documentation in `src/`
- Run individual functions for debugging:
  ```bash
  python3 -c "from molecular_graphs import MOLECULES; print(MOLECULES.keys())"
  ```

---

**Last Updated**: 2025-05-05
