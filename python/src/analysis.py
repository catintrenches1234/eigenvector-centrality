"""
Topological Indices and Analysis
=================================
Based on eigenvector centrality values (Rani & Balamurugan, 2025).
"""

import numpy as np


def build_adjacency_matrix(n, edges):
    """Build the n×n adjacency matrix for a carbon skeleton graph."""
    A = np.zeros((n, n), dtype=float)
    for i, j in edges:
        A[i][j] = 1.0
        A[j][i] = 1.0
    return A


def eigenvector_centrality_normalized(A, normalization="euclidean"):
    """
    Compute eigenvector centrality with various normalizations.
    
    Parameters:
        A              : adjacency matrix
        normalization  : 'euclidean', 'max', or 'sum'
    
    Returns:
        x_euclidean : raw eigenvector (L2 norm = 1)
        lambda1     : principal eigenvalue
        nc          : node-centrality (normalized by specified method)
    """
    eigenvalues, eigenvectors = np.linalg.eigh(A)
    idx = np.argmax(eigenvalues)
    lambda1 = eigenvalues[idx]
    x_raw = eigenvectors[:, idx]
    
    # Ensure positive
    if x_raw[0] < 0:
        x_raw = -x_raw
    
    # Euclidean normalization (L2)
    x_euclidean = x_raw / np.linalg.norm(x_raw)
    
    if normalization == "euclidean":
        # Scale so star center = 1 (Ruhnau 2000)
        nc = np.sqrt(2) * x_euclidean
    elif normalization == "max":
        nc = x_euclidean / np.max(x_euclidean)
    elif normalization == "sum":
        nc = x_euclidean / np.sum(x_euclidean)
    else:
        nc = x_euclidean
    
    return x_euclidean, lambda1, nc


def topological_indices(edges, x):
    """
    Compute seven eigenvector centrality-based topological indices
    (Rani & Balamurugan, 2025).
    
    Parameters:
        edges : list of (i,j) tuples (C-C bonds)
        x     : eigenvector centrality values (0-indexed)
    
    Returns dict of index values.
    """
    M1x = sum(x[i] + x[j] for i, j in edges)
    M2x = sum(x[i] * x[j] for i, j in edges)
    Fx  = sum(x[i]**2 + x[j]**2 for i, j in edges)
    Hx  = sum(2 / (x[i] + x[j]) for i, j in edges)
    ISx = sum((x[i] * x[j]) / (x[i] + x[j]) for i, j in edges)
    RRx = sum(np.sqrt(x[i] * x[j]) for i, j in edges)
    ABCx = sum(np.sqrt(abs((x[i] + x[j] - 2) / (x[i] * x[j]))) for i, j in edges)
    
    return {
        "M1x":  M1x,
        "M2x":  M2x,
        "Fx":   Fx,
        "Hx":   Hx,
        "ISx":  ISx,
        "RRx":  RRx,
        "ABCx": ABCx,
    }


def node_degrees(n, edges):
    """Return degree sequence for the carbon skeleton."""
    deg = [0] * n
    for i, j in edges:
        deg[i] += 1
        deg[j] += 1
    return deg


class MoleculeAnalysis:
    """Complete analysis of a single molecule."""
    
    def __init__(self, name, n, edges, formula=None, note=None, boiling_point=None):
        """
        Parameters:
            name           : molecule name
            n              : number of atoms
            edges          : list of (i,j) bonds
            formula        : molecular formula (optional)
            note           : structural note (optional)
            boiling_point  : known boiling point in °C (optional)
        """
        self.name = name
        self.n = n
        self.edges = edges
        self.formula = formula or ""
        self.note = note or ""
        self.boiling_point = boiling_point
        
        self.A = None
        self.lambda1 = None
        self.x_euclidean = None
        self.nc_euclidean = None
        self.degrees = None
        self.indices = None
        
    def analyze(self):
        """Run full analysis."""
        self.A = build_adjacency_matrix(self.n, self.edges)
        self.x_euclidean, self.lambda1, self.nc_euclidean = \
            eigenvector_centrality_normalized(self.A, normalization="euclidean")
        self.indices = topological_indices(self.edges, self.x_euclidean)
        self.degrees = node_degrees(self.n, self.edges)
        
        return self
    
    def get_result_dict(self):
        """Export as dictionary."""
        return {
            "name": self.name,
            "formula": self.formula,
            "note": self.note,
            "n_atoms": self.n,
            "n_bonds": len(self.edges),
            "principal_eigenvalue": self.lambda1,
            "degrees": self.degrees,
            "eigenvector_centrality_euclidean": self.x_euclidean,
            "node_centrality_scaled": self.nc_euclidean,
            "most_central_atom": int(np.argmax(self.x_euclidean)),
            "topological_indices": self.indices,
            "boiling_point_C": self.boiling_point,
        }


def format_result(result):
    """Format a result dict for pretty printing."""
    sep = "=" * 72
    lines = []
    lines.append(sep)
    lines.append(f"  {result['name']}  ({result['formula']})")
    lines.append(f"  {result['note']}")
    lines.append(f"  Atoms: {result['n_atoms']}  |  Bonds: {result['n_bonds']}  |  "
                 f"λ₁ = {result['principal_eigenvalue']:.4f}")
    if result['boiling_point_C'] is not None:
        lines.append(f"  Boiling point: {result['boiling_point_C']} °C")
    lines.append("")
    
    lines.append("  Atom   Degree   Eigenvector Centrality (L2)   Node-Centrality (×√2)")
    for i, (deg, xv, nc) in enumerate(zip(
        result['degrees'],
        result['eigenvector_centrality_euclidean'],
        result['node_centrality_scaled']
    )):
        marker = " ← most central" if i == result['most_central_atom'] else ""
        lines.append(f"   C{i+1}     {deg}        {xv:.6f}                     {nc:.6f}{marker}")
    
    lines.append("")
    lines.append("  Topological Indices:")
    ti = result['topological_indices']
    lines.append(f"    M1x  = {ti['M1x']:.6f}   (first: Σ(xi+xj) over bonds)")
    lines.append(f"    M2x  = {ti['M2x']:.6f}   (second: Σ(xi·xj))")
    lines.append(f"    Fx   = {ti['Fx']:.6f}   (forgotten: Σ(xi²+xj²))")
    lines.append(f"    Hx   = {ti['Hx']:.6f}   (harmonic: Σ 2/(xi+xj))")
    lines.append(f"    ISx  = {ti['ISx']:.6f}   (inverse sum: Σ xi·xj/(xi+xj))")
    lines.append(f"    RRx  = {ti['RRx']:.6f}   (reciprocal Randic: Σ √(xi·xj))")
    lines.append(f"    ABCx = {ti['ABCx']:.6f}   (atom-bond connectivity)")
    lines.append("")
    
    return "\n".join(lines)
