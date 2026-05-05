"""
Eigenvector Centrality Algorithms
==================================
Direct eigenvalue decomposition and Power Iteration method.

Based on:
  - Ruhnau (2000): Eigenvector-centrality as node-centrality
  - Rani & Balamurugan (2025): Novel topological indices
"""

import numpy as np
import time


class EigenvectorCentralityDirect:
    """
    Direct eigenvalue decomposition using numpy.linalg.eigh.
    
    For exact baseline computation of eigenvector centrality.
    Suitable for small to medium graphs (n < 5000).
    
    Time Complexity: O(n³) for dense adjacency matrix
    Space Complexity: O(n²) for storing the full matrix and eigenvectors
    """
    
    def __init__(self, A):
        """
        Parameters:
            A : adjacency matrix (n×n numpy array)
        """
        self.A = A
        self.n = A.shape[0]
        self.lambda1 = None
        self.x_raw = None
        self.x_normalized = None
        self.iterations = 1
        self.convergence_error = 0.0
        self.computation_time = 0.0
        
    def compute(self):
        """
        Compute eigenvector centrality via direct eigenvalue decomposition.
        
        Algorithm:
          1. Compute full eigendecomposition A = Q Λ Q^T
          2. Extract principal eigenvector (largest eigenvalue)
          3. Ensure positive entries
          4. Normalize to L2 norm = 1
        
        Returns:
            x_norm : normalized principal eigenvector (L2 norm = 1)
        """
        start_time = time.time()
        
        # Full eigendecomposition
        eigenvalues, eigenvectors = np.linalg.eigh(self.A)
        
        # Extract principal eigenvector (largest eigenvalue)
        idx = np.argmax(eigenvalues)
        self.lambda1 = eigenvalues[idx]
        self.x_raw = eigenvectors[:, idx]
        
        # Ensure all entries positive (principal eigenvector property)
        if self.x_raw[0] < 0:
            self.x_raw = -self.x_raw
        
        # L2 normalization: ||x|| = 1
        self.x_normalized = self.x_raw / np.linalg.norm(self.x_raw)
        
        self.convergence_error = 0.0
        self.computation_time = time.time() - start_time
        
        return self.x_normalized
    
    def get_metadata(self):
        """Return algorithm metadata for reporting."""
        return {
            "method": "Direct Eigenvalue Decomposition",
            "iterations": self.iterations,
            "convergence_error": self.convergence_error,
            "computation_time": self.computation_time,
            "principal_eigenvalue": self.lambda1,
        }



class EigenvectorCentralityPowerIteration:
    """
    Power Iteration method for computing the dominant eigenvector.
    
    Iteratively: v_{k+1} = A @ v_k / ||A @ v_k||
    
    Suitable for large, sparse graphs due to O(m) cost per iteration
    where m = number of edges. Memory efficient for sparse matrices.
    
    Time Complexity: O(k·m) where k = iterations, m = edges
    Space Complexity: O(n) for storing vectors (sparse-friendly)
    """
    
    def __init__(self, A, max_iterations=100, tolerance=1e-6, random_seed=42):
        """
        Parameters:
            A                : adjacency matrix (n×n numpy array, can be sparse)
            max_iterations   : maximum iterations (default 100)
            tolerance        : convergence threshold for L2 residual (default 1e-6)
            random_seed      : seed for initial vector reproducibility (default 42)
        """
        self.A = A
        self.n = A.shape[0]
        self.max_iterations = max_iterations
        self.tolerance = tolerance
        self.random_seed = random_seed
        
        self.lambda1 = None
        self.x_normalized = None
        self.iterations = 0
        self.convergence_error = None
        self.convergence_history = []
        self.computation_time = 0.0
        
    def compute(self):
        """
        Run Power Iteration algorithm.
        
        Algorithm:
          1. Initialize v with random unit vector
          2. Repeat until convergence:
               a. Compute v_new = A @ v
               b. Estimate λ₁ = v^T A v
               c. Normalize: v_new = v_new / ||v_new||
               d. Check convergence: ||v_new - v|| < ε
          3. Return normalized eigenvector
        
        Convergence is measured by L2 residual: ||v_{k+1} - v_k||
        
        Returns:
            x : normalized principal eigenvector (L2 norm = 1)
        """
        start_time = time.time()
        
        # Initialize with random vector
        np.random.seed(self.random_seed)
        v = np.random.randn(self.n)
        v = v / np.linalg.norm(v)
        
        # Power iteration loop
        for k in range(self.max_iterations):
            # Matrix-vector product: v_new = A @ v
            v_new = self.A @ v
            
            # Compute Rayleigh quotient (approximation of λ₁)
            # λ₁ ≈ v^T A v (for normalized v)
            lambda_est = np.dot(v, v_new)
            
            # Normalize: v_new = v_new / ||v_new||
            v_new_norm = np.linalg.norm(v_new)
            v_new = v_new / v_new_norm
            
            # Convergence check: L2 residual ||v_new - v||
            residual = np.linalg.norm(v_new - v)
            self.convergence_history.append(residual)
            
            if residual < self.tolerance:
                # Converged
                self.iterations = k + 1
                self.lambda1 = np.dot(v_new, self.A @ v_new)
                self.x_normalized = v_new
                self.convergence_error = residual
                self.computation_time = time.time() - start_time
                return self.x_normalized
            
            v = v_new
        
        # Max iterations reached (did not converge)
        self.iterations = self.max_iterations
        v_new = self.A @ v
        self.lambda1 = np.dot(v, v_new)
        self.x_normalized = v
        self.convergence_error = self.convergence_history[-1]
        self.computation_time = time.time() - start_time
        
        return self.x_normalized
    
    def get_metadata(self):
        """Return algorithm metadata for reporting."""
        return {
            "method": "Power Iteration",
            "iterations": self.iterations,
            "convergence_error": self.convergence_error,
            "computation_time": self.computation_time,
            "principal_eigenvalue": self.lambda1,
        }


def compare_methods(A, tolerance=1e-6, max_iterations=100):
    """
    Compare Direct and Power Iteration methods on the same adjacency matrix.
    
    Parameters:
        A                : adjacency matrix
        tolerance        : convergence tolerance for Power Iteration
        max_iterations   : max iterations for Power Iteration
    
    Returns dict with results and comparison metrics.
    """
    # Direct method
    direct = EigenvectorCentralityDirect(A)
    x_direct = direct.compute()
    
    # Power Iteration
    power = EigenvectorCentralityPowerIteration(
        A, max_iterations=max_iterations, tolerance=tolerance
    )
    x_power = power.compute()
    
    # Compare solutions
    # Note: eigenvectors are defined up to sign, so check both +v and -v
    l2_diff_positive = np.linalg.norm(x_direct - x_power)
    l2_diff_negative = np.linalg.norm(x_direct + x_power)
    l2_diff = min(l2_diff_positive, l2_diff_negative)  # Use minimum (correct sign)
    correlation = np.dot(x_direct, x_power)
    speedup = direct.get_metadata()["computation_time"] / max(1e-9, power.get_metadata()["computation_time"])
    
    return {
        "direct": {
            "solution": x_direct,
            "metadata": direct.get_metadata(),
        },
        "power_iteration": {
            "solution": x_power,
            "metadata": power.get_metadata(),
        },
        "comparison": {
            "L2_difference": l2_diff,
            "correlation": correlation,
            "speedup": speedup,
        },
    }
