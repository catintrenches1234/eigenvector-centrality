#ifndef EIGENVECTOR_CENTRALITY_HPP_
#define EIGENVECTOR_CENTRALITY_HPP_

#include <cstddef>
#include <vector>

class Graph;

/*
 * Power Iteration Method
 * Computes eigenvector centrality iteratively.
 */
std::vector<double> eigenvector_centrality_power_iteration(
    const Graph& g,
    std::size_t max_iterations = 100,
    double tolerance = 1e-6
);

/*
 * Direct Eigenvalue Decomposition
 * Computes eigenvector centrality using full matrix eigendecomposition.
 */
std::vector<double> eigenvector_centrality_direct(
    const Graph& g
);

#endif // EIGENVECTOR_CENTRALITY_HPP_
