#ifndef EIGENVECTOR_CENTRALITY_HPP_
#define EIGENVECTOR_CENTRALITY_HPP_

#include <cstddef>
#include <vector>

class Graph;

std::vector<double> eigenvector_centrality_power_iteration(
    const Graph& g,
    std::size_t max_iterations = 100,
    double tolerance = 1e-6
);

std::vector<double> eigenvector_centrality_direct(
    const Graph& g
);

#endif // EIGENVECTOR_CENTRALITY_HPP_
