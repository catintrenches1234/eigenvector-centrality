#ifndef EIGENVECTOR_CENTRALITY_HPP_
#define EIGENVECTOR_CENTRALITY_HPP_

#include <vector>

class Graph;

std::vector<double> eigenvector_centrality(
    const Graph& g,
    std::size_t max_iterations = 100,
    double tolerance = 1e-6
);

#endif // EIGENVECTOR_CENTRALITY_HPP_
