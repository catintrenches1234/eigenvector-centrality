#include "eigenvector_centrality.hpp"
#include "graph.hpp"

#include <cmath>
#include <stdexcept>

std::vector<double> eigenvector_centrality(
    const Graph& g,
    std::size_t max_iterations,
    double tolerance
) {
    const std::size_t n = g.num_vertices();
    if (n == 0) {
        return {};
    }

    std::vector<double> x(n, 1.0 / std::sqrt(static_cast<double>(n)));
    std::vector<double> x_new(n, 0.0);

    for (std::size_t iter = 0; iter < max_iterations; ++iter) {
        std::fill(x_new.begin(), x_new.end(), 0.0);

        for (std::size_t u = 0; u < n; ++u) {
            for (auto v : g.neighbours(u)) {
                x_new[u] += x[v];
            }
        }

        double norm = 0.0;
        for (double val : x_new) {
            norm += val * val;
        }
        norm = std::sqrt(norm);

        if (norm == 0.0) {
            throw std::runtime_error("zero norm encountered");
        }

        for (double& val : x_new) {
            val /= norm;
        }

        double diff = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            diff += std::abs(x_new[i] - x[i]);
        }

        if (diff < tolerance) {
            return x_new;
        }

        x = x_new;
    }

    return x;
}
