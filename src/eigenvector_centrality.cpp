#include "eigenvector_centrality.hpp"
#include "graph.hpp"

#include <eigen3/Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

std::vector<double> eigenvector_centrality_power_iteration(
    const Graph& g,
    std::size_t max_iterations,
    double tolerance
) {
    const std::size_t n = g.num_vertices();

    if (n == 0) {
        return {};
    }

    std::vector<double> x(
        n,
        1.0 / std::sqrt(static_cast<double>(n))
    );

    std::vector<double> x_new(n, 0.0);

    for (std::size_t iter = 0; iter < max_iterations; ++iter) {
        std::fill(x_new.begin(), x_new.end(), 0.0);

        for (std::size_t u = 0; u < n; ++u) {
            for (auto v : g.neighbours(u)) {
                x_new[u] += x[v];
            }
        }

        double l2_norm = 0.0;

        for (double val : x_new) {
            l2_norm += val * val;
        }

        l2_norm = std::sqrt(l2_norm);

        if (l2_norm == 0.0) {
            throw std::runtime_error(
                "zero norm encountered in power iteration"
            );
        }

        for (double& val : x_new) {
            val /= l2_norm;
        }

        double diff = 0.0;

        for (std::size_t i = 0; i < n; ++i) {
            diff += std::abs(x_new[i] - x[i]);
        }

        if (diff < tolerance) {
            double l1_norm = 0.0;

            for (double val : x_new) {
                l1_norm += std::abs(val);
            }

            if (l1_norm == 0.0) {
                throw std::runtime_error(
                    "zero L1 norm in final normalization"
                );
            }

            for (double& val : x_new) {
                val /= l1_norm;
            }

            return x_new;
        }

        x = x_new;
    }

    double l1_norm = 0.0;

    for (double val : x) {
        l1_norm += std::abs(val);
    }

    if (l1_norm == 0.0) {
        throw std::runtime_error(
            "zero L1 norm after max iterations"
        );
    }

    for (double& val : x) {
        val /= l1_norm;
    }

    return x;
}


std::vector<double> eigenvector_centrality_direct(
    const Graph& g
) {
    const std::size_t n = g.num_vertices();

    if (n == 0) {
        return {};
    }

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);

    for (std::size_t u = 0; u < n; ++u) {
        for (auto v : g.neighbours(u)) {
            A(u, v) = 1.0;
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(A);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error(
            "eigendecomposition failed"
        );
    }

    Eigen::VectorXd x =
        solver.eigenvectors().col(n - 1);

    std::vector<double> result(n);

    double l1_norm = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        result[i] = std::abs(x(i));
        l1_norm += result[i];
    }

    if (l1_norm == 0.0) {
        throw std::runtime_error(
            "zero norm encountered in direct method"
        );
    }

    for (double& val : result) {
        val /= l1_norm;
    }

    return result;
}
