#include "graph.hpp"
#include "eigenvector_centrality.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "usage: "
                  << argv[0]
                  << " <input_file> <output_file> <method>\n";

        std::cerr << "method options: power | direct\n";
        return 1;
    }

    const std::string input_path = argv[1];
    const std::string output_path = argv[2];
    const std::string method = argv[3];

    try {
        Graph g = load_graph_from_file(input_path);

        std::vector<double> centrality;

        auto start =
            std::chrono::high_resolution_clock::now();

        if (method == "power") {
            centrality =
                eigenvector_centrality_power_iteration(g);
        }
        else if (method == "direct") {
            centrality =
                eigenvector_centrality_direct(g);
        }
        else {
            throw std::runtime_error(
                "invalid method. Use: power or direct"
            );
        }

        auto end =
            std::chrono::high_resolution_clock::now();

        double runtime_ms =
            std::chrono::duration<double, std::milli>(
                end - start
            ).count();

        std::ofstream out(output_path);

        if (!out) {
            throw std::runtime_error(
                "failed to open output file"
            );
        }

        out << "# Method: " << method << "\n";
        out << "# Vertices: " << g.num_vertices() << "\n";
        out << "# Edges: " << g.num_edges() << "\n";
        out << "# Runtime(ms): " << runtime_ms << "\n";

        for (std::size_t i = 0; i < centrality.size(); ++i) {
            out << i << " "
                << centrality[i] << "\n";
        }

        std::cout << "completed: "
                  << method
                  << " | runtime = "
                  << runtime_ms
                  << " ms\n";
    }
    catch (const std::exception& e) {
        std::cerr << "error: "
                  << e.what()
                  << "\n";
        return 1;
    }

    return 0;
}
