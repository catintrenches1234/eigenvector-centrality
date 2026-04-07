#include "graph.hpp"
#include "eigenvector_centrality.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <input_file> <output_file>\n";
        return 1;
    }

    const std::string input_path  = argv[1];
    const std::string output_path = argv[2];

    try {
        Graph g = load_graph_from_file(input_path);

        auto centrality = eigenvector_centrality(g);

        std::ofstream out(output_path);
        if (!out) {
            throw std::runtime_error("failed to open output file");
        }

        for (std::size_t i = 0; i < centrality.size(); ++i) {
            out << i << " " << centrality[i] << "\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
