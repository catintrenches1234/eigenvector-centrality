#include "graph.hpp"

#include <fstream>
#include <stdexcept>

Graph::Graph(std::size_t n) : n_(n), m_(0), adj_(n) {}

std::size_t Graph::num_vertices() const noexcept {
    return n_;
}

std::size_t Graph::num_edges() const noexcept {
    return m_;
}

const std::vector<Graph::vertex_id>& Graph::neighbours(vertex_id u) const noexcept {
    return adj_[u];
}

void Graph::add_edge(vertex_id u, vertex_id v) {
    if (u >= n_ || v >= n_) {
        throw std::out_of_range("vertex_id out of range");
    }
    adj_[u].push_back(v);
    adj_[v].push_back(u);
    ++m_;
}

Graph load_graph_from_file(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("failed to open file");
    }

    std::size_t n, m;
    file >> n >> m;

    Graph g(n);

    std::size_t u, v;
    while (file >> u >> v) {
        g.add_edge(u, v);
    }

    return g;
}
