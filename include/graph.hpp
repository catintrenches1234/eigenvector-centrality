#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <cstddef>
#include <string>
#include <vector>

class Graph {
public:
    using vertex_id = std::size_t;

    explicit Graph(std::size_t n);

    std::size_t num_vertices() const noexcept;
    std::size_t num_edges() const noexcept;

    const std::vector<vertex_id>& neighbors(vertex_id u) const;

private:
    void add_edge(vertex_id u, vertex_id v);

    std::size_t n_ = 0;
    std::size_t m_ = 0;
    std::vector<std::vector<vertex_id>> adj_;

    friend Graph load_graph_from_file(const std::string& path);
};

Graph load_graph_from_file(const std::string& path);

#endif // GRAPH_HPP_
