import os
import random
import networkx as nx

OUTPUT_DIR = "data/input"

GRAPH_SIZES = [100, 200, 400, 600, 800, 1000]
GRAPHS_PER_SIZE = 5

EDGE_PROBABILITY = 0.01
SEED = 7498


def save_graph(graph, filepath):
    """
    Save graph in format required by C++ loader:

    n m
    u v
    u v
    ...
    """

    n = graph.number_of_nodes()
    m = graph.number_of_edges()

    with open(filepath, "w") as f:
        f.write(f"{n} {m}\n")

        for u, v in graph.edges():
            f.write(f"{u} {v}\n")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    random.seed(SEED)

    for n in GRAPH_SIZES:
        for i in range(1, GRAPHS_PER_SIZE + 1):
            graph_seed = SEED + n + i

            G = nx.erdos_renyi_graph(
                n=n,
                p=EDGE_PROBABILITY,
                seed=graph_seed
            )

            filename = f"graph_n{n}_{i}.txt"
            filepath = os.path.join(
                OUTPUT_DIR,
                filename
            )

            save_graph(G, filepath)

            print(
                f"Generated: {filename} "
                f"(nodes={G.number_of_nodes()}, "
                f"edges={G.number_of_edges()})"
            )


if __name__ == "__main__":
    main()
