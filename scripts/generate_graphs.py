import os
import random

OUTPUT_DIR = "data/input"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def generate_random_graph(n, p):
    edges = []
    for u in range(n):
        for v in range(u + 1, n):
            if random.random() < p:
                edges.append((u, v))
    return edges


def write_graph(path, n, edges):
    with open(path, "w") as f:
        f.write(f"{n} {len(edges)}\n")
        for u, v in edges:
            f.write(f"{u} {v}\n")


def main():
    random.seed(42)

    configs = [
        (100, 0.05),
        (500, 0.01),
        (1000, 0.005),
        (2000, 0.002),
    ]

    for n, p in configs:
        edges = generate_random_graph(n, p)
        filename = f"graph_n{n}_p{p}.txt"
        path = os.path.join(OUTPUT_DIR, filename)
        write_graph(path, n, edges)
        print(f"generated {path}")


if __name__ == "__main__":
    main()
