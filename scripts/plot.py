import os
import matplotlib.pyplot as plt

INPUT_DIR = "data/output"
FIG_DIR = "report/figures"

os.makedirs(FIG_DIR, exist_ok=True)


def read_centrality(path):
    values = []
    with open(path, "r") as f:
        for line in f:
            _, val = line.strip().split()
            values.append(float(val))
    return values


def plot_distribution(values, title, output_path):
    plt.figure()
    plt.hist(values, bins=50)
    plt.title(title)
    plt.xlabel("Eigenvector Centrality")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def main():
    for file in os.listdir(INPUT_DIR):
        if not file.endswith(".out"):
            continue

        input_path = os.path.join(INPUT_DIR, file)
        values = read_centrality(input_path)

        title = f"Centrality Distribution ({file})"
        output_file = file.replace(".out", ".png")
        output_path = os.path.join(FIG_DIR, output_file)

        plot_distribution(values, title, output_path)
        print(f"saved {output_path}")


if __name__ == "__main__":
    main()
