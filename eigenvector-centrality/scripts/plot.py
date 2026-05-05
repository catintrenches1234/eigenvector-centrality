import os
import math
import matplotlib.pyplot as plt

INPUT_DIR = "data/output"
FIG_DIR = "results/figures"

os.makedirs(FIG_DIR, exist_ok=True)


def extract_metadata(filepath):
    method = None
    vertices = None
    edges = None
    runtime = None
    values = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith("#"):
                if "Method:" in line:
                    method = line.split(":")[1].strip()

                elif "Vertices:" in line:
                    vertices = int(line.split(":")[1].strip())

                elif "Edges:" in line:
                    edges = int(line.split(":")[1].strip())

                elif "Runtime(ms):" in line:
                    runtime = float(line.split(":")[1].strip())

            else:
                parts = line.split()
                if len(parts) == 2:
                    values.append(float(parts[1]))

    return {
        "method": method,
        "vertices": vertices,
        "edges": edges,
        "runtime": runtime,
        "values": values,
        "filename": os.path.basename(filepath)
    }


def collect_results():
    results = []

    for file in os.listdir(INPUT_DIR):
        if not file.endswith(".out"):
            continue

        path = os.path.join(INPUT_DIR, file)
        results.append(extract_metadata(path))

    return results


def average(values):
    if not values:
        return None
    return sum(values) / len(values)


def compute_percentage_relative_error(power, direct):
    if len(power) != len(direct):
        return None

    numerator = 0.0
    denominator = 0.0

    for x, y in zip(power, direct):
        numerator += (x - y) ** 2
        denominator += y ** 2

    numerator = math.sqrt(numerator)
    denominator = math.sqrt(denominator)

    if denominator == 0.0:
        return None

    return (numerator / denominator) * 100.0


def plot_runtime_comparison(results):
    power_grouped = {}
    direct_grouped = {}

    for item in results:
        n = item["vertices"]

        if item["method"] == "power":
            power_grouped.setdefault(n, []).append(
                item["runtime"]
            )

        elif item["method"] == "direct":
            direct_grouped.setdefault(n, []).append(
                item["runtime"]
            )

    x_values = sorted(
        set(power_grouped.keys()) |
        set(direct_grouped.keys())
    )

    power_avg = []
    direct_avg = []

    for n in x_values:
        power_avg.append(
            average(power_grouped.get(n, []))
        )

        direct_avg.append(
            average(direct_grouped.get(n, []))
        )

    plt.figure()

    plt.plot(
        x_values,
        power_avg,
        marker="o",
        label="Power Iteration"
    )

    plt.plot(
        x_values,
        direct_avg,
        marker="o",
        label="Direct Eigendecomposition"
    )

    plt.ylim(bottom=0)

    plt.xlabel("Number of Vertices")
    plt.ylabel("Average Runtime (ms)")
    plt.title("Runtime Comparison")
    plt.legend()
    plt.tight_layout()

    output_path = os.path.join(
        FIG_DIR,
        "runtime_comparison.png"
    )

    plt.savefig(output_path)
    plt.close()

    print(f"Saved: {output_path}")


def plot_accuracy_comparison(results):
    grouped = {}

    for item in results:
        name = item["filename"]

        if "_power.out" in name:
            key = name.replace("_power.out", "")
            grouped.setdefault(key, {})["power"] = item

        elif "_direct.out" in name:
            key = name.replace("_direct.out", "")
            grouped.setdefault(key, {})["direct"] = item

    error_by_size = {}

    for key, methods in grouped.items():
        if "power" not in methods or "direct" not in methods:
            continue

        power_vals = methods["power"]["values"]
        direct_vals = methods["direct"]["values"]

        error = compute_percentage_relative_error(
            power_vals,
            direct_vals
        )

        if error is None:
            continue

        n = methods["power"]["vertices"]

        error_by_size.setdefault(n, []).append(error)

    x_values = sorted(error_by_size.keys())
    avg_errors = []

    for n in x_values:
        avg_errors.append(
            average(error_by_size[n])
        )

    plt.figure()

    plt.plot(
        x_values,
        avg_errors,
        marker="o"
    )

    plt.ylim(bottom=0)

    plt.xlabel("Number of Vertices")
    plt.ylabel("Average Percentage Error (%)")
    plt.title("Accuracy Comparison (Power vs Direct)")
    plt.tight_layout()

    output_path = os.path.join(
        FIG_DIR,
        "accuracy_comparison.png"
    )

    plt.savefig(output_path)
    plt.close()

    print(f"Saved: {output_path}")


def main():
    results = collect_results()

    if not results:
        print("No output files found.")
        return

    plot_runtime_comparison(results)
    plot_accuracy_comparison(results)

    print("All plots generated successfully.")


if __name__ == "__main__":
    main()
