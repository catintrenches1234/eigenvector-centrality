import os
import math
import csv

INPUT_DIR = "data/output"
OUTPUT_DIR = "report/tables"

os.makedirs(OUTPUT_DIR, exist_ok=True)


def extract_metadata(filepath):
    method = None
    vertices = None
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

                elif "Runtime(ms):" in line:
                    runtime = float(line.split(":")[1].strip())

            else:
                parts = line.split()
                if len(parts) == 2:
                    values.append(float(parts[1]))

    return {
        "method": method,
        "vertices": vertices,
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


def export_runtime_table(results):
    power_grouped = {}
    direct_grouped = {}

    for item in results:
        n = item["vertices"]

        if item["method"] == "power":
            power_grouped.setdefault(n, []).append(item["runtime"])

        elif item["method"] == "direct":
            direct_grouped.setdefault(n, []).append(item["runtime"])

    x_values = sorted(
        set(power_grouped.keys()) |
        set(direct_grouped.keys())
    )

    output_path = os.path.join(OUTPUT_DIR, "runtime_table.csv")

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Vertices", "Power(ms)", "Direct(ms)"])

        for n in x_values:
            power_avg = average(power_grouped.get(n, []))
            direct_avg = average(direct_grouped.get(n, []))
            writer.writerow([n, power_avg, direct_avg])

    print(f"Saved: {output_path}")


def export_accuracy_table(results):
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

        error = compute_percentage_relative_error(
            methods["power"]["values"],
            methods["direct"]["values"]
        )

        if error is None:
            continue

        n = methods["power"]["vertices"]
        error_by_size.setdefault(n, []).append(error)

    output_path = os.path.join(OUTPUT_DIR, "accuracy_table.csv")

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Vertices", "RelativeError(%)"])

        for n in sorted(error_by_size.keys()):
            avg_error = average(error_by_size[n])
            writer.writerow([n, avg_error])

    print(f"Saved: {output_path}")


def main():
    results = collect_results()

    if not results:
        print("No output files found.")
        return

    export_runtime_table(results)
    export_accuracy_table(results)

    print("CSV tables generated successfully.")


if __name__ == "__main__":
    main()
