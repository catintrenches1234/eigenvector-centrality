#!/usr/bin/env python3
"""
Main Eigenvector Centrality Project
===================================
Comprehensive analysis with Direct and Power Iteration methods,
benchmarking, and topological indices for molecular graphs.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import json
import os
from pathlib import Path
import time

from molecular_graphs import MOLECULES, BOILING_POINTS
from algorithms import EigenvectorCentralityDirect, EigenvectorCentralityPowerIteration, compare_methods
from analysis import MoleculeAnalysis, format_result, build_adjacency_matrix


def ensure_output_dirs():
    """Create output directories if they don't exist."""
    Path("../plots").mkdir(exist_ok=True)
    Path("../data").mkdir(exist_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def run_molecular_analysis():
    """Analyze all molecules."""
    print("\n" + "="*80)
    print("EIGENVECTOR CENTRALITY — MOLECULAR GRAPH ANALYSIS")
    print("="*80)
    print("Based on Ruhnau (2000) and Rani & Balamurugan (2025)\n")
    
    results = []
    for mol_name, mol_data in MOLECULES.items():
        mol = MoleculeAnalysis(
            name=mol_name,
            n=mol_data["n"],
            edges=mol_data["edges"],
            formula=mol_data.get("formula", ""),
            note=mol_data.get("note", ""),
            boiling_point=BOILING_POINTS.get(mol_name)
        )
        mol.analyze()
        result = mol.get_result_dict()
        results.append(result)
        print(format_result(result))
    
    print("="*80)
    return results


def benchmark_algorithms(num_nodes_range=None):
    """
    Benchmark Direct vs Power Iteration on graphs of varying sizes.
    
    NOTE: For small graphs (n < 100), direct eigendecomposition is typically faster
    due to lower overhead. Power Iteration excels on large sparse graphs (n > 1000).
    
    Creates linear chain + random edges (sparse, but not pathological).
    """
    if num_nodes_range is None:
        num_nodes_range = [10, 20, 50, 100]
    
    print("\n" + "="*80)
    print("ALGORITHM BENCHMARKING")
    print("="*80)
    print("Note: For small graphs, Direct method is typically faster.")
    print("Power Iteration advantage appears at n > 1000 with sparse matrices.\n")
    
    results_bench = {
        "num_nodes": [],
        "num_edges": [],
        "direct_time": [],
        "power_time": [],
        "speedup": [],
        "l2_difference": [],
    }
    
    for n in num_nodes_range:
        # Create well-connected sparse graph (good spectral properties for Power Iteration)
        # Linear chain + additional random edges for better conditioning
        np.random.seed(42)
        edges = [(i, i+1) for i in range(n-1)]  # Linear chain base
        
        # Add edges to improve spectral conditioning
        for _ in range(n//3):
            i, j = np.random.choice(n, 2, replace=False)
            if i > j:
                i, j = j, i
            if (i, j) not in edges and i != j:
                edges.append((i, j))
        
        A = build_adjacency_matrix(n, edges)
        
        # Run comparison with reasonable tolerance
        try:
            comp = compare_methods(A, tolerance=1e-6, max_iterations=500)
            
            direct_time = comp["direct"]["metadata"]["computation_time"]
            power_time = comp["power_iteration"]["metadata"]["computation_time"]
            speedup = power_time / direct_time if direct_time > 0 else 0  # Inverted for clarity
            l2_diff = comp["comparison"]["L2_difference"]
            
            results_bench["num_nodes"].append(n)
            results_bench["num_edges"].append(len(edges))
            results_bench["direct_time"].append(direct_time)
            results_bench["power_time"].append(power_time)
            results_bench["speedup"].append(speedup)
            results_bench["l2_difference"].append(l2_diff)
            
            # Show which method is faster
            faster = "Direct" if direct_time < power_time else "Power Iter"
            print(f"n={n:4d}, m={len(edges):3d} | "
                  f"Direct: {direct_time*1000:8.4f}ms | "
                  f"Power: {power_time*1000:8.4f}ms | "
                  f"Faster: {faster:10s} | "
                  f"L2 err: {l2_diff:.2e}")
        except Exception as e:
            print(f"n={n:4d}: ERROR - {e}")
    
    print("="*80)
    return results_bench


def convergence_analysis(mol_data):
    """Analyze Power Iteration convergence for selected molecules."""
    print("\n" + "="*80)
    print("CONVERGENCE ANALYSIS — Power Iteration Method")
    print("="*80 + "\n")
    
    convergence_data = {}
    
    for mol_name in ["n-butane", "n-pentane", "n-hexane"]:
        mol_info = mol_data[mol_name]
        n = mol_info["n"]
        edges = mol_info["edges"]
        A = build_adjacency_matrix(n, edges)
        
        pi = EigenvectorCentralityPowerIteration(
            A, max_iterations=200, tolerance=1e-10
        )
        pi.compute()
        
        convergence_data[mol_name] = {
            "history": pi.convergence_history,
            "iterations": pi.iterations,
            "final_error": pi.convergence_error,
        }
        
        print(f"{mol_name}: converged in {pi.iterations} iterations")
        print(f"  Final residual: {pi.convergence_error:.2e}")
        print()
    
    print("="*80)
    return convergence_data


def comparison_table(results):
    """Print comparison table of all molecules."""
    print("\n" + "="*90)
    print("COMPARISON TABLE — Topological Indices")
    print("="*90)
    header = f"{'Molecule':<35} {'λ₁':>6} {'M1x':>8} {'M2x':>8} {'Fx':>8} {'RRx':>8} {'BP(°C)':>8}"
    print(header)
    print("-"*90)
    for r in results:
        ti = r['topological_indices']
        bp = f"{r['boiling_point_C']:.1f}" if r['boiling_point_C'] is not None else "—"
        print(f"{r['name']:<35} {r['principal_eigenvalue']:>6.4f} "
              f"{ti['M1x']:>8.4f} {ti['M2x']:>8.4f} {ti['Fx']:>8.4f} "
              f"{ti['RRx']:>8.4f} {bp:>8}")
    print("="*90)


def verify_octane(results):
    """Verify n-octane against published values."""
    published = {
        "M1x": 5.0240, "M2x": 0.9395, "Fx": 1.9476,
        "Hx": 20.9806, "ISx": 1.2247, "RRx": 2.4801, "ABCx": 24.7841
    }
    
    octane = next(r for r in results if r['name'] == 'n-octane')
    computed = octane['topological_indices']
    
    print("\n" + "="*60)
    print("VERIFICATION: n-Octane vs. Rani & Balamurugan (2025)")
    print("="*60)
    print(f"{'Index':<10} {'Published':>12} {'Computed':>12} {'Match':>8}")
    print("-"*60)
    
    all_match = True
    for idx in published:
        pub = published[idx]
        comp = computed[idx]
        match = abs(pub - comp) < 0.001
        all_match = all_match and match
        marker = "✓" if match else "✗"
        print(f"{idx:<10} {pub:>12.4f} {comp:>12.4f} {marker:>8}")
    print("="*60)
    print(f"Verification: {'PASSED' if all_match else 'FAILED'}")
    print("="*60)


# ─────────────────────────────────────────────────────────────────────────────
# PLOTTING FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def plot_centrality_profiles(results):
    """Plot eigenvector centrality values for butane and pentane isomers."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle('Eigenvector Centrality Profiles for Alkane Isomers', fontsize=14, fontweight='bold')
    
    molecules_to_plot = [
        "n-butane", "isobutane (2-methylpropane)",
        "n-pentane", "isopentane (2-methylbutane)", "neopentane (2,2-dimethylpropane)",
        "n-hexane"
    ]
    
    axes = axes.flatten()
    for idx, mol_name in enumerate(molecules_to_plot):
        result = next(r for r in results if r['name'] == mol_name)
        x = result['eigenvector_centrality_euclidean']
        atoms = [f"C{i+1}" for i in range(len(x))]
        
        ax = axes[idx]
        colors = ['red' if i == result['most_central_atom'] else 'steelblue' for i in range(len(x))]
        ax.bar(atoms, x, color=colors, alpha=0.7, edgecolor='black')
        ax.set_ylabel('Centrality (L2)', fontsize=10)
        ax.set_title(mol_name, fontsize=11, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        ax.set_ylim(0, max(x) * 1.1)
    
    # Remove extra subplot
    fig.delaxes(axes[5])
    
    plt.tight_layout()
    plt.savefig('../plots/01_centrality_profiles.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: plots/01_centrality_profiles.png")


def plot_topological_indices(results):
    """Plot topological indices for all molecules."""
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    fig.suptitle('Topological Indices Based on Eigenvector Centrality', fontsize=14, fontweight='bold')
    
    indices = ["M1x", "M2x", "Fx", "Hx", "ISx", "RRx", "ABCx"]
    axes = axes.flatten()
    
    mol_names = [r['name'] for r in results]
    
    for idx_num, idx_name in enumerate(indices):
        values = [r['topological_indices'][idx_name] for r in results]
        ax = axes[idx_num]
        colors_list = ['#ff9999' if 'isobutane' in n or 'neopentane' in n or '2,2-dimethylbutane' in n 
                       else '#99ccff' for n in mol_names]
        ax.bar(range(len(values)), values, color=colors_list, alpha=0.7, edgecolor='black')
        ax.set_ylabel(idx_name, fontsize=11, fontweight='bold')
        ax.set_xticks(range(len(mol_names)))
        ax.set_xticklabels([n.split('(')[0].strip()[:12] for n in mol_names], rotation=45, ha='right', fontsize=8)
        ax.grid(axis='y', alpha=0.3)
        ax.set_title(f'{idx_name}', fontsize=11, fontweight='bold')
    
    # Remove extra subplot
    fig.delaxes(axes[7])
    
    plt.tight_layout()
    plt.savefig('../plots/02_topological_indices.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: plots/02_topological_indices.png")


def plot_benchmark_results(bench_data):
    """Plot algorithm benchmarking results."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Algorithm Performance Comparison', fontsize=14, fontweight='bold')
    
    # Timing comparison
    ax = axes[0]
    ax.loglog(bench_data['num_nodes'], bench_data['direct_time'], 'o-', label='Direct', linewidth=2, markersize=8)
    ax.loglog(bench_data['num_nodes'], bench_data['power_time'], 's-', label='Power Iteration', linewidth=2, markersize=8)
    ax.set_xlabel('Number of Nodes (n)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Computation Time (seconds)', fontsize=11, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('Execution Time vs. Graph Size', fontsize=11, fontweight='bold')
    
    # Speedup
    ax = axes[1]
    ax.semilogx(bench_data['num_nodes'], bench_data['speedup'], 'D-', color='green', linewidth=2, markersize=8)
    ax.set_xlabel('Number of Nodes (n)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Speedup Factor (Direct / Power Iter)', fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_title('Speedup of Power Iteration over Direct Method', fontsize=11, fontweight='bold')
    ax.axhline(y=1, color='r', linestyle='--', alpha=0.5, label='Equal performance')
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    plt.savefig('../plots/03_algorithm_benchmark.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: plots/03_algorithm_benchmark.png")


def plot_convergence(convergence_data):
    """Plot Power Iteration convergence curves."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for mol_name, data in convergence_data.items():
        history = data['history']
        ax.semilogy(history, marker='o', linewidth=2, markersize=5, label=mol_name, alpha=0.8)
    
    ax.set_xlabel('Iteration', fontsize=11, fontweight='bold')
    ax.set_ylabel('L2 Residual ||v_{k+1} - v_k||', fontsize=11, fontweight='bold')
    ax.set_title('Power Iteration Convergence: Residual vs. Iteration', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.axhline(y=1e-6, color='r', linestyle='--', alpha=0.5, label='Tolerance (1e-6)')
    
    plt.tight_layout()
    plt.savefig('../plots/04_convergence_curves.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: plots/04_convergence_curves.png")


def plot_boiling_point_correlation(results):
    """Correlate topological indices with boiling point."""
    # Filter molecules with known boiling points
    valid_results = [r for r in results if r['boiling_point_C'] is not None]
    
    if len(valid_results) < 2:
        print("⚠ Not enough molecules with boiling points for correlation plot")
        return
    
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle('Correlation: Topological Indices vs. Boiling Point', fontsize=14, fontweight='bold')
    
    indices_to_plot = ['M1x', 'M2x', 'RRx', 'Fx', 'ISx', 'Hx']
    axes = axes.flatten()
    
    for idx_num, idx_name in enumerate(indices_to_plot):
        x_vals = [r['topological_indices'][idx_name] for r in valid_results]
        y_vals = [r['boiling_point_C'] for r in valid_results]
        names = [r['name'].split('(')[0].strip() for r in valid_results]
        
        ax = axes[idx_num]
        ax.scatter(x_vals, y_vals, s=100, alpha=0.6, edgecolors='black', linewidth=1.5)
        
        # Add labels
        for i, name in enumerate(names):
            ax.annotate(name[:10], (x_vals[i], y_vals[i]), fontsize=8, alpha=0.7)
        
        # Fit line
        if len(x_vals) > 1:
            z = np.polyfit(x_vals, y_vals, 1)
            p = np.poly1d(z)
            x_line = np.linspace(min(x_vals), max(x_vals), 100)
            ax.plot(x_line, p(x_line), "r--", alpha=0.5, linewidth=2)
            
            # Correlation coefficient
            corr = np.corrcoef(x_vals, y_vals)[0, 1]
            ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes,
                   fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_xlabel(idx_name, fontsize=11, fontweight='bold')
        ax.set_ylabel('Boiling Point (°C)', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('../plots/05_boiling_point_correlation.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: plots/05_boiling_point_correlation.png")


def generate_all_plots(results, bench_data, convergence_data):
    """Generate all plots."""
    print("\n" + "="*80)
    print("GENERATING PLOTS")
    print("="*80 + "\n")
    
    ensure_output_dirs()
    plot_centrality_profiles(results)
    plot_topological_indices(results)
    plot_benchmark_results(bench_data)
    plot_convergence(convergence_data)
    plot_boiling_point_correlation(results)
    
    print("\n" + "="*80)


# ─────────────────────────────────────────────────────────────────────────────
# DATA EXPORT
# ─────────────────────────────────────────────────────────────────────────────

def save_results_json(results):
    """Save analysis results to JSON."""
    ensure_output_dirs()
    
    # Convert numpy arrays to lists
    export_data = []
    for r in results:
        r_copy = r.copy()
        r_copy['eigenvector_centrality_euclidean'] = [float(v) for v in r['eigenvector_centrality_euclidean']]
        r_copy['node_centrality_scaled'] = [float(v) for v in r['node_centrality_scaled']]
        r_copy['topological_indices'] = {k: float(v) for k, v in r['topological_indices'].items()}
        export_data.append(r_copy)
    
    with open('../data/analysis_results.json', 'w') as f:
        json.dump(export_data, f, indent=2)
    print("✓ Saved: data/analysis_results.json")


def save_benchmark_csv(bench_data):
    """Save benchmark data to CSV."""
    ensure_output_dirs()
    
    with open('../data/benchmark_results.csv', 'w') as f:
        f.write("NumNodes,NumEdges,DirectTime_s,PowerIterTime_s,Speedup,L2Difference\n")
        for i in range(len(bench_data['num_nodes'])):
            f.write(f"{bench_data['num_nodes'][i]},")
            f.write(f"{bench_data['num_edges'][i]},")
            f.write(f"{bench_data['direct_time'][i]:.6e},")
            f.write(f"{bench_data['power_time'][i]:.6e},")
            f.write(f"{bench_data['speedup'][i]:.4f},")
            f.write(f"{bench_data['l2_difference'][i]:.6e}\n")
    
    print("✓ Saved: data/benchmark_results.csv")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN EXECUTION
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Run analyses
    results = run_molecular_analysis()
    bench_data = benchmark_algorithms([10, 20, 50, 100, 200, 500])
    convergence_data = convergence_analysis(MOLECULES)
    
    # Print summaries
    comparison_table(results)
    verify_octane(results)
    
    # Generate plots
    generate_all_plots(results, bench_data, convergence_data)
    
    # Export data
    print("\n" + "="*80)
    print("EXPORTING DATA")
    print("="*80 + "\n")
    save_results_json(results)
    save_benchmark_csv(bench_data)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print("\nGenerated files:")
    print("  Plots:")
    print("    • plots/01_centrality_profiles.png")
    print("    • plots/02_topological_indices.png")
    print("    • plots/03_algorithm_benchmark.png")
    print("    • plots/04_convergence_curves.png")
    print("    • plots/05_boiling_point_correlation.png")
    print("  Data:")
    print("    • data/analysis_results.json")
    print("    • data/benchmark_results.csv")
    print("="*80 + "\n")
