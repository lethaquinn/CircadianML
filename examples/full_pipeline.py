#!/usr/bin/env python3
"""
Full pipeline example demonstrating all analysis capabilities.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from circadian_analysis.utils.helpers import generate_synthetic_data
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.rhythm_detection.cosinor import detect_rhythms_cosinor
from circadian_analysis.network_analysis.coexpression import build_coexpression_network, compute_network_metrics
from circadian_analysis.visualization.dashboard import create_analysis_dashboard


def main():
    print("=" * 70)
    print("Circadian Transcriptomics Analysis - Full Pipeline")
    print("=" * 70)

    # 1. Generate data
    print("\n[Step 1] Generating synthetic dataset...")
    data = generate_synthetic_data(
        n_genes=500,
        n_timepoints=24,
        n_rhythmic=150,
        period=24.0,
        random_state=42
    )
    print(f"   Dataset: {data.shape[0]} genes × {data.shape[1]} timepoints")

    # 2. Rhythm detection with multiple methods
    print("\n[Step 2] Running rhythm detection...")

    print("   a) JTK_CYCLE...")
    jtk_results = detect_rhythms(data, period=24.0)
    jtk_rhythmic = jtk_results[jtk_results['BH.Q'] < 0.05]
    print(f"      → {len(jtk_rhythmic)} rhythmic genes detected")

    print("   b) Cosinor regression...")
    cosinor_results = detect_rhythms_cosinor(data, period=24.0)
    cosinor_rhythmic = cosinor_results[cosinor_results['Qvalue'] < 0.05]
    print(f"      → {len(cosinor_rhythmic)} rhythmic genes detected")

    # 3. Network analysis
    print("\n[Step 3] Building co-expression network...")
    rhythmic_data = data.loc[jtk_rhythmic['GeneID'].values]

    adj_matrix, edge_list = build_coexpression_network(
        rhythmic_data,
        threshold=0.8,
        method='pearson'
    )
    print(f"   Network: {len(adj_matrix)} nodes, {len(edge_list)} edges")

    if len(edge_list) > 0:
        network_metrics = compute_network_metrics(adj_matrix)
        print(f"   Mean degree: {network_metrics['Degree'].mean():.2f}")

    # 4. Create comprehensive dashboard
    print("\n[Step 4] Creating analysis dashboard...")
    os.makedirs("results/full_pipeline", exist_ok=True)

    create_analysis_dashboard(
        data,
        jtk_results,
        save_path="results/full_pipeline/analysis_dashboard.png"
    )

    # 5. Save all results
    print("\n[Step 5] Saving results...")
    jtk_results.to_csv("results/full_pipeline/jtk_results.csv", index=False)
    cosinor_results.to_csv("results/full_pipeline/cosinor_results.csv", index=False)

    if len(edge_list) > 0:
        edge_list.to_csv("results/full_pipeline/network_edges.csv", index=False)
        network_metrics.to_csv("results/full_pipeline/network_metrics.csv", index=False)

    # 6. Summary report
    print("\n" + "=" * 70)
    print("ANALYSIS SUMMARY")
    print("=" * 70)
    print(f"\nDataset:")
    print(f"  • Total genes: {len(data)}")
    print(f"  • Timepoints: {data.shape[1]}")

    print(f"\nRhythm Detection:")
    print(f"  • JTK_CYCLE: {len(jtk_rhythmic)} rhythmic genes ({100*len(jtk_rhythmic)/len(data):.1f}%)")
    print(f"  • Cosinor: {len(cosinor_rhythmic)} rhythmic genes ({100*len(cosinor_rhythmic)/len(data):.1f}%)")

    if len(jtk_rhythmic) > 0:
        mean_amp = jtk_rhythmic['AMP'].mean()
        mean_phase = jtk_rhythmic['LAG'].mean()
        print(f"\nRhythmic Gene Characteristics:")
        print(f"  • Mean amplitude: {mean_amp:.3f}")
        print(f"  • Mean phase: {mean_phase:.2f} hours")

    if len(edge_list) > 0:
        print(f"\nNetwork Analysis:")
        print(f"  • Nodes: {len(adj_matrix)}")
        print(f"  • Edges: {len(edge_list)}")
        print(f"  • Mean degree: {network_metrics['Degree'].mean():.2f}")

    print(f"\nOutput files saved to: results/full_pipeline/")
    print("=" * 70)


if __name__ == "__main__":
    main()
