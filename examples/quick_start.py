#!/usr/bin/env python3
"""
Quick start example for circadian transcriptomics analysis.
"""

import sys
import os

# Add src to path if running directly
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from circadian_analysis.utils.helpers import generate_synthetic_data
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.visualization.phase_plots import plot_rhythmic_genes, plot_phase_distribution


def main():
    print("=" * 60)
    print("Circadian Transcriptomics Analysis - Quick Start")
    print("=" * 60)

    # 1. Generate synthetic data
    print("\n1. Generating synthetic circadian expression data...")
    data = generate_synthetic_data(
        n_genes=200,
        n_timepoints=24,
        n_rhythmic=50,
        period=24.0,
        random_state=42
    )
    print(f"   Generated: {data.shape[0]} genes × {data.shape[1]} timepoints")

    # 2. Detect rhythms
    print("\n2. Running JTK_CYCLE rhythm detection...")
    results = detect_rhythms(data, period=24.0, alpha=0.05)

    # Filter significant genes
    rhythmic_genes = results[results['BH.Q'] < 0.05]
    print(f"   Detected {len(rhythmic_genes)} rhythmic genes (q < 0.05)")

    # 3. Visualize results
    print("\n3. Creating visualizations...")
    os.makedirs("results/quick_start", exist_ok=True)

    # Plot top rhythmic genes
    plot_rhythmic_genes(
        data,
        rhythmic_genes,
        top_n=6,
        save_path="results/quick_start/top_rhythmic_genes.png"
    )

    # Plot phase distribution
    plot_phase_distribution(
        rhythmic_genes,
        phase_col='LAG',
        save_path="results/quick_start/phase_distribution.png"
    )

    # 4. Save results
    print("\n4. Saving results...")
    results.to_csv("results/quick_start/jtk_results.csv", index=False)
    print(f"   Results saved to results/quick_start/jtk_results.csv")

    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("Check results/quick_start/ for output files")
    print("=" * 60)


if __name__ == "__main__":
    main()
