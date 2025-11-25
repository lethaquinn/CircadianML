"""
Circadian Transcriptomics Analysis Pipeline
Advanced computational framework for discovering circadian rhythms and biomarkers in transcriptomic data
"""

__version__ = "0.1.0"
__author__ = "Lethaquinn"
__email__ = "2679066373@qq.com"

from .rhythm_detection import jtk_cycle, lomb_scargle, cosinor
from .phase_prediction import ml_models, ensemble
from .visualization import phase_plots, dashboard
from .utils import helpers

# Demo function for quick start
def demo_analysis():
    """
    Run a demo analysis with synthetic data.
    Generates circadian expression data, runs rhythm detection and phase prediction,
    and saves figures to results/demo/.
    """
    import numpy as np
    import os
    from .utils.helpers import generate_synthetic_data
    from .rhythm_detection.jtk_cycle import detect_rhythms
    from .visualization.phase_plots import plot_rhythmic_genes

    print("🧬 Circadian Transcriptomics Analysis Demo")
    print("=" * 50)

    # Generate synthetic data
    print("\n1. Generating synthetic circadian expression data...")
    data = generate_synthetic_data(n_genes=100, n_timepoints=24, n_rhythmic=30)
    print(f"   ✓ Created dataset: {data.shape[0]} genes × {data.shape[1]} timepoints")

    # Detect rhythms
    print("\n2. Running JTK_CYCLE rhythm detection...")
    results = detect_rhythms(data)
    rhythmic_genes = results[results['BH.Q'] < 0.05]
    print(f"   ✓ Detected {len(rhythmic_genes)} rhythmic genes (q < 0.05)")

    # Create output directory
    os.makedirs("results/demo", exist_ok=True)

    # Plot results
    print("\n3. Generating visualizations...")
    plot_rhythmic_genes(data, rhythmic_genes, save_path="results/demo/rhythmic_genes.png")
    print(f"   ✓ Saved figure to results/demo/rhythmic_genes.png")

    print("\n✅ Demo completed successfully!")
    print("   Check results/demo/ for output figures.")

    return data, results

__all__ = [
    "demo_analysis",
    "jtk_cycle",
    "lomb_scargle",
    "cosinor",
    "ml_models",
    "ensemble",
    "phase_plots",
    "dashboard",
    "helpers",
]
