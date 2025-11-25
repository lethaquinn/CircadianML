"""
Dashboard creation for comprehensive circadian analysis results.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional
from .phase_plots import plot_rhythmic_genes, plot_phase_distribution


def create_analysis_dashboard(
    data: pd.DataFrame,
    rhythmic_results: pd.DataFrame,
    save_path: Optional[str] = None
):
    """
    Create comprehensive analysis dashboard.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix
    rhythmic_results : pd.DataFrame
        Rhythm detection results
    save_path : str, optional
        Path to save figure
    """
    fig = plt.figure(figsize=(16, 12))

    # Summary statistics
    n_total = len(data)
    n_rhythmic = len(rhythmic_results[rhythmic_results['BH.Q'] < 0.05]) if 'BH.Q' in rhythmic_results.columns else 0

    # Title
    fig.suptitle(f'Circadian Analysis Dashboard\n'
                f'{n_rhythmic} / {n_total} rhythmic genes ({100*n_rhythmic/n_total:.1f}%)',
                fontsize=16, fontweight='bold')

    # 1. Top rhythmic genes (2×3 grid)
    for i in range(6):
        ax = plt.subplot(4, 3, i+1)

        if i < len(rhythmic_results):
            gene_info = rhythmic_results.iloc[i]
            gene_id = gene_info['GeneID']
            expression = data.loc[gene_id].values
            time = np.arange(len(expression))

            ax.plot(time, expression, 'o-', linewidth=2, markersize=5)
            ax.set_xlabel('Time (h)', fontsize=9)
            ax.set_ylabel('Expression', fontsize=9)

            if 'BH.Q' in gene_info:
                title = f"{gene_id}\nq={gene_info['BH.Q']:.2e}"
            else:
                title = gene_id

            ax.set_title(title, fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3)

    # 2. Phase distribution
    if 'LAG' in rhythmic_results.columns:
        phase_col = 'LAG'
    elif 'Acrophase' in rhythmic_results.columns:
        phase_col = 'Acrophase'
    else:
        phase_col = None

    if phase_col:
        phases = rhythmic_results[phase_col].values

        # Linear histogram
        ax = plt.subplot(4, 3, 7)
        ax.hist(phases, bins=24, color='steelblue', edgecolor='black', alpha=0.7)
        ax.set_xlabel('Phase (hours)', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title('Phase Distribution', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

        # Circular plot
        ax = plt.subplot(4, 3, 8, projection='polar')
        theta = 2 * np.pi * phases / 24
        ax.hist(theta, bins=24, color='steelblue', edgecolor='black', alpha=0.7)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_title('Circular Phase', fontsize=11, fontweight='bold', pad=20)

    # 3. Amplitude distribution
    if 'AMP' in rhythmic_results.columns:
        ax = plt.subplot(4, 3, 9)
        amplitudes = rhythmic_results['AMP'].values
        ax.hist(amplitudes, bins=30, color='coral', edgecolor='black', alpha=0.7)
        ax.set_xlabel('Amplitude', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title('Amplitude Distribution', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

    # 4. Q-value distribution
    if 'BH.Q' in rhythmic_results.columns:
        ax = plt.subplot(4, 3, 10)
        qvalues = rhythmic_results['BH.Q'].values
        ax.hist(np.log10(qvalues + 1e-300), bins=30, color='green', edgecolor='black', alpha=0.7)
        ax.set_xlabel('log10(Q-value)', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title('Q-value Distribution', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.axvline(np.log10(0.05), color='red', linestyle='--', label='q=0.05')
        ax.legend()

    # 5. Period distribution (if available)
    if 'PER' in rhythmic_results.columns:
        ax = plt.subplot(4, 3, 11)
        periods = rhythmic_results['PER'].values
        ax.hist(periods, bins=20, color='purple', edgecolor='black', alpha=0.7)
        ax.set_xlabel('Period (hours)', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title('Period Distribution', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

    # 6. Summary text
    ax = plt.subplot(4, 3, 12)
    ax.axis('off')

    summary_text = f"""
    === Summary Statistics ===

    Total genes: {n_total:,}
    Rhythmic genes: {n_rhythmic:,}
    Percentage: {100*n_rhythmic/n_total:.2f}%

    """

    if phase_col:
        summary_text += f"Mean phase: {np.mean(phases):.2f} h\n"
    if 'AMP' in rhythmic_results.columns:
        summary_text += f"Mean amplitude: {np.mean(amplitudes):.3f}\n"

    ax.text(0.1, 0.5, summary_text, fontsize=11, family='monospace',
           verticalalignment='center')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"   Dashboard saved to {save_path}")
    else:
        plt.show()

    plt.close()


def create_comparison_dashboard(
    results_dict: dict,
    save_path: Optional[str] = None
):
    """
    Create dashboard comparing multiple analysis methods.

    Parameters
    ----------
    results_dict : dict
        Dictionary of method_name -> results DataFrame
    save_path : str, optional
        Path to save figure
    """
    n_methods = len(results_dict)

    fig, axes = plt.subplots(2, n_methods, figsize=(5*n_methods, 10))

    if n_methods == 1:
        axes = axes.reshape(2, 1)

    for i, (method, results) in enumerate(results_dict.items()):
        # Number of rhythmic genes
        if 'BH.Q' in results.columns:
            n_rhythmic = len(results[results['BH.Q'] < 0.05])
        elif 'Qvalue' in results.columns:
            n_rhythmic = len(results[results['Qvalue'] < 0.05])
        else:
            n_rhythmic = len(results)

        # Phase distribution
        ax = axes[0, i]
        if 'LAG' in results.columns:
            phases = results['LAG'].values
        elif 'Acrophase' in results.columns:
            phases = results['Acrophase'].values
        else:
            phases = None

        if phases is not None:
            ax.hist(phases, bins=24, color='steelblue', edgecolor='black', alpha=0.7)
            ax.set_xlabel('Phase (hours)', fontsize=10)
            ax.set_ylabel('Count', fontsize=10)
            ax.set_title(f'{method}\n{n_rhythmic} rhythmic genes', fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3)

        # Q-value distribution
        ax = axes[1, i]
        if 'BH.Q' in results.columns:
            qvalues = results['BH.Q'].values
        elif 'Qvalue' in results.columns:
            qvalues = results['Qvalue'].values
        elif 'Pvalue' in results.columns:
            qvalues = results['Pvalue'].values
        else:
            qvalues = None

        if qvalues is not None:
            ax.hist(np.log10(qvalues + 1e-300), bins=30, color='green', edgecolor='black', alpha=0.7)
            ax.set_xlabel('log10(Q-value)', fontsize=10)
            ax.set_ylabel('Count', fontsize=10)
            ax.set_title(f'{method} - Significance', fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.axvline(np.log10(0.05), color='red', linestyle='--')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"   Comparison dashboard saved to {save_path}")
    else:
        plt.show()

    plt.close()
