"""
Visualization functions for circadian phase analysis.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional


def plot_rhythmic_genes(
    data: pd.DataFrame,
    rhythmic_results: pd.DataFrame,
    top_n: int = 6,
    save_path: Optional[str] = None
):
    """
    Plot expression profiles of top rhythmic genes.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes × timepoints)
    rhythmic_results : pd.DataFrame
        Results from rhythm detection
    top_n : int
        Number of top genes to plot
    save_path : str, optional
        Path to save figure
    """
    # Select top rhythmic genes
    top_genes = rhythmic_results.head(top_n)

    # Create subplots
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()

    for idx, (_, gene_info) in enumerate(top_genes.iterrows()):
        if idx >= top_n:
            break

        gene_id = gene_info['GeneID']
        expression = data.loc[gene_id].values
        time = np.arange(len(expression))

        # Plot
        axes[idx].plot(time, expression, 'o-', linewidth=2, markersize=6)
        axes[idx].set_xlabel('Time (hours)', fontsize=10)
        axes[idx].set_ylabel('Expression', fontsize=10)

        # Add title with statistics
        if 'BH.Q' in gene_info:
            title = f"{gene_id}\nq = {gene_info['BH.Q']:.3e}"
        elif 'Qvalue' in gene_info:
            title = f"{gene_id}\nq = {gene_info['Qvalue']:.3e}"
        else:
            title = gene_id

        axes[idx].set_title(title, fontsize=11, fontweight='bold')
        axes[idx].grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"   Figure saved to {save_path}")
    else:
        plt.show()

    plt.close()


def plot_phase_distribution(
    rhythmic_results: pd.DataFrame,
    phase_col: str = 'LAG',
    bins: int = 24,
    save_path: Optional[str] = None
):
    """
    Plot distribution of phases.

    Parameters
    ----------
    rhythmic_results : pd.DataFrame
        Results with phase information
    phase_col : str
        Column name for phase values
    bins : int
        Number of bins for histogram
    save_path : str, optional
        Path to save figure
    """
    if phase_col not in rhythmic_results.columns:
        # Try alternative column names
        if 'Acrophase' in rhythmic_results.columns:
            phase_col = 'Acrophase'
        else:
            print(f"Warning: No phase column found")
            return

    phases = rhythmic_results[phase_col].values

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Linear histogram
    ax1.hist(phases, bins=bins, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Phase (hours)', fontsize=12)
    ax1.set_ylabel('Number of genes', fontsize=12)
    ax1.set_title('Phase Distribution', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # Circular plot
    theta = 2 * np.pi * phases / 24
    ax2 = plt.subplot(122, projection='polar')
    ax2.hist(theta, bins=bins, color='steelblue', edgecolor='black', alpha=0.7)
    ax2.set_theta_zero_location('N')
    ax2.set_theta_direction(-1)
    ax2.set_title('Circular Phase Distribution', fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"   Figure saved to {save_path}")
    else:
        plt.show()

    plt.close()


def plot_phase_wheel(
    rhythmic_results: pd.DataFrame,
    phase_col: str = 'LAG',
    amplitude_col: Optional[str] = 'AMP',
    top_n: int = 50,
    save_path: Optional[str] = None
):
    """
    Create a phase wheel plot showing genes arranged by phase.

    Parameters
    ----------
    rhythmic_results : pd.DataFrame
        Results with phase information
    phase_col : str
        Column name for phase
    amplitude_col : str, optional
        Column name for amplitude (for point size)
    top_n : int
        Number of top genes to show
    save_path : str, optional
        Path to save figure
    """
    # Auto-detect phase column
    if phase_col not in rhythmic_results.columns:
        if 'Acrophase' in rhythmic_results.columns:
            phase_col = 'Acrophase'

    # Select top genes
    top_genes = rhythmic_results.head(top_n)

    phases = top_genes[phase_col].values
    theta = 2 * np.pi * phases / 24

    # Amplitude for sizing
    if amplitude_col and amplitude_col in top_genes.columns:
        amplitudes = top_genes[amplitude_col].values
        sizes = 100 * (amplitudes / amplitudes.max())
    else:
        sizes = 50

    # Create polar plot
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(111, projection='polar')

    # Plot points
    ax.scatter(theta, np.ones_like(theta), s=sizes, alpha=0.6, c=phases, cmap='twilight')

    # Formatting
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1.2)
    ax.set_yticks([])
    ax.set_xticks(np.linspace(0, 2*np.pi, 24, endpoint=False))
    ax.set_xticklabels([f'{i}h' for i in range(24)])
    ax.set_title(f'Phase Wheel - Top {top_n} Genes', fontsize=16, fontweight='bold', pad=20)

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"   Figure saved to {save_path}")
    else:
        plt.show()

    plt.close()


def plot_amplitude_vs_phase(
    rhythmic_results: pd.DataFrame,
    phase_col: str = 'LAG',
    amplitude_col: str = 'AMP',
    save_path: Optional[str] = None
):
    """
    Plot amplitude vs phase scatter plot.

    Parameters
    ----------
    rhythmic_results : pd.DataFrame
        Results with phase and amplitude
    phase_col : str
        Phase column name
    amplitude_col : str
        Amplitude column name
    save_path : str, optional
        Path to save figure
    """
    # Auto-detect columns
    if phase_col not in rhythmic_results.columns and 'Acrophase' in rhythmic_results.columns:
        phase_col = 'Acrophase'
    if amplitude_col not in rhythmic_results.columns and 'Amplitude' in rhythmic_results.columns:
        amplitude_col = 'Amplitude'

    phases = rhythmic_results[phase_col].values
    amplitudes = rhythmic_results[amplitude_col].values

    fig, ax = plt.subplots(figsize=(10, 6))

    scatter = ax.scatter(phases, amplitudes, alpha=0.6, c=phases, cmap='twilight', s=30)
    ax.set_xlabel('Phase (hours)', fontsize=12)
    ax.set_ylabel('Amplitude', fontsize=12)
    ax.set_title('Amplitude vs Phase', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    plt.colorbar(scatter, ax=ax, label='Phase (hours)')

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"   Figure saved to {save_path}")
    else:
        plt.show()

    plt.close()


def plot_heatmap(
    data: pd.DataFrame,
    rhythmic_results: Optional[pd.DataFrame] = None,
    top_n: int = 50,
    cluster: bool = True,
    save_path: Optional[str] = None
):
    """
    Plot heatmap of rhythmic gene expression.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes × timepoints)
    rhythmic_results : pd.DataFrame, optional
        Rhythmic gene results
    top_n : int
        Number of top genes to show
    cluster : bool
        Whether to cluster genes
    save_path : str, optional
        Path to save figure
    """
    # Select genes
    if rhythmic_results is not None:
        top_genes = rhythmic_results.head(top_n)['GeneID'].values
        plot_data = data.loc[top_genes]
    else:
        plot_data = data.head(top_n)

    # Normalize
    plot_data = (plot_data.T - plot_data.mean(axis=1)) / plot_data.std(axis=1)
    plot_data = plot_data.T

    # Plot
    fig, ax = plt.subplots(figsize=(12, 8))

    if cluster:
        sns.clustermap(plot_data, cmap='RdBu_r', center=0, figsize=(12, 10),
                      cbar_kws={'label': 'Z-score'}, yticklabels=False)
    else:
        sns.heatmap(plot_data, cmap='RdBu_r', center=0, ax=ax,
                   cbar_kws={'label': 'Z-score'}, yticklabels=False)
        ax.set_xlabel('Time', fontsize=12)
        ax.set_ylabel('Genes', fontsize=12)
        ax.set_title(f'Heatmap of Top {top_n} Rhythmic Genes', fontsize=14, fontweight='bold')

    if save_path and not cluster:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"   Figure saved to {save_path}")
    elif not save_path:
        plt.show()

    plt.close()
