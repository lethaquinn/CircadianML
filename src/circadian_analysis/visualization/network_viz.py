"""
Network visualization functions.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional, Dict


def plot_network_simple(
    edge_list: pd.DataFrame,
    node_positions: Optional[Dict] = None,
    save_path: Optional[str] = None
):
    """
    Simple network visualization.

    Parameters
    ----------
    edge_list : pd.DataFrame
        Edge list with 'source' and 'target' columns
    node_positions : dict, optional
        Custom node positions
    save_path : str, optional
        Path to save figure
    """
    # Extract nodes
    nodes = set(edge_list['source'].tolist() + edge_list['target'].tolist())

    # Create positions if not provided
    if node_positions is None:
        # Circular layout
        n = len(nodes)
        angles = np.linspace(0, 2*np.pi, n, endpoint=False)
        node_positions = {}
        for i, node in enumerate(sorted(nodes)):
            node_positions[node] = (np.cos(angles[i]), np.sin(angles[i]))

    # Plot
    fig, ax = plt.subplots(figsize=(12, 12))

    # Draw edges
    for _, edge in edge_list.iterrows():
        source = edge['source']
        target = edge['target']

        if source in node_positions and target in node_positions:
            x1, y1 = node_positions[source]
            x2, y2 = node_positions[target]

            ax.plot([x1, x2], [y1, y2], 'gray', alpha=0.3, linewidth=0.5)

    # Draw nodes
    for node, (x, y) in node_positions.items():
        ax.plot(x, y, 'o', markersize=8, color='steelblue')

    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Co-expression Network', fontsize=16, fontweight='bold')

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"   Figure saved to {save_path}")
    else:
        plt.show()

    plt.close()


def plot_network_circular(
    edge_list: pd.DataFrame,
    phase_info: Optional[pd.DataFrame] = None,
    save_path: Optional[str] = None
):
    """
    Circular network layout based on phase.

    Parameters
    ----------
    edge_list : pd.DataFrame
        Edge list
    phase_info : pd.DataFrame, optional
        Phase information for nodes
    save_path : str, optional
        Path to save figure
    """
    # Extract unique nodes
    nodes = set(edge_list['source'].tolist() + edge_list['target'].tolist())

    # Position nodes based on phase if available
    if phase_info is not None and 'LAG' in phase_info.columns:
        node_positions = {}
        for node in nodes:
            if node in phase_info['GeneID'].values:
                phase = phase_info[phase_info['GeneID'] == node]['LAG'].values[0]
                angle = 2 * np.pi * phase / 24
                node_positions[node] = (np.cos(angle), np.sin(angle))
    else:
        # Default circular layout
        n = len(nodes)
        angles = np.linspace(0, 2*np.pi, n, endpoint=False)
        node_positions = {}
        for i, node in enumerate(sorted(nodes)):
            node_positions[node] = (np.cos(angles[i]), np.sin(angles[i]))

    # Plot
    plot_network_simple(edge_list, node_positions, save_path)


def plot_degree_distribution(
    network_metrics: pd.DataFrame,
    save_path: Optional[str] = None
):
    """
    Plot degree distribution of network.

    Parameters
    ----------
    network_metrics : pd.DataFrame
        Network metrics with 'Degree' column
    save_path : str, optional
        Path to save figure
    """
    degrees = network_metrics['Degree'].values

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Linear scale
    ax1.hist(degrees, bins=30, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Degree', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Degree Distribution', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # Log-log scale
    degree_counts = {}
    for d in degrees:
        degree_counts[d] = degree_counts.get(d, 0) + 1

    degrees_sorted = sorted(degree_counts.keys())
    counts = [degree_counts[d] for d in degrees_sorted]

    ax2.loglog(degrees_sorted, counts, 'o', markersize=8, color='steelblue', alpha=0.7)
    ax2.set_xlabel('Degree (log)', fontsize=12)
    ax2.set_ylabel('Frequency (log)', fontsize=12)
    ax2.set_title('Degree Distribution (log-log)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"   Figure saved to {save_path}")
    else:
        plt.show()

    plt.close()


def plot_community_sizes(
    community_df: pd.DataFrame,
    save_path: Optional[str] = None
):
    """
    Plot community size distribution.

    Parameters
    ----------
    community_df : pd.DataFrame
        Community annotations with 'Community' column
    save_path : str, optional
        Path to save figure
    """
    comm_sizes = community_df['Community'].value_counts().sort_values(ascending=False)

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.bar(range(len(comm_sizes)), comm_sizes.values, color='steelblue', alpha=0.7)
    ax.set_xlabel('Community ID', fontsize=12)
    ax.set_ylabel('Number of genes', fontsize=12)
    ax.set_title('Community Size Distribution', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"   Figure saved to {save_path}")
    else:
        plt.show()

    plt.close()
