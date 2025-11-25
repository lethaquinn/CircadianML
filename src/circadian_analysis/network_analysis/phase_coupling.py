"""
Phase coupling analysis for circadian networks.
"""

import numpy as np
import pandas as pd
from typing import Tuple


def compute_phase_difference(
    phase1: np.ndarray,
    phase2: np.ndarray,
    period: float = 24.0
) -> np.ndarray:
    """
    Compute circular phase difference.

    Parameters
    ----------
    phase1 : np.ndarray
        First phase array
    phase2 : np.ndarray
        Second phase array
    period : float
        Period length

    Returns
    -------
    np.ndarray
        Phase differences in [-period/2, period/2]
    """
    diff = phase1 - phase2

    # Wrap to [-period/2, period/2]
    diff = ((diff + period/2) % period) - period/2

    return diff


def phase_coherence(
    phases: np.ndarray,
    period: float = 24.0
) -> float:
    """
    Compute phase coherence (measure of synchronization).

    Parameters
    ----------
    phases : np.ndarray
        Array of phase values
    period : float
        Period length

    Returns
    -------
    float
        Phase coherence in [0, 1]
    """
    # Convert to circular coordinates
    angles = 2 * np.pi * phases / period

    # Compute mean resultant length
    mean_cos = np.mean(np.cos(angles))
    mean_sin = np.mean(np.sin(angles))

    coherence = np.sqrt(mean_cos**2 + mean_sin**2)

    return coherence


def build_phase_coupling_network(
    rhythmic_results: pd.DataFrame,
    coupling_threshold: float = 3.0,
    period: float = 24.0
) -> pd.DataFrame:
    """
    Build phase coupling network.

    Parameters
    ----------
    rhythmic_results : pd.DataFrame
        Results with phase information (requires 'LAG' or 'Acrophase' column)
    coupling_threshold : float
        Maximum phase difference for coupling (in hours)
    period : float
        Period length

    Returns
    -------
    pd.DataFrame
        Edge list of coupled genes
    """
    # Extract phases
    if 'LAG' in rhythmic_results.columns:
        phase_col = 'LAG'
    elif 'Acrophase' in rhythmic_results.columns:
        phase_col = 'Acrophase'
    else:
        raise ValueError("No phase column found in results")

    genes = rhythmic_results['GeneID'].values
    phases = rhythmic_results[phase_col].values

    edges = []

    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            # Compute phase difference
            phase_diff = compute_phase_difference(
                np.array([phases[i]]),
                np.array([phases[j]]),
                period
            )[0]

            # Check if coupled
            if abs(phase_diff) <= coupling_threshold:
                edges.append({
                    'source': genes[i],
                    'target': genes[j],
                    'phase_diff': phase_diff,
                    'phase_source': phases[i],
                    'phase_target': phases[j]
                })

    return pd.DataFrame(edges)


def identify_phase_clusters(
    rhythmic_results: pd.DataFrame,
    n_clusters: int = 4,
    period: float = 24.0
) -> pd.DataFrame:
    """
    Identify phase clusters (groups of genes with similar phases).

    Parameters
    ----------
    rhythmic_results : pd.DataFrame
        Results with phase information
    n_clusters : int
        Number of phase clusters
    period : float
        Period length

    Returns
    -------
    pd.DataFrame
        Results with cluster assignments
    """
    # Extract phases
    if 'LAG' in rhythmic_results.columns:
        phase_col = 'LAG'
    elif 'Acrophase' in rhythmic_results.columns:
        phase_col = 'Acrophase'
    else:
        raise ValueError("No phase column found")

    phases = rhythmic_results[phase_col].values

    # Simple k-means clustering on circular data
    # Initialize cluster centers evenly
    cluster_centers = np.linspace(0, period, n_clusters, endpoint=False)

    # Iterative assignment
    for iteration in range(10):
        # Assign to nearest cluster
        labels = []
        for phase in phases:
            diffs = compute_phase_difference(
                np.array([phase] * n_clusters),
                cluster_centers,
                period
            )
            labels.append(np.argmin(np.abs(diffs)))

        labels = np.array(labels)

        # Update cluster centers (circular mean)
        new_centers = []
        for k in range(n_clusters):
            cluster_phases = phases[labels == k]
            if len(cluster_phases) > 0:
                # Circular mean
                angles = 2 * np.pi * cluster_phases / period
                mean_angle = np.arctan2(np.mean(np.sin(angles)), np.mean(np.cos(angles)))
                mean_phase = (mean_angle * period / (2 * np.pi)) % period
                new_centers.append(mean_phase)
            else:
                new_centers.append(cluster_centers[k])

        cluster_centers = np.array(new_centers)

    # Add cluster labels to results
    results = rhythmic_results.copy()
    results['PhaseCluster'] = labels

    return results
