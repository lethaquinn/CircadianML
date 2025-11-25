"""
Community detection in circadian networks.
"""

import numpy as np
import pandas as pd
from typing import Dict, List


def greedy_modularity_communities(
    adj_matrix: pd.DataFrame,
    max_iterations: int = 100
) -> Dict[str, int]:
    """
    Detect communities using greedy modularity optimization.

    Parameters
    ----------
    adj_matrix : pd.DataFrame
        Adjacency matrix
    max_iterations : int
        Maximum number of iterations

    Returns
    -------
    dict
        Gene to community mapping
    """
    genes = adj_matrix.index.tolist()
    n = len(genes)

    # Initialize: each node in its own community
    communities = {gene: i for i, gene in enumerate(genes)}

    # Compute total edges
    m = adj_matrix.sum().sum() / 2

    if m == 0:
        return communities

    # Iterative merging
    improved = True
    iteration = 0

    while improved and iteration < max_iterations:
        improved = False
        iteration += 1

        # Try merging communities
        unique_comms = sorted(set(communities.values()))

        for i, comm1 in enumerate(unique_comms):
            for comm2 in unique_comms[i+1:]:
                # Compute modularity change
                delta_q = compute_modularity_delta(
                    adj_matrix, communities, comm1, comm2, m
                )

                if delta_q > 0:
                    # Merge communities
                    for gene, comm in communities.items():
                        if comm == comm2:
                            communities[gene] = comm1
                    improved = True
                    break

            if improved:
                break

    # Renumber communities
    unique_comms = sorted(set(communities.values()))
    comm_map = {old: new for new, old in enumerate(unique_comms)}
    communities = {gene: comm_map[comm] for gene, comm in communities.items()}

    return communities


def compute_modularity_delta(
    adj_matrix: pd.DataFrame,
    communities: Dict[str, int],
    comm1: int,
    comm2: int,
    m: float
) -> float:
    """
    Compute change in modularity from merging two communities.

    Parameters
    ----------
    adj_matrix : pd.DataFrame
        Adjacency matrix
    communities : dict
        Current community assignments
    comm1, comm2 : int
        Communities to merge
    m : float
        Total number of edges

    Returns
    -------
    float
        Change in modularity
    """
    # Get genes in each community
    genes1 = [g for g, c in communities.items() if c == comm1]
    genes2 = [g for g, c in communities.items() if c == comm2]

    # Compute edges between communities
    e_cross = 0
    for g1 in genes1:
        for g2 in genes2:
            e_cross += adj_matrix.loc[g1, g2]

    # Compute degree sums
    k1 = adj_matrix.loc[genes1].sum().sum()
    k2 = adj_matrix.loc[genes2].sum().sum()

    # Modularity change
    delta_q = (e_cross / m) - (k1 * k2) / (2 * m**2)

    return delta_q


def label_propagation_communities(
    adj_matrix: pd.DataFrame,
    max_iterations: int = 100
) -> Dict[str, int]:
    """
    Detect communities using label propagation algorithm.

    Parameters
    ----------
    adj_matrix : pd.DataFrame
        Adjacency matrix
    max_iterations : int
        Maximum number of iterations

    Returns
    -------
    dict
        Gene to community mapping
    """
    genes = adj_matrix.index.tolist()

    # Initialize: each node has unique label
    labels = {gene: i for i, gene in enumerate(genes)}

    # Iterative label propagation
    for iteration in range(max_iterations):
        old_labels = labels.copy()

        # Random order
        np.random.shuffle(genes)

        for gene in genes:
            # Get neighbors
            neighbors = adj_matrix.loc[gene][adj_matrix.loc[gene] > 0].index.tolist()

            if len(neighbors) == 0:
                continue

            # Count neighbor labels
            neighbor_labels = [labels[n] for n in neighbors]
            label_counts = {}
            for label in neighbor_labels:
                label_counts[label] = label_counts.get(label, 0) + 1

            # Assign most common label
            max_count = max(label_counts.values())
            most_common = [l for l, c in label_counts.items() if c == max_count]
            labels[gene] = np.random.choice(most_common)

        # Check convergence
        if labels == old_labels:
            break

    # Renumber communities
    unique_labels = sorted(set(labels.values()))
    label_map = {old: new for new, old in enumerate(unique_labels)}
    labels = {gene: label_map[label] for gene, label in labels.items()}

    return labels


def compute_modularity(
    adj_matrix: pd.DataFrame,
    communities: Dict[str, int]
) -> float:
    """
    Compute modularity of community partition.

    Parameters
    ----------
    adj_matrix : pd.DataFrame
        Adjacency matrix
    communities : dict
        Community assignments

    Returns
    -------
    float
        Modularity score
    """
    genes = adj_matrix.index.tolist()
    m = adj_matrix.sum().sum() / 2

    if m == 0:
        return 0.0

    Q = 0.0

    for i, gene1 in enumerate(genes):
        for gene2 in genes[i:]:
            if communities[gene1] == communities[gene2]:
                A_ij = adj_matrix.loc[gene1, gene2]
                k_i = adj_matrix.loc[gene1].sum()
                k_j = adj_matrix.loc[gene2].sum()

                Q += (A_ij - k_i * k_j / (2 * m))

    Q = Q / (2 * m)

    return Q


def annotate_communities(
    communities: Dict[str, int]
) -> pd.DataFrame:
    """
    Create community annotation dataframe.

    Parameters
    ----------
    communities : dict
        Gene to community mapping

    Returns
    -------
    pd.DataFrame
        Community annotations
    """
    data = []
    for gene, comm in communities.items():
        data.append({'Gene': gene, 'Community': comm})

    df = pd.DataFrame(data)

    # Compute community sizes
    comm_sizes = df['Community'].value_counts().to_dict()
    df['CommunitySize'] = df['Community'].map(comm_sizes)

    return df.sort_values('Community')
