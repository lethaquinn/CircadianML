"""
Co-expression network analysis for circadian genes.
"""

import numpy as np
import pandas as pd
from scipy import stats
from typing import Optional, Tuple


def compute_correlation_matrix(
    data: pd.DataFrame,
    method: str = 'pearson'
) -> pd.DataFrame:
    """
    Compute pairwise correlation matrix.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes × timepoints)
    method : str
        Correlation method: 'pearson', 'spearman', 'kendall'

    Returns
    -------
    pd.DataFrame
        Correlation matrix (genes × genes)
    """
    if method == 'pearson':
        corr_matrix = data.T.corr(method='pearson')
    elif method == 'spearman':
        corr_matrix = data.T.corr(method='spearman')
    elif method == 'kendall':
        corr_matrix = data.T.corr(method='kendall')
    else:
        raise ValueError(f"Unknown correlation method: {method}")

    return corr_matrix


def build_coexpression_network(
    data: pd.DataFrame,
    threshold: float = 0.7,
    method: str = 'pearson'
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build co-expression network from expression data.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes × timepoints)
    threshold : float
        Correlation threshold for edge creation
    method : str
        Correlation method

    Returns
    -------
    tuple
        (adjacency_matrix, edge_list)
    """
    # Compute correlation matrix
    corr_matrix = compute_correlation_matrix(data, method)

    # Create adjacency matrix based on threshold
    adj_matrix = (corr_matrix.abs() >= threshold).astype(int)

    # Remove self-loops
    np.fill_diagonal(adj_matrix.values, 0)

    # Create edge list
    edges = []
    genes = corr_matrix.index.tolist()

    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            if adj_matrix.iloc[i, j] == 1:
                edges.append({
                    'source': genes[i],
                    'target': genes[j],
                    'weight': corr_matrix.iloc[i, j]
                })

    edge_list = pd.DataFrame(edges)

    return adj_matrix, edge_list


def compute_network_metrics(adj_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Compute basic network metrics for each node.

    Parameters
    ----------
    adj_matrix : pd.DataFrame
        Adjacency matrix

    Returns
    -------
    pd.DataFrame
        Network metrics (degree, clustering coefficient, etc.)
    """
    genes = adj_matrix.index.tolist()
    metrics = []

    for gene in genes:
        # Degree
        degree = adj_matrix.loc[gene].sum()

        # Clustering coefficient
        neighbors = adj_matrix.loc[gene][adj_matrix.loc[gene] == 1].index.tolist()
        if len(neighbors) > 1:
            # Count edges among neighbors
            subgraph = adj_matrix.loc[neighbors, neighbors]
            actual_edges = subgraph.sum().sum() / 2
            possible_edges = len(neighbors) * (len(neighbors) - 1) / 2
            clustering = actual_edges / possible_edges if possible_edges > 0 else 0
        else:
            clustering = 0

        metrics.append({
            'Gene': gene,
            'Degree': degree,
            'Clustering': clustering
        })

    return pd.DataFrame(metrics)


def identify_hub_genes(
    adj_matrix: pd.DataFrame,
    top_n: int = 10
) -> pd.DataFrame:
    """
    Identify hub genes (highly connected genes).

    Parameters
    ----------
    adj_matrix : pd.DataFrame
        Adjacency matrix
    top_n : int
        Number of top hubs to return

    Returns
    -------
    pd.DataFrame
        Top hub genes
    """
    metrics = compute_network_metrics(adj_matrix)
    hubs = metrics.sort_values('Degree', ascending=False).head(top_n)

    return hubs
