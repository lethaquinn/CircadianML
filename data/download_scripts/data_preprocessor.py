"""
Data preprocessing utilities for circadian transcriptomics.
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple


def filter_low_expression(
    data: pd.DataFrame,
    threshold: float = 1.0,
    min_samples: int = 3
) -> pd.DataFrame:
    """
    Filter out genes with low expression.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes × samples)
    threshold : float
        Minimum expression threshold
    min_samples : int
        Minimum number of samples above threshold

    Returns
    -------
    pd.DataFrame
        Filtered expression matrix
    """
    # Count samples above threshold for each gene
    n_above_threshold = (data > threshold).sum(axis=1)

    # Keep genes with sufficient samples above threshold
    filtered_data = data[n_above_threshold >= min_samples]

    print(f"Filtered from {len(data)} to {len(filtered_data)} genes")

    return filtered_data


def normalize_by_quantile(
    data: pd.DataFrame
) -> pd.DataFrame:
    """
    Quantile normalization across samples.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes × samples)

    Returns
    -------
    pd.DataFrame
        Normalized expression matrix
    """
    # Rank genes within each sample
    ranks = data.rank(method='average')

    # Compute mean across samples at each rank
    sorted_data = np.sort(data.values, axis=0)
    mean_values = np.mean(sorted_data, axis=1)

    # Map ranks back to mean values
    normalized = ranks.apply(lambda x: mean_values[x.astype(int) - 1], axis=0)

    print("Quantile normalization complete")

    return normalized


def organize_by_timepoints(
    data: pd.DataFrame,
    sample_info: pd.DataFrame,
    time_column: str = 'Time'
) -> pd.DataFrame:
    """
    Organize samples by timepoints.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes × samples)
    sample_info : pd.DataFrame
        Sample metadata with time information
    time_column : str
        Column name for time values

    Returns
    -------
    pd.DataFrame
        Expression matrix organized by time
    """
    # Sort samples by time
    sorted_samples = sample_info.sort_values(time_column).index

    # Reorder columns
    organized_data = data[sorted_samples]

    # Rename columns with time labels
    time_labels = [f"ZT{int(t):02d}" for t in sample_info.loc[sorted_samples, time_column]]
    organized_data.columns = time_labels

    print(f"Organized {len(organized_data.columns)} samples by timepoint")

    return organized_data


def handle_replicates(
    data: pd.DataFrame,
    replicate_info: pd.DataFrame,
    time_column: str = 'Time',
    method: str = 'mean'
) -> pd.DataFrame:
    """
    Combine technical/biological replicates.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix
    replicate_info : pd.DataFrame
        Information about replicates
    time_column : str
        Column with time information
    method : str
        Aggregation method: 'mean', 'median'

    Returns
    -------
    pd.DataFrame
        Expression matrix with replicates combined
    """
    combined_data = []

    for time_point in replicate_info[time_column].unique():
        # Get samples at this timepoint
        samples = replicate_info[replicate_info[time_column] == time_point].index

        # Aggregate
        if method == 'mean':
            combined = data[samples].mean(axis=1)
        elif method == 'median':
            combined = data[samples].median(axis=1)
        else:
            raise ValueError(f"Unknown method: {method}")

        combined_data.append(combined)

    # Create combined DataFrame
    result = pd.DataFrame(combined_data).T
    result.columns = [f"ZT{int(t):02d}" for t in sorted(replicate_info[time_column].unique())]

    print(f"Combined replicates: {len(data.columns)} samples -> {len(result.columns)} timepoints")

    return result


def preprocess_pipeline(
    data: pd.DataFrame,
    filter_threshold: float = 1.0,
    normalize: bool = True,
    log_transform: bool = True
) -> pd.DataFrame:
    """
    Complete preprocessing pipeline.

    Parameters
    ----------
    data : pd.DataFrame
        Raw expression matrix
    filter_threshold : float
        Expression threshold for filtering
    normalize : bool
        Whether to perform quantile normalization
    log_transform : bool
        Whether to log2 transform

    Returns
    -------
    pd.DataFrame
        Preprocessed expression matrix
    """
    print("=" * 50)
    print("Starting preprocessing pipeline")
    print("=" * 50)

    # Filter low expression
    print("\n1. Filtering low expression genes...")
    processed = filter_low_expression(data, threshold=filter_threshold)

    # Log transform
    if log_transform:
        print("\n2. Log2 transformation...")
        processed = np.log2(processed + 1)

    # Normalize
    if normalize:
        print("\n3. Quantile normalization...")
        processed = normalize_by_quantile(processed)

    print("\n" + "=" * 50)
    print("Preprocessing complete!")
    print(f"Final data shape: {processed.shape}")
    print("=" * 50)

    return processed
