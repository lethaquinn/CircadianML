

import numpy as np
import pandas as pd
from typing import Optional, Tuple


def filter_low_expression(
    data: pd.DataFrame,
    threshold: float = 1.0,
    min_samples: int = 3
) -> pd.DataFrame:

   
    n_above_threshold = (data > threshold).sum(axis=1)

    filtered_data = data[n_above_threshold >= min_samples]

    print(f"Filtered from {len(data)} to {len(filtered_data)} genes")

    return filtered_data


def normalize_by_quantile(
    data: pd.DataFrame
) -> pd.DataFrame:

    ranks = data.rank(method='average')
    sorted_data = np.sort(data.values, axis=0)
    mean_values = np.mean(sorted_data, axis=1)

    normalized = ranks.apply(lambda x: mean_values[x.astype(int) - 1], axis=0)

    print("Quantile normalization complete")

    return normalized


def organize_by_timepoints(
    data: pd.DataFrame,
    sample_info: pd.DataFrame,
    time_column: str = 'Time'
) -> pd.DataFrame:


    sorted_samples = sample_info.sort_values(time_column).index
    organized_data = data[sorted_samples]


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
    
    combined_data = []

    for time_point in replicate_info[time_column].unique():

        samples = replicate_info[replicate_info[time_column] == time_point].index


        if method == 'mean':
            combined = data[samples].mean(axis=1)
        elif method == 'median':
            combined = data[samples].median(axis=1)
        else:
            raise ValueError(f"Unknown method: {method}")

        combined_data.append(combined)

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
    
    print("=" * 50)
    print("Starting preprocessing pipeline")
    print("=" * 50)

    # Filter low expression
    print("\n1. Filtering low expression genes...")
    processed = filter_low_expression(data, threshold=filter_threshold)


    if log_transform:
        print("\n2. Log2 transformation...")
        processed = np.log2(processed + 1)

    if normalize:
        print("\n3. Quantile normalization...")
        processed = normalize_by_quantile(processed)

    print("\n" + "=" * 50)
    print("Preprocessing complete!")
    print(f"Final data shape: {processed.shape}")
    print("=" * 50)

    return processed
