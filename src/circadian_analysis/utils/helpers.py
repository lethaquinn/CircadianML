"""
Helper functions for circadian transcriptomics analysis.
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, Union


def generate_synthetic_data(
    n_genes: int = 100,
    n_timepoints: int = 24,
    n_rhythmic: int = 30,
    period: float = 24.0,
    noise_level: float = 0.3,
    random_state: Optional[int] = 42
) -> pd.DataFrame:
    """
    Generate synthetic circadian expression data for testing.

    Parameters
    ----------
    n_genes : int
        Number of genes to generate
    n_timepoints : int
        Number of timepoints
    n_rhythmic : int
        Number of rhythmic genes
    period : float
        Period length in hours
    noise_level : float
        Standard deviation of Gaussian noise
    random_state : int, optional
        Random seed for reproducibility

    Returns
    -------
    pd.DataFrame
        Expression matrix with genes as rows and timepoints as columns
    """
    if random_state is not None:
        np.random.seed(random_state)

    # Time vector
    time = np.arange(n_timepoints)

    # Initialize expression matrix
    data = np.zeros((n_genes, n_timepoints))

    # Generate rhythmic genes
    for i in range(n_rhythmic):
        # Random amplitude and phase
        amplitude = np.random.uniform(2, 5)
        phase = np.random.uniform(0, 2 * np.pi)
        baseline = np.random.uniform(5, 10)

        # Cosine wave
        signal = baseline + amplitude * np.cos(2 * np.pi * time / period - phase)

        # Add noise
        noise = np.random.normal(0, noise_level, n_timepoints)
        data[i, :] = signal + noise

    # Generate non-rhythmic genes
    for i in range(n_rhythmic, n_genes):
        baseline = np.random.uniform(5, 10)
        noise = np.random.normal(0, noise_level * 3, n_timepoints)
        data[i, :] = baseline + noise

    # Create DataFrame
    gene_names = [f"Gene_{i+1:03d}" for i in range(n_genes)]
    time_names = [f"ZT{t:02d}" for t in time]
    df = pd.DataFrame(data, index=gene_names, columns=time_names)

    return df


def normalize_expression(
    data: pd.DataFrame,
    method: str = "zscore"
) -> pd.DataFrame:
    """
    Normalize expression data.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes × timepoints)
    method : str
        Normalization method: 'zscore', 'minmax', or 'log2'

    Returns
    -------
    pd.DataFrame
        Normalized expression matrix
    """
    if method == "zscore":
        # Z-score normalization per gene
        normalized = (data - data.mean(axis=1).values.reshape(-1, 1)) / data.std(axis=1).values.reshape(-1, 1)
    elif method == "minmax":
        # Min-max normalization per gene
        min_val = data.min(axis=1).values.reshape(-1, 1)
        max_val = data.max(axis=1).values.reshape(-1, 1)
        normalized = (data - min_val) / (max_val - min_val)
    elif method == "log2":
        # Log2 transformation
        normalized = np.log2(data + 1)
    else:
        raise ValueError(f"Unknown normalization method: {method}")

    return normalized


def extract_time_from_labels(labels: list) -> np.ndarray:
    """
    Extract numeric time values from time labels.

    Parameters
    ----------
    labels : list
        Time labels (e.g., ['ZT00', 'ZT01', ..., 'ZT23'])

    Returns
    -------
    np.ndarray
        Numeric time values
    """
    time = []
    for label in labels:
        # Try to extract numbers from label
        import re
        numbers = re.findall(r'\d+', str(label))
        if numbers:
            time.append(float(numbers[0]))
        else:
            time.append(0)

    return np.array(time)


def adjust_pvalues(pvalues: np.ndarray, method: str = "fdr_bh") -> np.ndarray:
    """
    Adjust p-values for multiple testing.

    Parameters
    ----------
    pvalues : np.ndarray
        Array of p-values
    method : str
        Correction method: 'fdr_bh' (Benjamini-Hochberg) or 'bonferroni'

    Returns
    -------
    np.ndarray
        Adjusted p-values
    """
    n = len(pvalues)
    if n == 0:
        return pvalues

    if method == "fdr_bh":
        # Benjamini-Hochberg FDR correction
        sorted_indices = np.argsort(pvalues)
        sorted_pvalues = pvalues[sorted_indices]

        adjusted = np.zeros(n)
        for i in range(n):
            adjusted[sorted_indices[i]] = min(sorted_pvalues[i] * n / (i + 1), 1.0)

        # Enforce monotonicity
        for i in range(n - 1, 0, -1):
            if adjusted[sorted_indices[i - 1]] > adjusted[sorted_indices[i]]:
                adjusted[sorted_indices[i - 1]] = adjusted[sorted_indices[i]]

    elif method == "bonferroni":
        # Bonferroni correction
        adjusted = np.minimum(pvalues * n, 1.0)
    else:
        raise ValueError(f"Unknown correction method: {method}")

    return adjusted


def compute_phase(expression: np.ndarray, time: np.ndarray, period: float = 24.0) -> float:
    """
    Compute phase of a circadian signal using cosinor regression.

    Parameters
    ----------
    expression : np.ndarray
        Expression values
    time : np.ndarray
        Time values
    period : float
        Period length

    Returns
    -------
    float
        Phase in hours
    """
    # Cosinor regression: y = M + A*cos(2π*t/T - φ)
    # Rewrite as: y = M + β*cos(2π*t/T) + γ*sin(2π*t/T)
    # where β = A*cos(φ), γ = A*sin(φ)
    # Then φ = atan2(γ, β)

    omega = 2 * np.pi / period
    cos_term = np.cos(omega * time)
    sin_term = np.sin(omega * time)

    # Simple linear regression
    X = np.column_stack([np.ones_like(time), cos_term, sin_term])
    try:
        coeffs = np.linalg.lstsq(X, expression, rcond=None)[0]
        beta = coeffs[1]
        gamma = coeffs[2]

        # Phase in radians
        phase_rad = np.arctan2(gamma, beta)

        # Convert to hours
        phase_hours = (phase_rad / omega) % period

        return phase_hours
    except:
        return 0.0


def create_circular_time(time: np.ndarray, period: float = 24.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert linear time to circular coordinates for circular statistics.

    Parameters
    ----------
    time : np.ndarray
        Linear time values
    period : float
        Period length

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Cosine and sine components
    """
    angle = 2 * np.pi * time / period
    return np.cos(angle), np.sin(angle)
