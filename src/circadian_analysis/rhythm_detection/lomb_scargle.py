"""
Lomb-Scargle periodogram for circadian rhythm detection.
"""

import numpy as np
import pandas as pd
from scipy import signal
from typing import Optional, Tuple


def lomb_scargle_periodogram(
    expression: np.ndarray,
    time: np.ndarray,
    frequencies: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute Lomb-Scargle periodogram for a single time series.

    Parameters
    ----------
    expression : np.ndarray
        Expression values
    time : np.ndarray
        Time points
    frequencies : np.ndarray, optional
        Frequencies to evaluate

    Returns
    -------
    tuple
        (frequencies, power)
    """
    if frequencies is None:
        # Default: test periods from 18 to 30 hours
        periods = np.linspace(18, 30, 100)
        frequencies = 1.0 / periods

    # Compute Lomb-Scargle periodogram
    power = signal.lombscargle(time, expression, 2 * np.pi * frequencies, normalize=True)

    return frequencies, power


def detect_rhythms_ls(
    data: pd.DataFrame,
    time: Optional[np.ndarray] = None,
    period_range: Tuple[float, float] = (18.0, 30.0),
    n_periods: int = 100
) -> pd.DataFrame:
    """
    Detect rhythms using Lomb-Scargle periodogram.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes × timepoints)
    time : np.ndarray, optional
        Time values (extracted from column names if not provided)
    period_range : tuple
        Range of periods to test (min, max)
    n_periods : int
        Number of periods to test

    Returns
    -------
    pd.DataFrame
        Results with dominant period and power for each gene
    """
    from ..utils.helpers import extract_time_from_labels

    # Extract time if not provided
    if time is None:
        time = extract_time_from_labels(data.columns)

    # Frequency grid
    periods = np.linspace(period_range[0], period_range[1], n_periods)
    frequencies = 1.0 / periods

    results = []

    for idx, gene in enumerate(data.index):
        expression = data.loc[gene].values

        # Compute periodogram
        freqs, power = lomb_scargle_periodogram(expression, time, frequencies)

        # Find dominant period
        max_idx = np.argmax(power)
        dominant_period = 1.0 / freqs[max_idx]
        max_power = power[max_idx]

        results.append({
            'GeneID': gene,
            'Period': dominant_period,
            'Power': max_power,
            'Amplitude': np.std(expression)
        })

        if (idx + 1) % 100 == 0 or idx == len(data) - 1:
            print(f"   Processed {idx + 1}/{len(data)} genes", end='\r')

    print()

    return pd.DataFrame(results)


def compute_significance(
    power: float,
    n_timepoints: int,
    n_independent: Optional[int] = None
) -> float:
    """
    Estimate significance of Lomb-Scargle power.

    Parameters
    ----------
    power : float
        Normalized Lomb-Scargle power
    n_timepoints : int
        Number of time points
    n_independent : int, optional
        Number of independent frequencies tested

    Returns
    -------
    float
        Approximate p-value
    """
    if n_independent is None:
        n_independent = n_timepoints // 2

    # Approximate false alarm probability
    # P(Power > z) ≈ 1 - (1 - e^(-z))^N
    prob = 1.0 - (1.0 - np.exp(-power)) ** n_independent

    return prob
