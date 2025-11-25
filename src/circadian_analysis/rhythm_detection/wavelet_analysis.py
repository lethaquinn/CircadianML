"""
Wavelet analysis for circadian rhythm detection.
"""

import numpy as np
import pandas as pd
from scipy import signal
from typing import Optional, Tuple, Dict


def continuous_wavelet_transform(
    expression: np.ndarray,
    time: np.ndarray,
    scales: Optional[np.ndarray] = None,
    wavelet: str = 'morlet'
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute continuous wavelet transform.

    Parameters
    ----------
    expression : np.ndarray
        Expression values
    time : np.ndarray
        Time values
    scales : np.ndarray, optional
        Wavelet scales to use
    wavelet : str
        Wavelet type ('morlet', 'ricker', etc.)

    Returns
    -------
    tuple
        (scales, frequencies, coefficients)
    """
    # Default scales corresponding to periods 18-30 hours
    if scales is None:
        periods = np.linspace(18, 30, 50)
        # Approximate scale for Morlet wavelet: scale ≈ period / (4π)
        scales = periods / (4 * np.pi)

    # Choose wavelet function
    if wavelet == 'morlet':
        wavelet_func = signal.morlet2
    elif wavelet == 'ricker':
        wavelet_func = signal.ricker
    else:
        raise ValueError(f"Unknown wavelet: {wavelet}")

    # Compute CWT
    coeffs = []
    for scale in scales:
        if wavelet == 'morlet':
            wavelet_data = wavelet_func(len(expression), scale)
        else:
            wavelet_data = wavelet_func(len(expression), scale)

        # Convolve
        coef = signal.convolve(expression, wavelet_data, mode='same')
        coeffs.append(coef)

    coeffs = np.array(coeffs)

    # Convert scales to frequencies
    # For Morlet wavelet with w0=6: f = 1 / (4π * scale)
    frequencies = 1.0 / (4 * np.pi * scales)

    return scales, frequencies, coeffs


def detect_dominant_period_wavelet(
    expression: np.ndarray,
    time: np.ndarray,
    period_range: Tuple[float, float] = (18.0, 30.0)
) -> Dict[str, float]:
    """
    Detect dominant period using wavelet analysis.

    Parameters
    ----------
    expression : np.ndarray
        Expression values
    time : np.ndarray
        Time values
    period_range : tuple
        Range of periods to search

    Returns
    -------
    dict
        Dominant period and power
    """
    # Generate scales
    periods = np.linspace(period_range[0], period_range[1], 50)
    scales = periods / (4 * np.pi)

    # Compute CWT
    scales, freqs, coeffs = continuous_wavelet_transform(expression, time, scales)

    # Compute time-averaged power spectrum
    power = np.mean(np.abs(coeffs)**2, axis=1)

    # Find dominant period
    max_idx = np.argmax(power)
    dominant_period = 1.0 / freqs[max_idx]
    max_power = power[max_idx]

    return {
        'period': dominant_period,
        'power': max_power,
        'power_spectrum': power,
        'periods': 1.0 / freqs
    }


def wavelet_ridge_analysis(
    expression: np.ndarray,
    time: np.ndarray,
    period_range: Tuple[float, float] = (18.0, 30.0)
) -> pd.DataFrame:
    """
    Perform wavelet ridge analysis to track time-varying periods.

    Parameters
    ----------
    expression : np.ndarray
        Expression values
    time : np.ndarray
        Time values
    period_range : tuple
        Range of periods to search

    Returns
    -------
    pd.DataFrame
        Time-resolved period and amplitude
    """
    periods = np.linspace(period_range[0], period_range[1], 50)
    scales = periods / (4 * np.pi)

    scales, freqs, coeffs = continuous_wavelet_transform(expression, time, scales)

    # For each time point, find the scale with maximum power
    power = np.abs(coeffs)**2
    max_scale_idx = np.argmax(power, axis=0)

    # Extract ridge
    ridge_periods = 1.0 / freqs[max_scale_idx]
    ridge_power = power[max_scale_idx, np.arange(len(time))]

    results = pd.DataFrame({
        'Time': time,
        'Period': ridge_periods,
        'Power': ridge_power
    })

    return results


def detect_rhythms_wavelet(
    data: pd.DataFrame,
    time: Optional[np.ndarray] = None,
    period_range: Tuple[float, float] = (18.0, 30.0)
) -> pd.DataFrame:
    """
    Detect rhythms using wavelet analysis.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes × timepoints)
    time : np.ndarray, optional
        Time values
    period_range : tuple
        Period range to search

    Returns
    -------
    pd.DataFrame
        Results with dominant period and power
    """
    from ..utils.helpers import extract_time_from_labels

    if time is None:
        time = extract_time_from_labels(data.columns)

    results = []

    for idx, gene in enumerate(data.index):
        expression = data.loc[gene].values

        # Detect dominant period
        result = detect_dominant_period_wavelet(expression, time, period_range)

        results.append({
            'GeneID': gene,
            'Period': result['period'],
            'Power': result['power']
        })

        if (idx + 1) % 100 == 0 or idx == len(data) - 1:
            print(f"   Processed {idx + 1}/{len(data)} genes", end='\r')

    print()

    return pd.DataFrame(results)
