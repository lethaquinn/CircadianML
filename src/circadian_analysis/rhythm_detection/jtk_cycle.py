"""
JTK_CYCLE algorithm for circadian rhythm detection.
Implementation of the Jonckheere-Terpstra-Kendall test for detecting rhythmic signals.

Reference:
Hughes ME, Hogenesch JB, Kornacker K (2010). JTK_CYCLE: an efficient nonparametric
algorithm for detecting rhythmic components in genome-scale data sets.
J Biol Rhythms 25(5):372-80.
"""

import numpy as np
import pandas as pd
from scipy import stats
from typing import Optional, Tuple
from ..utils.helpers import adjust_pvalues, extract_time_from_labels


def detect_rhythms(
    data: pd.DataFrame,
    period: float = 24.0,
    lag_range: Optional[Tuple[int, int]] = None,
    alpha: float = 0.05
) -> pd.DataFrame:
    """
    Detect circadian rhythms using JTK_CYCLE algorithm.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix with genes as rows and timepoints as columns
    period : float
        Expected period length
    lag_range : tuple, optional
        Range of phases to test (default: all phases)
    alpha : float
        Significance threshold for FDR correction

    Returns
    -------
    pd.DataFrame
        Results table with columns:
        - ADJ.P: Adjusted p-value (Benjamini-Hochberg)
        - BH.Q: FDR q-value
        - PER: Period
        - LAG: Phase (lag)
        - AMP: Amplitude
    """
    n_genes = data.shape[0]
    n_timepoints = data.shape[1]

    # Default lag range: test all possible phases
    if lag_range is None:
        lag_range = (0, n_timepoints - 1)

    results = []

    for idx, gene in enumerate(data.index):
        expression = data.loc[gene].values

        # Run JTK test
        pval, best_lag, best_tau = jtk_test(expression, period, lag_range)

        # Compute amplitude
        amplitude = (np.max(expression) - np.min(expression)) / 2.0

        results.append({
            'GeneID': gene,
            'P': pval,
            'TAU': best_tau,
            'PER': period,
            'LAG': best_lag,
            'AMP': amplitude
        })

        # Progress indicator
        if (idx + 1) % 100 == 0 or idx == n_genes - 1:
            print(f"   Processed {idx + 1}/{n_genes} genes", end='\r')

    print()  # New line after progress

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # FDR correction
    results_df['BH.Q'] = adjust_pvalues(results_df['P'].values, method='fdr_bh')
    results_df['ADJ.P'] = results_df['BH.Q']

    # Sort by adjusted p-value
    results_df = results_df.sort_values('BH.Q')

    return results_df


def jtk_test(
    expression: np.ndarray,
    period: float,
    lag_range: Tuple[int, int]
) -> Tuple[float, float, float]:
    """
    Perform JTK test on a single gene.

    Parameters
    ----------
    expression : np.ndarray
        Expression values
    period : float
        Expected period
    lag_range : tuple
        Range of phases to test

    Returns
    -------
    tuple
        (p-value, best_lag, best_tau)
    """
    n_timepoints = len(expression)

    # Generate reference waveforms for each phase
    best_tau = -np.inf
    best_lag = 0
    best_pval = 1.0

    for lag in range(lag_range[0], lag_range[1] + 1):
        # Generate cosine reference waveform
        time = np.arange(n_timepoints)
        reference = np.cos(2 * np.pi * (time - lag) / period)

        # Compute Kendall's tau
        tau, pval = stats.kendalltau(expression, reference)

        # Keep track of best correlation
        if tau > best_tau:
            best_tau = tau
            best_lag = lag
            best_pval = pval

    return best_pval, best_lag, best_tau


def jtk_cycle_comprehensive(
    data: pd.DataFrame,
    periods: list = [24.0],
    alpha: float = 0.05
) -> pd.DataFrame:
    """
    Run JTK_CYCLE with multiple period lengths.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix
    periods : list
        List of periods to test
    alpha : float
        Significance threshold

    Returns
    -------
    pd.DataFrame
        Results for best period per gene
    """
    all_results = []

    for period in periods:
        print(f"\n Testing period: {period}h")
        results = detect_rhythms(data, period=period, alpha=alpha)
        results['TestPeriod'] = period
        all_results.append(results)

    # Combine results
    combined = pd.concat(all_results, ignore_index=True)

    # For each gene, keep result with lowest p-value
    best_results = combined.loc[combined.groupby('GeneID')['BH.Q'].idxmin()]

    return best_results.sort_values('BH.Q')


def filter_rhythmic_genes(results: pd.DataFrame, qvalue_cutoff: float = 0.05) -> pd.DataFrame:
    """
    Filter for significantly rhythmic genes.

    Parameters
    ----------
    results : pd.DataFrame
        JTK_CYCLE results
    qvalue_cutoff : float
        Q-value threshold

    Returns
    -------
    pd.DataFrame
        Filtered results
    """
    return results[results['BH.Q'] < qvalue_cutoff].copy()
