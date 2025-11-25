"""
Cosinor regression for circadian rhythm analysis.
"""

import numpy as np
import pandas as pd
from scipy import stats
from typing import Dict, Tuple, Optional


def cosinor_regression(
    expression: np.ndarray,
    time: np.ndarray,
    period: float = 24.0
) -> Dict[str, float]:
    """
    Fit cosinor model to expression data.

    Model: y = M + A*cos(2π*t/T - φ) + ε
    Rewritten as: y = M + β*cos(ω*t) + γ*sin(ω*t) + ε
    where ω = 2π/T, β = A*cos(φ), γ = A*sin(φ)

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
    dict
        Fitted parameters:
        - mesor: M (mean level)
        - amplitude: A
        - acrophase: φ (phase in hours)
        - rsquared: R²
        - pvalue: p-value for rhythmicity
    """
    omega = 2 * np.pi / period

    # Design matrix
    cos_term = np.cos(omega * time)
    sin_term = np.sin(omega * time)
    X = np.column_stack([np.ones_like(time), cos_term, sin_term])

    # Least squares fit
    coeffs, residuals, rank, s = np.linalg.lstsq(X, expression, rcond=None)

    mesor = coeffs[0]
    beta = coeffs[1]
    gamma = coeffs[2]

    # Compute amplitude and acrophase
    amplitude = np.sqrt(beta**2 + gamma**2)
    acrophase_rad = np.arctan2(gamma, beta)
    acrophase_hours = (acrophase_rad / omega) % period

    # Compute R-squared
    y_pred = X @ coeffs
    ss_res = np.sum((expression - y_pred)**2)
    ss_tot = np.sum((expression - np.mean(expression))**2)
    rsquared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    # F-test for significance
    n = len(expression)
    k = 2  # Number of predictors (cos, sin)
    f_stat = (rsquared / k) / ((1 - rsquared) / (n - k - 1)) if (n - k - 1) > 0 else 0
    pvalue = 1 - stats.f.cdf(f_stat, k, n - k - 1) if f_stat > 0 else 1.0

    return {
        'mesor': mesor,
        'amplitude': amplitude,
        'acrophase': acrophase_hours,
        'beta': beta,
        'gamma': gamma,
        'rsquared': rsquared,
        'pvalue': pvalue,
        'fitted': y_pred
    }


def detect_rhythms_cosinor(
    data: pd.DataFrame,
    time: Optional[np.ndarray] = None,
    period: float = 24.0,
    alpha: float = 0.05
) -> pd.DataFrame:
    """
    Detect rhythms using cosinor regression.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes × timepoints)
    time : np.ndarray, optional
        Time values
    period : float
        Expected period
    alpha : float
        Significance threshold

    Returns
    -------
    pd.DataFrame
        Results table
    """
    from ..utils.helpers import extract_time_from_labels, adjust_pvalues

    if time is None:
        time = extract_time_from_labels(data.columns)

    results = []

    for idx, gene in enumerate(data.index):
        expression = data.loc[gene].values

        # Fit cosinor
        fit_result = cosinor_regression(expression, time, period)

        results.append({
            'GeneID': gene,
            'MESOR': fit_result['mesor'],
            'Amplitude': fit_result['amplitude'],
            'Acrophase': fit_result['acrophase'],
            'Rsquared': fit_result['rsquared'],
            'Pvalue': fit_result['pvalue']
        })

        if (idx + 1) % 100 == 0 or idx == len(data) - 1:
            print(f"   Processed {idx + 1}/{len(data)} genes", end='\r')

    print()

    results_df = pd.DataFrame(results)

    # FDR correction
    results_df['Qvalue'] = adjust_pvalues(results_df['Pvalue'].values, method='fdr_bh')

    # Sort by q-value
    results_df = results_df.sort_values('Qvalue')

    return results_df


def multicomponent_cosinor(
    expression: np.ndarray,
    time: np.ndarray,
    periods: list = [24.0, 12.0]
) -> Dict[str, float]:
    """
    Fit multi-component cosinor model with multiple harmonics.

    Parameters
    ----------
    expression : np.ndarray
        Expression values
    time : np.ndarray
        Time values
    periods : list
        List of period components to fit

    Returns
    -------
    dict
        Fitted parameters for each component
    """
    # Design matrix with multiple harmonics
    X_terms = [np.ones_like(time)]

    for period in periods:
        omega = 2 * np.pi / period
        X_terms.append(np.cos(omega * time))
        X_terms.append(np.sin(omega * time))

    X = np.column_stack(X_terms)

    # Least squares fit
    coeffs = np.linalg.lstsq(X, expression, rcond=None)[0]

    # Extract parameters
    results = {'mesor': coeffs[0]}

    for i, period in enumerate(periods):
        beta = coeffs[1 + 2*i]
        gamma = coeffs[2 + 2*i]

        amplitude = np.sqrt(beta**2 + gamma**2)
        acrophase_rad = np.arctan2(gamma, beta)
        acrophase_hours = (acrophase_rad * period / (2 * np.pi)) % period

        results[f'amplitude_{period}'] = amplitude
        results[f'acrophase_{period}'] = acrophase_hours

    # R-squared
    y_pred = X @ coeffs
    ss_res = np.sum((expression - y_pred)**2)
    ss_tot = np.sum((expression - np.mean(expression))**2)
    results['rsquared'] = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    return results
