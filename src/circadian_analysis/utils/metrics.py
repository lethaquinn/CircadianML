"""
Evaluation metrics for circadian analysis.
"""

import numpy as np
from typing import Dict, Optional
from scipy import stats


def compute_amplitude(expression: np.ndarray) -> float:
    """
    Compute amplitude of expression signal.

    Parameters
    ----------
    expression : np.ndarray
        Expression values

    Returns
    -------
    float
        Amplitude (half of peak-to-trough difference)
    """
    return (np.max(expression) - np.min(expression)) / 2.0


def compute_mesor(expression: np.ndarray) -> float:
    """
    Compute MESOR (Midline Estimating Statistic of Rhythm).

    Parameters
    ----------
    expression : np.ndarray
        Expression values

    Returns
    -------
    float
        MESOR (mean expression level)
    """
    return np.mean(expression)


def compute_relative_amplitude(expression: np.ndarray) -> float:
    """
    Compute relative amplitude (amplitude / MESOR).

    Parameters
    ----------
    expression : np.ndarray
        Expression values

    Returns
    -------
    float
        Relative amplitude
    """
    amplitude = compute_amplitude(expression)
    mesor = compute_mesor(expression)

    if mesor > 0:
        return amplitude / mesor
    else:
        return 0.0


def compute_rsquared(observed: np.ndarray, predicted: np.ndarray) -> float:
    """
    Compute R-squared (coefficient of determination).

    Parameters
    ----------
    observed : np.ndarray
        Observed values
    predicted : np.ndarray
        Predicted values

    Returns
    -------
    float
        R-squared value
    """
    ss_res = np.sum((observed - predicted) ** 2)
    ss_tot = np.sum((observed - np.mean(observed)) ** 2)

    if ss_tot > 0:
        return 1 - (ss_res / ss_tot)
    else:
        return 0.0


def compute_rmse(observed: np.ndarray, predicted: np.ndarray) -> float:
    """
    Compute Root Mean Squared Error.

    Parameters
    ----------
    observed : np.ndarray
        Observed values
    predicted : np.ndarray
        Predicted values

    Returns
    -------
    float
        RMSE value
    """
    return np.sqrt(np.mean((observed - predicted) ** 2))


def compute_mae(observed: np.ndarray, predicted: np.ndarray) -> float:
    """
    Compute Mean Absolute Error.

    Parameters
    ----------
    observed : np.ndarray
        Observed values
    predicted : np.ndarray
        Predicted values

    Returns
    -------
    float
        MAE value
    """
    return np.mean(np.abs(observed - predicted))


def compute_phase_error(true_phase: np.ndarray, pred_phase: np.ndarray, period: float = 24.0) -> float:
    """
    Compute circular phase error.

    Parameters
    ----------
    true_phase : np.ndarray
        True phase values
    pred_phase : np.ndarray
        Predicted phase values
    period : float
        Period length

    Returns
    -------
    float
        Mean absolute phase error
    """
    # Convert to radians
    true_rad = 2 * np.pi * true_phase / period
    pred_rad = 2 * np.pi * pred_phase / period

    # Compute circular distance
    diff = np.angle(np.exp(1j * (true_rad - pred_rad)))

    # Convert back to hours
    diff_hours = np.abs(diff) * period / (2 * np.pi)

    return np.mean(diff_hours)


def compute_circular_correlation(phase1: np.ndarray, phase2: np.ndarray, period: float = 24.0) -> float:
    """
    Compute circular-circular correlation coefficient.

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
    float
        Circular correlation coefficient
    """
    # Convert to radians
    rad1 = 2 * np.pi * phase1 / period
    rad2 = 2 * np.pi * phase2 / period

    # Compute circular correlation
    sin1 = np.sin(rad1 - np.arctan2(np.sum(np.sin(rad1)), np.sum(np.cos(rad1))))
    sin2 = np.sin(rad2 - np.arctan2(np.sum(np.sin(rad2)), np.sum(np.cos(rad2))))

    numerator = np.sum(sin1 * sin2)
    denominator = np.sqrt(np.sum(sin1**2) * np.sum(sin2**2))

    if denominator > 0:
        return numerator / denominator
    else:
        return 0.0


def evaluate_rhythm_detection(
    true_labels: np.ndarray,
    pred_labels: np.ndarray,
    pred_scores: Optional[np.ndarray] = None
) -> Dict[str, float]:
    """
    Evaluate rhythm detection performance.

    Parameters
    ----------
    true_labels : np.ndarray
        True binary labels (1=rhythmic, 0=non-rhythmic)
    pred_labels : np.ndarray
        Predicted binary labels
    pred_scores : np.ndarray, optional
        Prediction scores for ROC/PR curves

    Returns
    -------
    Dict[str, float]
        Dictionary of evaluation metrics
    """
    # Compute confusion matrix components
    tp = np.sum((true_labels == 1) & (pred_labels == 1))
    tn = np.sum((true_labels == 0) & (pred_labels == 0))
    fp = np.sum((true_labels == 0) & (pred_labels == 1))
    fn = np.sum((true_labels == 1) & (pred_labels == 0))

    # Compute metrics
    accuracy = (tp + tn) / len(true_labels) if len(true_labels) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

    metrics = {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1_score": f1,
        "specificity": specificity,
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
    }

    # Add AUC if scores are provided
    if pred_scores is not None:
        try:
            from sklearn.metrics import roc_auc_score
            auc = roc_auc_score(true_labels, pred_scores)
            metrics["auc_roc"] = auc
        except:
            pass

    return metrics


def compute_kendall_tau(x: np.ndarray, y: np.ndarray) -> tuple:
    """
    Compute Kendall's tau correlation coefficient.

    Parameters
    ----------
    x : np.ndarray
        First array
    y : np.ndarray
        Second array

    Returns
    -------
    tuple
        (tau, p-value)
    """
    return stats.kendalltau(x, y)
