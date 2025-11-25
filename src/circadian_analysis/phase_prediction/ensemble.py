"""
Ensemble methods for phase prediction.
"""

import numpy as np
import pandas as pd
from typing import List, Tuple, Dict
from .ml_models import train_phase_predictor, predict_phase


class PhaseEnsemble:
    """
    Ensemble of phase prediction models.
    """

    def __init__(self, model_types: List[str] = None):
        """
        Initialize ensemble.

        Parameters
        ----------
        model_types : list
            List of model types to include in ensemble
        """
        if model_types is None:
            model_types = ['random_forest', 'gradient_boosting']

        self.model_types = model_types
        self.models = []
        self.scalers = []
        self.weights = None

    def fit(self, X: np.ndarray, y: np.ndarray, **model_params):
        """
        Train all models in the ensemble.

        Parameters
        ----------
        X : np.ndarray
            Training features
        y : np.ndarray
            Training targets
        **model_params
            Model parameters
        """
        self.models = []
        self.scalers = []

        for model_type in self.model_types:
            print(f"   Training {model_type}...")
            model, scaler = train_phase_predictor(X, y, model_type, **model_params)
            self.models.append(model)
            self.scalers.append(scaler)

        # Equal weights by default
        self.weights = np.ones(len(self.models)) / len(self.models)

    def predict(self, X: np.ndarray, method: str = 'weighted_mean') -> np.ndarray:
        """
        Make ensemble predictions.

        Parameters
        ----------
        X : np.ndarray
            Input features
        method : str
            Ensemble method: 'weighted_mean', 'median', 'circular_mean'

        Returns
        -------
        np.ndarray
            Ensemble predictions
        """
        # Collect predictions from all models
        predictions = []
        for model, scaler in zip(self.models, self.scalers):
            pred = predict_phase(model, scaler, X)
            predictions.append(pred)

        predictions = np.array(predictions)

        # Combine predictions
        if method == 'weighted_mean':
            # Weighted average
            ensemble_pred = np.average(predictions, axis=0, weights=self.weights)
        elif method == 'median':
            # Median
            ensemble_pred = np.median(predictions, axis=0)
        elif method == 'circular_mean':
            # Circular mean for phase values
            ensemble_pred = circular_mean(predictions, axis=0)
        else:
            raise ValueError(f"Unknown ensemble method: {method}")

        # Ensure in [0, 24) range
        ensemble_pred = ensemble_pred % 24

        return ensemble_pred

    def optimize_weights(self, X_val: np.ndarray, y_val: np.ndarray):
        """
        Optimize ensemble weights on validation set.

        Parameters
        ----------
        X_val : np.ndarray
            Validation features
        y_val : np.ndarray
            Validation targets
        """
        from ..utils.metrics import compute_phase_error

        # Get predictions from each model
        predictions = []
        for model, scaler in zip(self.models, self.scalers):
            pred = predict_phase(model, scaler, X_val)
            predictions.append(pred)

        predictions = np.array(predictions)

        # Grid search for optimal weights
        best_error = float('inf')
        best_weights = None

        # Simple grid search (can be improved with optimization)
        for w1 in np.linspace(0, 1, 11):
            for w2 in np.linspace(0, 1 - w1, 11):
                weights = np.array([w1, w2])
                if len(self.models) > 2:
                    weights = np.append(weights, 1 - np.sum(weights))

                weights = weights / np.sum(weights)  # Normalize

                # Compute weighted prediction
                ensemble_pred = np.average(predictions, axis=0, weights=weights)
                ensemble_pred = ensemble_pred % 24

                # Compute error
                error = compute_phase_error(y_val, ensemble_pred)

                if error < best_error:
                    best_error = error
                    best_weights = weights

        self.weights = best_weights
        print(f"   Optimized weights: {self.weights}")
        print(f"   Validation error: {best_error:.2f} hours")


def circular_mean(phases: np.ndarray, axis: int = 0, period: float = 24.0) -> np.ndarray:
    """
    Compute circular mean of phase values.

    Parameters
    ----------
    phases : np.ndarray
        Phase values
    axis : int
        Axis along which to compute mean
    period : float
        Period length

    Returns
    -------
    np.ndarray
        Circular mean
    """
    # Convert to angles
    angles = 2 * np.pi * phases / period

    # Compute mean of complex exponentials
    mean_cos = np.mean(np.cos(angles), axis=axis)
    mean_sin = np.mean(np.sin(angles), axis=axis)

    # Convert back to phase
    mean_angle = np.arctan2(mean_sin, mean_cos)
    mean_phase = (mean_angle * period / (2 * np.pi)) % period

    return mean_phase


def stacking_ensemble(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    base_models: List[str] = None,
    meta_model: str = 'random_forest'
) -> np.ndarray:
    """
    Stacking ensemble for phase prediction.

    Parameters
    ----------
    X_train : np.ndarray
        Training features
    y_train : np.ndarray
        Training targets
    X_test : np.ndarray
        Test features
    base_models : list
        List of base model types
    meta_model : str
        Meta-learner model type

    Returns
    -------
    np.ndarray
        Predictions
    """
    if base_models is None:
        base_models = ['random_forest', 'gradient_boosting']

    # Train base models and get predictions
    base_predictions_train = []
    base_predictions_test = []

    for model_type in base_models:
        print(f"   Training base model: {model_type}")
        model, scaler = train_phase_predictor(X_train, y_train, model_type)

        pred_train = predict_phase(model, scaler, X_train)
        pred_test = predict_phase(model, scaler, X_test)

        base_predictions_train.append(pred_train)
        base_predictions_test.append(pred_test)

    # Create meta-features
    X_meta_train = np.column_stack(base_predictions_train)
    X_meta_test = np.column_stack(base_predictions_test)

    # Train meta-learner
    print(f"   Training meta-learner: {meta_model}")
    meta_model_obj, meta_scaler = train_phase_predictor(X_meta_train, y_train, meta_model)

    # Final predictions
    predictions = predict_phase(meta_model_obj, meta_scaler, X_meta_test)

    return predictions
