"""
Machine learning models for circadian phase prediction.
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, KFold
from sklearn.preprocessing import StandardScaler
from typing import Dict, Optional, Tuple


def train_phase_predictor(
    X: np.ndarray,
    y: np.ndarray,
    model_type: str = 'random_forest',
    **model_params
) -> Tuple[object, object]:
    """
    Train a phase prediction model.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix (n_samples × n_features)
    y : np.ndarray
        Phase labels (in hours)
    model_type : str
        Model type: 'random_forest', 'gradient_boosting'
    **model_params
        Additional parameters for the model

    Returns
    -------
    tuple
        (trained_model, scaler)
    """
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Initialize model
    if model_type == 'random_forest':
        model = RandomForestRegressor(
            n_estimators=model_params.get('n_estimators', 100),
            max_depth=model_params.get('max_depth', None),
            random_state=model_params.get('random_state', 42),
            n_jobs=model_params.get('n_jobs', -1)
        )
    elif model_type == 'gradient_boosting':
        model = GradientBoostingRegressor(
            n_estimators=model_params.get('n_estimators', 100),
            max_depth=model_params.get('max_depth', 3),
            learning_rate=model_params.get('learning_rate', 0.1),
            random_state=model_params.get('random_state', 42)
        )
    else:
        raise ValueError(f"Unknown model type: {model_type}")

    # Train model
    model.fit(X_scaled, y)

    return model, scaler


def predict_phase(
    model: object,
    scaler: object,
    X: np.ndarray
) -> np.ndarray:
    """
    Predict circadian phase.

    Parameters
    ----------
    model : object
        Trained model
    scaler : object
        Fitted scaler
    X : np.ndarray
        Feature matrix

    Returns
    -------
    np.ndarray
        Predicted phases (in hours)
    """
    X_scaled = scaler.transform(X)
    predictions = model.predict(X_scaled)

    # Ensure predictions are in [0, 24) range
    predictions = predictions % 24

    return predictions


def evaluate_phase_predictor(
    X: np.ndarray,
    y: np.ndarray,
    model_type: str = 'random_forest',
    cv: int = 5,
    **model_params
) -> Dict[str, float]:
    """
    Evaluate phase prediction model using cross-validation.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix
    y : np.ndarray
        True phases
    model_type : str
        Model type
    cv : int
        Number of cross-validation folds
    **model_params
        Model parameters

    Returns
    -------
    dict
        Evaluation metrics
    """
    from ..utils.metrics import compute_phase_error, compute_mae

    # Perform cross-validation
    kfold = KFold(n_splits=cv, shuffle=True, random_state=42)

    mae_scores = []
    phase_errors = []

    for train_idx, test_idx in kfold.split(X):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Train model
        model, scaler = train_phase_predictor(X_train, y_train, model_type, **model_params)

        # Predict
        y_pred = predict_phase(model, scaler, X_test)

        # Compute metrics
        mae = compute_mae(y_test, y_pred)
        phase_error = compute_phase_error(y_test, y_pred)

        mae_scores.append(mae)
        phase_errors.append(phase_error)

    return {
        'mean_mae': np.mean(mae_scores),
        'std_mae': np.std(mae_scores),
        'mean_phase_error': np.mean(phase_errors),
        'std_phase_error': np.std(phase_errors)
    }


def extract_gene_features(
    data: pd.DataFrame,
    rhythmic_results: Optional[pd.DataFrame] = None
) -> Tuple[np.ndarray, list]:
    """
    Extract features from gene expression data for phase prediction.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes × timepoints)
    rhythmic_results : pd.DataFrame, optional
        Results from rhythm detection (includes amplitude, phase, etc.)

    Returns
    -------
    tuple
        (feature_matrix, feature_names)
    """
    features = []
    feature_names = []

    # Raw expression values
    for col in data.columns:
        features.append(data[col].values)
        feature_names.append(f'expr_{col}')

    # Statistical features
    features.append(data.mean(axis=1).values)
    feature_names.append('mean_expr')

    features.append(data.std(axis=1).values)
    feature_names.append('std_expr')

    features.append(data.max(axis=1).values)
    feature_names.append('max_expr')

    features.append(data.min(axis=1).values)
    feature_names.append('min_expr')

    # If rhythmic results available, add rhythm features
    if rhythmic_results is not None:
        if 'AMP' in rhythmic_results.columns:
            amp_map = dict(zip(rhythmic_results['GeneID'], rhythmic_results['AMP']))
            features.append(data.index.map(amp_map).fillna(0).values)
            feature_names.append('amplitude')

        if 'TAU' in rhythmic_results.columns:
            tau_map = dict(zip(rhythmic_results['GeneID'], rhythmic_results['TAU']))
            features.append(data.index.map(tau_map).fillna(0).values)
            feature_names.append('kendall_tau')

    # Stack features
    X = np.column_stack(features)

    return X, feature_names


def get_feature_importance(
    model: object,
    feature_names: list,
    top_n: int = 20
) -> pd.DataFrame:
    """
    Get feature importance from trained model.

    Parameters
    ----------
    model : object
        Trained model
    feature_names : list
        List of feature names
    top_n : int
        Number of top features to return

    Returns
    -------
    pd.DataFrame
        Feature importance ranking
    """
    if hasattr(model, 'feature_importances_'):
        importance = model.feature_importances_
    else:
        return pd.DataFrame()

    importance_df = pd.DataFrame({
        'Feature': feature_names,
        'Importance': importance
    })

    importance_df = importance_df.sort_values('Importance', ascending=False)

    return importance_df.head(top_n)
