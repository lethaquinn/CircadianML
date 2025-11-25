"""
Tests for phase prediction modules.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import pandas as pd
from circadian_analysis.utils.helpers import generate_synthetic_data
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.phase_prediction.ml_models import (
    train_phase_predictor, predict_phase, extract_gene_features
)


def test_feature_extraction():
    """Test feature extraction from expression data."""
    data = generate_synthetic_data(n_genes=50, n_timepoints=24, n_rhythmic=20)

    features, feature_names = extract_gene_features(data)

    assert features.shape[0] == 50, "Number of genes mismatch"
    assert features.shape[1] > 0, "No features extracted"
    assert len(feature_names) == features.shape[1], "Feature names mismatch"


def test_phase_prediction():
    """Test phase prediction model training and prediction."""
    # Generate data
    data = generate_synthetic_data(n_genes=100, n_timepoints=24, n_rhythmic=50, random_state=42)

    # Detect rhythms to get true phases
    results = detect_rhythms(data, period=24.0)
    rhythmic = results[results['BH.Q'] < 0.05].head(50)

    # Extract features
    rhythmic_data = data.loc[rhythmic['GeneID'].values]
    X, _ = extract_gene_features(rhythmic_data)
    y = rhythmic['LAG'].values

    # Split data
    n_train = int(0.8 * len(X))
    X_train, X_test = X[:n_train], X[n_train:]
    y_train, y_test = y[:n_train], y[n_train:]

    # Train model
    model, scaler = train_phase_predictor(X_train, y_train, model_type='random_forest')

    # Predict
    predictions = predict_phase(model, scaler, X_test)

    assert len(predictions) == len(y_test), "Prediction length mismatch"
    assert all(predictions >= 0) and all(predictions < 24), "Predictions out of range"


def test_phase_predictor_types():
    """Test different model types."""
    data = generate_synthetic_data(n_genes=50, n_timepoints=24, n_rhythmic=30, random_state=42)
    X, _ = extract_gene_features(data)
    y = np.random.uniform(0, 24, len(X))

    for model_type in ['random_forest', 'gradient_boosting']:
        model, scaler = train_phase_predictor(
            X, y, model_type=model_type, n_estimators=10
        )
        predictions = predict_phase(model, scaler, X)

        assert len(predictions) == len(y), f"Prediction length mismatch for {model_type}"


if __name__ == "__main__":
    print("Running phase prediction tests...")

    test_feature_extraction()
    print("✓ test_feature_extraction")

    test_phase_prediction()
    print("✓ test_phase_prediction")

    test_phase_predictor_types()
    print("✓ test_phase_predictor_types")

    print("\nAll tests passed! ✓")
