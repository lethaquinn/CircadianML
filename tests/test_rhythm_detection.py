"""
Tests for rhythm detection modules.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import pandas as pd
import pytest
from circadian_analysis.utils.helpers import generate_synthetic_data
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms, jtk_test
from circadian_analysis.rhythm_detection.cosinor import cosinor_regression, detect_rhythms_cosinor


def test_synthetic_data_generation():
    """Test synthetic data generation."""
    data = generate_synthetic_data(n_genes=100, n_timepoints=24, n_rhythmic=30)

    assert data.shape == (100, 24), "Data shape mismatch"
    assert len(data.index) == 100, "Number of genes mismatch"
    assert len(data.columns) == 24, "Number of timepoints mismatch"


def test_jtk_cycle_detection():
    """Test JTK_CYCLE rhythm detection."""
    # Generate data with known rhythmic genes
    data = generate_synthetic_data(n_genes=50, n_timepoints=24, n_rhythmic=20, random_state=42)

    # Run JTK_CYCLE
    results = detect_rhythms(data, period=24.0)

    assert len(results) == 50, "Results length mismatch"
    assert 'GeneID' in results.columns, "Missing GeneID column"
    assert 'BH.Q' in results.columns, "Missing BH.Q column"
    assert 'LAG' in results.columns, "Missing LAG column"

    # Check that we detect some rhythmic genes
    rhythmic = results[results['BH.Q'] < 0.05]
    assert len(rhythmic) > 0, "Should detect at least some rhythmic genes"


def test_jtk_test():
    """Test single gene JTK test."""
    # Create a perfect cosine wave
    time = np.arange(24)
    expression = np.cos(2 * np.pi * time / 24)

    pval, lag, tau = jtk_test(expression, period=24.0, lag_range=(0, 23))

    assert pval < 0.05, "Should detect significant rhythm"
    assert tau > 0, "Tau should be positive"


def test_cosinor_regression():
    """Test cosinor regression."""
    # Create cosine wave
    time = np.arange(24)
    true_amplitude = 2.0
    true_phase = 6.0
    expression = 10 + true_amplitude * np.cos(2 * np.pi * (time - true_phase) / 24)

    # Fit cosinor
    result = cosinor_regression(expression, time, period=24.0)

    assert 'amplitude' in result, "Missing amplitude"
    assert 'acrophase' in result, "Missing acrophase"
    assert 'pvalue' in result, "Missing pvalue"

    # Check amplitude is close to true value
    assert abs(result['amplitude'] - true_amplitude) < 0.5, "Amplitude estimate error"


def test_cosinor_detection():
    """Test cosinor-based rhythm detection."""
    data = generate_synthetic_data(n_genes=50, n_timepoints=24, n_rhythmic=20, random_state=42)

    results = detect_rhythms_cosinor(data, period=24.0)

    assert len(results) == 50, "Results length mismatch"
    assert 'Amplitude' in results.columns, "Missing Amplitude column"
    assert 'Acrophase' in results.columns, "Missing Acrophase column"
    assert 'Qvalue' in results.columns, "Missing Qvalue column"


def test_no_rhythm_detection():
    """Test that non-rhythmic data is correctly identified."""
    # Generate random noise
    np.random.seed(42)
    data = pd.DataFrame(
        np.random.randn(20, 24),
        index=[f"Gene_{i}" for i in range(20)],
        columns=[f"ZT{i:02d}" for i in range(24)]
    )

    results = detect_rhythms(data, period=24.0)

    # Most genes should not be significant
    rhythmic = results[results['BH.Q'] < 0.05]
    assert len(rhythmic) < 5, "Should not detect many rhythms in noise"


if __name__ == "__main__":
    # Run tests
    print("Running rhythm detection tests...")

    test_synthetic_data_generation()
    print("✓ test_synthetic_data_generation")

    test_jtk_cycle_detection()
    print("✓ test_jtk_cycle_detection")

    test_jtk_test()
    print("✓ test_jtk_test")

    test_cosinor_regression()
    print("✓ test_cosinor_regression")

    test_cosinor_detection()
    print("✓ test_cosinor_detection")

    test_no_rhythm_detection()
    print("✓ test_no_rhythm_detection")

    print("\nAll tests passed! ✓")
