"""Tests for Dorris-Gray module."""

import numpy as np
import pytest

from igcsea.analysis.dorris_gray import calculate_rtlnv, fit_dorris_gray, prepare_alkane_data


def test_calculate_rtlnv_scalar():
    """Test RTlnV calculation with scalar inputs."""
    result = calculate_rtlnv(10.0, 300.0)
    assert np.isscalar(result) or isinstance(result, np.ndarray)
    assert np.asarray(result) > 5000  # Reasonable range


def test_calculate_rtlnv_array():
    """Test RTlnV calculation with array inputs."""
    volumes = np.array([5.0, 10.0, 15.0])
    result = calculate_rtlnv(volumes, 300.0)
    assert isinstance(result, np.ndarray)
    assert len(result) == 3
    assert np.all(result > 0)


def test_calculate_rtlnv_zero_volume():
    """Test that zero or negative volumes return NaN."""
    result = calculate_rtlnv(0.0, 300.0)
    assert np.isnan(result)

    result = calculate_rtlnv(-5.0, 300.0)
    assert np.isnan(result)


def test_prepare_alkane_data(sample_igc_result):
    """Test alkane data preparation."""
    alkane_data = prepare_alkane_data(sample_igc_result)

    assert "Carbon Number" in alkane_data.columns
    assert "RTlnVg" in alkane_data.columns
    assert not alkane_data.empty

    # Check that only alkanes are present
    assert alkane_data["Carbon Number"].notna().all()


def test_prepare_alkane_data_with_filter(sample_igc_result):
    """Test alkane data preparation with user-specified alkane filter."""
    # Get all available alkanes first
    all_alkanes = prepare_alkane_data(sample_igc_result)
    all_alkane_names = set(all_alkanes["Solvent"].unique())

    # Test with subset of alkanes
    selected_alkanes = ["HEPTANE", "OCTANE", "NONANE"]
    # Only use alkanes that are actually in the data
    selected_alkanes = [a for a in selected_alkanes if a in all_alkane_names]

    if selected_alkanes:  # Only run if we have matching alkanes
        filtered_data = prepare_alkane_data(sample_igc_result, alkanes=selected_alkanes)

        # Check that only selected alkanes are present
        assert set(filtered_data["Solvent"].unique()) == set(selected_alkanes)
        assert len(filtered_data) <= len(all_alkanes)

        # Check structure is still correct
        assert "Carbon Number" in filtered_data.columns
        assert "RTlnVg" in filtered_data.columns
        assert filtered_data["Carbon Number"].notna().all()


def test_fit_dorris_gray(sample_igc_result):
    """Test Dorris-Gray linear fitting."""
    alkane_data = prepare_alkane_data(sample_igc_result)
    fits = fit_dorris_gray(alkane_data)

    assert "slope" in fits.columns
    assert "intercept" in fits.columns
    assert "r_squared" in fits.columns
    assert not fits.empty
