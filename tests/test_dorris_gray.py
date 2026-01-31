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


def test_fit_dorris_gray(sample_igc_result):
    """Test Dorris-Gray linear fitting."""
    alkane_data = prepare_alkane_data(sample_igc_result)
    fits = fit_dorris_gray(alkane_data)

    assert "slope" in fits.columns
    assert "intercept" in fits.columns
    assert "r_squared" in fits.columns
    assert not fits.empty
