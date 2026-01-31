"""Tests for surface energy module."""

import numpy as np
import pytest

from igcsea.analysis.surface_energy import (
    calculate_surface_energy_profile,
    exp_asym,
    fit_exponential_decay,
)
from igcsea.core.models import SurfaceEnergyProfile


def test_exp_asym():
    """Test exponential asymptotic model."""
    x = np.array([0.0, 0.1, 0.2])
    result = exp_asym(x, c=40, a=6, b=10)

    assert isinstance(result, np.ndarray)
    assert len(result) == 3
    assert result[0] > result[1] > result[2]  # Decreasing


def test_fit_exponential_decay():
    """Test exponential decay fitting."""
    coverage = np.array([0.005, 0.01, 0.025, 0.05, 0.10, 0.14])
    energy = np.array([46.5, 44.8, 43.3, 42.3, 40.7, 40.5])

    params = fit_exponential_decay(coverage, energy)

    assert "c" in params
    assert "a" in params
    assert "b" in params
    assert params["c"] > 0  # Asymptote should be positive
    assert params["a"] > 0  # Initial excess should be positive


def test_fit_exponential_decay_insufficient_points():
    """Test that fitting with too few points raises error."""
    coverage = np.array([0.005, 0.01])
    energy = np.array([46.5, 44.8])

    with pytest.raises(ValueError, match="at least 3 data points"):
        fit_exponential_decay(coverage, energy)


def test_calculate_surface_energy_profile(sample_igc_result):
    """Test complete surface energy profile calculation."""
    profile = calculate_surface_energy_profile(sample_igc_result)

    assert isinstance(profile, SurfaceEnergyProfile)
    assert len(profile.coverage) > 0
    assert len(profile.yd) == len(profile.coverage)
    assert len(profile.yab) == len(profile.coverage)
    assert len(profile.yt) == len(profile.coverage)


def test_surface_energy_profile_with_fits(sample_igc_result):
    """Test that fit parameters are included when requested."""
    profile = calculate_surface_energy_profile(sample_igc_result, fit_exponential=True)

    assert profile.fit_params is not None
    assert "yd" in profile.fit_params
    assert "yab" in profile.fit_params
    assert "yt" in profile.fit_params


def test_surface_energy_profile_without_fits(sample_igc_result):
    """Test that fit parameters are None when not requested."""
    profile = calculate_surface_energy_profile(sample_igc_result, fit_exponential=False)

    assert profile.fit_params is None


def test_surface_energy_yt_equals_yd_plus_yab(sample_igc_result):
    """Test that yt = yd + yab."""
    profile = calculate_surface_energy_profile(sample_igc_result)

    np.testing.assert_allclose(profile.yt, profile.yd + profile.yab, rtol=1e-10)
