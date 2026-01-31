"""Tests for interpolation module."""

import numpy as np
import pandas as pd
import pytest

from igcsea.analysis.interpolation import interpolate_to_standard_coverages


def test_interpolate_to_standard_coverages(sample_igc_result):
    """Test basic interpolation functionality."""
    result = interpolate_to_standard_coverages(sample_igc_result)

    assert isinstance(result, pd.DataFrame)
    assert not result.empty
    assert "Solvent" in result.columns
    assert "Target Fractional Surface Coverage" in result.columns
    assert "Sp. Ret Volume (Com) [ml/g] (interp)" in result.columns


def test_interpolate_has_temperature_data(sample_igc_result):
    """Test that temperature data is included."""
    result = interpolate_to_standard_coverages(sample_igc_result)

    assert "Closest Column Temperature [Kelvin]" in result.columns
    # At least some temperature values should be non-null
    assert result["Closest Column Temperature [Kelvin]"].notna().any()


def test_interpolate_with_custom_coverages(sample_igc_result):
    """Test interpolation with custom target coverages."""
    custom_coverages = np.array([0.01, 0.05, 0.10])
    result = interpolate_to_standard_coverages(sample_igc_result, target_coverages=custom_coverages)

    unique_targets = result["Target Fractional Surface Coverage"].unique()
    assert len(unique_targets) == 3
    assert 0.01 in unique_targets
    assert 0.05 in unique_targets
    assert 0.10 in unique_targets
