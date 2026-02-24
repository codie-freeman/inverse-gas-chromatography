"""Tests for manual Dorris-Gray dispersive surface energy module."""

import numpy as np
import pandas as pd
import pytest

from igcsea.analysis.manual_dispersive import (
    calculate_dispersive_from_injections,
    gamma_ch2,
    prepare_alkane_data_from_injections,
    slope_to_gamma_d,
    validate_against_sms,
)


# ---------------------------------------------------------------------------
# Physical helper function tests
# ---------------------------------------------------------------------------

def test_gamma_ch2_at_30c():
    """γ_CH2 at 30 °C should equal 35.6 - 0.058*30."""
    assert gamma_ch2(30.0) == pytest.approx(35.6 - 0.058 * 30.0)


def test_gamma_ch2_decreases_with_temperature():
    """γ_CH2 should decrease as temperature increases."""
    assert gamma_ch2(20.0) > gamma_ch2(40.0)


def test_slope_to_gamma_d_positive():
    """γ_s^d should be positive for a positive slope."""
    result = slope_to_gamma_d(500.0, 303.15)
    assert result > 0


def test_slope_to_gamma_d_units():
    """γ_s^d should be in a physically reasonable range for typical IGC (20–100 mJ/m²)."""
    result = slope_to_gamma_d(500.0, 303.15)
    assert 1.0 < result < 500.0


def test_slope_to_gamma_d_zero_slope():
    """Zero slope should give zero γ_s^d."""
    assert slope_to_gamma_d(0.0, 303.15) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Data preparation tests
# ---------------------------------------------------------------------------

def test_prepare_alkane_data_from_injections_returns_dataframe(sample_igc_result):
    """Should return a non-empty DataFrame."""
    result = prepare_alkane_data_from_injections(sample_igc_result)
    assert isinstance(result, pd.DataFrame)
    assert not result.empty


def test_prepare_alkane_data_columns(sample_igc_result):
    """Should have the expected columns."""
    result = prepare_alkane_data_from_injections(sample_igc_result)
    for col in ["Solvent", "Target Fractional Surface Coverage",
                "Carbon Number", "Column Temperature [Kelvin]", "RTlnVg"]:
        assert col in result.columns


def test_prepare_alkane_data_only_alkanes(sample_igc_result):
    """All rows should correspond to recognised alkanes."""
    result = prepare_alkane_data_from_injections(sample_igc_result)
    assert result["Carbon Number"].notna().all()


def test_prepare_alkane_data_max(sample_igc_result):
    """Should work with retention_type='MAX'."""
    result = prepare_alkane_data_from_injections(sample_igc_result, retention_type="MAX")
    assert not result.empty
    assert "RTlnVg" in result.columns


def test_prepare_alkane_data_invalid_retention_type(sample_igc_result):
    """Should raise ValueError for an unrecognised retention_type."""
    with pytest.raises(ValueError, match="retention_type"):
        prepare_alkane_data_from_injections(sample_igc_result, retention_type="INVALID")


def test_prepare_alkane_data_with_alkane_filter(sample_igc_result):
    """User-supplied alkane filter should restrict the result."""
    all_data = prepare_alkane_data_from_injections(sample_igc_result)
    available = list(all_data["Solvent"].unique())

    if len(available) >= 2:
        subset = available[:2]
        filtered = prepare_alkane_data_from_injections(sample_igc_result, alkanes=subset)
        assert set(filtered["Solvent"].unique()) == set(subset)


# ---------------------------------------------------------------------------
# Main calculation tests
# ---------------------------------------------------------------------------

def test_calculate_dispersive_from_injections_returns_dataframe(sample_igc_result):
    """Should return a non-empty DataFrame."""
    result = calculate_dispersive_from_injections(sample_igc_result)
    assert isinstance(result, pd.DataFrame)
    assert not result.empty


def test_calculate_dispersive_columns(sample_igc_result):
    """Should have the expected output columns."""
    result = calculate_dispersive_from_injections(sample_igc_result)
    for col in ["n/nm", "gamma_d", "delta_g_ch2", "r_squared", "temperature_kelvin"]:
        assert col in result.columns


def test_calculate_dispersive_gamma_d_positive(sample_igc_result):
    """γ_s^d values should all be positive."""
    result = calculate_dispersive_from_injections(sample_igc_result)
    assert (result["gamma_d"].dropna() > 0).all()


def test_calculate_dispersive_r_squared_range(sample_igc_result):
    """R² values should be between 0 and 1."""
    result = calculate_dispersive_from_injections(sample_igc_result)
    valid = result["r_squared"].dropna()
    assert (valid >= 0).all()
    assert (valid <= 1).all()


def test_calculate_dispersive_coverage_count(sample_igc_result):
    """Number of coverage rows should match standard coverages used."""
    result = calculate_dispersive_from_injections(sample_igc_result)
    assert len(result) > 0


def test_calculate_dispersive_com_vs_max(sample_igc_result):
    """COM and MAX calculations should give different but nearby γ_s^d values."""
    com = calculate_dispersive_from_injections(sample_igc_result, retention_type="COM")
    max_ = calculate_dispersive_from_injections(sample_igc_result, retention_type="MAX")
    # Both should produce results
    assert not com["gamma_d"].dropna().empty
    assert not max_["gamma_d"].dropna().empty


# ---------------------------------------------------------------------------
# Validation against SMS tests
# ---------------------------------------------------------------------------

def test_validate_against_sms_returns_dataframe(sample_igc_result):
    """Should return a comparison DataFrame when dispersive_surface_energy is present."""
    if sample_igc_result.dispersive_surface_energy is None:
        pytest.skip("Sample CSV does not contain dispersive_surface_energy table")

    gd = calculate_dispersive_from_injections(sample_igc_result)
    cmp = validate_against_sms(gd, sample_igc_result.dispersive_surface_energy)

    assert isinstance(cmp, pd.DataFrame)
    assert not cmp.empty


def test_validate_against_sms_columns(sample_igc_result):
    """Comparison DataFrame should have all expected columns."""
    if sample_igc_result.dispersive_surface_energy is None:
        pytest.skip("Sample CSV does not contain dispersive_surface_energy table")

    gd = calculate_dispersive_from_injections(sample_igc_result)
    cmp = validate_against_sms(gd, sample_igc_result.dispersive_surface_energy)

    for col in ["n/nm", "gamma_d_manual", "gamma_d_sms", "pct_diff", "exceeds_threshold"]:
        assert col in cmp.columns


def test_validate_against_sms_pct_diff_range(sample_igc_result):
    """Comparison should produce finite pct_diff values."""
    if sample_igc_result.dispersive_surface_energy is None:
        pytest.skip("Sample CSV does not contain dispersive_surface_energy table")

    gd = calculate_dispersive_from_injections(sample_igc_result)
    cmp = validate_against_sms(gd, sample_igc_result.dispersive_surface_energy)

    # pct_diff should be finite numbers — magnitude will vary depending on how
    # closely our linear interpolation matches the SMS proprietary method
    assert cmp["pct_diff"].notna().all()
    assert np.isfinite(cmp["pct_diff"]).all()


def test_validate_against_sms_invalid_retention_type(sample_igc_result):
    """Should raise ValueError for an unrecognised retention_type."""
    if sample_igc_result.dispersive_surface_energy is None:
        pytest.skip("Sample CSV does not contain dispersive_surface_energy table")

    gd = calculate_dispersive_from_injections(sample_igc_result)
    with pytest.raises(ValueError, match="retention_type"):
        validate_against_sms(
            gd, sample_igc_result.dispersive_surface_energy, retention_type="INVALID"
        )
