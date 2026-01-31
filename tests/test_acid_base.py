"""Tests for acid-base module."""

import pandas as pd

from igcsea.analysis.acid_base import (
    calculate_acid_base_params,
    calculate_yab_della_volpe,
    extract_probe_data,
)
from igcsea.core.models import AcidBaseParams


def test_extract_probe_data(sample_igc_result):
    """Test probe data extraction."""
    ea = extract_probe_data(sample_igc_result.free_energy, "ETHYL ACETATE")

    assert isinstance(ea, pd.DataFrame)
    assert "n/nm" in ea.columns
    assert "En. (Pol Com)" in ea.columns
    assert not ea.empty


def test_calculate_yab_della_volpe(sample_igc_result):
    """Test YAB calculation using Della Volpe method."""
    ea = extract_probe_data(sample_igc_result.free_energy, "ETHYL ACETATE")
    dcm = extract_probe_data(sample_igc_result.free_energy, "DICHLOROMETHANE")

    yab_df = calculate_yab_della_volpe(ea, dcm)

    assert isinstance(yab_df, pd.DataFrame)
    assert "n/nm" in yab_df.columns
    assert "yab" in yab_df.columns
    assert "ys+" in yab_df.columns
    assert "ys-" in yab_df.columns
    assert not yab_df.empty


def test_yab_values_positive(sample_igc_result):
    """Test that YAB values are positive."""
    ea = extract_probe_data(sample_igc_result.free_energy, "ETHYL ACETATE")
    dcm = extract_probe_data(sample_igc_result.free_energy, "DICHLOROMETHANE")
    yab_df = calculate_yab_della_volpe(ea, dcm)

    assert (yab_df["yab"] > 0).all()
    assert (yab_df["ys+"] > 0).all()
    assert (yab_df["ys-"] > 0).all()


def test_calculate_acid_base_params(sample_igc_result):
    """Test complete acid-base parameter calculation."""
    params = calculate_acid_base_params(sample_igc_result)

    assert isinstance(params, AcidBaseParams)
    assert len(params.coverage) > 0
    assert len(params.yab) == len(params.coverage)
    assert len(params.ys_plus) == len(params.coverage)
    assert len(params.ys_minus) == len(params.coverage)
