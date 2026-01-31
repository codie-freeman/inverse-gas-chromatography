"""Tests for IGC-SEA CSV parsing."""

from pathlib import Path

import pandas as pd
import pytest

from igcsea.core.models import IGCResult
from igcsea.parsing import parse_igc_csv


def test_parse_igc_csv_returns_igc_result(sample_igc_result):
    """Test that parse_igc_csv returns an IGCResult instance."""
    assert isinstance(sample_igc_result, IGCResult)


def test_parse_igc_csv_has_all_tables(sample_igc_result):
    """Test that all three required tables are present."""
    assert isinstance(sample_igc_result.free_energy, pd.DataFrame)
    assert isinstance(sample_igc_result.dispersive_surface_energy, pd.DataFrame)
    assert isinstance(sample_igc_result.injection_items, pd.DataFrame)


def test_parse_igc_csv_tables_not_empty(sample_igc_result):
    """Test that parsed tables contain data."""
    assert not sample_igc_result.free_energy.empty
    assert not sample_igc_result.dispersive_surface_energy.empty
    assert not sample_igc_result.injection_items.empty


def test_parse_igc_csv_source_path(sample_csv_path, sample_igc_result):
    """Test that source_path is correctly stored."""
    assert sample_igc_result.source_path == sample_csv_path


def test_free_energy_has_expected_columns(sample_igc_result):
    """Test that free_energy DataFrame has expected columns."""
    expected_columns = ["n/nm", "Solvent Name"]
    for col in expected_columns:
        assert col in sample_igc_result.free_energy.columns


def test_dispersive_surface_energy_has_coverage(sample_igc_result):
    """Test that dispersive_surface_energy has coverage column."""
    assert "n/nm" in sample_igc_result.dispersive_surface_energy.columns


def test_injection_items_has_expected_columns(sample_igc_result):
    """Test that injection_items has expected columns."""
    expected_columns = ["Solvent", "Actual Fractional Surface Coverage"]
    for col in expected_columns:
        assert col in sample_igc_result.injection_items.columns


def test_parse_igc_csv_numeric_conversion(sample_igc_result):
    """Test that numeric columns are properly converted."""
    # Coverage should be numeric
    assert pd.api.types.is_numeric_dtype(sample_igc_result.free_energy["n/nm"])

    # Solvent Name should remain as object/string
    assert pd.api.types.is_object_dtype(sample_igc_result.free_energy["Solvent Name"])


def test_parse_igc_csv_file_not_found():
    """Test that FileNotFoundError is raised for missing file."""
    with pytest.raises(FileNotFoundError):
        parse_igc_csv("nonexistent_file.csv")


def test_igc_result_as_dict(sample_igc_result):
    """Test IGCResult.as_dict() method."""
    result_dict = sample_igc_result.as_dict()

    assert isinstance(result_dict, dict)
    assert "free_energy" in result_dict
    assert "dispersive_surface_energy" in result_dict
    assert "injection_items" in result_dict


def test_igc_result_to_csv_dir(sample_igc_result, tmp_path):
    """Test IGCResult.to_csv_dir() method."""
    output_dir = tmp_path / "output"

    sample_igc_result.to_csv_dir(output_dir)

    # Check that files were created
    assert (output_dir / "free_energy.csv").exists()
    assert (output_dir / "dispersive_surface_energy.csv").exists()
    assert (output_dir / "injection_items.csv").exists()

    # Verify files can be read back
    df = pd.read_csv(output_dir / "free_energy.csv")
    assert not df.empty


def test_parse_igc_csv_with_string_path(sample_csv_path):
    """Test that parse_igc_csv works with string paths."""
    result = parse_igc_csv(str(sample_csv_path))
    assert isinstance(result, IGCResult)


def test_parse_igc_csv_with_path_object(sample_csv_path):
    """Test that parse_igc_csv works with Path objects."""
    result = parse_igc_csv(sample_csv_path)
    assert isinstance(result, IGCResult)
