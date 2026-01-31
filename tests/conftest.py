"""Pytest fixtures for IGC-SEA tests."""

from pathlib import Path

import pytest

from igcsea.parsing import parse_igc_csv


@pytest.fixture
def data_dir() -> Path:
    """Return path to test data directory."""
    return Path(__file__).parent.parent / "data" / "examples"


@pytest.fixture
def sample_csv_path(data_dir) -> Path:
    """Return path to sample IGC-SEA CSV file."""
    return data_dir / "sample_igc_export.csv"


@pytest.fixture
def sample_igc_result(sample_csv_path):
    """Return parsed IGCResult from sample CSV."""
    return parse_igc_csv(sample_csv_path)
