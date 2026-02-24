"""CSV parsing utilities for IGC-SEA data exports.

This module provides functions to parse multi-table CSV exports from IGC-SEA
instruments (e.g., SMS Cirrus Plus) into structured DataFrames.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import pandas as pd

from igcsea.core.models import IGCResult


def _read_rows(path: Path) -> pd.DataFrame:
    """Read raw CSV rows with Windows-1252 encoding.

    Args:
        path: Path to CSV file.

    Returns:
        DataFrame with all columns as strings, empty strings preserved.
    """
    return pd.read_csv(path, encoding="cp1252", header=None, dtype=str, keep_default_na=False)


def _is_blank_row(row: pd.Series) -> bool:
    """Check if a row is entirely blank.

    Args:
        row: pandas Series representing a row.

    Returns:
        True if all values are empty strings after stripping whitespace.
    """
    return all((str(v).strip() == "" for v in row.tolist()))


def _maybe_to_numeric(s: pd.Series, *, threshold: float = 0.8) -> pd.Series:
    """Convert Series to numeric if most values are numeric.

    This preserves text columns (like 'Solvent Name') while converting
    numeric columns to proper numeric dtype.

    Args:
        s: pandas Series to potentially convert.
        threshold: Fraction of non-empty values that must be numeric to trigger
            conversion. Defaults to 0.8 (80%).

    Returns:
        Series converted to numeric if threshold met, otherwise original with
        empty strings replaced by pd.NA.
    """
    s2 = s.replace("", pd.NA)
    non_empty = s2.dropna()
    if non_empty.empty:
        return s2

    numeric = pd.to_numeric(non_empty, errors="coerce")
    ratio = float(numeric.notna().mean())
    if ratio >= threshold:
        out = pd.to_numeric(s2, errors="coerce")
        return out
    return s2


def extract_section_table(
    rows: pd.DataFrame,
    section_title: str,
    *,
    stop_titles: Optional[List[str]] = None,
    blank_run_to_stop: int = 2,
) -> pd.DataFrame:
    """Extract a single section table from multi-table CSV.

    IGC-SEA CSV exports contain multiple sections separated by title rows.
    This function extracts data for a single section.

    Args:
        rows: DataFrame containing all raw CSV rows.
        section_title: Title string marking the start of the section
            (e.g., "Free Energy", "Injection Items").
        stop_titles: List of section titles that mark the end of this section.
            If None, uses blank rows to detect end.
        blank_run_to_stop: Number of consecutive blank rows required to stop
            extraction. Defaults to 2.

    Returns:
        DataFrame with section data, header row as column names, numeric
        columns converted, blank rows/columns removed.

    Raises:
        ValueError: If section_title not found or no header row found.

    Examples:
        >>> rows = _read_rows("data.csv")
        >>> free_energy = extract_section_table(
        ...     rows, "Free Energy",
        ...     stop_titles=["Dispersive Surface Energy"]
        ... )
    """
    stop_titles = stop_titles or []
    col0 = rows[0].astype(str).str.strip()

    # Locate section title row
    matches = col0[col0 == section_title]
    if matches.empty:
        raise ValueError(f"Section title not found: {section_title!r}")
    start_idx = int(matches.index[0])

    # Header row = next non-blank row after title
    header_idx = None
    for i in range(start_idx + 1, len(rows)):
        if not _is_blank_row(rows.loc[i]):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"No header row found after section {section_title!r}")

    header = [str(x).strip() for x in rows.loc[header_idx].tolist()]
    while header and header[-1] == "":
        header.pop()

    data_start = header_idx + 1

    # Find end: stop at next title OR a run of blank rows
    end_idx = len(rows)
    blank_run = 0
    for i in range(data_start, len(rows)):
        if str(rows.loc[i, 0]).strip() in stop_titles:
            end_idx = i
            break

        if _is_blank_row(rows.loc[i]):
            blank_run += 1
            if blank_run >= blank_run_to_stop:
                end_idx = i - blank_run + 1  # End at first blank in the run
                break
        else:
            blank_run = 0

    data = rows.iloc[data_start:end_idx, : len(header)].copy()
    data.columns = header

    # Drop fully empty rows/cols
    data = data.replace("", pd.NA)
    data = data.dropna(axis=0, how="all").dropna(axis=1, how="all")

    # Numeric conversion
    for c in data.columns:
        data[c] = _maybe_to_numeric(data[c])

    return data


def parse_igc_csv(path: Path | str) -> IGCResult:
    """Parse an IGC-SEA CSV export into structured DataFrames.

    This function reads a multi-table CSV export from an IGC-SEA instrument
    (e.g., SMS Cirrus Plus) and extracts three main sections:
    - Free Energy: Free energy by coverage and solvent
    - Dispersive Surface Energy: Dispersive component values
    - Injection Items: Raw injection data and adsorption measurements

    Args:
        path: Path to the IGC-SEA CSV export file.

    Returns:
        IGCResult dataclass containing three DataFrames and source path.

    Raises:
        ValueError: If required sections are not found in the CSV.
        FileNotFoundError: If the file does not exist.

    Examples:
        >>> result = parse_igc_csv("data/sample.csv")
        >>> print(result.free_energy.head())
        >>> print(f"Found {len(result.injection_items)} injections")

    References:
        SMS Surface Measurement Systems - Cirrus Plus IGC-SEA
        https://www.surfacemeasurementsystems.com/
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"CSV file not found: {path}")

    rows = _read_rows(path)

    try:
        free_energy = extract_section_table(
            rows,
            "Free Energy",
            stop_titles=["Dispersive Surface Energy", "Injection Items"],
        )
    except ValueError:
        free_energy = None

    try:
        dispersive_surface_energy = extract_section_table(
            rows,
            "Dispersive Surface Energy",
            stop_titles=["Injection Items"],
        )
    except ValueError:
        dispersive_surface_energy = None

    injection_items = extract_section_table(
        rows,
        "Injection Items",
        stop_titles=[],
    )

    return IGCResult(
        injection_items=injection_items,
        source_path=path,
        free_energy=free_energy,
        dispersive_surface_energy=dispersive_surface_energy,
    )
