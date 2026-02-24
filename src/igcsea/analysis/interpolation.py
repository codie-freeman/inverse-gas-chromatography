"""Coverage interpolation and extrapolation utilities.

This module provides functions to interpolate/extrapolate retention volumes
from actual measured coverages to standard target coverages for analysis.
"""

from typing import Optional

import numpy as np
import pandas as pd

from igcsea.core.constants import STANDARD_COVERAGES
from igcsea.core.models import IGCResult


def interpolate_to_standard_coverages(
    igc_result: IGCResult,
    target_coverages: Optional[np.ndarray] = None,
    coverage_col: str = "Actual Fractional Surface Coverage",
    value_col: str = "Sp. Ret Volume (Com) [ml/g]",
) -> pd.DataFrame:
    """Interpolate retention volumes to standard coverage values.

    IGC experiments measure retention at various actual coverages, which don't
    necessarily match the desired analysis points. This function interpolates
    (or linearly extrapolates) retention volumes to standard target coverages
    for consistent analysis across solvents.

    Args:
        igc_result: Parsed IGC-SEA result containing injection data.
        target_coverages: Array of target coverage values (n/nm). If None,
            uses STANDARD_COVERAGES from constants.
        coverage_col: Column name for actual coverage values.
        value_col: Column name for retention volume values to interpolate.

    Returns:
        DataFrame with columns:
        - Solvent: Solvent name
        - Target Fractional Surface Coverage: Target coverage value
        - Sp. Ret Volume (Com) [ml/g] (interp): Interpolated retention volume
        - Closest Actual Fractional Surface Coverage: Closest measured coverage
        - Closest Sp. Ret Volume (Com) [ml/g]: Retention volume at closest coverage
        - Closest Column Temperature [Kelvin]: Temperature at closest coverage

    Examples:
        >>> result = parse_igc_csv("data.csv")
        >>> interpolated = interpolate_to_standard_coverages(result)
        >>> print(interpolated[interpolated["Solvent"] == "HEXANE"])

    References:
        Linear extrapolation beyond data range uses endpoint slopes to maintain
        physically reasonable behavior.
    """
    if target_coverages is None:
        target_coverages = np.array(STANDARD_COVERAGES, dtype=float)

    # Extract injection items and filter out methane
    inj = igc_result.injection_items
    data = inj.loc[inj["Solvent"].ne("Methane"), ["Solvent", coverage_col, value_col]].copy()

    # Ensure numeric
    data[coverage_col] = pd.to_numeric(data[coverage_col], errors="coerce")
    data[value_col] = pd.to_numeric(data[value_col], errors="coerce")

    # Drop rows with missing values
    data = data.dropna(subset=["Solvent", coverage_col, value_col])

    # Build temperature lookup table
    temp_lookup = _build_temperature_lookup(inj, coverage_col)

    # Apply interpolation per solvent (explicit loop avoids pandas 3.0
    # include_groups=False default, which drops the groupby key column)
    groups = [
        _interpolate_single_solvent(
            group_df,
            solvent_name=solvent_name,
            target_coverages=target_coverages,
            coverage_col=coverage_col,
            value_col=value_col,
            temp_lookup=temp_lookup,
        )
        for solvent_name, group_df in data.groupby("Solvent")
    ]
    combined = pd.concat(groups, ignore_index=True)

    combined = combined.sort_values(
        ["Target Fractional Surface Coverage", "Solvent"]
    ).reset_index(drop=True)

    return combined


def _build_temperature_lookup(
    injection_items: pd.DataFrame,
    coverage_col: str,
) -> pd.DataFrame:
    """Build lookup table mapping (Solvent, Coverage) -> Temperature.

    Args:
        injection_items: DataFrame with injection data.
        coverage_col: Column name for actual coverage values.

    Returns:
        DataFrame with Solvent, coverage, and temperature columns.
    """
    temp_lookup = (
        injection_items.loc[
            injection_items["Solvent"].ne("Methane"),
            ["Solvent", coverage_col, "Column Temperature [Kelvin]"],
        ]
        .assign(
            **{
                coverage_col: pd.to_numeric(
                    injection_items[coverage_col], errors="coerce"
                ),
                "Column Temperature [Kelvin]": pd.to_numeric(
                    injection_items["Column Temperature [Kelvin]"], errors="coerce"
                ),
            }
        )
        .dropna(subset=["Solvent", coverage_col, "Column Temperature [Kelvin]"])
        # If duplicates exist at same solvent+coverage, average them
        .groupby(["Solvent", coverage_col], as_index=False)["Column Temperature [Kelvin]"]
        .mean()
    )
    return temp_lookup


def _interpolate_single_solvent(
    df: pd.DataFrame,
    solvent_name: str,
    target_coverages: np.ndarray,
    coverage_col: str,
    value_col: str,
    temp_lookup: pd.DataFrame,
) -> pd.DataFrame:
    """Interpolate/extrapolate retention volumes for a single solvent.

    This function handles three cases:
    1. Linear interpolation within the measured range
    2. Linear extrapolation below the minimum using left endpoint slope
    3. Linear extrapolation above the maximum using right endpoint slope

    Args:
        df: DataFrame for single solvent with actual coverage and retention data.
        target_coverages: Array of target coverage values.
        coverage_col: Column name for actual coverage values.
        value_col: Column name for retention volume values.
        temp_lookup: Temperature lookup DataFrame.

    Returns:
        DataFrame with interpolated values and closest actual measurements.
    """
    solvent = solvent_name
    df = df.sort_values(coverage_col).copy()

    x = df[coverage_col].to_numpy(float)
    y = df[value_col].to_numpy(float)

    # Filter out non-finite values
    m = np.isfinite(x) & np.isfinite(y)
    xk, yk = x[m], y[m]

    out = pd.DataFrame({
        "Solvent": solvent,
        "Target Fractional Surface Coverage": target_coverages,
    })

    interp_col = f"{value_col} (interp)"
    closest_val_col = f"Closest {value_col}"

    # Need at least 2 points for interpolation
    if len(xk) < 2:
        out["Closest Actual Fractional Surface Coverage"] = np.nan
        out[closest_val_col] = np.nan
        out[interp_col] = np.nan
        out["Closest Column Temperature [Kelvin]"] = np.nan
        return out

    # Interpolation (numpy.interp clamps at boundaries by default)
    y_interp = np.interp(target_coverages, xk, yk)

    # Linear extrapolation on the left (below minimum coverage)
    left = target_coverages < xk.min()
    if left.any():
        x0, x1 = xk[0], xk[1]
        y0, y1 = yk[0], yk[1]
        slope = (y1 - y0) / (x1 - x0)
        y_interp[left] = y0 + slope * (target_coverages[left] - x0)

    # Linear extrapolation on the right (above maximum coverage)
    right = target_coverages > xk.max()
    if right.any():
        x0, x1 = xk[-2], xk[-1]
        y0, y1 = yk[-2], yk[-1]
        slope = (y1 - y0) / (x1 - x0)
        y_interp[right] = y1 + slope * (target_coverages[right] - x1)

    out[interp_col] = y_interp

    # Find closest actual coverage for each target (for temperature lookup)
    closest_actual = np.empty_like(target_coverages, dtype=float)
    closest_spret = np.empty_like(target_coverages, dtype=float)

    xmin, xmax = float(xk.min()), float(xk.max())

    for i, target in enumerate(target_coverages):
        if target <= xmin:
            j = 0  # Lowest actual obtained
        elif target >= xmax:
            j = len(xk) - 1  # Highest actual obtained
        else:
            # Find nearest actual coverage (ties broken toward smaller)
            diffs = np.abs(xk - target)
            j_candidates = np.where(diffs == diffs.min())[0]
            j = int(j_candidates[0])

        closest_actual[i] = xk[j]
        closest_spret[i] = yk[j]

    out["Closest Actual Fractional Surface Coverage"] = closest_actual
    out[closest_val_col] = closest_spret

    # Merge temperature data from lookup table
    tmp = pd.DataFrame({
        "Solvent": solvent,
        "Closest Actual Fractional Surface Coverage": closest_actual,
    }).merge(
        temp_lookup,
        left_on=["Solvent", "Closest Actual Fractional Surface Coverage"],
        right_on=["Solvent", coverage_col],
        how="left",
    )

    out["Closest Column Temperature [Kelvin]"] = tmp["Column Temperature [Kelvin]"].to_numpy()

    return out
