"""Dorris-Gray alkane analysis for IGC-SEA data.

This module implements the Dorris-Gray method for determining dispersive
surface energy from n-alkane retention volumes.
"""

from typing import Optional

import numpy as np
import pandas as pd

from igcsea.core.constants import ALKANE_CARBON_NUMBERS, R
from igcsea.core.models import IGCResult


def calculate_rtlnv(
    volume: np.ndarray | float,
    temperature: np.ndarray | float,
    gas_constant: float = R,
) -> np.ndarray | float:
    """Calculate RT·ln(V) for dispersive energy analysis.

    The Dorris-Gray parameter RT·ln(V) is used in determining dispersive
    surface energy from n-alkane retention volumes. This quantity is
    linearly related to the carbon number for alkanes.

    Args:
        volume: Specific retention volume in ml/g (scalar or array).
        temperature: Column temperature in Kelvin (scalar or array).
        gas_constant: Universal gas constant in J/(mol·K). Defaults to R.

    Returns:
        RT·ln(V) in J/mol. Returns NaN for non-positive volumes.

    Examples:
        >>> calculate_rtlnv(10.0, 300.0)
        5765.18...

        >>> volumes = np.array([5, 10, 15])
        >>> calculate_rtlnv(volumes, 300.0)
        array([4028.85..., 5765.18..., 6794.43...])

    References:
        Dorris, G. M., & Gray, D. G. (1980). Adsorption of n-alkanes at zero
        surface coverage on cellulose paper and wood fibers. Journal of Colloid
        and Interface Science, 77(2), 353-362.
    """
    V = np.asarray(volume, dtype=float)
    T = np.asarray(temperature, dtype=float)
    return np.where(V > 0, gas_constant * T * np.log(V), np.nan)


_RETENTION_VOL_COLS = {
    "COM": "Interpolated Retention Volume (Com)",
    "MAX": "Interpolated Retention Volume (Max)",
}


def prepare_alkane_data(
    igc_result: IGCResult,
    target_coverages: Optional[np.ndarray] = None,
    alkanes: Optional[list[str]] = None,
    retention_type: str = "COM",
) -> pd.DataFrame:
    """Prepare alkane probe data for Dorris-Gray analysis.

    This function extracts alkane data from the free_energy table, which contains
    the instrument's pre-calculated interpolated retention volumes. This ensures
    accurate RT·ln(V) calculations that match the instrument's output.

    Args:
        igc_result: Parsed IGC-SEA result.
        target_coverages: Unused (kept for API compatibility). The free_energy
            table already has standard coverage values.
        alkanes: Optional list of alkane names to include in analysis. If None,
            uses all alkanes found in the data. Use this to exclude alkanes with
            poor retention or to specify a custom range. Example:
            ["HEPTANE", "OCTANE", "NONANE", "DECANE"]
        retention_type: ``'COM'`` (centre of mass) or ``'MAX'`` (peak maximum).
            Selects which interpolated retention volume column to use from the
            free_energy table. Defaults to ``'COM'``.

    Returns:
        DataFrame with columns:
        - Solvent
        - Target Fractional Surface Coverage
        - Interpolated Retention Volume (Com/Max)
        - Column Temperature [Kelvin]
        - Carbon Number
        - RTlnVg (RT·ln(V) in J/mol)

    Raises:
        ValueError: If retention_type is not ``'COM'`` or ``'MAX'``.

    Examples:
        >>> result = parse_igc_csv("data.csv")
        >>> # Use all available alkanes (default)
        >>> alkane_data = prepare_alkane_data(result)
        >>>
        >>> # Use only specific alkanes (e.g., exclude hexane)
        >>> alkane_data = prepare_alkane_data(
        ...     result,
        ...     alkanes=["HEPTANE", "OCTANE", "NONANE", "DECANE"]
        ... )
    """
    retention_type = retention_type.upper()
    if retention_type not in _RETENTION_VOL_COLS:
        raise ValueError(f"retention_type must be 'COM' or 'MAX', got {retention_type!r}")

    if igc_result.free_energy is None:
        raise ValueError(
            "prepare_alkane_data() requires igc_result.free_energy to be present. "
            "For injection-items-only exports use prepare_alkane_data_from_injections() "
            "from igcsea.analysis.manual_dispersive instead."
        )

    vol_col = _RETENTION_VOL_COLS[retention_type]

    # Extract alkane data from free_energy table
    free = igc_result.free_energy.copy()

    # Rename columns to match expected output format
    alkane = free.rename(columns={
        "Solvent Name": "Solvent",
        "n/nm": "Target Fractional Surface Coverage",
    })

    # Filter for alkanes only
    alkane = alkane[alkane["Solvent"].map(ALKANE_CARBON_NUMBERS).notna()].copy()

    # Apply user-specified alkane filter if provided
    if alkanes is not None:
        alkane = alkane[alkane["Solvent"].isin(alkanes)].copy()

    alkane["Carbon Number"] = alkane["Solvent"].map(ALKANE_CARBON_NUMBERS)

    # Calculate RT·ln(V) using the pre-calculated interpolated retention volumes
    alkane["RTlnVg"] = calculate_rtlnv(
        alkane[vol_col],
        alkane["Column Temperature [Kelvin]"],
    )

    # Select and reorder columns for output
    alkane = alkane[[
        "Solvent",
        "Target Fractional Surface Coverage",
        vol_col,
        "Column Temperature [Kelvin]",
        "Carbon Number",
        "RTlnVg",
    ]].reset_index(drop=True)

    return alkane


def fit_dorris_gray(alkane_data: pd.DataFrame) -> pd.DataFrame:
    """Fit linear regression for RT·ln(V) vs carbon number at each coverage.

    For each coverage level, fits the model:
        RT·ln(V) = slope * carbon_number + intercept

    The slope is related to the dispersive surface energy per CH2 group.

    Args:
        alkane_data: DataFrame from prepare_alkane_data() with alkane data.

    Returns:
        DataFrame with columns:
        - Target Fractional Surface Coverage
        - slope: Linear fit slope (J/mol per carbon)
        - intercept: Linear fit intercept (J/mol)
        - r_squared: Goodness of fit (R²)

    Examples:
        >>> alkane_data = prepare_alkane_data(result)
        >>> fits = fit_dorris_gray(alkane_data)
        >>> print(fits)
    """
    results = []

    for coverage in alkane_data["Target Fractional Surface Coverage"].unique():
        subset = alkane_data[
            alkane_data["Target Fractional Surface Coverage"] == coverage
        ]

        # Drop rows with NaN in RTlnVg
        subset = subset.dropna(subset=["Carbon Number", "RTlnVg"])

        if len(subset) < 2:
            # Need at least 2 points for linear fit
            results.append({
                "Target Fractional Surface Coverage": coverage,
                "slope": np.nan,
                "intercept": np.nan,
                "r_squared": np.nan,
            })
            continue

        x = subset["Carbon Number"].to_numpy(float)
        y = subset["RTlnVg"].to_numpy(float)

        # Linear fit: y = mx + b
        slope, intercept = np.polyfit(x, y, 1)

        # Calculate R²
        y_pred = slope * x + intercept
        ss_res = np.sum((y - y_pred) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot != 0 else np.nan

        results.append({
            "Target Fractional Surface Coverage": coverage,
            "slope": slope,
            "intercept": intercept,
            "r_squared": r_squared,
        })

    return pd.DataFrame(results)
