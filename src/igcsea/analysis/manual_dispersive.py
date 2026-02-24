"""Manual Dorris-Gray dispersive surface energy from raw injection data.

This module calculates dispersive surface energy directly from injection_items
without relying on the SMS Cirrus Plus pre-calculated Free Energy or
Dispersive Surface Energy tables.  It works with both full CSV exports and
injection-items-only exports, enabling independent reproducible analysis and
optional validation against the instrument's output.

Workflow:
    1. Interpolate raw per-injection retention volumes to standard coverages.
    2. Compute RT·ln(V) for each alkane at each coverage.
    3. Fit RT·ln(V) vs carbon number (Dorris-Gray linear regression).
    4. Convert the ΔG_CH2 slope to γ_s^d using the temperature-dependent
       Dorris-Gray formula.
    5. Optionally compare results to the SMS DnG columns.

References:
    Dorris, G. M., & Gray, D. G. (1980). Adsorption of n-alkanes at zero
    surface coverage on cellulose paper and wood fibers. Journal of Colloid
    and Interface Science, 77(2), 353-362.
"""

from typing import List, Optional

import numpy as np
import pandas as pd

from igcsea.analysis.dorris_gray import calculate_rtlnv, fit_dorris_gray
from igcsea.analysis.interpolation import interpolate_to_standard_coverages
from igcsea.core.constants import A_CH2, ALKANE_CARBON_NUMBERS, NA
from igcsea.core.models import IGCResult

# Injection items column names keyed by retention type
_RETENTION_COLS = {
    "COM": "Sp. Ret Volume (Com) [ml/g]",
    "MAX": "Sp. Ret Volume (Max) [ml/g]",
}

# SMS dispersive surface energy column names keyed by retention type
_SMS_DNG_COLS = {
    "COM": "Disp. Surf. En. (mJ/m^2) - DnG & Com",
    "MAX": "Disp. Surf. En. (mJ/m^2) - DnG & Max",
}

_DEVIATION_THRESHOLD = 5.0  # Default % threshold for flagging deviations


# ---------------------------------------------------------------------------
# Physical helper functions
# ---------------------------------------------------------------------------

def gamma_ch2(temperature_celsius: float) -> float:
    """Temperature-dependent dispersive surface energy of a CH2 group (mJ/m²).

    Linear approximation:

        γ_CH2 = 35.6 − 0.058·t  [mJ/m²]

    Args:
        temperature_celsius: Column temperature in °C.

    Returns:
        γ_CH2 in mJ/m².

    Examples:
        >>> gamma_ch2(30.0)
        33.86
    """
    return 35.6 - 0.058 * temperature_celsius


def slope_to_gamma_d(delta_g_ch2: float, temperature_kelvin: float) -> float:
    """Convert a Dorris-Gray slope (ΔG_CH2) to dispersive surface energy γ_s^d.

    Applies the Dorris-Gray formula:

        γ_s^d (mJ/m²) = 10⁶ · ΔG_CH2² / (4 · γ_CH2 · N_A² · a_CH2²)

    where γ_CH2 = (35.6 − 0.058·t) mJ/m² is the temperature-dependent
    surface energy of a CH2 group.

    Args:
        delta_g_ch2: ΔG_CH2 slope from RT·ln(V) vs carbon number fit (J/mol).
        temperature_kelvin: Column temperature in Kelvin.

    Returns:
        γ_s^d in mJ/m².

    Examples:
        >>> slope_to_gamma_d(500.0, 303.15)   # doctest: +ELLIPSIS
        ...
    """
    t_celsius = temperature_kelvin - 273.15
    y_ch2 = gamma_ch2(t_celsius)  # mJ/m²
    return 1e6 * delta_g_ch2 ** 2 / (4.0 * y_ch2 * NA ** 2 * A_CH2 ** 2)


# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

def prepare_alkane_data_from_injections(
    igc_result: IGCResult,
    retention_type: str = "COM",
    alkanes: Optional[List[str]] = None,
    target_coverages: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Prepare alkane data from raw injection items for Dorris-Gray analysis.

    Interpolates raw per-injection retention volumes to standard coverage
    points and computes RT·ln(V) for each alkane.  The resulting DataFrame
    has the same shape as ``prepare_alkane_data()`` from dorris_gray.py and
    can be fed directly into ``fit_dorris_gray()``.

    Args:
        igc_result: Parsed IGC result.  Only injection_items is required.
        retention_type: ``'COM'`` (centre of mass) or ``'MAX'`` (peak maximum).
            Defaults to ``'COM'``.
        alkanes: Optional list of alkane names to include in the fit, e.g.
            ``["HEPTANE", "OCTANE", "NONANE", "DECANE"]``.  If None, all
            alkanes in ALKANE_CARBON_NUMBERS are used.
        target_coverages: Coverage points to interpolate to.  If None, uses
            STANDARD_COVERAGES from constants.

    Returns:
        DataFrame with columns:
        - Solvent
        - Target Fractional Surface Coverage
        - Carbon Number
        - Column Temperature [Kelvin]
        - RTlnVg (J/mol)

    Raises:
        ValueError: If retention_type is not ``'COM'`` or ``'MAX'``.
    """
    retention_type = retention_type.upper()
    if retention_type not in _RETENTION_COLS:
        raise ValueError(f"retention_type must be 'COM' or 'MAX', got {retention_type!r}")

    value_col = _RETENTION_COLS[retention_type]
    interp_col = f"{value_col} (interp)"

    interpolated = interpolate_to_standard_coverages(
        igc_result,
        target_coverages=target_coverages,
        value_col=value_col,
    )

    # Filter for recognised alkanes
    alkane_interp = interpolated[
        interpolated["Solvent"].isin(ALKANE_CARBON_NUMBERS.keys())
    ].copy()

    if alkanes is not None:
        alkane_interp = alkane_interp[alkane_interp["Solvent"].isin(alkanes)].copy()

    alkane_interp["Carbon Number"] = alkane_interp["Solvent"].map(ALKANE_CARBON_NUMBERS)

    alkane_interp["RTlnVg"] = calculate_rtlnv(
        alkane_interp[interp_col],
        alkane_interp["Closest Column Temperature [Kelvin]"],
    )

    return (
        alkane_interp[[
            "Solvent",
            "Target Fractional Surface Coverage",
            "Carbon Number",
            "Closest Column Temperature [Kelvin]",
            "RTlnVg",
        ]]
        .rename(columns={"Closest Column Temperature [Kelvin]": "Column Temperature [Kelvin]"})
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# Main calculation
# ---------------------------------------------------------------------------

def calculate_dispersive_from_injections(
    igc_result: IGCResult,
    retention_type: str = "COM",
    alkanes: Optional[List[str]] = None,
    target_coverages: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Calculate dispersive surface energy from raw injection data.

    Performs the full Dorris-Gray pipeline using only injection_items:
    interpolation → RT·ln(V) → linear fit → γ_s^d conversion.  Works with
    both full CSV exports and injection-items-only exports.

    Args:
        igc_result: Parsed IGC result.  Only injection_items is required.
        retention_type: ``'COM'`` (centre of mass) or ``'MAX'`` (peak maximum).
        alkanes: Optional list of alkane names to restrict the fit.
        target_coverages: Coverage points.  Defaults to STANDARD_COVERAGES.

    Returns:
        DataFrame with columns:

        - ``n/nm``: Surface coverage
        - ``gamma_d``: Dispersive surface energy (mJ/m²)
        - ``delta_g_ch2``: ΔG_CH2 slope from linear fit (J/mol per CH2)
        - ``r_squared``: Goodness of fit (R²)
        - ``temperature_kelvin``: Mean column temperature at this coverage (K)

    Examples:
        >>> result = parse_igc_csv("data/sample.csv")
        >>> gd = calculate_dispersive_from_injections(result)
        >>> print(gd[["n/nm", "gamma_d", "r_squared"]])

    References:
        Dorris, G. M., & Gray, D. G. (1980). Journal of Colloid and Interface
        Science, 77(2), 353-362.
    """
    alkane_data = prepare_alkane_data_from_injections(
        igc_result,
        retention_type=retention_type,
        alkanes=alkanes,
        target_coverages=target_coverages,
    )

    fits = fit_dorris_gray(alkane_data)

    # Mean temperature per coverage for the γ_CH2(t) calculation
    mean_temps = (
        alkane_data
        .groupby("Target Fractional Surface Coverage")["Column Temperature [Kelvin]"]
        .mean()
        .rename("temperature_kelvin")
    )

    result = (
        fits
        .rename(columns={
            "Target Fractional Surface Coverage": "n/nm",
            "slope": "delta_g_ch2",
        })
        .merge(mean_temps, left_on="n/nm", right_index=True)
    )

    result["gamma_d"] = result.apply(
        lambda row: slope_to_gamma_d(row["delta_g_ch2"], row["temperature_kelvin"])
        if pd.notna(row["delta_g_ch2"]) and pd.notna(row["temperature_kelvin"])
        else np.nan,
        axis=1,
    )

    return result[["n/nm", "gamma_d", "delta_g_ch2", "r_squared", "temperature_kelvin"]]


# ---------------------------------------------------------------------------
# Validation against SMS output
# ---------------------------------------------------------------------------

def validate_against_sms(
    manual_df: pd.DataFrame,
    dispersive_surface_energy: pd.DataFrame,
    retention_type: str = "COM",
    threshold: float = _DEVIATION_THRESHOLD,
) -> pd.DataFrame:
    """Compare manual Dorris-Gray results to SMS pre-calculated DnG values.

    Args:
        manual_df: Output DataFrame from ``calculate_dispersive_from_injections()``.
        dispersive_surface_energy: The ``IGCResult.dispersive_surface_energy``
            DataFrame (must not be None).
        retention_type: ``'COM'`` or ``'MAX'`` — selects the SMS DnG column to
            compare against.
        threshold: Percentage deviation above which rows are flagged.
            Defaults to 5.0.

    Returns:
        DataFrame with columns:

        - ``n/nm``: Surface coverage
        - ``gamma_d_manual``: Manually calculated γ_s^d (mJ/m²)
        - ``gamma_d_sms``: SMS DnG γ_s^d (mJ/m²)
        - ``pct_diff``: (manual − sms) / sms × 100
        - ``exceeds_threshold``: True if |pct_diff| > threshold

    Raises:
        ValueError: If retention_type is not ``'COM'`` or ``'MAX'``.
        KeyError: If the expected SMS column is not present.

    Examples:
        >>> result = parse_igc_csv("data/sample.csv")
        >>> gd = calculate_dispersive_from_injections(result)
        >>> if result.dispersive_surface_energy is not None:
        ...     cmp = validate_against_sms(gd, result.dispersive_surface_energy)
        ...     print(cmp)
    """
    retention_type = retention_type.upper()
    if retention_type not in _SMS_DNG_COLS:
        raise ValueError(f"retention_type must be 'COM' or 'MAX', got {retention_type!r}")

    sms_col = _SMS_DNG_COLS[retention_type]
    if sms_col not in dispersive_surface_energy.columns:
        raise KeyError(
            f"Expected column {sms_col!r} not found in dispersive_surface_energy. "
            f"Available: {list(dispersive_surface_energy.columns)}"
        )

    sms = dispersive_surface_energy[["n/nm", sms_col]].rename(
        columns={sms_col: "gamma_d_sms"}
    )
    comparison = (
        manual_df[["n/nm", "gamma_d"]]
        .rename(columns={"gamma_d": "gamma_d_manual"})
        .merge(sms, on="n/nm", how="inner")
    )

    comparison["pct_diff"] = (
        (comparison["gamma_d_manual"] - comparison["gamma_d_sms"])
        / comparison["gamma_d_sms"]
        * 100
    )
    comparison["exceeds_threshold"] = comparison["pct_diff"].abs() > threshold

    return comparison
