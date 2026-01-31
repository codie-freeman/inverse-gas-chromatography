"""Surface energy profile analysis with exponential fitting.

This module provides functions to calculate complete surface energy profiles
(yd, yab, yt) and fit exponential decay models to the data.
"""

from typing import Dict, Optional

import numpy as np
from scipy.optimize import curve_fit

from igcsea.analysis.acid_base import calculate_acid_base_params
from igcsea.core.models import IGCResult, SurfaceEnergyProfile


def exp_asym(x: np.ndarray, c: float, a: float, b: float) -> np.ndarray:
    """Exponential asymptotic model for surface energy vs coverage.

    Model: y(x) = c + a·exp(-b·x)

    Args:
        x: Surface coverage values (n/nm).
        c: Asymptotic value at high coverage (mJ/m²).
        a: Excess above asymptote at x=0 (mJ/m²).
        b: Decay rate constant (nm).

    Returns:
        Fitted surface energy values (mJ/m²).
    """
    return c + a * np.exp(-b * x)


def fit_exponential_decay(
    coverage: np.ndarray,
    energy: np.ndarray,
) -> Dict[str, float]:
    """Fit exponential decay model to surface energy data.

    Fits the model: y(x) = c + a·exp(-b·x)

    Args:
        coverage: Array of surface coverage values (n/nm).
        energy: Array of surface energy values (mJ/m²).

    Returns:
        Dictionary with keys 'c', 'a', 'b' containing fitted parameters.

    Raises:
        RuntimeError: If curve fitting fails to converge.

    Examples:
        >>> coverage = np.array([0.005, 0.01, 0.025, 0.05, 0.10, 0.14])
        >>> energy = np.array([46.5, 44.8, 43.3, 42.3, 40.7, 40.5])
        >>> params = fit_exponential_decay(coverage, energy)
        >>> print(f"Asymptote: {params['c']:.2f} mJ/m²")
    """
    # Convert to numpy and ensure float
    x = np.asarray(coverage, dtype=float)
    y = np.asarray(energy, dtype=float)

    # Keep only finite values
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]

    if len(x) < 3:
        raise ValueError("Need at least 3 data points for exponential fit")

    # Initial parameter guesses
    c0 = float(np.min(y))  # Asymptote ~ minimum value
    a0 = float(np.max(y) - np.min(y))  # Initial excess
    b0 = 10.0  # Moderate decay rate

    # Fit using non-linear least squares
    popt, _ = curve_fit(exp_asym, x, y, p0=[c0, a0, b0], maxfev=20000)

    return {"c": float(popt[0]), "a": float(popt[1]), "b": float(popt[2])}


def calculate_surface_energy_profile(
    igc_result: IGCResult,
    fit_exponential: bool = True,
) -> SurfaceEnergyProfile:
    """Calculate complete surface energy profile (yd, yab, yt).

    This function:
    1. Extracts dispersive component (yd) from dispersive_surface_energy table
    2. Calculates acid-base component (yab) using Della Volpe method
    3. Calculates total surface energy (yt = yd + yab)
    4. Optionally fits exponential decay models to each component

    Args:
        igc_result: Parsed IGC-SEA result.
        fit_exponential: If True, fit exponential models to yd, yab, yt.
            Defaults to True.

    Returns:
        SurfaceEnergyProfile with coverage, yd, yab, yt arrays and optional
        fit parameters.

    Examples:
        >>> result = parse_igc_csv("data.csv")
        >>> profile = calculate_surface_energy_profile(result)
        >>> print(f"YD at 0.005: {profile.yd[0]:.2f} mJ/m²")
        >>> print(f"Fit params: {profile.fit_params['yd']}")
    """
    # Extract dispersive component (yd)
    dse = igc_result.dispersive_surface_energy
    yd_df = dse[["n/nm", "Disp. Surf. En. (mJ/m^2) - DnG & Com"]].rename(
        columns={"Disp. Surf. En. (mJ/m^2) - DnG & Com": "yd"}
    ).copy()

    # Calculate acid-base component (yab)
    acid_base = calculate_acid_base_params(igc_result)
    yab_df = acid_base.to_dataframe()[["coverage", "yab"]].rename(
        columns={"coverage": "n/nm"}
    )

    # Merge yd and yab
    yt_df = (
        yd_df.merge(yab_df, on="n/nm", how="inner")
        .sort_values("n/nm")
        .reset_index(drop=True)
    )

    # Calculate total (yt = yd + yab)
    yt_df["yt"] = yt_df["yd"] + yt_df["yab"]

    # Extract arrays
    coverage = yt_df["n/nm"].to_numpy()
    yd = yt_df["yd"].to_numpy()
    yab = yt_df["yab"].to_numpy()
    yt = yt_df["yt"].to_numpy()

    # Fit exponential models if requested
    fit_params = None
    if fit_exponential:
        fit_params = {
            "yd": fit_exponential_decay(coverage, yd),
            "yab": fit_exponential_decay(coverage, yab),
            "yt": fit_exponential_decay(coverage, yt),
        }

    return SurfaceEnergyProfile(
        coverage=coverage,
        yd=yd,
        yab=yab,
        yt=yt,
        fit_params=fit_params,
    )
