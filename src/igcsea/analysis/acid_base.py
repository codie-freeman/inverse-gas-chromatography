"""Acid-base surface energy analysis using Della Volpe method.

This module provides functions to calculate acid-base surface energy components
from polar probe measurements.
"""

import numpy as np
import pandas as pd

from igcsea.core.constants import NA, PROBE_PARAMETERS
from igcsea.core.models import AcidBaseParams, IGCResult


def extract_probe_data(free_energy: pd.DataFrame, solvent_name: str) -> pd.DataFrame:
    """Extract free energy data for a specific probe solvent.

    Args:
        free_energy: Free energy DataFrame from IGCResult.
        solvent_name: Name of the solvent to extract (e.g., "ETHYL ACETATE").

    Returns:
        DataFrame with n/nm and En. (Pol Com) columns for the specified solvent.
    """
    probe_data = free_energy[free_energy["Solvent Name"] == solvent_name][
        ["n/nm", "En. (Pol Com)"]
    ].copy()
    return probe_data


def calculate_yab_della_volpe(
    ethyl_acetate_df: pd.DataFrame,
    dichloromethane_df: pd.DataFrame,
    probe_params: dict = None,
) -> pd.DataFrame:
    """Calculate acid-base component using Della Volpe method.

    This method uses two polar probe molecules (ethyl acetate as basic probe,
    dichloromethane as acidic probe) to determine the acidic (ys+) and basic (ys-)
    components of surface energy, and their geometric mean (yab).

    Args:
        ethyl_acetate_df: DataFrame with n/nm and En. (Pol Com) for ethyl acetate.
        dichloromethane_df: DataFrame with n/nm and En. (Pol Com) for dichloromethane.
        probe_params: Optional dict with probe parameters. If None, uses PROBE_PARAMETERS.

    Returns:
        DataFrame with columns:
        - n/nm: Surface coverage
        - yab: Acid-base component (mJ/m²), yab = 2·sqrt(ys+·ys-)
        - ys+: Acidic component from dichloromethane (mJ/m²)
        - ys-: Basic component from ethyl acetate (mJ/m²)

    Examples:
        >>> ea = extract_probe_data(result.free_energy, "ETHYL ACETATE")
        >>> dcm = extract_probe_data(result.free_energy, "DICHLOROMETHANE")
        >>> yab_df = calculate_yab_della_volpe(ea, dcm)

    References:
        Della Volpe, C., & Siboni, S. (1997). Some reflections on acid–base solid
        surface free energy theories. Journal of Colloid and Interface Science,
        195(1), 121-136.
    """
    if probe_params is None:
        probe_params = PROBE_PARAMETERS

    # Extract probe parameters (area in m², Lewis parameters in J/m²)
    A_ea = probe_params["ETHYL ACETATE"]["area"]
    LP_ea = probe_params["ETHYL ACETATE"]["lp"]
    A_dcm = probe_params["DICHLOROMETHANE"]["area"]
    LP_dcm = probe_params["DICHLOROMETHANE"]["lp"]

    # Rename columns for clarity
    ea = ethyl_acetate_df[["n/nm", "En. (Pol Com)"]].rename(
        columns={"En. (Pol Com)": "pol_ea"}
    ).copy()
    dcm = dichloromethane_df[["n/nm", "En. (Pol Com)"]].rename(
        columns={"En. (Pol Com)": "pol_dcm"}
    ).copy()

    # Merge on coverage
    merged = ea.merge(dcm, on="n/nm", how="inner").sort_values("n/nm")

    # Calculate ys- from ethyl acetate (basic probe)
    # deltaG = (En. (Pol Com) [kJ/mol] * 1000) / (Na * A_probe)
    # ys = (deltaG / -2)² / LP_probe
    deltaG_ea = (merged["pol_ea"] * 1000) / (NA * A_ea)
    ys_minus = (((deltaG_ea / -2) ** 2) / LP_ea) * 1000  # Convert to mJ/m²

    # Calculate ys+ from dichloromethane (acidic probe)
    deltaG_dcm = (merged["pol_dcm"] * 1000) / (NA * A_dcm)
    ys_plus = (((deltaG_dcm / -2) ** 2) / LP_dcm) * 1000  # Convert to mJ/m²

    # Calculate yab as geometric mean
    merged["ys-"] = ys_minus
    merged["ys+"] = ys_plus
    merged["yab"] = 2 * np.sqrt(merged["ys+"] * merged["ys-"])

    return merged[["n/nm", "yab", "ys+", "ys-"]]


def calculate_acid_base_params(igc_result: IGCResult) -> AcidBaseParams:
    """Calculate complete acid-base parameters from IGC result.

    This convenience function extracts probe data and calculates acid-base
    components, returning them as an AcidBaseParams dataclass.

    Args:
        igc_result: Parsed IGC-SEA result.

    Returns:
        AcidBaseParams with coverage, ys+, ys-, yab arrays.
        (Ka and Kb will be None - to be added in future).

    Examples:
        >>> result = parse_igc_csv("data.csv")
        >>> params = calculate_acid_base_params(result)
        >>> print(f"YAB at 0.005 n/nm: {params.yab[0]:.2f}")
    """
    # Extract probe data
    ea = extract_probe_data(igc_result.free_energy, "ETHYL ACETATE")
    dcm = extract_probe_data(igc_result.free_energy, "DICHLOROMETHANE")

    # Calculate yab
    yab_df = calculate_yab_della_volpe(ea, dcm)

    # Convert to AcidBaseParams
    return AcidBaseParams(
        coverage=yab_df["n/nm"].to_numpy(),
        ys_plus=yab_df["ys+"].to_numpy(),
        ys_minus=yab_df["ys-"].to_numpy(),
        yab=yab_df["yab"].to_numpy(),
        ka=None,  # To be implemented
        kb=None,  # To be implemented
    )
