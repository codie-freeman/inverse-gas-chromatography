"""Data models for IGC-SEA analysis results."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class IGCResult:
    """Results from parsing an IGC-SEA CSV export.

    This dataclass contains the three main tables extracted from an IGC instrument
    CSV export, along with metadata about the source file.

    Attributes:
        free_energy: DataFrame with free energy data by coverage and solvent.
            Columns include: n/nm, Solvent Name, En. (Pol Com), etc.
        dispersive_surface_energy: DataFrame with dispersive surface energy values.
            Columns include: n/nm, Disp. Surf. En. (mJ/m^2) - DnG & Com
        injection_items: DataFrame with raw injection data and adsorption measurements.
            28 columns including retention volumes, temperatures, adsorbed amounts.
        source_path: Path to the source CSV file that was parsed.

    Examples:
        >>> from igcsea import parse
        >>> result = parse("data/sample.csv")
        >>> print(result.free_energy.head())
        >>> print(f"Parsed from: {result.source_path}")
    """

    free_energy: pd.DataFrame
    dispersive_surface_energy: pd.DataFrame
    injection_items: pd.DataFrame
    source_path: Path

    def as_dict(self) -> Dict[str, pd.DataFrame]:
        """Convert to dictionary of DataFrames.

        Returns:
            Dictionary with keys 'free_energy', 'dispersive_surface_energy', 'injection_items'
        """
        return {
            "free_energy": self.free_energy,
            "dispersive_surface_energy": self.dispersive_surface_energy,
            "injection_items": self.injection_items,
        }

    def to_csv_dir(self, output_dir: Path) -> None:
        """Export all tables to separate CSV files in a directory.

        Args:
            output_dir: Directory to write CSV files to. Will be created if it doesn't exist.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.free_energy.to_csv(output_dir / "free_energy.csv", index=False)
        self.dispersive_surface_energy.to_csv(
            output_dir / "dispersive_surface_energy.csv", index=False
        )
        self.injection_items.to_csv(output_dir / "injection_items.csv", index=False)


@dataclass
class SurfaceEnergyProfile:
    """Surface energy components across coverage range with exponential fits.

    This dataclass represents the complete surface energy characterization including
    dispersive (yd), acid-base (yab), and total (yt) components, along with
    exponential fit parameters for each component.

    Attributes:
        coverage: Array of surface coverage values (n/nm).
        yd: Array of dispersive surface energy values (mJ/m²).
        yab: Array of acid-base surface energy values (mJ/m²).
        yt: Array of total surface energy values (mJ/m²), yt = yd + yab.
        fit_params: Dictionary containing exponential fit parameters for each component.
            Keys: 'yd', 'yab', 'yt'. Each value is a dict with keys 'c', 'a', 'b'
            for the model y(x) = c + a·exp(-b·x).

    Examples:
        >>> profile = analyze_surface_energy(result)
        >>> print(f"YD at 0.005 n/nm: {profile.yd[0]:.2f} mJ/m²")
        >>> print(f"YAB fit params: {profile.fit_params['yab']}")
    """

    coverage: np.ndarray
    yd: np.ndarray
    yab: np.ndarray
    yt: np.ndarray
    fit_params: Optional[Dict[str, Dict[str, float]]] = None

    def to_dataframe(self) -> pd.DataFrame:
        """Convert to pandas DataFrame.

        Returns:
            DataFrame with columns: coverage, yd, yab, yt
        """
        return pd.DataFrame({
            "coverage": self.coverage,
            "yd": self.yd,
            "yab": self.yab,
            "yt": self.yt,
        })

    def evaluate_fit(self, coverage_values: np.ndarray) -> Dict[str, np.ndarray]:
        """Evaluate exponential fits at given coverage values.

        Args:
            coverage_values: Array of coverage values to evaluate fits at.

        Returns:
            Dictionary with keys 'yd', 'yab', 'yt' containing fitted values.

        Raises:
            ValueError: If fit_params is None (no fits available).
        """
        if self.fit_params is None:
            raise ValueError("No fit parameters available. Run fitting first.")

        results = {}
        for component in ['yd', 'yab', 'yt']:
            params = self.fit_params[component]
            c, a, b = params['c'], params['a'], params['b']
            results[component] = c + a * np.exp(-b * coverage_values)

        return results


@dataclass
class AcidBaseParams:
    """Acid-base parameters across coverage range.

    This dataclass contains the complete acid-base characterization including
    individual acidic (ys+) and basic (ys-) components, their geometric mean (yab),
    and Gutmann acid/base numbers (Ka, Kb).

    Attributes:
        coverage: Array of surface coverage values (n/nm).
        ys_plus: Array of acidic component values (ys+) from dichloromethane probe.
        ys_minus: Array of basic component values (ys-) from ethyl acetate probe.
        yab: Array of acid-base component values (mJ/m²), yab = 2·sqrt(ys+·ys-).
        ka: Array of Gutmann acid numbers.
        kb: Array of Gutmann base numbers.

    Examples:
        >>> acid_base = analyze_acid_base(result)
        >>> print(f"Ka at 0.005 n/nm: {acid_base.ka[0]:.2f}")
        >>> print(f"Kb at 0.005 n/nm: {acid_base.kb[0]:.2f}")
    """

    coverage: np.ndarray
    ys_plus: np.ndarray
    ys_minus: np.ndarray
    yab: np.ndarray
    ka: Optional[np.ndarray] = None
    kb: Optional[np.ndarray] = None

    def to_dataframe(self) -> pd.DataFrame:
        """Convert to pandas DataFrame.

        Returns:
            DataFrame with columns: coverage, ys_plus, ys_minus, yab, ka, kb
        """
        data = {
            "coverage": self.coverage,
            "ys_plus": self.ys_plus,
            "ys_minus": self.ys_minus,
            "yab": self.yab,
        }

        if self.ka is not None:
            data["ka"] = self.ka
        if self.kb is not None:
            data["kb"] = self.kb

        return pd.DataFrame(data)
