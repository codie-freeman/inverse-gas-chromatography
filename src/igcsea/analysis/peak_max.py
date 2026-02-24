"""Peak MAX vs Centre of Mass (COM) retention analysis for IGC-SEA.

Provides tools to run Dorris-Gray dispersive surface energy calculations using
both Peak MAX and Peak COM (centre of mass) retention volumes, compare the two
methods, and recommend which is more appropriate for a given dataset.

References:
    Dorris, G. M., & Gray, D. G. (1980). Adsorption of n-alkanes at zero
    surface coverage on cellulose paper and wood fibers. Journal of Colloid
    and Interface Science, 77(2), 353-362.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
import pandas as pd

from igcsea.analysis.manual_dispersive import calculate_dispersive_from_injections
from igcsea.core.models import IGCResult

_DEFAULT_THRESHOLD = 10.0  # % deviation above which methods are considered to disagree


@dataclass
class PeakComparisonResult:
    """Comparison of COM and MAX Dorris-Gray dispersive surface energy results.

    Attributes:
        coverage: Array of surface coverage values (n/nm).
        gamma_d_com: Dispersive surface energy using COM retention (mJ/m²).
        gamma_d_max: Dispersive surface energy using MAX retention (mJ/m²).
        r_squared_com: Linear fit R² for COM at each coverage.
        r_squared_max: Linear fit R² for MAX at each coverage.
        pct_diff: Percentage difference (COM − MAX) / MAX × 100 at each coverage.
        exceeds_threshold: Boolean array; True where |pct_diff| > threshold.
        recommendation: Suggested retention method — ``'COM'``, ``'MAX'``, or
            ``'AGREE'`` (methods agree within threshold; COM used by convention).
        threshold: Percentage threshold used to flag deviations. Default 10.0.

    Examples:
        >>> analyzer = PeakMAXAnalyzer(result)
        >>> cmp = analyzer.compare_com_vs_max()
        >>> print(cmp.recommendation)
        >>> print(cmp.to_dataframe())
    """

    coverage: np.ndarray
    gamma_d_com: np.ndarray
    gamma_d_max: np.ndarray
    r_squared_com: np.ndarray
    r_squared_max: np.ndarray
    pct_diff: np.ndarray
    exceeds_threshold: np.ndarray
    recommendation: str
    threshold: float = field(default=_DEFAULT_THRESHOLD)

    def to_dataframe(self) -> pd.DataFrame:
        """Return comparison as a tidy DataFrame.

        Returns:
            DataFrame with columns: n/nm, gamma_d_com, gamma_d_max,
            r_squared_com, r_squared_max, pct_diff, exceeds_threshold.
        """
        return pd.DataFrame({
            "n/nm": self.coverage,
            "gamma_d_com": self.gamma_d_com,
            "gamma_d_max": self.gamma_d_max,
            "r_squared_com": self.r_squared_com,
            "r_squared_max": self.r_squared_max,
            "pct_diff": self.pct_diff,
            "exceeds_threshold": self.exceeds_threshold,
        })

    def summary(self) -> str:
        """Return a human-readable one-line summary of the comparison.

        Returns:
            String describing the recommendation, mean R² values, and how many
            coverage points exceeded the deviation threshold.
        """
        flagged = int(self.exceeds_threshold.sum())
        mean_com = float(np.nanmean(self.r_squared_com))
        mean_max = float(np.nanmean(self.r_squared_max))
        return (
            f"Recommendation: {self.recommendation} | "
            f"Mean R² — COM: {mean_com:.4f}, MAX: {mean_max:.4f} | "
            f"Points >{self.threshold:.0f}% deviation: {flagged}/{len(self.coverage)}"
        )


class PeakMAXAnalyzer:
    """Analyse and compare COM vs MAX retention methods for Dorris-Gray.

    All Dorris-Gray calculations are performed using raw injection_items data
    (via ``calculate_dispersive_from_injections``), so no pre-calculated SMS
    tables are required.

    Args:
        igc_result: Parsed IGC result.  Only ``injection_items`` is required.
        alkanes: Optional list of alkane names to restrict the fit, e.g.
            ``["HEPTANE", "OCTANE", "NONANE", "DECANE"]``.  If None, all
            recognised alkanes in the data are used.
        target_coverages: Coverage points to interpolate to.  If None, uses
            ``STANDARD_COVERAGES`` from constants.

    Examples:
        >>> from igcsea.parsing import parse_igc_csv
        >>> from igcsea.analysis.peak_max import PeakMAXAnalyzer
        >>>
        >>> result = parse_igc_csv("data/sample.csv")
        >>> analyzer = PeakMAXAnalyzer(result, alkanes=["HEPTANE", "OCTANE", "NONANE", "DECANE"])
        >>> cmp = analyzer.compare_com_vs_max()
        >>> print(cmp.summary())
        >>> print(cmp.to_dataframe())
    """

    def __init__(
        self,
        igc_result: IGCResult,
        alkanes: Optional[List[str]] = None,
        target_coverages: Optional[np.ndarray] = None,
    ) -> None:
        self._igc_result = igc_result
        self._alkanes = alkanes
        self._target_coverages = target_coverages

    def parse_peak_max_from_csv(self) -> pd.DataFrame:
        """Extract MAX retention volume data from ``injection_items``.

        Returns a filtered view of the injection data containing only the MAX
        retention volume column, limited to recognised alkane solvents.

        Returns:
            DataFrame with columns: Solvent, Actual Fractional Surface Coverage,
            Sp. Ret Volume (Max) [ml/g], Column Temperature [Kelvin].

        Raises:
            KeyError: If the expected MAX retention volume column is absent from
                ``injection_items``.
        """
        col = "Sp. Ret Volume (Max) [ml/g]"
        inj = self._igc_result.injection_items
        if col not in inj.columns:
            raise KeyError(
                f"Column {col!r} not found in injection_items. "
                f"Available columns: {list(inj.columns)}"
            )
        return inj[[
            "Solvent",
            "Actual Fractional Surface Coverage",
            col,
            "Column Temperature [Kelvin]",
        ]].copy()

    def calculate_dispersive_from_max(self) -> pd.DataFrame:
        """Run Dorris-Gray dispersive surface energy using MAX retention volumes.

        Returns:
            DataFrame with columns: n/nm, gamma_d, delta_g_ch2, r_squared,
            temperature_kelvin — identical schema to
            ``calculate_dispersive_from_injections(retention_type='MAX')``.
        """
        return calculate_dispersive_from_injections(
            self._igc_result,
            retention_type="MAX",
            alkanes=self._alkanes,
            target_coverages=self._target_coverages,
        )

    def compare_com_vs_max(self, threshold: float = _DEFAULT_THRESHOLD) -> PeakComparisonResult:
        """Calculate γ_s^d with both COM and MAX retention, then compare.

        Runs Dorris-Gray using centre-of-mass (COM) and peak maximum (MAX)
        retention volumes, computes percentage differences per coverage point,
        and recommends which retention method to use based on linear fit quality
        (R²) and the magnitude of deviation between the two methods.

        Recommendation rules (applied in order):

        1. If mean R² of MAX exceeds mean R² of COM by more than 0.005 → ``'MAX'``
        2. If mean R² of COM exceeds mean R² of MAX by more than 0.005 → ``'COM'``
        3. If any coverage point exceeds ``threshold``% deviation:
           → whichever method has the higher mean R²
        4. Both methods agree within threshold and similar R² → ``'AGREE'``

        Args:
            threshold: Percentage deviation above which a coverage point is
                flagged.  Defaults to 10.0.

        Returns:
            ``PeakComparisonResult`` with per-coverage arrays, a recommendation
            string, and a convenience ``to_dataframe()`` method.

        Examples:
            >>> cmp = PeakMAXAnalyzer(result).compare_com_vs_max()
            >>> print(cmp.summary())
            >>> flagged = cmp.to_dataframe()[cmp.to_dataframe()["exceeds_threshold"]]
        """
        com_df = calculate_dispersive_from_injections(
            self._igc_result,
            retention_type="COM",
            alkanes=self._alkanes,
            target_coverages=self._target_coverages,
        )
        max_df = calculate_dispersive_from_injections(
            self._igc_result,
            retention_type="MAX",
            alkanes=self._alkanes,
            target_coverages=self._target_coverages,
        )

        merged = (
            com_df[["n/nm", "gamma_d", "r_squared"]]
            .rename(columns={"gamma_d": "gamma_d_com", "r_squared": "r_squared_com"})
            .merge(
                max_df[["n/nm", "gamma_d", "r_squared"]].rename(
                    columns={"gamma_d": "gamma_d_max", "r_squared": "r_squared_max"}
                ),
                on="n/nm",
                how="inner",
            )
        )

        pct_diff = (
            (merged["gamma_d_com"] - merged["gamma_d_max"])
            / merged["gamma_d_max"]
            * 100
        )
        exceeds = pct_diff.abs() > threshold

        mean_r2_com = float(np.nanmean(merged["r_squared_com"].to_numpy(float)))
        mean_r2_max = float(np.nanmean(merged["r_squared_max"].to_numpy(float)))
        any_flagged = bool(exceeds.any())

        if mean_r2_max - mean_r2_com > 0.005:
            recommendation = "MAX"
        elif mean_r2_com - mean_r2_max > 0.005:
            recommendation = "COM"
        elif any_flagged:
            recommendation = "MAX" if mean_r2_max >= mean_r2_com else "COM"
        else:
            recommendation = "AGREE"

        return PeakComparisonResult(
            coverage=merged["n/nm"].to_numpy(),
            gamma_d_com=merged["gamma_d_com"].to_numpy(),
            gamma_d_max=merged["gamma_d_max"].to_numpy(),
            r_squared_com=merged["r_squared_com"].to_numpy(),
            r_squared_max=merged["r_squared_max"].to_numpy(),
            pct_diff=pct_diff.to_numpy(),
            exceeds_threshold=exceeds.to_numpy(),
            recommendation=recommendation,
            threshold=threshold,
        )
