"""Analysis modules for IGC-SEA data."""

from igcsea.analysis.acid_base import (
    calculate_acid_base_params,
    calculate_yab_della_volpe,
    extract_probe_data,
)
from igcsea.analysis.dorris_gray import (
    calculate_rtlnv,
    fit_dorris_gray,
    prepare_alkane_data,
)
from igcsea.analysis.interpolation import interpolate_to_standard_coverages
from igcsea.analysis.manual_dispersive import (
    calculate_dispersive_from_injections,
    gamma_ch2,
    prepare_alkane_data_from_injections,
    slope_to_gamma_d,
    validate_against_sms,
)
from igcsea.analysis.peak_max import PeakComparisonResult, PeakMAXAnalyzer
from igcsea.analysis.surface_energy import (
    calculate_surface_energy_profile,
    fit_exponential_decay,
)

__all__ = [
    # Acid-base (Della Volpe)
    "calculate_acid_base_params",
    "calculate_yab_della_volpe",
    "extract_probe_data",
    # Dorris-Gray — SMS free_energy path
    "calculate_rtlnv",
    "fit_dorris_gray",
    "prepare_alkane_data",
    # Coverage interpolation
    "interpolate_to_standard_coverages",
    # Dorris-Gray — manual injection-items path
    "calculate_dispersive_from_injections",
    "gamma_ch2",
    "prepare_alkane_data_from_injections",
    "slope_to_gamma_d",
    "validate_against_sms",
    # Peak COM vs MAX comparison
    "PeakComparisonResult",
    "PeakMAXAnalyzer",
    # Full surface energy profile
    "calculate_surface_energy_profile",
    "fit_exponential_decay",
]
