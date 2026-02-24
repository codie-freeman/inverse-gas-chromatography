"""IGC-SEA: Inverse Gas Chromatography Surface Energy Analysis toolkit.

This package provides tools for analyzing IGC-SEA (Inverse Gas Chromatography -
Surface Energy Analysis) data exports.

Core functionality:
- CSV parsing for multi-table IGC instrument exports
- Coverage interpolation/extrapolation for standardization
- Dorris-Gray alkane analysis for dispersive surface energy
- Acid-base surface energy components (Della Volpe method)
- Complete surface energy profiles with exponential fitting

Example usage:
    >>> from igcsea.parsing import parse_igc_csv
    >>> from igcsea.analysis import calculate_surface_energy_profile
    >>>
    >>> # Parse CSV export
    >>> result = parse_igc_csv("path/to/igc_export.csv")
    >>>
    >>> # Generate complete surface energy profile
    >>> profile = calculate_surface_energy_profile(result)
    >>>
    >>> # Access results
    >>> print(profile.yd)   # Dispersive component
    >>> print(profile.yab)  # Acid-base component
    >>> print(profile.yt)   # Total surface energy

See examples/test_package.py for a complete working example with visualization.
"""

__version__ = "0.1.0"

# Core models and constants
from igcsea.core import (
    AcidBaseParams,
    ALKANE_CARBON_NUMBERS,
    IGCResult,
    PROBE_PARAMETERS,
    R,
    STANDARD_COVERAGES,
    SurfaceEnergyProfile,
)

# High-level API — parsing
from igcsea.parsing import parse_igc_csv

# High-level API — analysis
from igcsea.analysis.acid_base import calculate_acid_base_params
from igcsea.analysis.manual_dispersive import (
    calculate_dispersive_from_injections,
    validate_against_sms,
)
from igcsea.analysis.peak_max import PeakComparisonResult, PeakMAXAnalyzer
from igcsea.analysis.surface_energy import calculate_surface_energy_profile

__all__ = [
    "__version__",
    # Models
    "IGCResult",
    "SurfaceEnergyProfile",
    "AcidBaseParams",
    # Constants
    "R",
    "ALKANE_CARBON_NUMBERS",
    "PROBE_PARAMETERS",
    "STANDARD_COVERAGES",
    # Parsing
    "parse_igc_csv",
    # Analysis — dispersive
    "calculate_dispersive_from_injections",
    "validate_against_sms",
    # Analysis — acid-base
    "calculate_acid_base_params",
    # Analysis — full profile
    "calculate_surface_energy_profile",
    # Analysis — COM vs MAX comparison
    "PeakMAXAnalyzer",
    "PeakComparisonResult",
]
