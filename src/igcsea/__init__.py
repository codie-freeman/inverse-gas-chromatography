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

# Import core models for easy access
from igcsea.core import (
    AcidBaseParams,
    ALKANE_CARBON_NUMBERS,
    IGCResult,
    PROBE_PARAMETERS,
    R,
    STANDARD_COVERAGES,
    SurfaceEnergyProfile,
)

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
    # High-level API functions will be added as we build them
]
