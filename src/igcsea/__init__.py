"""IGC-SEA: Inverse Gas Chromatography Surface Energy Analysis toolkit.

This package provides tools for analyzing IGC-SEA (Inverse Gas Chromatography -
Surface Energy Analysis) data exports, including:

- Dorris-Gray alkane analysis
- Polar probe free energy analysis
- Acid-base surface energy components (Della Volpe method)
- Gutmann acid/base parameters (Ka, Kb)
- Surface energy distributions vs area increment
- Multi-CSV batch processing
- Publication-quality plotting

Examples:
    >>> import igcsea
    >>> result = igcsea.parse("path/to/file.csv")
    >>> profile = igcsea.analyze_surface_energy(result)
    >>> fig = igcsea.plot_energy_profile(profile)
    >>> fig.savefig("energy_profile.png")
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
