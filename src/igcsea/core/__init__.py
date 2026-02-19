"""Core data models and constants for IGC-SEA analysis."""

from igcsea.core.constants import (
    A_CH2,
    ALKANE_CARBON_NUMBERS,
    NA,
    PROBE_PARAMETERS,
    R,
    STANDARD_COVERAGES,
)
from igcsea.core.models import AcidBaseParams, IGCResult, SurfaceEnergyProfile

__all__ = [
    # Models
    "IGCResult",
    "SurfaceEnergyProfile",
    "AcidBaseParams",
    # Constants
    "R",
    "NA",
    "A_CH2",
    "ALKANE_CARBON_NUMBERS",
    "PROBE_PARAMETERS",
    "STANDARD_COVERAGES",
]
