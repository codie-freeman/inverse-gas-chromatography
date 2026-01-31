"""Physical constants and probe parameters for IGC-SEA analysis."""

from typing import Dict

# Universal gas constant (J/(mol·K))
R = 8.31446261815324

# Avogadro's number
NA = 6.022e23

# Alkane carbon numbers for Dorris-Gray analysis
ALKANE_CARBON_NUMBERS = {
    "HEXANE": 6,
    "HEPTANE": 7,
    "OCTANE": 8,
    "NONANE": 9,
    "DECANE": 10,
    "UNDECANE": 11,
    "DODECANE": 12,
}

# Probe parameters for acid-base calculations
# Area (A_probe) in m²
# Lewis parameters (LP_probe) in J/m²
# Values from SMS Cirrus Plus software and literature
PROBE_PARAMETERS: Dict[str, Dict[str, float]] = {
    "ETHYL ACETATE": {
        "area": 3.3e-19,  # m²
        "lp": 0.47567,  # J/m² (Lewis parameter)
        "type": "basic",  # Basic probe (electron donor)
    },
    "DICHLOROMETHANE": {
        "area": 2.45e-19,  # m²
        "lp": 0.12458,  # J/m²
        "type": "acidic",  # Acidic probe (electron acceptor)
    },
}

# Standard target coverages for analysis (n/nm)
STANDARD_COVERAGES = [0.005, 0.01, 0.025, 0.05, 0.10, 0.14]
