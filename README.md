# IGC‑SEA Data Processing Toolkit

**Version 0.1.0 (Alpha)**

This repository contains a Python toolkit for processing **Inverse Gas Chromatography – Surface Energy Analyzer (IGC‑SEA)** data. It's both a practical tool for analysing my own IGC data and an educational project for deepening my understanding of surface chemistry and developing computational skills.

The toolkit takes messy instrument CSV exports and transforms them into analysed surface energy profiles using established methods from the IGC literature.

---

## Why this toolkit?

### The current workflow problem

In typical IGC-SEA research, the workflow looks like this:

1. **Data export from proprietary software** (e.g., SMS Cirrus Plus analysis software, version 1.5)
2. **Manual copy-paste into Excel** for plotting and further analysis
3. **Limited analysis capabilities** in Excel:
   - Cannot reproduce exponential distribution plots
   - Poor handling of non-linear thermodynamic models
   - No programmatic control over equations or fitting procedures
4. **Incomplete data access** — not all analysis parameters are freely available from the proprietary software

This manual workflow has several problems:
- **Limited validation** — difficult to verify calculations against known methods
- **Inflexible plotting** — Excel cannot easily produce publication-quality scientific figures
- **No automation** — each new dataset requires the same repetitive manual steps

### What this toolkit provides

**igcsea** addresses these limitations by providing:

✓ **Automated parsing** of multi-file CSV exports from IGC instruments
✓ **Complete parameter extraction** including:
  - Peak center-of-mass (COM) parameters
  - Della Volpe acid-base parameters (γₛ⁺, γₛ⁻, γᴬᴮ)
  - Dorris-Gray dispersive parameters (γᴰ)
  - Polarisation parameters and retention volumes

✓ **Independent testing** against proprietary software outputs for comparison purposes
✓ **Full control** over thermodynamic equations, fitting procedures, and analysis methods
✓ **Reproducible workflows** with transparent calculations and version-controlled code
✓ **Publication-ready outputs** with customizable plots and data exports
✓ **Foundation for expansion** — a platform that can grow to include:
  - Exponential distribution plots
  - BET-style surface energy distributions
  - Flow-rate normalization
  - Retention time corrections
  - Advanced thermodynamic modeling

By replacing manual Excel workflows with robust Python code, this toolkit improves **data integrity, transparency, and reproducibility** in IGC surface energy research.

---

## What this toolkit does

**igcsea** provides a complete workflow for IGC-SEA data analysis:

1. **Parse complex CSV exports** from IGC instruments (like SMS Cirrus Plus)
   - Extracts multiple tables from single-file exports
   - Handles encoding issues and numeric conversion automatically
   - Returns structured data models (not just raw DataFrames)

2. **Standardize coverage points** through interpolation/extrapolation
   - Maps retention volumes to consistent surface coverage values
   - Enables comparison across different temperature conditions

3. **Calculate dispersive surface energy** using the **Dorris-Gray method**
   - Analyzes n-alkane (C6-C12) retention data
   - Performs linear regression to determine dispersive component (γᴰ)
   - Provides goodness-of-fit metrics (R²)

4. **Determine acid-base components** using the **Della Volpe approach**
   - Uses polar probe molecules (ethyl acetate and dichloromethane)
   - Calculates acidic (γₛ⁺) and basic (γₛ⁻) surface parameters
   - Computes acid-base contribution (γᴬᴮ)

5. **Generate complete surface energy profiles**
   - Combines dispersive and acid-base components
   - Calculates total surface energy: **γᵀ = γᴰ + γᴬᴮ**
   - Fits exponential asymptotic models to coverage-energy relationships

6. **Create publication-quality visualizations**
   - Multi-panel Dorris-Gray plots with fit statistics
   - Surface energy component profiles with exponential fits

This is primarily a **learning project and personal research tool**, not production software. I'm refining it as I develop both my understanding of IGC theory and my Python skills.

---

## Quick start

### Installation

```bash
# Clone the repository
git clone https://github.com/codiefreeman/inverse-gas-chromatography.git
cd inverse-gas-chromatography

# Install in development mode
pip install -e ".[dev]"
```

### Basic usage

```python
from igcsea.parsing import parse_igc_csv
from igcsea.analysis import (
    interpolate_to_standard_coverages,
    prepare_alkane_data,
    fit_dorris_gray,
    calculate_acid_base_params,
    calculate_surface_energy_profile
)

# Parse the CSV export
result = parse_igc_csv("data/examples/sample_igc_export.csv")

# Standardize to target coverages
target_coverages = [0.005, 0.01, 0.025, 0.05, 0.10, 0.14]
interp_data = interpolate_to_standard_coverages(result, target_coverages)

# Perform Dorris-Gray analysis on alkanes
alkane_data = prepare_alkane_data(interp_data)
dorris_gray_fits = fit_dorris_gray(alkane_data)

# Calculate acid-base parameters (Della Volpe method)
acid_base = calculate_acid_base_params(interp_data)

# Generate complete surface energy profile with exponential fits
profile = calculate_surface_energy_profile(interp_data)

# Verify: yt = yd + yab
print(f"Total energy: {profile.yt}")
print(f"Sum of components: {profile.yd + profile.yab}")
```

See [examples/test_package.py](examples/test_package.py) for a complete working example with visualization.

### Running tests

```bash
pytest tests/
```

The test suite includes comprehensive unit tests for all analysis modules, using the sample CSV data as a fixture.

---

## Project structure

```
src/igcsea/
├── core/                  # Data models and constants
│   ├── models.py         # IGCResult, SurfaceEnergyProfile, AcidBaseParams
│   └── constants.py      # Physical constants, alkane data, probe parameters
├── parsing/              # CSV parsing
│   └── parser.py         # Multi-table CSV export parser
├── analysis/             # Analysis modules
│   ├── interpolation.py  # Coverage interpolation/extrapolation
│   ├── dorris_gray.py    # Dispersive energy (Dorris-Gray method)
│   ├── acid_base.py      # Acid-base components (Della Volpe method)
│   └── surface_energy.py # Complete profile with exponential fitting
├── plotting/             # Visualization (in development)
└── utils/                # Utilities (in development)
```

**Key modules:**

- **`igcsea.parsing.parse_igc_csv()`**: Parses multi-table IGC CSV exports into structured `IGCResult` objects
- **`igcsea.analysis.interpolate_to_standard_coverages()`**: Standardizes retention data to target coverage values
- **`igcsea.analysis.fit_dorris_gray()`**: Performs Dorris-Gray linear regression on alkane data
- **`igcsea.analysis.calculate_acid_base_params()`**: Computes γₛ⁺, γₛ⁻, and γᴬᴮ using Della Volpe method
- **`igcsea.analysis.calculate_surface_energy_profile()`**: Generates complete γᴰ, γᴬᴮ, γᵀ profiles with exponential fits

---

## Scientific methods implemented

This toolkit implements established methods from IGC surface energy literature:

**Dorris-Gray Method** (1980)
- Analyses n-alkane retention to determine dispersive surface energy
- Linear regression of RT·ln(Vₙ) vs. carbon number
- Slope relates to dispersive component γᴰ

**Della Volpe Acid-Base Analysis** (1991)
- Uses polar probes (ethyl acetate and dichloromethane) to separate acid-base contributions
- Calculates acidic (γₛ⁺) and basic (γₛ⁻) surface parameters
- Computes acid-base component: γᴬᴮ = 2√(γₛ⁺ × γₛ⁻)

**Exponential Asymptotic Fitting**
- Models energy-coverage relationships: γ(n) = c + a·exp(-b·n)
- Provides parameters for dispersive, acid-base, and total surface energy profiles

All calculations use standard physical constants (R = 8.314 J/(mol·K), Avogadro's number) and probe-specific parameters from the literature.

---

## Status and future development

This project is in **alpha** (v0.1.0). The core analysis pipeline is functional and **tested against SMS proprietary software outputs for comparison purposes**.

**Current focus:**
- Improving the plotting module with publication-ready visualizations
- Refining existing analysis methods and adding better error handling
- Expanding documentation and examples
- Fixing bugs and improving test coverage
- Continued testing against proprietary software for new features

**Potential future expansion:**

As this toolkit matures, it could expand to replicate more capabilities currently only available in proprietary IGC software:

- **Advanced thermodynamic models**
  - Exponential distribution plots for surface energy heterogeneity
  - BET-style surface energy distribution analysis
  - Non-linear model fitting for complex surface interactions

- **Data preprocessing and corrections**
  - Flow-rate normalization
  - Retention time corrections
  - Dead volume compensation
  - Temperature-dependent parameter adjustments

- **Batch processing and automation**
  - Multi-file CSV processing pipelines
  - Automated report generation
  - Comparison across sample sets

- **Extended parameter extraction**
  - Full Gutmann acid-base number implementation (Ka, Kb)
  - Work of adhesion calculations
  - Spreading coefficients

These expansions would further reduce dependence on proprietary software and improve reproducibility in IGC research.

**What this project is:**
- A learning space for scientific Python, data processing, and IGC surface chemistry
- A practical tool for my own IGC-SEA data analysis, tested against commercial software
- An educational resource for understanding IGC methods and computational analysis
- A foundation for building open, reproducible IGC analysis workflows

**What this project is not (yet):**
- Production-ready software for broad distribution
- A complete replacement for commercial IGC analysis software
- Generalised to handle all possible IGC instrument formats

I'm building this incrementally as both a useful research tool and a way to develop computational skills. Feedback and suggestions are welcome!

---

## Development with AI tools

This project has been developed with assistance from AI coding tools, which has been part of the learning process:

- **GitHub Copilot** helped create the initial CSV parser and boilerplate code
- **Claude Code** has supported ongoing development, refactoring, and implementation of analysis methods

All scientific decisions—which methods to implement, how to interpret IGC theory, what calculations to perform—are mine. The AI tools help translate those ideas into working Python code more efficiently, especially as I'm still developing my programming skills. I review and understand all generated code, and I'm using this project to learn proper software design patterns.

This transparent approach to AI-assisted development is intentional: I want to learn *how* to code well, not just produce code. The toolkit works and produces correct results, but I recognize it will benefit from continued refinement as my skills grow.

---

## License

MIT License - see LICENSE file for details.

## Dependencies

- Python ≥3.8
- NumPy ≥1.21
- Pandas ≥1.3
- Matplotlib ≥3.4
- SciPy ≥1.7

See [pyproject.toml](pyproject.toml) for complete dependency information.