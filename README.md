# IGC‑SEA Data Processing Toolkit (Early Draft)

This repository contains an early‑stage Python toolkit for processing **Inverse Gas Chromatography – Surface Energy Analysis (IGC‑SEA)** data.

Right now, the main purpose of this project is to **get from messy instrument CSV exports to usable data for analysis**, rather than to build a perfect, general‑purpose CSV parser. The parsing code works, but it’s something I plan to revisit and improve as my skills grow.

---

## Project goals

- **Parse IGC‑SEA CSV exports** that contain multiple tables in a single sheet (free energy, dispersive surface energy, injection items, etc.).
- **Extract and clean key tables** into structured `pandas` DataFrames.
- **Compute surface energy components**, including:
  - `yd` — dispersive component  
  - `yab` — acid–base component  
  - `yt` — total surface energy (`yt = yd + yab`)
- **Generate basic plots and profiles** of surface energy vs surface coverage.
- **Provide a foundation** for more advanced modelling and visualisation later.

This is primarily a **data analysis and processing project**, not a polished parsing library (yet).

---

## How the parser was created

The CSV parser in this project was developed using **GitHub Copilot** as an assistive coding tool.

- I defined the problem: multi‑table IGC‑SEA CSV exports with section titles like *“Free Energy”*, *“Dispersive Surface Energy”*, and *“Injection Items”*.
- My initial attempts at writing the parser from scratch had limited success.
- I then used GitHub Copilot to help generate the parsing logic, including:
  - detecting section titles
  - finding header rows
  - slicing out tables
  - handling blank rows
  - doing basic numeric conversion

The result is a parser that **mostly works for my current dataset and use case**, but:

- it is **not yet generalised** for all possible IGC‑SEA exports  
- it will likely need refactoring and hardening in the future  

I consider this parser a **first working draft**, created with heavy assistance from Copilot, which I plan to revisit when I’m more confident with Python, `pandas`, and software design.

---

## Current components

- **`parse_sugar_se`**  
  A function that reads a specific IGC‑SEA CSV export and returns a dictionary of `pandas` DataFrames, e.g.:
  - `free_energy`
  - `dispersive_surface_energy`
  - `injection_items`

- **Free energy processing**  
  From the `free_energy` table, I extract:
  - surface coverage (`n/nm`)
  - solvent name
  - polar free energy (`En. (Pol Com)`)

  These are then used to compute acid–base related quantities.

- **Surface energy components**  
  Using processed data, I compute:
  - `yd` — dispersive surface energy  
  - `yab` — acid–base surface energy  
  - `yt` — total surface energy (`yt = yd + yab`)

- **Basic plotting**  
  Simple plots of `yd`, `yab`, and `yt` vs `n/nm` are implemented to visualise how surface energy changes with coverage.

---

## Status and future work

This project is **very much a work in progress**. Planned improvements include:

- Refactoring the parser to be:
  - more robust to format changes  
  - less hard‑coded to a single file layout  
- Adding clearer validation and error messages.
- Building a more structured analysis pipeline for:
  - acid–base parameter extraction  
  - exponential distribution regression  
  - γᴰ, γᴬᴮ, γᵀ profiling
- Improving documentation and examples.

For now, the repository serves as:

- a **learning space** for scientific Python and data processing  
- a **practical tool** for analysing my own IGC‑SEA data  
- a **starting point** that I will refine as I become more competent

---

## Use of AI tools

This project was developed with the assistance of **GitHub Copilot** for code suggestions and boilerplate generation.

All domain decisions (how to interpret the IGC‑SEA export, which tables to extract, what to compute from them) were made by me. Copilot helped turn those ideas into working code more quickly, especially for the parser, which I intend to revisit and improve over time.
