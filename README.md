<h1 align="center">FunMaP: <b>Fun</b>ctional <b>Ma</b>gnetic <b>P</b>articles analysis pipeline for FePt thin films on spherical SiO₂ substrates</h1>

<p align="center">
  <a href="https://is.mpg.de/person/lsmith"><strong>Natalia Gonzalez-Vazquez</strong></a> ·
  <a href="https://hi.is.mpg.de/person/aschulz"><strong>Andrew K. Schulz</strong></a>
</p>

<p align="center">
  <a href="https://arxiv.org/abs/2504.07143">
    <img src="https://img.shields.io/badge/arXiv-Preprint-B31B1B.svg" alt="arXiv Preprint">
  </a>
  <a href="https://doi.org/10.17617/3.ROQPWZ">
    <img src="https://img.shields.io/badge/Data%20Repository-Edmond-005BBB.svg" alt="Edmond Repository">
  </a>
</p>

A thesis-linked (soon to be pre-printed) research repository for micromagnetic simulations, SQUID magnetometry analysis, and XRD visualisation of FePt thin films on spherical SiO₂ substrates.

---

## Overview

This repository contains three independent Python-based modules developed for the characterisation of **FePt thin films on spherical SiO₂ substrates**. Together they provide the full computational workflow for micromagnetic simulations, XRD visualisation, and SQUID magnetometry analysis presented in the associated thesis and manuscript.

- **Part 1:** Micromagnetic simulations of hemispherical FePt caps (ubermag/OOMMF)
- **Part 2:** Automated visualisation of X-ray diffraction (XRD) data
- **Part 3:** SQUID magnetometry analysis and SQUID–simulation comparison

Each module is self-contained but designed to function within a unified research workflow.

<!-- Preview images — uncomment once figures/ folder is populated
## Preview

### SQUID analysis
![SQUID analysis example](figures/squid_example.png)

### SQUID–simulation comparison
![SQUID–simulation overlay](figures/squid_oommf_overlay.png)

### XRD plotting
![XRD plot example](figures/xrd_example.png)
-->

---

## Repository Structure

```
FunMaP/
│
├── analysis/                              # Experimental data analysis
│   ├── SQUID_analysis_Caps.ipynb       # Batch averaging & background correction
│   ├── SQUID-OOMMF-analysis.ipynb                     # SQUID vs simulation overlay
│   ├── OOMMF-analysis.ipynb                 # Simulation hysteresis + SFD plots
│   └── XRDplot.ipynb                        # XRD masked/naked plotting
│
├── simulations/                           # Micromagnetic simulation scripts
│   ├── FePt_L10_A1_MultipleCaps_HystLoop_radial_merged.ipynb
│   └── FePt_L10_MultipleCaps_HystLoop_DiameterSweep.ipynb
│
├── sample_data/                           # Synthetic demo data (not real measurements)
│   ├── squid/
│   │   ├── sample_3um_measurement_A.dat
│   │   └── sample_3um_measurement_B.dat
│   ├── simulations/
│   │   ├── data_sphere_1.0um_synthetic.csv
│   │   ├── data_sphere_3.0um_synthetic.csv
│   │   ├── data_sphere_5.0um_synthetic.csv
│   │   ├── data_sphere_8.0um_synthetic.csv
│   │   └── data_sphere_10.0um_synthetic.csv
│   └── xrd/
│       └── sample_FePt_SiO2_on_Si.xy
│
├── environment.yml                        # Conda environment (ubermag_env)
├── .gitignore
└── README.md
```

---

## Notebook Guide

| Goal | Notebook | Main output |
|---|---|---|
| Run mixed-phase FePt cap simulations | `simulations/FePt_L10_A1_MultipleCaps_HystLoop_radial_merged.ipynb` | Hysteresis data for multiple diameters and phase fractions |
| Run pure L1₀ diameter-sweep simulations | `simulations/FePt_L10_MultipleCaps_HystLoop_DiameterSweep.ipynb` | Hysteresis data for pure L1₀ caps |
| Analyse SQUID batches | `analysis/SQUID_analysis_Caps.ipynb` | Averaged loops, statistics, corrected plots, CSV export |
| Compare SQUID and simulation results | `analysis/SQUID-OOMMF-analysis.ipynb` | Two-panel SQUID–simulation overlay |
| Analyse converted OOMMF simulation files | `analysis/OOMMF-analysis.ipynb` | Overlay plots, SFD analysis, per-file reports |
| Plot XRD files | `analysis/XRDplot.ipynb` | Publication-ready XRD PNG/SVG plots |

---

## Pipeline

```
SQUID measurements (.dat)
        │
        ▼
SQUID_analysis_Caps.ipynb
  → background correction (diamagnetic slope subtraction)
  → per-batch averaging + statistics (Hc, Mr/Ms, W_hyst)
  → Averaged_Data_{diam}um_{ts}.csv
        │
        └──────────────────────────────┐
                                       ▼
Simulations (ubermag/OOMMF)     SQUID-OOMMF-analysis.ipynb
        │                         → overlay: SQUID mean vs sims
        ▼                         → two-panel: raw A·m² / M/Msat
[Simulation scripts]
  → CONVERTED CSVs
    (B_ext, Mz/Ms, Moment_Am2)
        │
        ▼
OOMMF-analysis.ipynb
  → multi-file overlay plot
  → per-file hysteresis + SFD + report

XRD measurements (.xy)
        │
        ▼
XRDplot.ipynb
  → masked / naked mode
  → reference markers: FePt, Si, SiO₂, Fe–O phases
```

---

## Sample Data

This repository includes synthetic example files in `sample_data/` for testing the plotting and analysis workflow without requiring access to unpublished experimental data.

| Folder | Contents |
|---|---|
| `sample_data/squid/` | Demo SQUID-style `.dat` files (Quantum Design format) |
| `sample_data/simulations/` | Synthetic converted simulation CSVs |
| `sample_data/xrd/` | Example `.xy` diffractogram |

These files are intended only for workflow demonstration and code testing — they do not represent the original experimental datasets used in the thesis.

---

## Dependencies

All scripts run in a single conda environment based on [ubermag](https://ubermag.github.io/).

| Dependency | Version | Purpose |
|------------|---------|---------|
| Python | ≥ 3.9 | Core execution |
| Ubermag | Latest stable | Micromagnetic simulation interface |
| OOMMF | ≥ 1.2 | Micromagnetic solver backend |
| NumPy | ≥ 1.20 | Numerical operations |
| pandas | ≥ 1.3 | Data I/O and tabulation |
| Matplotlib | ≥ 3.5 | Visualisation |
| scipy | ≥ 1.7 | Signal processing |

```bash
git clone https://github.com/nagova/FunMaP.git
cd FunMaP
conda env create -f environment.yml
conda activate ubermag_env
jupyter notebook
```

Verify OOMMF is available:
```bash
oommf.tcl +version
```

---

## Quick Start

After installation, open Jupyter Notebook inside the `ubermag_env` environment and choose the notebook matching your task:

- **Run micromagnetic hysteresis simulations**
  - `simulations/FePt_L10_A1_MultipleCaps_HystLoop_radial_merged.ipynb`
  - `simulations/FePt_L10_MultipleCaps_HystLoop_DiameterSweep.ipynb`
- **Analyse SQUID measurements** → `analysis/SQUID_analysis_Caps.ipynb`
- **Compare SQUID data with simulations** → `analysis/SQUID-OOMMF-analysis.ipynb`
- **Plot XRD diffractograms** → `analysis/XRDplot.ipynb`
- **Analyse converted simulation files** → `analysis/OOMMF-analysis.ipynb`

Recommended order for first-time users:
1. Create and activate the conda environment
2. Verify OOMMF is available (`oommf.tcl +version`)
3. Open Jupyter Notebook
4. Start from the notebook matching your workflow goal

---

## Part 1 — Micromagnetic Simulations

<details>
<summary><strong>Click to expand</strong></summary>

### Scientific Scope

This module models the intrinsic magnetic switching behaviour of hemispherical FePt shells deposited on diamagnetic SiO₂ spheres. The SiO₂ substrate is treated as magnetically inactive (Ms = 0) and is excluded from the computational mesh.

The simulations explicitly account for:
- Curvature-dependent demagnetising effects
- Competition between exchange and magnetocrystalline anisotropy energies
- Phase composition (ordered L1₀ vs mixed L1₀/A1 FePt)
- Size-dependent discretisation for numerical convergence

### Simulation Scripts

**`FePt_L10_A1_MultipleCaps_HystLoop_radial_merged.ipynb`**  
Hysteresis loops for diameters 1, 3, 5, 8, 10, 20 µm (60 nm cap). Mixed L1₀/A1 phase with radial anisotropy. Runtime prompts select temperature (MinDriver at 0 K / TimeDriver at T > 0 K). Supports checkpoint-based resume.

**`FePt_L10_MultipleCaps_HystLoop_DiameterSweep.ipynb`**  
Same diameter sweep but pure L1₀ phase — no A1 soft fraction. Runtime prompts select temperature and anisotropy mode (Radial or Uniaxial_Vertical). Supports checkpoint-based resume.

### Key Parameters

> **Note on notation:** Throughout this repository, `L10` in file names and code refers to the **L1₀** (L1-zero) ordered intermetallic phase of FePt.

| Parameter | Value | Units |
|-----------|-------|-------|
| Ms (saturation magnetisation) | 1 × 10⁶ | A/m |
| A (exchange stiffness) | 1 × 10⁻¹¹ | J/m |
| Ku (L1₀ hard phase) | 6.6 × 10⁶ | J/m³ |
| Ku (A1 soft phase) | 1 × 10⁴ | J/m³ |
| B_max (field sweep range) | ±18 | T |
| Cap thickness | 60 | nm |
| Sphere diameters | 1, 3, 5, 8, 10, 20 | µm |

### Modelling Assumptions
- Single isolated hemispherical FePt cap; no inter-particle dipolar coupling
- Zero-temperature approximation for static driver (no thermal activation)
- Open boundary conditions; self-consistent magnetostatic field within the cap
- Homogeneous material parameters within each defined phase

</details>

---

## Part 2 — XRD Data Visualisation

<details>
<summary><strong>Click to expand</strong></summary>

### Scientific Purpose

Automated visualisation of Cu Kα XRD diffractograms of FePt films on SiO₂ sphere monolayers, enabling rapid structural phase assessment and L1₀ ordering evaluation.

### Script: `XRDplot.ipynb`

- GUI file selector → one or more `.xy` files
- Terminal prompt: **masked** (apply angular masks) or **naked** (raw data)
- log₁₀ intensity axis by default (switchable to linear)
- Colour-coded vertical reference markers with (hkl) + 2θ labels in an external legend:
  - FePt L1₀ — orange
  - Si substrate — gray
  - SiO₂ amorphous hump — blue
  - Fe–O phases (magnetite / hematite / maghemite) — green (togglable)
- Angular masks configurable via `MASKS` list at the top of the script
- **Outputs:** PNG (600 dpi) + SVG per file → `xrd_outputs/`

</details>

---

## Part 3 — SQUID Magnetometry Analysis

<details>
<summary><strong>Click to expand</strong></summary>

### Scientific Purpose

This module bridges experimental SQUID magnetometry data and micromagnetic simulation results, enabling direct comparison between ensemble-averaged experimental hysteresis loops and single-cap simulated loops.

### Script: `SQUID_analysis_Caps.ipynb`

Processes raw SQUID `.dat` files (Quantum Design format) for a single particle batch.

- GUI file selector → ≥2 `.dat` files per batch
- Prompts: sphere diameter [µm], cap thickness [nm], substrate area [mm²]
- Optional diamagnetic background subtraction (slope fitted to ±95% saturation tails independently, then averaged)
- Averaging by branch (descending and ascending separately) → proper closed mean loop
- Extracts: Hc [T], Mr/Ms, hysteresis loss W_hyst [J/kg]
- **Outputs:** PNG + SVG two-panel plot, `.txt` statistical report, averaged CSV

**Plot layout:**
- Top panel: raw emu — individual curves + mean ± 1 SD
- Bottom panel: background-corrected M/Msat — individual curves + mean ± 1 SD

**Averaged CSV columns:** `Field_T`, `M_norm_desc`, `M_norm_asc`, `M_emu_desc`, `M_emu_asc`

---

### Script: `SQUID-OOMMF-analysis.ipynb`

Two-panel overlay of averaged SQUID data and simulation results.

- Step 1: select averaged SQUID CSV (output of `SQUID_analysis_Caps`)
- Step 2: select one or more CONVERTED simulation CSVs
- **Outputs:** PNG + SVG → `Comparison_Results_{ts}/`

**Plot layout:**
- Top panel: physical moment [A·m²] — SQUID (emu → A·m² converted) + sim `Moment_Am2`
- Bottom panel: normalised M/Msat — SQUID `M_norm` + sim `Mz/Ms`

---

### Script: `OOMMF-analysis.ipynb`

Multi-file simulation overlay and per-file SFD analysis.

- GUI file selector → any number of CONVERTED simulation CSVs
- Terminal prompt: sphere diameter + cap thickness per file
- **Outputs per session:**
  - `Overlay_{ts}.png/.svg` — all files on one two-panel figure
  - Per-file `{name}_{ts}.png/.svg` + `_report.txt` (Hc, Mr/Ms, SFD FWHM)

</details>

---

## Citation

If you use this repository, please cite:

```bibtex
@misc{gonzalezvazquez_funmap_2026,
  title  = {FunMaP: Customizable simulations for Janus particles' magnetic properties
            with associated visualizations},
  author = {Gonzalez-Vazquez, Natalia and Schulz, Andrew K.},
  year   = {2026},
  note   = {In preparation},
}
```

For thesis citation or exact reproducibility, users are encouraged to reference a tagged repository release when available.

---

## Known Limitations

This repository was developed as a research workflow accompanying a thesis and associated manuscript. The notebooks are intended for transparent scientific analysis rather than as a general-purpose software package.

- Micromagnetic models treat the SiO₂ substrate as magnetically inactive (Ms = 0)
- Simulations represent isolated hemispherical caps — no inter-particle dipolar coupling
- Static simulations use a zero-temperature approximation unless thermal drivers are explicitly selected
- Material properties are homogeneous within each defined phase region
- Sample files included in the repository are synthetic demonstration data only

These assumptions should be considered when comparing simulation outputs with ensemble-averaged experimental measurements.

---

## License

This project is licensed under the **GNU General Public License v3.0**.  
See the `LICENSE` file for details.

---

## Acknowledgements

We acknowledge the International Max Planck Research School for Intelligent Systems (IMPRS-IS) for support.  
We thank J. Burns and J.-C. Passy for assistance in repository preparation.  
We thank N. Rokhmanova for the ARIADNE repository inspiration.  
We thank G. Richter for scientific feedback and discussion.

---

## Contact

Maintained by:
- **Natalia Gonzalez-Vazquez** — https://github.com/nagova
- **Andrew K. Schulz** — https://github.com/Aschulz94

If you find this repository useful, consider giving it a ⭐
