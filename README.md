<h1 align="center">FunMaP: Functional Magnetic Particles analysis pipeline for FePt thin films on spherical SiO₂ substrates</h1>

<p align="center">
  <a href="https://is.mpg.de/person/lsmith"><strong>Natalia Gonzalez-Vazquez*</strong></a> ·
  <a href="https://hi.is.mpg.de/person/aschulz"><strong>Andrew K. Schulz*</strong></a>
</p>

<p align="center"><strong>*</strong> denotes equal contributions to this repository.</p>

<p align="center">
  <a href="https://arxiv.org/abs/2504.07143">
    <img src="https://img.shields.io/badge/arXiv-Preprint-B31B1B.svg" alt="arXiv Preprint">
  </a>
  <a href="https://doi.org/10.17617/3.ROQPWZ">
    <img src="https://img.shields.io/badge/Data%20Repository-Edmond-005BBB.svg" alt="Edmond Repository">
  </a>
</p>

---

## Overview

This repository contains three independent Python-based modules developed for the characterisation of **FePt thin films on spherical SiO₂ substrates**. Together they provide the full computational workflow for micromagnetic simulations, XRD visualisation, and SQUID magnetometry analysis presented in the associated thesis and manuscript.

- **Part 1:** Micromagnetic simulations of hemispherical FePt caps (ubermag/OOMMF)
- **Part 2:** Automated visualisation of X-ray diffraction (XRD) data
- **Part 3:** SQUID magnetometry analysis and SQUID–simulation comparison

Each module is self-contained but designed to function within a unified research workflow.

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

## Installation

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

## Part 1 — Micromagnetic Simulations

### Scientific Scope

This module models the intrinsic magnetic switching behaviour of hemispherical FePt shells deposited on diamagnetic SiO₂ spheres. The SiO₂ substrate is treated as magnetically inactive (Ms = 0) and is excluded from the computational mesh.

The simulations explicitly account for:
- Curvature-dependent demagnetising effects
- Competition between exchange and magnetocrystalline anisotropy energies
- Phase composition (ordered L1₀ vs mixed L1₀/A1 FePt)
- Size-dependent discretisation for numerical convergence

### Simulation Scripts

**`FePt_L10_A1_MultipleCaps_HystLoop_radial_merged.ipynb`**  
Hysteresis loops for multiple sphere diameters. A temperature prompt at runtime selects:
- `T = 0` → static energy minimisation (MinDriver, 0 K)
- `T > 0` → dynamic LLG integration with thermal field (TimeDriver, T K)

Mixed L1₀/A1 phase with radial anisotropy. Supports checkpoint-based resume.

Cap thickness sweep (10–80 nm) on a 1 µm sphere. Pure L1₀ phase, both radial and uniaxial-vertical anisotropy modes, volume-corrected output.

### Key Parameters

| Parameter | Value |
|-----------|-------|
| Ms | 1 × 10⁶ A/m |
| A (exchange) | 1 × 10⁻¹¹ J/m |
| Ku (L1₀) | 6.6 × 10⁶ J/m³ |
| Ku (A1) | 1 × 10⁴ J/m³ |
| B_max | ±18 T |
| Cap thickness (default) | 60 nm |
| Sphere diameters | 1, 3, 5, 8, 10 µm |

### Modelling Assumptions
- Single isolated hemispherical FePt cap; no inter-particle dipolar coupling
- Zero-temperature approximation for static driver (no thermal activation)
- Open boundary conditions; self-consistent magnetostatic field within the cap
- Homogeneous material parameters within each defined phase

---

## Part 2 — XRD Data Visualisation

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

---

## Part 3 — SQUID Magnetometry Analysis

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
