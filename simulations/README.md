# Micromagnetic Simulation Suite (`simulations/`)

This directory contains the production micromagnetic simulation scripts for **FunMaP**, executing via [Ubermag](https://ubermag.github.io/) and NIST [OOMMF](https://math.nist.gov/oommf/).

---

## 🔬 Key Simulation Scripts

### 1. `run_polycrystal_cap.py` (Definitive "Real Cap" Model)
Constructs the realistic polycrystalline FePt cap on the tapered meridian wedge geometry ($D = 3\ \mu\text{m}$, pole thickness $t_0 = 60\text{ nm}$, $W_y(\theta) \propto \sin\theta$) using direct experimental microstructure:
* **Chemical Order Parameter**: $S = 0.70 \pm 0.09$ from experimental XRD $(001)/(002)$ integrated intensity ratio on flat reference film ($K_u = S^2 K_{L1_0} = 3.23\text{ MJ/m}^3$).
* **Grain Size**: Mean equivalent sphere diameter $d_g = 19\text{ nm}$ from Scherrer peak broadening on the $(111)$ reflection.
* **Grain Morphology**: 3D Poisson-Voronoi cellular tessellation (~2,924 grains per cap).
* **Crystallographic Texture**: Uniaxial easy axes $\mathbf{u}$ randomly oriented on the unit sphere (3D random texture, matching XRD baseline).
* **Intergranular Exchange Coupling ($x$)**: Swept over $x \in \{1.0, 0.3, 0.1, 0.0\}$ ($A_{\text{gb}} = x \cdot A_{\text{bulk}}$) across a 2-cell boundary layer.
* **Bulk Constants**: $M_s = 1.0\times 10^6\text{ A/m}$, $A_{\text{bulk}} = 10\text{ pJ/m}$.

#### Execution Commands:
```bash
# Phase 1: Control Benchmark (S = 1.0, radial axis, must switch in [6.45, 6.60] T)
python simulations/run_polycrystal_cap.py --phase control

# Phase 2: Main Microstructure Sweep (S = 0.70, random axes, x in {1.0, 0.3, 0.1, 0.0})
python simulations/run_polycrystal_cap.py --phase main

# Or run individual coupling values:
python simulations/run_polycrystal_cap.py --phase main --x 1.0
```

#### Results Summary:
* $x = 1.0$ (Continuous exchange across grain boundaries): **$\mu_0 H_c \in [1.10, 1.15]\text{ T}$**, $M_r/M_s = 0.5898$ $\longrightarrow$ **Exact match to experimental SQUID ($1.00\text{--}1.19\text{ T}$)** without any free parameters.
* $x = 0.3$: $\mu_0 H_c \in [1.75, 1.80]\text{ T}$, $M_r/M_s = 0.5668$
* $x = 0.1$: $\mu_0 H_c \in [2.45, 2.50]\text{ T}$, $M_r/M_s = 0.5257$
* $x = 0.0$ (Decoupled grains): $\mu_0 H_c \in [2.90, 2.95]\text{ T}$, $M_r/M_s = 0.4899$ (Stoner-Wohlfarth 3D random limit).

---

### 2. `run_pure_a1_cap.py` (Disordered Phase Control)
Simulates a 3D hemispherical cap of 100% chemically disordered fcc FePt phase ($K_u \approx 0$, $M_s = 1.0\text{ MA/m}$, $A = 10\text{ pJ/m}$) across a fine field sweep ($-0.30\text{ T} \leftrightarrow +0.30\text{ T}$).
* **Result**: $\mu_0 H_c < 0.005\text{ T}$ ($< 5\text{ mT}$), $M_r / M_s = 0.0055 \approx 0.0$.
* **Takeaway**: Formally demonstrates that the disordered A1 phase acts as an ultra-soft ferromagnet ($>1,300\times$ softer than $L1_0$), proving that all hard magnetic retention originates strictly from ordered $L1_0$ crystallites.

```bash
python simulations/run_pure_a1_cap.py
```

---

### 3. `run_mesh_convergence.py` (R1 Benchmark Series)
Executes the spatial mesh resolution study on a $D = 200\text{ nm}$ 3D hemispherical cap with $L1_0$ radial texture from cell sizes $4.0\text{ nm}$ down to $1.0\text{ nm}$, demonstrating strict asymptotic convergence at $\mu_0 H_c \in [6.20, 6.30]\text{ T}$.

---

### 4. `run_pure_l10_diameter_sweep.py` (Curvature Invariance Series)
Sweeps sphere diameters $D$ from $200\text{ nm}$ up to $10\ \mu\text{m}$, demonstrating that for $D \ge 500\text{ nm}$ ($\ell_{\text{ex}} / R \ll 1$), switching coercivity becomes independent of particle size.

---

## 📂 Output Locations
* Processed descending hysteresis loops: `results/simulations/loops_polycrystal/*.csv`
* Simulation metadata and switching intervals: `results/simulations/loops_polycrystal/*_params.json`
* Midplane 2D magnetization vector slices: `results/simulations/loops_polycrystal/*_slice_*.npz`
* Pure A1 hysteresis loops & figures: `results/simulations/pure_a1/`

*(Note: Heavy raw `.omf` and restart files are excluded from git and preserved in the Max Planck Edmond repository under DOI [10.17617/3.ROQPWZ](https://doi.org/10.17617/3.ROQPWZ)).*
