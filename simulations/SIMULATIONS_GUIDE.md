# FunMaP Simulation Guide

This guide is for the micromagnetic simulation notebooks only. It explains what each notebook does, when to use it, what parameters matter, and how to interpret the generated outputs.

The original production notebooks are kept intact. If you need extra figures, spatial snapshots, or exploratory visualizations, use a separate notebook rather than modifying the working simulation notebooks.

---

## What These Simulations Do

The simulation notebooks model FePt hemispherical caps on spherical SiO2 particles using Ubermag/OOMMF. The SiO2 sphere is treated as magnetically inactive; only the FePt cap is included in the magnetic mesh.

The basic workflow is:

```text
Choose cap geometry -> Build magnetic mesh -> Sweep external field -> Save Mz/Ms loop -> Analyze or visualize outputs
```

The main simulation output is a hysteresis loop:

```text
B_ext (T), Mz/Ms
```

where `Mz/Ms` is the average z-component of the magnetization normalized over the magnetic FePt material.

---

## Simulation Notebooks

| File | Purpose | Use when |
|------|---------|----------|
| `FePt_L10_MultipleCaps_HystLoop_DiameterSweep.ipynb` | Pure L10 FePt cap simulations across multiple sphere diameters | You want baseline hard-phase FePt behavior |
| `FePt_L10_A1_MultipleCaps_HystLoop_radial_merged.ipynb` | Mixed L10/A1 FePt cap simulations across multiple sphere diameters | You want to study soft-phase A1 contributions |
| `FePt_real_magnetization_snapshots.ipynb` | Real spatial magnetization snapshots from `system.m` along one selected loop | You want XZ/XY state maps, selected-state exports, or an interactive loop viewer |

---

## Recommended Starting Point

For standard hysteresis-loop generation, start with:

```text
FePt_L10_A1_MultipleCaps_HystLoop_radial_merged.ipynb
```

This notebook is closest to the mixed-phase Janus-particle picture used in the manuscript workflow. It simulates FePt caps with a hard L10 phase and a soft A1 phase.

For clean reference behavior, use:

```text
FePt_L10_MultipleCaps_HystLoop_DiameterSweep.ipynb
```

This notebook assigns all magnetic cells the L10 anisotropy constant.

For figure-making or presentation snapshots, use:

```text
FePt_real_magnetization_snapshots.ipynb
```

This notebook reruns one selected cap with a lighter field schedule and saves actual spatial magnetization states from `system.m`. It can also render selected states and build an interactive Plotly/HTML viewer from a saved result folder.

---

## Physical Model

The FePt cap is represented as a hemispherical shell:

```text
R_inner = sphere radius
R_outer = sphere radius + cap thickness
```

Only cells satisfying the shell condition are magnetic:

```text
R_inner <= r <= R_outer and z >= 0
```

Cells outside the cap have:

```text
Ms = 0
```

The energy terms are:

| Term | Meaning |
|------|---------|
| Exchange | Penalizes rapid spatial variation of magnetization |
| Uniaxial anisotropy | Encodes L10 or A1 easy-axis behavior |
| Demagnetization | Includes magnetostatic self-interaction |
| Zeeman | Applies the external magnetic field sweep |

At `T = 0 K`, the notebooks use `MinDriver`. At `T > 0 K`, they use `TimeDriver` with LLG dynamics and a stochastic thermal field.

---

## Key Parameters

| Parameter | Typical value | Meaning |
|-----------|---------------|---------|
| `diameters` | `[1, 3, 5, 8, 10, 20] um` | Sphere diameters to simulate |
| `cap_thickness` | `60 nm` | FePt film/cap thickness |
| `Ms` | `1.0e6 A/m` | Saturation magnetization |
| `A` | `1.0e-11 J/m` | Exchange stiffness |
| `Ku_hard` | `6.6e6 J/m^3` | L10 anisotropy constant |
| `Ku_soft` | `1.0e4 J/m^3` | A1 soft-phase anisotropy constant |
| `SOFT_FRACTION_A1` | commonly `0.15` or `0.50` | Fraction of magnetic cells assigned A1 behavior |
| `B_max` | `18 T` | Maximum applied field magnitude |
| `B_tilt` | `0.01 T` | Small transverse field used to avoid perfectly symmetric saddle states |

The field sweep is:

```text
+B_max -> -B_max -> +B_max
```

with 121 points per branch by default in the production sweep notebooks. The snapshot notebook uses 41 points per branch by default, producing 82 saved states for lighter visualization datasets.

---

## Anisotropy Modes

Two anisotropy modes are supported in the simulation logic.

| Mode | Description |
|------|-------------|
| `Radial` | The easy axis follows the local outward surface normal of the spherical cap |
| `Uniaxial_Vertical` | The easy axis is fixed along the global z direction |

`Radial` is the default for the mixed L10/A1 notebook and is usually the best starting point for curved FePt caps.

---

## Mesh Resolution

The notebooks use an adaptive mesh strategy:

```text
d <= 1.5 um: approximately 18 nm cells
d > 1.5 um: fixed n = 110 cells per axis
```

This keeps large-diameter simulations computationally manageable. For very large particles, the cell size can exceed the FePt exchange length, so the results should be interpreted as qualitative or comparative rather than atomically resolved micromagnetic detail.

Always check the printed volume error. A smaller volume error means the discretized shell better matches the analytical cap volume.

---

## Standard Outputs

The production notebooks save timestamped output folders containing files such as:

| Output | Description |
|--------|-------------|
| `data_*.csv` | Per-cap hysteresis data |
| `MASTER_DATA_*.csv` | Combined hysteresis data for all completed runs |
| `SUMMARY_PLOT_*.svg` | Summary plot of simulated hysteresis loops |

The CSV files are the inputs expected by the analysis notebooks in `analysis/`.

---

## Real Magnetization Snapshots

The standard loop CSV contains only averaged magnetization. It cannot reconstruct the spatial magnetization pattern inside the cap.

Use:

```text
FePt_real_magnetization_snapshots.ipynb
```

when you need real magnetization-state images along the loop.

This notebook saves:

| Output | Description |
|--------|-------------|
| `*_hysteresis.csv` | Loop for the selected cap |
| `*_snapshot_manifest.csv` | Table linking each saved state to its field, averaged magnetization, and NPZ file |
| `state_*.npz` | Saved XZ and XY `m_z/M_s` projection arrays for each field step |
| selected-state `.png/.svg` | Publication/presentation exports generated from one chosen state file |
| `interactive_snapshot_viewer.html` | Interactive Plotly viewer generated from a folder of saved states |

By default, the snapshot notebook saves 82 states:

```text
41 descending-branch states
41 ascending-branch states
```

The snapshot maps are generated from `system.m`, so they are real micromagnetic states from the simulation rather than stylized reconstructions from the averaged loop. Empty/white pixels in the XZ or XY maps represent projected lines of sight with no FePt magnetic material.

---

## Suggested Workflow

1. Run one production simulation notebook to generate the main hysteresis loop data.
2. Analyze the resulting CSV files with the notebooks in `analysis/`.
3. If you need spatial state images, run `FePt_real_magnetization_snapshots.ipynb` for one representative cap.
4. Use the final selected-state cell to export a specific state as PNG/SVG.
5. Use the folder-based interactive cell to build an HTML loop viewer from saved `state_*.npz` files.
6. Keep exploratory visualization notebooks separate from the validated production notebooks.

---

## Interpreting the Loop

Useful landmarks along the loop:

| Region | Meaning |
|--------|---------|
| Positive saturation | Most magnetic moments aligned with +z |
| Descending remanence | Field is near zero after coming down from positive saturation |
| Descending coercivity | Magnetization crosses near zero while sweeping toward negative field |
| Negative saturation | Most magnetic moments aligned with -z |
| Ascending remanence | Field is near zero after coming up from negative saturation |
| Ascending coercivity | Magnetization crosses near zero while sweeping back toward positive field |

For mixed L10/A1 caps, the soft A1 regions may switch earlier than the hard L10 regions. This can create intermediate spatial states even when the averaged `Mz/Ms` value looks simple.

---

## Practical Notes

- OOMMF simulations can take a long time, especially for large caps and thermal runs.
- Run a single representative case before launching a full diameter sweep.
- Keep the timestamped output folders; they preserve the parameters and outputs from each run.
- Do not rely on the loop CSV alone for spatial interpretation.
- Use the `.npz` snapshot files if you want to redesign figures later without rerunning the simulation.

---

## Troubleshooting

**OOMMF is not found**

Make sure the Ubermag/OOMMF environment is installed and active before opening Jupyter.

**The notebook runs but no GUI appears**

These simulation notebooks are not GUI tools. Run the cells directly in Jupyter Lab or Jupyter Notebook.

**The simulation is very slow**

Start with a smaller diameter or a single selected cap. Large caps use many mesh cells and can be computationally expensive.

**The loop has missing or repeated-looking values**

The production notebooks carry forward the last known value if OOMMF crashes at a field step. Check the terminal/notebook output for warnings.

**The interactive viewer has an empty XY panel**

Older snapshot result folders may only contain XZ arrays. Re-run the current `FePt_real_magnetization_snapshots.ipynb` to save both XZ and XY projections.

**The HTML viewer looks blank during animation**

Re-export the HTML with the current notebook version. The viewer should update full heatmap frames for every state.

---

## Dependency Reminder

The simulation environment needs:

```text
ubermag
discretisedfield
micromagneticmodel
oommfc
numpy
pandas
matplotlib
```

See the repository-level `environment.yml` for the recommended environment.
