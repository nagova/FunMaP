# Fig. 5 (v2) — simulation panels, FunMaP

```
python scripts/fig5_v2.py          # -> output/Fig5_simulation_v2.{svg,pdf,png}
python scripts/sim_intervals.py    # -> output/sim_switching_intervals.{csv,tex}  (SI table)
```
Requires: numpy, pandas, matplotlib (Helvetica/Arial or Liberation Sans for SVG text). `--data` / `--out` override paths.

## Inputs (data/)
| Path | Used for |
|---|---|
| R1_Benchmark/descending_branch_B_vs_Mz_cell1p0nm.csv | (a) loop, R1 3D, D = 200 nm, Δx = 1.0 nm |
| R1_Benchmark/midplane_slices_y0/*.npz | (b–d) y = 0 slices at 0, −6.20, −6.30 T |
| sim_loops/L10_pure_wedge_D*_cell1p3nm.csv | (e) pure L1₀ wedges, D = 0.2–10 µm |
| sim_loops/R1_hemi_d200nm_cell_1p0nm.csv | (e) 3D reference point |
| sim_loops/M1_wedge_D3um_cell1p3nm_fA1_*.csv | (e, f) soft-fraction sweep, D = 3 µm |
| sim_loops/M6_snapshot_rerun_D10um_cell1p3nm_15pctA1.csv | (e) 15 % soft phase at D = 10 µm |
| sim_loops/R1_hemi_*, R2b_wedge_* | SI table only (mesh check, wedge validation) |
| squid/fig4_extracted_values.csv | (e, g) measured H_c (output of fig4_squid.py, fixed-baseline pipeline) |

## Rules encoded
- Simulated H_c = switching interval, never an interpolated zero crossing.
- (f): ⟨K_u⟩/K_L1₀ = 1 − f_A1 (soft cells below exchange length → exchange-averaged K).
  XRD band: S = 0.70 ± 0.09, K_u ∝ S² → [0.37, 0.62]; H_c at S² read from a linear fit through
  the interval midpoints (printed to stdout: 2.1–3.8 T, centre 2.9 T).
- Excluded on purpose (do not add back): M2, M3, M4, M5 (300 K), IP_vs_OOP, Sharrock 1.127 T.
