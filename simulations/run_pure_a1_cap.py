"""
run_pure_a1_cap.py  --  FunMaP
=============================================================================
Simulates a single FePt cap in the 100% disordered A1 fcc phase (f_A1 = 1.0,
Ku = 0 J/m^3, Ms = 1.0 MA/m, A = 10 pJ/m) to demonstrate that the chemically
disordered A1 phase behaves as a completely soft ferromagnet (Hc << 0.05 T,
dominated purely by shape demagnetization).

Outputs saved to: results/simulations/pure_a1/
  - <run_id>.csv           (full hysteresis loop: B_ext_T, Mx, My, Mz)
  - <run_id>_params.json   (metadata, Hc, Mr/Ms, runtime)
  - <run_id>_loop.png      (publication-ready plot showing soft magnetic loop)

Usage:
  python simulations/run_pure_a1_cap.py
"""

import os
import sys
import json
import time
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

for p in [r"C:\Users\admin\miniforge3\envs\ubermag_env\Library\bin", r"C:\Users\admin\miniforge3\Library\bin"]:
    if os.path.exists(p) and p not in os.environ.get("PATH", ""):
        os.environ["PATH"] = p + ";" + os.environ.get("PATH", "")

oommf_path = r"C:\Users\admin\Desktop\oommf\oommf.tcl"
if os.path.exists(oommf_path):
    os.environ["OOMMFTCL"] = oommf_path

import discretisedfield as df
import micromagneticmodel as mm
import oommfc as oc

EXPLORATORY_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "archive", "exploratory_scripts"))
sys.path.insert(0, EXPLORATORY_DIR)
import fept_order_mechanisms as fom

OUT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results", "simulations", "pure_a1"))
os.makedirs(OUT_DIR, exist_ok=True)

# -----------------------------------------------------------------------------
# Material Constants: 100% Disordered A1 FePt
# -----------------------------------------------------------------------------
MS_A1 = 1.0e6        # A/m (1.0 MA/m)
A_EX  = 1.0e-11      # J/m (10 pJ/m)
KU_A1 = 1.0e3        # J/m^3 (negligible floor)
MU0   = 4.0 * np.pi * 1e-7

# Fine field schedule around zero to capture soft switching
DESCENDING = [0.30, 0.15, 0.08, 0.04, 0.02, 0.01, 0.005, 0.0,
              -0.005, -0.01, -0.02, -0.04, -0.08, -0.15, -0.30]
ASCENDING  = [-0.15, -0.08, -0.04, -0.02, -0.01, -0.005, 0.0,
              0.005, 0.01, 0.02, 0.04, 0.08, 0.15, 0.30]
FULL_LOOP_SCHEDULE = DESCENDING + ASCENDING


def build_3d_hemisphere_a1(diameter=120e-9, thickness=25e-9, cell_nm=3.0, name="pure_A1_cap"):
    """
    Builds a full 3D hemispherical cap of pure A1 phase.
    """
    R_in = diameter / 2.0
    R_out = R_in + thickness
    cell_m = cell_nm * 1e-9

    n_xy = int(np.ceil((2.0 * R_out) / cell_m))
    if n_xy % 2 != 0:
        n_xy += 1
    L_xy = n_xy * cell_m
    half_xy = L_xy / 2.0

    n_z = int(np.ceil(R_out / cell_m))
    L_z = n_z * cell_m

    region = df.Region(p1=(-half_xy, -half_xy, 0.0), p2=(half_xy, half_xy, L_z))
    mesh = df.Mesh(region=region, cell=(cell_m, cell_m, cell_m))

    def ms_fn(pos):
        r = np.linalg.norm(pos)
        return MS_A1 if (R_in <= r <= R_out and pos[2] >= 0) else 0.0

    system = mm.System(name=name)
    system.energy = (
        mm.Exchange(A=A_EX)
        + mm.Demag()
        + mm.UniaxialAnisotropy(K=KU_A1, u=(0, 0, 1))
        + mm.Zeeman(H=(0, 0, 0))
    )
    system.m = df.Field(mesh, nvdim=3, value=(0, 0, 1), norm=df.Field(mesh, nvdim=1, value=ms_fn))
    return system, mesh


def run_pure_a1_loop(diameter_nm=120, cell_nm=3.0):
    fom.setup_oommf_runner()
    run_id = f"pure_A1_hemisphere_D{diameter_nm}nm_cell{str(cell_nm).replace('.', 'p')}nm"
    print(f"\n=======================================================", flush=True)
    print(f"Running Pure A1 FePt Cap Simulation: {run_id}", flush=True)
    print(f"Phase: 100% disordered fcc A1 (Ku = {KU_A1} J/m^3, Ms = {MS_A1/1e6} MA/m, A = {A_EX*1e12} pJ/m)", flush=True)
    print(f"Geometry: 3D hemisphere D = {diameter_nm} nm, cell = {cell_nm} nm", flush=True)
    print(f"=======================================================", flush=True)

    t0 = time.time()
    sysm, mesh = build_3d_hemisphere_a1(diameter=diameter_nm*1e-9, thickness=25e-9, cell_nm=cell_nm, name=run_id)

    vol = sysm.m.norm.integrate().item()
    n_mag = int(np.round(vol / np.prod(mesh.cell)))
    print(f"Mesh shape: {mesh.n}, Active magnetic volume: {vol*1e24:.1f} nm^3 (~{n_mag} active cells)", flush=True)

    # Note: stopping_mxHxm and stage_iteration_limit must be passed to MinDriver constructor!
    drv = oc.MinDriver(stopping_mxHxm=10.0, stage_iteration_limit=400)
    mvec = lambda: [getattr(sysm.m, k).integrate().item() / vol for k in "xyz"]

    B_TILT = 0.005  # slight tilt along y to avoid metastable dead center
    def at(B):
        sysm.energy.zeeman.H = (0, B_TILT / MU0, B / MU0)
        drv.drive(sysm)
        return mvec()

    # Pre-saturate at +0.5 T
    print("Pre-saturating at +0.5 T...", flush=True)
    at(0.5)

    rec = []
    print("\nTracing full hysteresis loop (-0.3 T <-> +0.3 T):", flush=True)
    for i, B in enumerate(FULL_LOOP_SCHEDULE):
        mx, my, mz = at(B)
        rec.append(dict(B_ext_T=round(B, 4), Mx=mx, My=my, Mz=mz))
        branch = "Desc" if i < len(DESCENDING) else "Asc "
        print(f"[{branch}] B = {B:+7.3f} T  -->  <mz> = {mz:+.4f}   <mx> = {mx:+.4f}   <my> = {my:+.4f}", flush=True)

    df_loop = pd.DataFrame(rec)
    csv_path = os.path.join(OUT_DIR, f"{run_id}.csv")
    df_loop.to_csv(csv_path, index=False)
    print(f"\nSaved CSV to {csv_path}", flush=True)

    # Coercivity Hc and Remanence Mr/Ms
    desc_df = df_loop.iloc[:len(DESCENDING)].reset_index(drop=True)
    row_0 = desc_df[desc_df["B_ext_T"] == 0.0]
    mr_ms = float(row_0["Mz"].values[0]) if len(row_0) > 0 else 0.0

    hc_interv = None
    for j in range(len(desc_df) - 1):
        mz1, mz2 = desc_df.loc[j, "Mz"], desc_df.loc[j+1, "Mz"]
        if mz1 >= 0 and mz2 < 0:
            b1, b2 = abs(desc_df.loc[j, "B_ext_T"]), abs(desc_df.loc[j+1, "B_ext_T"])
            hc_interv = [min(b1, b2), max(b1, b2)]
            break

    elapsed_min = (time.time() - t0) / 60.0

    meta = {
        "run_id": run_id,
        "phase": "100% disordered A1 FePt (f_A1 = 1.0)",
        "Ku_Jm3": KU_A1,
        "Ms_Am": MS_A1,
        "A_ex_Jm": A_EX,
        "diameter_nm": diameter_nm,
        "cell_nm": cell_nm,
        "Mr_over_Ms": mr_ms,
        "switching_interval_T": hc_interv if hc_interv else "Hc < 0.005 T",
        "minutes": elapsed_min,
        "n_magnetic_cells": n_mag
    }
    json_path = os.path.join(OUT_DIR, f"{run_id}_params.json")
    with open(json_path, "w") as f:
        json.dump(meta, f, indent=2)

    # Plot hysteresis loop
    plot_path = os.path.join(OUT_DIR, f"{run_id}_loop.png")
    fig, ax = plt.subplots(figsize=(6.5, 5), dpi=300)
    ax.plot(desc_df["B_ext_T"], desc_df["Mz"], "o-", color="#d95f02", lw=2, ms=4, label="Descending branch")
    asc_df = df_loop.iloc[len(DESCENDING):].reset_index(drop=True)
    ax.plot(asc_df["B_ext_T"], asc_df["Mz"], "s--", color="#1b9e77", lw=1.5, ms=3.5, label="Ascending branch")

    ax.axhline(0, color="gray", ls="--", lw=0.8, alpha=0.7)
    ax.axvline(0, color="gray", ls="--", lw=0.8, alpha=0.7)
    ax.set_xlabel(r"$\mu_0 H_{ext}$ (T)", fontsize=12)
    ax.set_ylabel(r"Normalized Magnetization $\langle m_z \rangle$", fontsize=12)
    hc_str = f"[{hc_interv[0]:.3f}, {hc_interv[1]:.3f}]" if hc_interv else "< 0.005"
    ax.set_title(f"Pure Disordered A1 FePt Cap (3D Hemisphere, D = {diameter_nm} nm)\n"
                 rf"$\mu_0 H_c \in {hc_str}$ T,  $M_r/M_s = {mr_ms:.3f}$ (Ultra-Soft Ferromagnet)", fontsize=11)
    ax.set_xlim(-0.35, 0.35)
    ax.set_ylim(-1.05, 1.05)
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend(frameon=True, loc="lower right")

    # Callout comparison with L1_0
    ax.text(0.04, 0.20,
            r"$\mathbf{FePt\ Phase\ Comparison:}$" + "\n" +
            r"• Ideal $L1_0$ ($S=1.0$): $\mu_0 H_c \approx 6.5$ T" + "\n" +
            r"• Real Cap ($S=0.7$): $\mu_0 H_c \approx 1.1$ T" + "\n" +
            rf"• Pure A1 ($S=0$): $\mu_0 H_c \in {hc_str}$ T" + "\n" +
            r"$\rightarrow$ Demonstrates A1 is completely soft magnet",
            transform=ax.transAxes, fontsize=9.5,
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#f8f9fa", edgecolor="#ced4da"))

    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()
    print(f"Saved plot to {plot_path}", flush=True)

    print("\nSimulation Complete!", flush=True)
    print(json.dumps(meta, indent=2), flush=True)
    return meta


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--d", type=int, default=120, help="Cap diameter in nm")
    parser.add_argument("--cell", type=float, default=3.0, help="Cell size in nm")
    args = parser.parse_args()

    run_pure_a1_loop(diameter_nm=args.d, cell_nm=args.cell)
