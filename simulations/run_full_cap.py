"""
run_full_cap.py
Full Hemispherical Cap Simulations (Corrected Geometry).

Corrects the previous octant/quarter-wedge geometry:
- Full hemisphere: box 420 x 420 x 210 nm, z >= 0 (no artificial symmetry cuts)
- Mask: Rin <= r <= Rout and z >= 0, with Rin = 150 nm, Rout = 210 nm (t = 60 nm)
- Resolutions:
    * cell = 2.0 nm (baseline mesh)
    * cell = 1.2 nm (strictly resolves L_anis = 1.231 nm, ~3 GB RAM)
- Evaluates pure L1_0 (soft_fraction = 0) and mixed phase (soft_fraction = 0.26)
- Measures:
    * Coercivity Hc
    * Remanence Mr/Ms (testing for Mr/Ms > 0.50 and approaching experimental ~0.67)
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Ensure tclsh is in PATH for Windows OOMMF runner
for p in [r"C:\Users\admin\miniforge3\Library\bin", r"C:\Users\admin\miniforge3\envs\ubermag_env\Library\bin"]:
    if os.path.exists(p) and p not in os.environ.get("PATH", ""):
        os.environ["PATH"] = p + ";" + os.environ.get("PATH", "")

import discretisedfield as df
import micromagneticmodel as mm
import oommfc as oc

sys.path.append(os.path.dirname(__file__))
import fept_micromagnetics as fm

MU0   = 4e-7 * np.pi
MS    = 1.0e6      # A/m
A_EX  = 1.0e-11    # J/m
KU_L10 = 6.6e6     # J/m^3
KU_A1  = 1.0e4     # J/m^3


def build_full_hemisphere_cap(
    name,
    diameter=300e-9,
    cell_size=1.2e-9,
    soft_fraction=0.0,
    seed=42
):
    """
    Builds the FULL hemispherical cap (no artificial quarter/octant cut faces).
    Box: 420 x 420 x 210 nm.
    """
    np.random.seed(seed)
    R_inner = diameter / 2.0
    R_outer = R_inner + fm.T_CAP

    # Snap box dimensions to integer multiples of cell_size
    n_xy = int(np.ceil((2.0 * R_outer) / cell_size))
    if n_xy % 2 != 0:
        n_xy += 1
    L_xy = n_xy * cell_size
    half_xy = L_xy / 2.0

    n_z = int(np.ceil(R_outer / cell_size))
    L_z = n_z * cell_size

    p1 = (-half_xy, -half_xy, 0.0)
    p2 = (half_xy, half_xy, L_z)
    region = df.Region(p1=p1, p2=p2)
    cell = (cell_size, cell_size, cell_size)
    mesh = df.Mesh(region=region, cell=cell)

    def ms_mask(pos):
        r = np.linalg.norm(pos)
        return MS if (R_inner <= r <= R_outer and pos[2] >= 0) else 0.0

    def ku_distribution(pos):
        r = np.linalg.norm(pos)
        if R_inner <= r <= R_outer and pos[2] >= 0:
            if soft_fraction > 0.0:
                return KU_A1 if np.random.random() < soft_fraction else KU_L10
            return KU_L10
        return 0.0

    def u_radial(pos):
        r = np.linalg.norm(pos)
        return (pos[0] / r, pos[1] / r, pos[2] / r) if r != 0 else (0, 0, 1)

    system = mm.System(name=name)
    system.energy = (
        mm.Exchange(A=A_EX)
        + mm.Demag()
        + mm.UniaxialAnisotropy(
            K=df.Field(mesh, nvdim=1, value=ku_distribution),
            u=df.Field(mesh, nvdim=3, value=u_radial)
        )
        + mm.Zeeman(H=(0, 0, 0))
    )

    system.m = df.Field(mesh, nvdim=3, value=(0, 0, 1),
                        norm=df.Field(mesh, nvdim=1, value=ms_mask))

    return system, mesh


def run_full_cap_study():
    fm.setup_oommf_runner()
    outdir = os.path.join("results", "full_hemisphere_cap")
    os.makedirs(outdir, exist_ok=True)
    
    runs = [
        # (name, cell_nm, soft_frac)
        ("full_cap_300nm_cell_2p0nm_pureL10", 2.0, 0.00),
        ("full_cap_300nm_cell_1p2nm_pureL10", 1.2, 0.00),
        ("full_cap_300nm_cell_1p2nm_26pctA1", 1.2, 0.26),
    ]
    
    results = []
    print("=" * 75)
    print("FULL HEMISPHERICAL CAP STUDY (d = 300 nm, t = 60 nm)")
    print("=" * 75)
    
    for run_name, cell_nm, soft_frac in runs:
        print(f"\n--- Running: {run_name} (cell = {cell_nm} nm, soft = {soft_frac*100:.0f}%) ---")
        t0 = time.time()
        system, mesh = build_full_hemisphere_cap(
            name=run_name,
            diameter=300e-9,
            cell_size=cell_nm * 1e-9,
            soft_fraction=soft_frac
        )
        
        driver = oc.MinDriver()
        # 1. Saturation at +16 T
        system.energy.zeeman.H = (0, 0.01 / MU0, 16.0 / MU0)
        driver.drive(system, stopping_mxHxm=10.0)
        
        # 2. Monotonic descending sweep to find Hc
        hc, mr_ms = fm.find_switching_field_monotonic(
            system,
            b_start_T=2.0 if soft_frac > 0 else 4.5,
            b_step_coarse_T=0.20,
            b_max_T=7.5,
            tol_fine_T=0.02,
            stopping_mxhxm=10.0
        )
        
        elapsed = time.time() - t0
        print(f"  Result: Hc = {hc:.3f} T, Mr/Ms = {mr_ms:.4f} (elapsed: {elapsed:.1f} s)")
        
        results.append({
            "run_name": run_name,
            "cell_nm": cell_nm,
            "soft_fraction": soft_frac,
            "total_cells": mesh.n[0] * mesh.n[1] * mesh.n[2],
            "Hc_T": hc,
            "Mr_Ms": mr_ms,
            "elapsed_s": elapsed
        })
        
        df_res = pd.DataFrame(results)
        df_res.to_csv(os.path.join(outdir, "full_cap_results.csv"), index=False)
        
    print("\n" + "=" * 75)
    print("FULL CAP STUDY COMPLETE")
    print(df_res.to_string(index=False))
    print("=" * 75)
    return df_res


if __name__ == "__main__":
    run_full_cap_study()
