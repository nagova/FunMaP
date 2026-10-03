"""
run_mesh_convergence.py
Gate 5 Mesh Convergence Study for FePt Cap (d = 300 nm).

Evaluates coercivity Hc as a function of mesh cell size:
  cell = [4.0, 3.0, 2.0, 1.5, 1.0, 0.6] nm
Demonstrates that coercivity plateaus as cell size resolves the
anisotropy length L_anis = 1.231 nm, converting referee objection
into a rigorous methods contribution for the manuscript.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Ensure local imports work
sys.path.append(os.path.dirname(__file__))
import fept_micromagnetics as fm

def run_convergence_study():
    fm.setup_oommf_runner()
    diameter = 300e-9  # 300 nm cap
    cells_nm = [4.0, 3.0, 2.0, 1.5, 1.0, 0.6]
    
    outdir = os.path.join("results", "mesh_convergence_300nm")
    os.makedirs(outdir, exist_ok=True)
    csv_path = os.path.join(outdir, "mesh_convergence_results.csv")
    
    results = []
    
    print("=" * 75)
    print("GATE 5: MESH CONVERGENCE STUDY (d = 300 nm, Pure L1_0 FePt Cap)")
    print(f"Target Anisotropy Length L_anis = {fm.L_ANIS*1e9:.3f} nm | L_wall = {fm.L_WALL*1e9:.3f} nm")
    print("=" * 75)

    for cell_nm in cells_nm:
        cell_m = cell_nm * 1e-9
        name = f"conv_300nm_cell_{str(cell_nm).replace('.', 'p')}nm"
        print(f"\n--- Testing cell size: {cell_nm:.2f} nm ---")
        t0 = time.time()
        
        # Build quarter-wedge system for efficiency
        system, mesh, mask = fm.build_cap_system(
            name=name,
            diameter=diameter,
            cell_size=cell_m,
            mode="Radial",
            soft_fraction=0.0,
            quarter_wedge=True,
            strict_mesh=False,
            check_gate1=False, # We benchmark coarse cells to show convergence
            include_demag=True
        )
        
        n_mag_cells = int(mask.sum())
        print(f"  Mesh: {mesh.n[0]}x{mesh.n[1]}x{mesh.n[2]} ({mesh.n[0]*mesh.n[1]*mesh.n[2]:,d} total cells, {n_mag_cells:,d} magnetic cells)")
        
        # Bisection to determine Hc to <= 0.02 T
        hc, mr = fm.find_switching_field_bisection(
            system,
            b_low_T=4.7,
            b_high_T=5.5,
            tol_T=0.02,
            stopping_mxhxm=10.0
        )
        
        elapsed = time.time() - t0
        print(f"  Result: cell = {cell_nm:.2f} nm -> Hc = {hc:.3f} T, Mr/Ms = {mr:.3f} (elapsed: {elapsed:.1f} s)")
        
        results.append({
            "cell_nm": cell_nm,
            "cell_m": cell_m,
            "total_cells": mesh.n[0] * mesh.n[1] * mesh.n[2],
            "magnetic_cells": n_mag_cells,
            "Hc_T": hc,
            "Mr_Ms": mr,
            "elapsed_s": elapsed
        })
        
        # Save intermediate CSV
        df_res = pd.DataFrame(results)
        df_res.to_csv(csv_path, index=False)
        
    print("\n" + "=" * 75)
    print("CONVERGENCE STUDY COMPLETE")
    print(df_res.to_string(index=False))
    print("=" * 75)
    
    # Generate publication figures
    plot_convergence_figure(df_res, outdir)
    return df_res


def plot_convergence_figure(df_res, outdir):
    """Plots Hc vs. cell size showing the convergence plateau."""
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    
    cells = df_res["cell_nm"].values
    hcs = df_res["Hc_T"].values
    
    ax.plot(cells, hcs, "o-", color="#1f77b4", linewidth=2.2, markersize=8, label=r"Cap $d = 300\,$nm ($t = 60\,$nm)")
    
    # Add vertical guidelines for governing physical lengths
    ax.axvline(fm.L_ANIS * 1e9, color="crimson", linestyle="--", linewidth=1.5,
               label=rf"$L_\mathrm{{anis}} = \sqrt{{A/K_u}} = {fm.L_ANIS*1e9:.2f}\,$nm (Governing)")
    ax.axvline(fm.CELL_MAX_STRICT * 1e9, color="darkorange", linestyle=":", linewidth=1.5,
               label=rf"Strict limit $= L_\mathrm{{anis}}/2 = {fm.CELL_MAX_STRICT*1e9:.2f}\,$nm")
    ax.axvline(fm.L_WALL * 1e9, color="gray", linestyle="-.", linewidth=1.2,
               label=rf"Wall width $\pi L_\mathrm{{anis}} = {fm.L_WALL*1e9:.2f}\,$nm")
    
    # Stoner-Wohlfarth non-interacting radial ceiling
    ax.axhline(fm.RADIAL_HC_CEILING_T, color="purple", linestyle="--", alpha=0.7,
               label=rf"SW Non-interacting ceiling $= {fm.RADIAL_HC_CEILING_T:.2f}\,$T")
    
    ax.set_xlabel("Mesh Cell Size (nm)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Coercivity $\mu_0 H_c$ (T)", fontsize=12, fontweight="bold")
    ax.set_title("Mesh Convergence Study: Resolving FePt Exchange Length", fontsize=13, fontweight="bold", pad=12)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(frameon=True, facecolor="white", edgecolor="none", fontsize=9, loc="upper right")
    
    plt.tight_layout()
    png_path = os.path.join(outdir, "mesh_convergence_300nm.png")
    svg_path = os.path.join(outdir, "mesh_convergence_300nm.svg")
    fig.savefig(png_path, dpi=300)
    fig.savefig(svg_path)
    plt.close(fig)
    print(f"Saved publication figures:\n  -> {png_path}\n  -> {svg_path}")


if __name__ == "__main__":
    run_convergence_study()
