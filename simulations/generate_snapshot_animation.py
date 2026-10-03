"""
generate_snapshot_animation.py
Generates 3-panel magnetization state animation (Loop + XZ slice + XY top view)
for FePt Janus caps with graded order S(theta) = S_pole * (cos theta)^p.
Matches the publication format of FunMaP.
"""

import os
import sys
import glob
import time
import numpy as np
import pandas as pd
import imageio

# Ensure conda Library\bin DLLs are in PATH
for p in [r"C:\Users\admin\miniforge3\Library\bin", r"C:\Users\admin\miniforge3\envs\ubermag_env\Library\bin"]:
    if os.path.exists(p) and p not in os.environ.get("PATH", ""):
        os.environ["PATH"] = p + ";" + os.environ.get("PATH", "")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import fept_order_mechanisms as fom
import discretisedfield as df
import micromagneticmodel as mm
import oommfc as oc


def run_snapshot_simulation(
    D_sphere=200e-9,
    cell_nm=1.3,
    Wy0=80e-9,
    S_pole=0.90,
    p_order=0.5,
    b_max=3.0,
    outdir=None
):
    """Runs descending loop saving magnetization arrays at each field step."""
    d_label = f"{int(round(D_sphere*1e9))}nm" if D_sphere < 1e-6 else f"{D_sphere*1e6:g}um"
    if outdir is None:
        outdir = os.path.join(fom.RESULTS_DIR, f"snapshots_D{d_label}_Sp90_p0p5")
    os.makedirs(outdir, exist_ok=True)

    run_id = f"snap_D{d_label}_cell{str(cell_nm).replace('.', 'p')}nm"
    print(f"Building system for {run_id} (D = {d_label})...")

    system, mesh = fom.build_wedge_system(
        D_sphere=D_sphere,
        cell_nm=cell_nm,
        Wy0=Wy0,
        S_pole=S_pole,
        p=p_order,
        name=run_id
    )

    fom.setup_oommf_runner()
    driver = oc.MinDriver()
    b_tilt = 0.01

    m_norm_int = system.m.norm.integrate().item()
    if m_norm_int == 0:
        raise ValueError("System contains zero magnetic volume!")

    # Field schedule focused on switching
    # For D = 200 nm, switching is around -2.0 to -2.3 T
    b_schedule = [
        3.0, 1.5, 0.0,
        -0.5, -1.0, -1.3, -1.5, -1.7, -1.8, -1.9,
        -2.0, -2.1, -2.2, -2.3, -2.4, -2.5, -3.0
    ]
    b_schedule = sorted(list(set(b_schedule)), reverse=True)

    print(f"Executing descending sweep across {len(b_schedule)} field steps...")

    loop_records = []
    state_files = []

    # Pre-saturate
    print(f"Pre-saturating at +{b_max:.1f} T...")
    system.energy.zeeman.H = (0, b_tilt / fom.MU0, b_max / fom.MU0)
    driver.drive(system, stopping_mxHxm=10.0)

    R = D_sphere / 2.0
    cell_m = cell_nm * 1e-9

    for idx, b_ext in enumerate(b_schedule):
        system.energy.zeeman.H = (0, b_tilt / fom.MU0, b_ext / fom.MU0)
        driver.drive(system, stopping_mxHxm=10.0)

        mx = system.m.x.integrate().item() / m_norm_int
        my = system.m.y.integrate().item() / m_norm_int
        mz = system.m.z.integrate().item() / m_norm_int
        print(f"  Step {idx:02d}/{len(b_schedule)-1} | B = {b_ext:+6.2f} T | Mz/Ms = {mz:+6.3f}")

        loop_records.append({"step": idx, "B_ext_T": b_ext, "Mx": mx, "My": my, "Mz": mz})

        # Extract 2D central slice m_z(x, z)
        m_arr = system.m.array # shape (nx, ny, nz, 3)
        ny = m_arr.shape[1]
        iy_mid = ny // 2
        # Extract XZ slice along meridian
        mz_slice = m_arr[:, iy_mid, :, 2] / fom.MS

        # Save slice and metadata
        state_npz = os.path.join(outdir, f"state_{idx:03d}_B{b_ext:+.3f}T.npz")
        np.savez_compressed(
            state_npz,
            step=idx,
            B_ext_T=b_ext,
            Mz_over_Ms=mz,
            mz_slice=mz_slice,
            D_sphere=D_sphere,
            cell_m=cell_m,
            t0=fom.T0_CAP,
            theta_max=fom.THETA_MAX
        )
        state_files.append(state_npz)

    df_loop = pd.DataFrame(loop_records)
    loop_csv = os.path.join(outdir, "hysteresis.csv")
    df_loop.to_csv(loop_csv, index=False)
    print(f"Hysteresis loop saved: {loop_csv}")

    return outdir, loop_csv, state_files


def render_animation(outdir, loop_csv=None, gif_name="snapshot_animation.gif", fps=3):
    """Renders 3-panel frames and compiles into an animated GIF."""
    if loop_csv is None:
        loop_csv = os.path.join(outdir, "hysteresis.csv")
    df_loop = pd.read_csv(loop_csv)

    state_files = sorted(glob.glob(os.path.join(outdir, "state_*.npz")))
    if not state_files:
        raise FileNotFoundError(f"No state files found in {outdir}")

    frames_dir = os.path.join(outdir, "frames")
    os.makedirs(frames_dir, exist_ok=True)

    frame_images = []

    # Get geometry parameters from first state
    first_data = np.load(state_files[0])
    D = float(first_data["D_sphere"])
    R = D / 2.0
    t0 = float(first_data["t0"])
    theta_max = float(first_data["theta_max"])

    use_nm = D < 1e-6
    scale = 1e9 if use_nm else 1e6
    unit = "nm" if use_nm else r"$\mu$m"

    R_scaled = R * scale
    t0_scaled = t0 * scale

    for i, p_state in enumerate(state_files):
        data = np.load(p_state)
        b_ext = float(data["B_ext_T"])
        mz_avg = float(data["Mz_over_Ms"])
        step = int(data["step"])
        mz_slice = data["mz_slice"] # (nx, nz)

        nx, nz = mz_slice.shape
        theta_arr = np.linspace(0, theta_max, nx)
        r_norm = np.linspace(0, 1, nz)

        TH, RN = np.meshgrid(theta_arr, r_norm, indexing='ij')
        R_grid = R_scaled + RN * (t0_scaled * np.cos(TH))

        X = R_grid * np.sin(TH)
        Z = R_grid * np.cos(TH)

        # Build figure
        fig = plt.figure(figsize=(12, 3.8))
        gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.2, 1.2], wspace=0.35)

        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax3 = fig.add_subplot(gs[2])

        # 1. Hysteresis Loop
        # Descending branch (solid)
        df_desc = df_loop[df_loop.get("branch", "descending") == "descending"] if "branch" in df_loop.columns else df_loop
        ax1.plot(df_desc["B_ext_T"], df_desc["Mz"], 'k-', alpha=0.6, lw=1.8, label="descending")

        # Ascending branch (dashed)
        if "branch" in df_loop.columns and "ascending" in df_loop["branch"].values:
            df_asc = df_loop[df_loop["branch"] == "ascending"]
        else:
            # Physical inversion symmetry: M_asc(B) = -M_desc(-B)
            df_asc = pd.DataFrame({
                "B_ext_T": -df_desc["B_ext_T"],
                "Mz": -df_desc["Mz"]
            }).sort_values("B_ext_T")
        ax1.plot(df_asc["B_ext_T"], df_asc["Mz"], 'k--', alpha=0.6, lw=1.8, label="ascending")

        ax1.scatter([b_ext], [mz_avg], color="#d7191c", s=90, zorder=5)
        ax1.axhline(0, color="gray", lw=0.6, ls="--")
        ax1.axvline(0, color="gray", lw=0.6, ls="--")
        b_lim = max(abs(df_loop["B_ext_T"].min()), abs(df_loop["B_ext_T"].max())) + 0.5
        ax1.set_xlim(-b_lim, b_lim)
        ax1.set_ylim(-1.08, 1.08)
        ax1.set_xlabel(r"$B_{\mathrm{ext}}$ (T)", fontsize=10)
        ax1.set_ylabel(r"$M_z / M_s$", fontsize=10)
        ax1.set_title("Hysteresis loop", fontsize=11, fontweight="bold")
        ax1.grid(True, ls=":", alpha=0.5)
        ax1.legend(loc="lower right", fontsize=8)

        # 2. XZ Central Slice
        ax2.pcolormesh(X, Z, mz_slice, cmap="RdBu_r", vmin=-1.0, vmax=1.0, shading="auto")
        ax2.pcolormesh(-X, Z, mz_slice, cmap="RdBu_r", vmin=-1.0, vmax=1.0, shading="auto")
        ax2.set_aspect("equal")
        ax2.set_xlabel(f"x ({unit})", fontsize=10)
        ax2.set_ylabel(f"z ({unit})", fontsize=10)
        ax2.set_title("XZ central slice", fontsize=11, fontweight="bold")
        ax2.set_xlim(- (R_scaled + t0_scaled) * 1.05, (R_scaled + t0_scaled) * 1.05)
        ax2.set_ylim(0, (R_scaled + t0_scaled) * 1.05)

        # 3. XY Top View
        phi = np.linspace(0, 2*np.pi, 80)
        TH_top, PHI = np.meshgrid(theta_arr, phi, indexing='ij')
        R_mid = (R_scaled + 0.5 * t0_scaled * np.cos(TH_top)) * np.sin(TH_top)
        X_top = R_mid * np.cos(PHI)
        Y_top = R_mid * np.sin(PHI)
        mz_th_avg = np.nanmean(mz_slice, axis=1)
        mz_top = np.repeat(mz_th_avg[:, np.newaxis], len(phi), axis=1)

        c2 = ax3.pcolormesh(X_top, Y_top, mz_top, cmap="RdBu_r", vmin=-1.0, vmax=1.0, shading="auto")
        ax3.set_aspect("equal")
        ax3.set_xlabel(f"x ({unit})", fontsize=10)
        ax3.set_ylabel(f"y ({unit})", fontsize=10)
        ax3.set_title("XY top view", fontsize=11, fontweight="bold")
        limit_xy = (R_scaled + t0_scaled * np.cos(theta_max)) * np.sin(theta_max) * 1.15
        ax3.set_xlim(-limit_xy, limit_xy)
        ax3.set_ylim(-limit_xy, limit_xy)

        cbar = fig.colorbar(c2, ax=ax3, fraction=0.046, pad=0.04)
        cbar.set_label(r"$m_z / M_s$", fontsize=10)

        d_str = f"{D*1e9:.0f} nm" if use_nm else f"{D*1e6:g} " + r"$\mu$m"
        fig.suptitle(f"FePt cap, {d_str} sphere — graded order $S(\\theta)$, 0 K", fontsize=12, fontweight="bold", y=1.03)

        caption = f"state {step:03d}/{len(state_files):03d}   descending   B = {b_ext:+6.2f} T   Mz/Ms = {mz_avg:+6.3f}"
        fig.text(0.5, -0.06, caption, ha="center", fontsize=10, family="monospace")

        frame_path = os.path.join(frames_dir, f"frame_{step:03d}.png")
        fig.savefig(frame_path, dpi=180, bbox_inches="tight")
        plt.close(fig)

        frame_images.append(imageio.imread(frame_path))

    # Save animated GIF
    gif_path = os.path.join(outdir, gif_name)
    imageio.mimsave(gif_path, frame_images, fps=fps, loop=0)
    print(f"Animated GIF saved: {gif_path}")

    return gif_path


if __name__ == "__main__":
    t0 = time.time()
    outdir, loop_csv, state_files = run_snapshot_simulation(
        D_sphere=200e-9,
        cell_nm=1.3,
        Wy0=80e-9,
        S_pole=0.90,
        p_order=0.5
    )
    gif_path = render_animation(outdir, loop_csv)
    print(f"Done in {time.time()-t0:.1f} s!")
