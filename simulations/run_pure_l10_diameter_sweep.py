"""
run_pure_l10_diameter_sweep.py
Executes the Pure L1_0 baseline across all particle diameters:
D = 200 nm, 500 nm, 1.0 um, 3.0 um, 5.0 um, 8.0 um, 10.0 um.

Configuration:
- Geometry: Tapered wedge (open y-BC) with realistic line-of-sight thickness t(theta) = t0 * cos(theta)
- Material: 100% pure ordered L1_0 phase (f_A1 = 0.0, S = 1.0, Ku = 6.59 MJ/m^3, Ms = 1.0 MA/m)
- Mesh: Sub-exchange resolution cell_nm = 1.3 nm (< l_ex = 2.1 nm)
- Field solver: Monotonic descending with optimized switching schedule around 6.0 - 6.8 T.
- Output: Logs each run to results/pure_l10_diameter_summary.csv and results/loops/
"""

import os
import sys
import time
import psutil
import numpy as np
import pandas as pd

# Ensure tclsh is in PATH and OOMMFTCL is set for Windows OOMMF runner
oommf_path = r"C:\Users\admin\Desktop\oommf\oommf.tcl"
if os.path.exists(oommf_path):
    os.environ["OOMMFTCL"] = oommf_path

for p in [r"C:\Users\admin\miniforge3\envs\ubermag_env\Library\bin", r"C:\Users\admin\miniforge3\Library\bin"]:
    if os.path.exists(p) and p not in os.environ.get("PATH", ""):
        os.environ["PATH"] = p + ";" + os.environ.get("PATH", "")

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

import fept_order_mechanisms as fom
import discretisedfield as df
import micromagneticmodel as mm
import oommfc as oc

DIAMETERS = [200e-9, 500e-9, 1.0e-6, 3.0e-6, 5.0e-6, 8.0e-6, 10.0e-6]
CELL_NM = 1.3
SUMMARY_CSV = os.path.join(fom.RESULTS_DIR, "pure_l10_diameter_summary.csv")


def run_pure_l10_sweep(diameters=None):
    if diameters is None:
        diameters = DIAMETERS

    fom.setup_oommf_runner()
    print("=" * 80)
    print("PURE L1_0 DIAMETER INVARIANCE SWEEP (t(theta) = t0*cos(theta), sub-l_ex cell = 1.3 nm)")
    print(f"Target Diameters: {[f'{d*1e6:g}um' if d>=1e-6 else f'{d*1e9:g}nm' for d in diameters]}")
    print("=" * 80)

    # Load existing summary if present
    if os.path.exists(SUMMARY_CSV):
        df_summary = pd.read_csv(SUMMARY_CSV)
    else:
        df_summary = pd.DataFrame(columns=[
            "run_id", "geometry", "D", "cell_nm", "f_A1", "mu0_Hc", "Mr_Ms", "M_perp_Ms", "RAM_GB", "wall_s"
        ])

    for d_val in diameters:
        d_um = d_val * 1e6
        d_label = f"{int(round(d_val*1e9))}nm" if d_val < 1e-6 else f"{d_um:g}um"
        run_id = f"L10_pure_wedge_D{d_label}_cell1p3nm"
        loop_file = os.path.join(fom.LOOPS_DIR, f"{run_id}.csv")

        # Check if already completed in summary CSV
        if run_id in df_summary["run_id"].values and os.path.exists(loop_file):
            row = df_summary[df_summary["run_id"] == run_id].iloc[0]
            print(f">> Skipping {run_id} (already completed): Hc = {row['mu0_Hc']:.4f} T, Mr = {row['Mr_Ms']:.4f}")
            continue

        # Special check: D = 3 um already completed in M1 sweep
        if np.isclose(d_val, 3.0e-6):
            m1_3um_file = os.path.join(fom.LOOPS_DIR, "M1_wedge_D3um_cell1p3nm_fA1_00pct.csv")
            if os.path.exists(m1_3um_file) and os.path.exists(fom.MASTER_CSV):
                df_m = pd.read_csv(fom.MASTER_CSV)
                m1_row = df_m[df_m["run_id"] == "M1_wedge_D3um_cell1p3nm_fA1_00pct"]
                if len(m1_row) > 0:
                    r = m1_row.iloc[0]
                    print(f">> Adopting existing D=3um pure L10 result from M1: Hc = {r['mu0_Hc']:.4f} T, Mr = {r['Mr_Ms']:.4f}")
                    # Copy loop file
                    df_loop_exist = pd.read_csv(m1_3um_file)
                    df_loop_exist.to_csv(loop_file, index=False)
                    new_entry = {
                        "run_id": run_id,
                        "geometry": "tapered_wedge_open_y",
                        "D": d_val,
                        "cell_nm": CELL_NM,
                        "f_A1": 0.0,
                        "mu0_Hc": r["mu0_Hc"],
                        "Mr_Ms": r["Mr_Ms"],
                        "M_perp_Ms": r["M_perp_Ms"],
                        "RAM_GB": r["RAM_GB"],
                        "wall_s": r["wall_s"]
                    }
                    df_summary = pd.concat([df_summary, pd.DataFrame([new_entry])], ignore_index=True)
                    df_summary.to_csv(SUMMARY_CSV, index=False)
                    continue

        print(f"\n" + "-" * 60)
        print(f"Starting {run_id} (D = {d_label})")
        print("-" * 60)

        t0 = time.time()
        ram_init = psutil.virtual_memory().used / 1e9

        # Set lateral width Wy0 to keep cell count tractable
        wy0 = 80.0e-9 if d_val <= 3.0e-6 else 50.0e-9

        system, mesh = fom.build_wedge_system(
            D_sphere=d_val,
            cell_nm=CELL_NM,
            Wy0=wy0,
            f_A1=0.0,  # 100% pure L1_0
            name=run_id
        )

        n_cells = mesh.n[0] * mesh.n[1] * mesh.n[2]
        print(f"  Mesh: {mesh.n[0]} x {mesh.n[1]} x {mesh.n[2]} = {n_cells:,} bounding cells")

        m_norm_int = system.m.norm.integrate().item()
        if m_norm_int == 0:
            raise ValueError(f"System {run_id} contains zero magnetic volume!")

        def get_m_vec():
            mx = system.m.x.integrate().item() / m_norm_int
            my = system.m.y.integrate().item() / m_norm_int
            mz = system.m.z.integrate().item() / m_norm_int
            return mx, my, mz

        driver = oc.MinDriver()
        b_tilt = 0.01

        # 1. Pre-saturation at +12.0 T
        print("  [1/3] Pre-saturating at +12.0 T...")
        system.energy.zeeman.H = (0, b_tilt / fom.MU0, 12.0 / fom.MU0)
        driver.drive(system, stopping_mxHxm=10.0)

        # 2. Relax at 0.0 T (Remanence)
        print("  [2/3] Relaxing at 0.0 T (Remanence)...")
        system.energy.zeeman.H = (0, b_tilt / fom.MU0, 0.0)
        driver.drive(system, stopping_mxHxm=10.0)
        mx_0, my_0, mz_0 = get_m_vec()
        mr_ms = mz_0
        m_perp_ms = np.sqrt(mx_0**2 + my_0**2)
        print(f"  Remanence: Mz/Ms = {mr_ms:.4f}, |M_perp|/Ms = {m_perp_ms:.4f}")

        loop_records = [
            {"B_ext_T": 12.0, "Mx": 0.0, "My": 0.0, "Mz": 1.0},
            {"B_ext_T": 0.0,  "Mx": mx_0, "My": my_0, "Mz": mz_0},
        ]

        # 3. Optimized descending field schedule around switching threshold (~6.0 - 6.8 T)
        # Approach fields with moderate spacing, then dense spacing in the critical bracket
        b_schedule = [
            -1.0, -3.0, -5.0, -5.5, -5.8, -6.0, -6.1, -6.2, -6.3, -6.4, -6.5, -6.6, -6.8, -7.0
        ]

        print("  [3/3] Descending field sweep...")
        switched = False
        hc = None
        ram_peak = psutil.virtual_memory().used / 1e9

        for b_ext in b_schedule:
            system.energy.zeeman.H = (0, b_tilt / fom.MU0, b_ext / fom.MU0)
            driver.drive(system, stopping_mxHxm=10.0)
            mx, my, mz = get_m_vec()
            loop_records.append({"B_ext_T": b_ext, "Mx": mx, "My": my, "Mz": mz})
            ram_curr = psutil.virtual_memory().used / 1e9
            ram_peak = max(ram_peak, ram_curr)
            print(f"    B = {b_ext:6.2f} T -> Mz/Ms = {mz:+.4f} (RAM: {ram_curr:.1f} GB)")

            if mz < 0.0 and not switched:
                switched = True
                prev = loop_records[-2]
                b_prev, mz_prev = prev["B_ext_T"], prev["Mz"]
                hc_interp = b_prev + (0.0 - mz_prev) * (b_ext - b_prev) / (mz - mz_prev)
                hc = abs(hc_interp)
                print(f"  ==> SWITCHED at mu0*Hc = {hc:.4f} T!")
                # Take one post-switching verification point if available
                next_b = b_ext - 0.2
                system.energy.zeeman.H = (0, b_tilt / fom.MU0, next_b / fom.MU0)
                driver.drive(system, stopping_mxHxm=10.0)
                mx_p, my_p, mz_p = get_m_vec()
                loop_records.append({"B_ext_T": next_b, "Mx": mx_p, "My": my_p, "Mz": mz_p})
                print(f"    Post-switch B = {next_b:6.2f} T -> Mz/Ms = {mz_p:+.4f}")
                break

        if hc is None:
            # Fallback interpolation
            mzs = [r["Mz"] for r in loop_records]
            bs = [r["B_ext_T"] for r in loop_records]
            idx = np.where(np.diff(np.sign(mzs)))[0]
            if len(idx) > 0:
                i = idx[0]
                hc = abs(bs[i] + (0.0 - mzs[i]) * (bs[i+1] - bs[i]) / (mzs[i+1] - mzs[i]))
            else:
                hc = abs(bs[-1])

        wall_s = time.time() - t0
        print(f"  Completed {run_id} in {wall_s:.1f} s ({wall_s/60:.2f} min). Hc = {hc:.4f} T, Mr = {mr_ms:.4f}")

        # Save loop CSV
        loop_df = pd.DataFrame(loop_records)
        loop_df.to_csv(loop_file, index=False)

        # Update summary CSV
        new_entry = {
            "run_id": run_id,
            "geometry": "tapered_wedge_open_y",
            "D": d_val,
            "cell_nm": CELL_NM,
            "f_A1": 0.0,
            "mu0_Hc": hc,
            "Mr_Ms": mr_ms,
            "M_perp_Ms": m_perp_ms,
            "RAM_GB": round(ram_peak, 2),
            "wall_s": round(wall_s, 1)
        }
        df_summary = pd.concat([df_summary, pd.DataFrame([new_entry])], ignore_index=True)
        df_summary.to_csv(SUMMARY_CSV, index=False)

    print("\n" + "=" * 80)
    print("ALL REQUESTED PURE L1_0 RUNS PROCESSED!")
    print(df_summary.to_string(index=False))
    print("=" * 80)
    return df_summary


if __name__ == "__main__":
    run_pure_l10_sweep()
