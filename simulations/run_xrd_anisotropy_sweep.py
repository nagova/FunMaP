"""
run_xrd_anisotropy_sweep.py
Simulation of D = 3.0 um FePt Janus cap tapered wedge (R2b geometry)
using the exact constructor from Model M1 (fept_order_mechanisms.build_wedge_system).

Physical Framework:
- Uniform anisotropy: K = S^2 * Ku0 with Ku0 = 6.59e6 J/m^3 for S in {1.00, 0.61, 0.70, 0.79}.
  * S = 1.00 serves as the strict control benchmark (must reproduce M1 f_A1 = 0% in [6.45, 6.60] T).
  * S = 0.61, 0.70, 0.79 correspond to the experimental XRD order parameter S = 0.70 +/- 0.09.
- Constant saturation magnetization: Ms = 1.0e6 A/m (held fixed to isolate the anisotropy effect).
- Constant exchange stiffness: A = 1.0e-11 J/m (strictly uniform).
- Sub-exchange cell size: cell = 1.3 nm (< l_ex = 3.99 nm).
- Transverse width: Wy0 = 80.0 nm (identical to M1, yielding ~4.98M cells and ~12 GB RAM).
- Pre-saturation field: +12.0 T (required because for S = 0.79, mu0*H_K = 2K/Ms ~ 8.2 T).
- Transverse tilt: By = 0.01 T (breaks numerical symmetry, identical to M1).
- Stopping criterion: stopping_mxHxm = 10.0 A/m (identical to M1).
- Switching window for S < 1.0: fine steps dH = 0.05 T between -1.50 T and -5.00 T.
- Output: writes CSV files to results/simulations/loops/ with columns B_ext_T,Mx,My,Mz.
"""

import os
import sys
import time
import numpy as np
import pandas as pd

# 1. Environment & OOMMF Runner Setup
for p in [r"C:\Users\admin\miniforge3\envs\ubermag_env\Library\bin", r"C:\Users\admin\miniforge3\Library\bin"]:
    if os.path.exists(p) and p not in os.environ.get("PATH", ""):
        os.environ["PATH"] = p + ";" + os.environ.get("PATH", "")

# Add exploratory_scripts to import the exact M1 constructor
EXPLORATORY_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "archive", "exploratory_scripts"))
if EXPLORATORY_DIR not in sys.path:
    sys.path.insert(0, EXPLORATORY_DIR)

import fept_order_mechanisms as fom
import micromagneticmodel as mm
import oommfc as oc

D_SPHERE = 3.0e-6     # 3.0 um particle diameter
CELL_NM = 1.3         # 1.3 nm sub-exchange cell
WY0 = 80.0e-9         # 80 nm max width (identical to M1)
B_SAT = 12.0          # 12 T pre-saturation (identical to M1)
B_TILT = 0.01         # 0.01 T y-tilt (identical to M1)
STOPPING = 10.0       # stopping_mxHxm = 10 A/m (identical to M1)

OUT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "MPI-IS", "results", "simulations", "loops"))
os.makedirs(OUT_DIR, exist_ok=True)


def run_targeted_sweep(system, S_val, b_schedule=None, b_fine_start=-1.50, b_fine_end=-5.00, dH=0.05):
    """
    Executes pre-saturation at +12 T, remanence at 0 T, and descending sweep.
    If b_schedule is provided (e.g. for S = 1.0 control), it executes those exact fields.
    Otherwise, it executes coarse steps down to b_fine_start and fine steps dH down to b_fine_end.
    """
    assert b_fine_start < 0, f"b_fine_start must be negative, got {b_fine_start}"
    assert b_fine_end < b_fine_start, f"b_fine_end ({b_fine_end}) must be more negative than b_fine_start ({b_fine_start})"

    driver = oc.MinDriver()
    m_norm_int = system.m.norm.integrate().item()
    if m_norm_int == 0:
        raise ValueError("System contains zero magnetic volume!")

    def get_m_vec():
        mx = system.m.x.integrate().item() / m_norm_int
        my = system.m.y.integrate().item() / m_norm_int
        mz = system.m.z.integrate().item() / m_norm_int
        return mx, my, mz

    print(f"\n{'='*70}")
    print(f"RUNNING: S = {S_val:.2f} (Ku = {(S_val**2)*fom.KU_L10*1e-6:.3f} MJ/m^3)")
    print(f"{'='*70}")

    # 1. Pre-saturation at +12.0 T
    print(f"  Pre-saturating at +{B_SAT:.1f} T...")
    system.energy.zeeman.H = (0, B_TILT / fom.MU0, B_SAT / fom.MU0)
    driver.drive(system, stopping_mxHxm=STOPPING)

    # 2. Relax at 0.0 T for remanence
    print("  Relaxing at 0.0 T (Remanence)...")
    system.energy.zeeman.H = (0, B_TILT / fom.MU0, 0.0)
    driver.drive(system, stopping_mxHxm=STOPPING)
    mx_0, my_0, mz_0 = get_m_vec()
    mr_ms = mz_0
    print(f"  Remanence Mr/Ms = {mr_ms:.4f}")

    records = [
        {"B_ext_T": round(B_SAT, 3), "Mx": 0.0, "My": 0.0, "Mz": 1.0},
        {"B_ext_T": 0.0, "Mx": mx_0, "My": my_0, "Mz": mz_0}
    ]

    # Build field list
    if b_schedule is None:
        field_list = []
        # Coarse steps
        b_c = -0.5
        while b_c > b_fine_start:
            field_list.append(round(b_c, 3))
            b_c -= 0.5
        # Fine steps
        for b_f in np.arange(b_fine_start, b_fine_end - 0.001, -dH):
            field_list.append(round(float(b_f), 3))
    else:
        field_list = [round(float(b), 3) for b in b_schedule]

    h_pre, h_post = None, None
    switched = False

    for b_ext in field_list:
        if switched:
            break
        system.energy.zeeman.H = (0, B_TILT / fom.MU0, b_ext / fom.MU0)
        driver.drive(system, stopping_mxHxm=STOPPING)
        mx, my, mz = get_m_vec()
        records.append({"B_ext_T": round(b_ext, 3), "Mx": mx, "My": my, "Mz": mz})
        print(f"    B = {b_ext:6.3f} T -> Mz/Ms = {mz:+.4f}")

        if mz < 0.0 and not switched:
            switched = True
            h_post = abs(b_ext)
            h_pre = abs(records[-2]["B_ext_T"])
            print(f"  >>> SWITCHED! Interval: [{min(h_pre, h_post):.3f}, {max(h_pre, h_post):.3f}] T")
            # Extra point past switching
            b_post_extra = round(b_ext - 0.20, 3)
            system.energy.zeeman.H = (0, B_TILT / fom.MU0, b_post_extra / fom.MU0)
            driver.drive(system, stopping_mxHxm=STOPPING)
            mx_e, my_e, mz_e = get_m_vec()
            records.append({"B_ext_T": b_post_extra, "Mx": mx_e, "My": my_e, "Mz": mz_e})
            break

    interval = [min(h_pre, h_post), max(h_pre, h_post)] if switched else None
    res = {
        "S": S_val,
        "Ku_MJm3": (S_val**2) * fom.KU_L10 * 1e-6,
        "Mr_Ms": mr_ms,
        "interval": interval,
        "interval_str": f"[{interval[0]:.2f}, {interval[1]:.2f}]" if interval else f"no switch down to {abs(field_list[-1]):.2f} T"
    }
    return res, records


def main():
    fom.setup_oommf_runner()
    t_start = time.time()

    # Field schedule for M1 control S = 1.0 (exact descending list of M1 f_A1 = 0%)
    m1_0pct_schedule = [-0.5, -1.5, -2.5, -3.5, -4.5, -5.1, -5.25, -5.4, -5.55, -5.7, -5.85, -6.0, -6.15, -6.3, -6.45, -6.6]

    targets = [
        {"S": 1.00, "schedule": m1_0pct_schedule},
        {"S": 0.61, "schedule": None, "start": -1.50, "end": -5.00},
        {"S": 0.70, "schedule": None, "start": -1.50, "end": -5.00},
        {"S": 0.79, "schedule": None, "start": -1.50, "end": -5.00},
    ]

    summary = []
    for cfg in targets:
        S = cfg["S"]
        run_id = f"Ku_S2_wedge_D3um_cell1p3nm_S{int(round(S*100)):03d}"
        
        # Build using the EXACT Model M1 constructor: p=0.0 means uniform S(theta) = S
        system, mesh = fom.build_wedge_system(
            D_sphere=D_SPHERE,
            cell_nm=CELL_NM,
            Wy0=WY0,
            S_pole=S,
            p=0.0,
            name=run_id
        )

        if cfg["schedule"] is not None:
            res, recs = run_targeted_sweep(system, S, b_schedule=cfg["schedule"])
        else:
            res, recs = run_targeted_sweep(system, S, b_fine_start=cfg["start"], b_fine_end=cfg["end"], dH=0.05)

        # Write CSV with exact columns matching sim_intervals.py
        out_csv = os.path.join(OUT_DIR, f"{run_id}.csv")
        df_out = pd.DataFrame(recs)[["B_ext_T", "Mx", "My", "Mz"]]
        df_out.to_csv(out_csv, index=False)
        print(f"Saved loop to: {out_csv}")
        print(f"Result: {res['interval_str']}")
        summary.append(res)

        # Strict validation on S = 1.00 control run before spending compute on sweeps
        if S == 1.00:
            ctrl_interval = res.get("interval")
            ctrl_mr = res.get("Mr_Ms")
            if ctrl_interval != [6.45, 6.60] or not (0.64 <= ctrl_mr <= 0.67):
                raise RuntimeError(
                    f"CRITICAL: Control benchmark S = 1.00 failed! "
                    f"Expected interval [6.45, 6.60] T and Mr/Ms ~ 0.656, "
                    f"but obtained {res['interval_str']} and Mr/Ms = {ctrl_mr:.4f}. "
                    f"Aborting execution to prevent invalid sweep calculations."
                )
            print(">>> Control run S = 1.00 PASSED benchmark ([6.45, 6.60] T, Mr/Ms ~ 0.656). Continuing sweep.")

    print("\n" + "="*70)
    print("ALL RUNS COMPLETE")
    print(pd.DataFrame(summary)[["S", "Ku_MJm3", "Mr_Ms", "interval_str"]].to_string(index=False))
    print(f"Total elapsed: {(time.time() - t_start)/60:.1f} minutes")
    print("="*70)


if __name__ == "__main__":
    main()
