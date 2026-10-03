"""
validate_micromagnetics.py
Acceptance tests for the corrected FePt Janus-cap micromagnetic pipeline.

These are PHYSICS tests, not unit tests. A simulation that fails any of them is
wrong regardless of whether it runs to completion. Run every one BEFORE
producing any production result.

Material parameters are those of Table S4 of the manuscript.
"""
import os
import sys
import numpy as np
from scipy import ndimage

# Ensure tclsh is in PATH for Windows OOMMF runner
for p in [r"C:\Users\admin\miniforge3\Library\bin", r"C:\Users\admin\miniforge3\envs\ubermag_env\Library\bin"]:
    if os.path.exists(p) and p not in os.environ.get("PATH", ""):
        os.environ["PATH"] = p + ";" + os.environ.get("PATH", "")

# ----------------------------------------------------------------------------
# Material constants (Table S4). Do not change without changing the targets.
# ----------------------------------------------------------------------------
MU0   = 4e-7 * np.pi
MS    = 1.0e6      # A/m
A_EX  = 1.0e-11    # J/m
KU_L10 = 6.6e6     # J/m^3   hard, ordered
KU_A1  = 1.0e4     # J/m^3   soft, disordered
T_CAP  = 60e-9     # m
ALPHA  = 0.5       # Gilbert damping (quasi-static relaxation)

HK_AM = 2 * KU_L10 / (MU0 * MS)     # A/m
HK_T  = 2 * KU_L10 / MS             # T   = 13.200 T

L_ANIS = np.sqrt(A_EX / KU_L10)               # 1.231 nm  <-- governs the mesh
L_WALL = np.pi * L_ANIS                       # 3.867 nm
L_MAGSTAT = np.sqrt(2 * A_EX / (MU0 * MS**2)) # 3.989 nm
L_MIN = min(L_ANIS, L_MAGSTAT)                # 1.231 nm

CELL_MAX_STRICT = L_MIN / 2     # 0.615 nm
CELL_MAX_MARGINAL = L_MIN       # 1.231 nm


# ----------------------------------------------------------------------------
# GATE 1 — mesh resolution. Refuse to run if this fails.
# ----------------------------------------------------------------------------
def assert_cell_size_ok(cell_m, strict=True):
    limit = CELL_MAX_STRICT if strict else CELL_MAX_MARGINAL
    if cell_m > limit:
        raise ValueError(
            f"cell size {cell_m*1e9:.3f} nm exceeds {limit*1e9:.3f} nm "
            f"(= sqrt(A/Ku){'/2' if strict else ''}). Exchange is not resolved; "
            f"the simulation will return a non-interacting Stoner-Wohlfarth "
            f"average, not micromagnetics."
        )
    n_through = T_CAP / cell_m
    if n_through < 20:
        raise ValueError(
            f"only {n_through:.1f} cells through the 60 nm cap thickness; "
            f"need >= 20 to resolve a wall crossing the film."
        )
    return True


# ----------------------------------------------------------------------------
# GATE 2 — the discretized body must be CONNECTED.
# This is the check that would have caught the 8/10/20 um runs.
# ----------------------------------------------------------------------------
def assert_mask_connected(mask, min_mean_neighbours=5.0):
    """mask: 3D boolean array of magnetic cells."""
    struct = ndimage.generate_binary_structure(3, 1)   # face connectivity
    lab, ncomp = ndimage.label(mask, structure=struct)
    if ncomp == 0 or mask.sum() == 0:
        raise ValueError("Mask is empty! No magnetic cells found.")
    sizes = np.bincount(lab.ravel())[1:]
    largest = sizes.max() / mask.sum()

    kern = np.zeros((3, 3, 3), dtype=np.int8)
    kern[1,1,0] = kern[1,1,2] = 1
    kern[1,0,1] = kern[1,2,1] = 1
    kern[0,1,1] = kern[2,1,1] = 1
    nb = ndimage.convolve(mask.astype(np.int8), kern, mode="constant")
    mean_nb = nb[mask].mean()

    if ncomp != 1 or largest < 0.999:
        raise ValueError(
            f"magnetic body is in {ncomp} disconnected pieces "
            f"(largest holds {largest*100:.1f}% of cells). Exchange cannot "
            f"couple them; each fragment reverses independently."
        )
    if mean_nb < min_mean_neighbours:
        raise ValueError(
            f"mean face-neighbour count {mean_nb:.2f} < {min_mean_neighbours}. "
            f"The shell is too thin relative to the cell size; it is a chain of "
            f"weakly-linked cells, not a continuous film."
        )
    return dict(ncomp=ncomp, largest_frac=largest, mean_neighbours=mean_nb)


# ----------------------------------------------------------------------------
# GATE 3 — solver convergence.
# ----------------------------------------------------------------------------
def assert_solver_tolerance(stopping_mxhxm_Am, dH_T):
    rel = stopping_mxhxm_Am / HK_AM
    if rel > 1e-5:
        raise ValueError(
            f"stopping mxHxm = {stopping_mxhxm_Am:.2e} A/m is {rel:.2e} of Hk. "
            f"Use <= 1e1 A/m (1e-6 of Hk). The manuscript's 2e4 A/m is 1.9e-3 "
            f"of Hk -- states are not converged at the switching field."
        )
    if dH_T > 0.05:
        raise ValueError(
            f"field step {dH_T} T gives sigma_Hc = {dH_T/2:.3f} T. Use <= 0.02 T "
            f"near switching, or bisect the transition."
        )
    return True


# ----------------------------------------------------------------------------
# GATE 4 — analytic limits the solver MUST reproduce.
# ----------------------------------------------------------------------------
def stoner_wohlfarth_loop(psis, weights, h_over_hk):
    """Quasi-static descending branch of an ensemble of NON-interacting
    uniaxial particles. psis = easy-axis angles to the field (rad)."""
    th_grid = np.linspace(-np.pi, np.pi, 3601)
    mz = np.zeros_like(h_over_hk, dtype=float)
    for psi, w in zip(psis, weights):
        th = 0.0
        for i, h in enumerate(h_over_hk):
            E = 0.5*np.sin(th_grid - psi)**2 - h*np.cos(th_grid)
            j = int(np.argmin(np.abs(th_grid - th)))
            while True:
                nbrs = [k for k in (j-1, j+1) if 0 <= k < len(th_grid)]
                k = min(nbrs, key=lambda k: E[k])
                if E[k] < E[j]: j = k
                else: break
            th = th_grid[j]
            mz[i] += w*np.cos(th)
    return mz / np.sum(weights)


def _metrics(h, m):
    mr = np.interp(0.0, h[::-1], m[::-1])
    i = np.where(np.diff(np.sign(m)))[0][0]
    hc = -np.interp(0.0, [m[i+1], m[i]], [h[i+1], h[i]])
    return mr / np.max(np.abs(m)), hc


# Reference values, in tesla, for Ku = 6.6e6 J/m^3, Ms = 1.0e6 A/m.
ANALYTIC_TARGETS = {
    "uniform_easy_axis_parallel_H": dict(Hc_T=13.200, Mr_Ms=1.000, tol_Hc=0.30),
    "easy_axis_45deg":              dict(Hc_T=6.600,  Mr_Ms=0.707, tol_Hc=0.20),
    "random_3D_noninteracting":     dict(Hc_T=6.320,  Mr_Ms=0.500, tol_Hc=0.20),
    "radial_hemisphere_noninteracting": dict(Hc_T=6.340, Mr_Ms=0.494, tol_Hc=0.20),
}

# HARD CEILING: any radial-easy-axis configuration whose coercivity EXCEEDS this
# without exchange is reporting a bug, not a result.
RADIAL_HC_CEILING_T = 0.479 * HK_T      # 6.32 T


def check_against_analytic(config_name, hc_T, mr_ms):
    tgt = ANALYTIC_TARGETS[config_name]
    ok_hc = abs(hc_T - tgt["Hc_T"]) <= tgt["tol_Hc"]
    ok_mr = abs(mr_ms - tgt["Mr_Ms"]) <= 0.03
    return dict(config=config_name, Hc_ok=ok_hc, Mr_ok=ok_mr,
                Hc_got=hc_T, Hc_expected=tgt["Hc_T"],
                Mr_got=mr_ms, Mr_expected=tgt["Mr_Ms"])


def flag_impossible_radial(hc_T, exchange_resolved=False):
    """The manuscript's 0% A1 radial result (12.17 T) trips this."""
    if not exchange_resolved and hc_T > RADIAL_HC_CEILING_T + 0.2:
        raise ValueError(
            f"radial-anisotropy Hc = {hc_T:.2f} T exceeds the non-interacting "
            f"ceiling of {RADIAL_HC_CEILING_T:.2f} T (0.48*Hk). With exchange "
            f"unresolved this is physically unreachable -- the run almost "
            f"certainly used a UNIFORM easy axis, not a radial one."
        )


def configure_oommf():
    """Helper to ensure ubermag oommfc runner uses the local OOMMF installation."""
    import oommfc as oc
    for tcl_candidate in [r"C:\Users\admin\Desktop\oommf\oommf.tcl"]:
        if os.path.exists(tcl_candidate):
            runner = oc.oommf.TclOOMMFRunner(oommf_tcl=tcl_candidate)
            oc.runner.runner = runner
            return runner
    return oc.runner.runner


def run_gate4_tests():
    """Runs the full Gate 4 analytical test suite."""
    print("=" * 70)
    print("GATE 4: Analytic Reproductions (Stoner-Wohlfarth Non-Interacting)")
    print("=" * 70)
    h = np.linspace(2.0, -2.0, 1601)

    # 1. Uniform || H
    m_uni = stoner_wohlfarth_loop([1e-4], [1.0], h)
    mr_u, hc_u = _metrics(h, m_uni)
    res_uni = check_against_analytic("uniform_easy_axis_parallel_H", hc_u * HK_T, mr_u)
    print(f"1. Uniform || H       : Hc = {res_uni['Hc_got']:6.2f} T (exp {res_uni['Hc_expected']:.2f}), "
          f"Mr/Ms = {res_uni['Mr_got']:.3f} | PASS: {res_uni['Hc_ok'] and res_uni['Mr_ok']}")

    # 2. Uniform 45 deg
    m_45 = stoner_wohlfarth_loop([np.pi/4], [1.0], h)
    mr_45, hc_45 = _metrics(h, m_45)
    res_45 = check_against_analytic("easy_axis_45deg", hc_45 * HK_T, mr_45)
    print(f"2. Uniform 45 deg     : Hc = {res_45['Hc_got']:6.2f} T (exp {res_45['Hc_expected']:.2f}), "
          f"Mr/Ms = {res_45['Mr_got']:.3f} | PASS: {res_45['Hc_ok'] and res_45['Mr_ok']}")

    # 3. 3D Random
    u = np.linspace(0.001, 0.999, 120)
    m_rand = stoner_wohlfarth_loop(np.arccos(u), np.ones_like(u), h)
    mr_rand, hc_rand = _metrics(h, m_rand)
    res_rand = check_against_analytic("random_3D_noninteracting", hc_rand * HK_T, mr_rand)
    print(f"3. 3D Random          : Hc = {res_rand['Hc_got']:6.2f} T (exp {res_rand['Hc_expected']:.2f}), "
          f"Mr/Ms = {res_rand['Mr_got']:.3f} | PASS: {res_rand['Hc_ok'] and res_rand['Mr_ok']}")

    # 4. Radial Hemisphere
    tt = np.linspace(0.001, np.pi/2 - 0.001, 120)
    m_rad = stoner_wohlfarth_loop(tt, np.sin(tt), h)
    mr_rad, hc_rad = _metrics(h, m_rad)
    res_rad = check_against_analytic("radial_hemisphere_noninteracting", hc_rad * HK_T, mr_rad)
    print(f"4. Radial Hemisphere  : Hc = {res_rad['Hc_got']:6.2f} T (exp {res_rad['Hc_expected']:.2f}), "
          f"Mr/Ms = {res_rad['Mr_got']:.3f} | PASS: {res_rad['Hc_ok'] and res_rad['Mr_ok']}")

    print("\nCeiling Check on original manuscript claim (12.17 T, Mr/Ms 1.00):")
    try:
        flag_impossible_radial(12.17, exchange_resolved=False)
    except ValueError as e:
        print(f"  -> CORRECTLY FLAGGED & REJECTED:\n     {e}")

    return all([res_uni['Hc_ok'], res_uni['Mr_ok'],
                res_45['Hc_ok'], res_45['Mr_ok'],
                res_rand['Hc_ok'], res_rand['Mr_ok'],
                res_rad['Hc_ok'], res_rad['Mr_ok']])


if __name__ == "__main__":
    print(f"PHYSICAL GOVERNING QUANTITIES:")
    print(f"  Anisotropy field mu0*Hk = {HK_T:.3f} T")
    print(f"  Anisotropy length L_anis = sqrt(A/Ku) = {L_ANIS*1e9:.3f} nm  <-- GOVERNING")
    print(f"  Domain wall width L_wall = pi*L_anis = {L_WALL*1e9:.3f} nm")
    print(f"  Magnetostatic exchange length = {L_MAGSTAT*1e9:.3f} nm")
    print(f"  Strict mesh limit = {CELL_MAX_STRICT*1e9:.3f} nm\n")

    # Run Gate 1 test
    print("GATE 1 Test: Checking cell size assertions...")
    try:
        assert_cell_size_ok(18e-9)
    except ValueError as e:
        print(f"  -> Successfully caught coarse 18 nm cell: {e}")
    assert_cell_size_ok(0.6e-9, strict=True)
    print("  -> Passed valid 0.6 nm cell assertion.\n")

    # Run Gate 2 test on sample masks
    print("GATE 2 Test: Checking mask connectivity assertions...")
    bad_mask = np.zeros((20, 20, 20), dtype=bool)
    bad_mask[2, 2, 2] = bad_mask[10, 10, 10] = True
    try:
        assert_mask_connected(bad_mask)
    except ValueError as e:
        print(f"  -> Successfully caught disconnected mask: {e}")
    good_mask = np.zeros((10, 10, 10), dtype=bool)
    good_mask[2:8, 2:8, 2:8] = True
    res_conn = assert_mask_connected(good_mask)
    print(f"  -> Passed continuous solid mask: {res_conn}\n")

    # Run Gate 3 test
    print("GATE 3 Test: Checking solver tolerance assertions...")
    try:
        assert_solver_tolerance(2e4, 0.3)
    except ValueError as e:
        print(f"  -> Successfully caught sloppy tolerance (2e4 A/m, 0.3 T): {e}")
    assert_solver_tolerance(10.0, 0.02)
    print("  -> Passed tight tolerance (10 A/m, 0.02 T).\n")

    # Run Gate 4 tests
    all_g4_pass = run_gate4_tests()
    print(f"\nOVERALL STATUS: {'ALL GATES PASSED' if all_g4_pass else 'FAILURES DETECTED'}")
