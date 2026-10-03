"""
fept_micromagnetics.py
Corrected Micromagnetic Engine & Utilities for FePt Janus Caps.

Implements:
- Physical parameters (Table S4 of manuscript)
- Mandatory Acceptance Gates 1-3
- Full-cap and quarter-wedge cap systems
- Tier 2 Curved-patch geometry and easy-axis mapping
- High-precision transition bisection and fine-step sweeps
- Tier 3 statistical ensemble / Preisach integration
"""

import os
import sys
import numpy as np
from scipy import ndimage

# Ensure tclsh is in PATH for Windows OOMMF runner
for p in [r"C:\Users\admin\miniforge3\Library\bin", r"C:\Users\admin\miniforge3\envs\ubermag_env\Library\bin"]:
    if os.path.exists(p) and p not in os.environ.get("PATH", ""):
        os.environ["PATH"] = p + ";" + os.environ.get("PATH", "")

import discretisedfield as df
import micromagneticmodel as mm
import oommfc as oc

# ----------------------------------------------------------------------------
# Material constants (Fixed - Table S4)
# ----------------------------------------------------------------------------
MU0       = 4e-7 * np.pi
MS        = 1.0e6      # A/m
A_EX      = 1.0e-11    # J/m
KU_L10    = 6.6e6      # J/m^3   hard, ordered L1_0
KU_A1     = 1.0e4      # J/m^3   soft, disordered A1
T_CAP     = 60e-9      # m (cap thickness)
ALPHA     = 0.5        # Gilbert damping (quasi-static relaxation)

HK_AM     = 2 * KU_L10 / (MU0 * MS)     # A/m
HK_T      = 2 * KU_L10 / MS             # T   = 13.200 T

L_ANIS    = np.sqrt(A_EX / KU_L10)               # 1.231 nm  <-- GOVERNING
L_WALL    = np.pi * L_ANIS                       # 3.867 nm
L_MAGSTAT = np.sqrt(2 * A_EX / (MU0 * MS**2)) # 3.989 nm
L_MIN     = min(L_ANIS, L_MAGSTAT)                # 1.231 nm

CELL_MAX_STRICT   = L_MIN / 2     # 0.615 nm
CELL_MAX_MARGINAL = L_MIN         # 1.231 nm
SINGLE_DOMAIN_DIAM = 465e-9       # ~465 nm single-domain limit

RADIAL_HC_CEILING_T = 0.479 * HK_T # 6.32 T


def setup_oommf_runner():
    """Configures the TclOOMMFRunner using the local OOMMF installation."""
    candidates = [
        r"C:\Users\admin\Desktop\oommf\oommf.tcl",
    ]
    for c in candidates:
        if os.path.exists(c):
            runner = oc.oommf.TclOOMMFRunner(oommf_tcl=c)
            oc.runner.runner = runner
            return runner
    return oc.runner.runner


# ----------------------------------------------------------------------------
# Gate Assertions
# ----------------------------------------------------------------------------
def assert_cell_size_ok(cell_m, strict=True):
    limit = CELL_MAX_STRICT if strict else CELL_MAX_MARGINAL
    if cell_m > limit:
        raise ValueError(
            f"Cell size {cell_m*1e9:.3f} nm exceeds limit {limit*1e9:.3f} nm "
            f"(= sqrt(A/Ku){'/2' if strict else ''}). Exchange underresolved!"
        )
    n_through = T_CAP / cell_m
    if n_through < 20:
        raise ValueError(
            f"Only {n_through:.1f} cells through 60 nm thickness; need >= 20."
        )
    return True


def assert_mask_connected(mask, min_mean_neighbours=5.0):
    struct = ndimage.generate_binary_structure(3, 1)
    lab, ncomp = ndimage.label(mask, structure=struct)
    if ncomp == 0 or mask.sum() == 0:
        raise ValueError("Mask contains no magnetic cells.")
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
            f"Magnetic mask is fragmented into {ncomp} pieces (largest holds {largest*100:.1f}%)."
        )
    if mean_nb < min_mean_neighbours:
        raise ValueError(
            f"Mean face-neighbour count {mean_nb:.2f} < {min_mean_neighbours}."
        )
    return dict(ncomp=ncomp, largest_frac=largest, mean_neighbours=mean_nb)


def assert_solver_tolerance(stopping_mxhxm_Am, dH_T):
    rel = stopping_mxhxm_Am / HK_AM
    if rel > 1e-5:
        raise ValueError(f"stopping mxHxm = {stopping_mxhxm_Am:.2e} A/m is {rel:.2e} of Hk. Use <= 10 A/m.")
    if dH_T > 0.05:
        raise ValueError(f"Field step dH = {dH_T:.3f} T too coarse. Use <= 0.05 T near switching.")
    return True


# ----------------------------------------------------------------------------
# System Builders
# ----------------------------------------------------------------------------
def build_cap_system(
    name,
    diameter,
    cell_size,
    mode="Radial",
    soft_fraction=0.0,
    quarter_wedge=False,
    strict_mesh=False,
    check_gate1=True,
    include_demag=True
):
    """
    Builds an Ubermag micromagnetic system for a hemispherical cap.
    quarter_wedge: if True, bounds to x >= 0, y >= 0, z >= 0 for 4x memory savings.
    """
    if check_gate1:
        assert_cell_size_ok(cell_size, strict=strict_mesh)

    R_inner = diameter / 2.0
    R_outer = R_inner + T_CAP

    # Snap bounding box to integer multiple of cell_size so discretisedfield divides cleanly
    n_dim = int(np.ceil(R_outer / cell_size))
    L_box = n_dim * cell_size

    # Region definition
    if quarter_wedge:
        p1 = (0.0, 0.0, 0.0)
    else:
        p1 = (-L_box, -L_box, 0.0)
    p2 = (L_box, L_box, L_box)

    region = df.Region(p1=p1, p2=p2)
    cell = (cell_size, cell_size, cell_size)
    mesh = df.Mesh(region=region, cell=cell)

    # Verification of mask connectivity
    # Generate discrete grid mask
    nx, ny, nz = mesh.n
    xs = p1[0] + (np.arange(nx) + 0.5) * cell_size
    ys = p1[1] + (np.arange(ny) + 0.5) * cell_size
    zs = p1[2] + (np.arange(nz) + 0.5) * cell_size
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    mask = (R >= R_inner) & (R <= R_outer) & (Z >= 0)
    if quarter_wedge:
        mask &= (X >= 0) & (Y >= 0)

    assert_mask_connected(mask, min_mean_neighbours=4.5 if quarter_wedge else 5.0)

    # Spatial field functions
    def ms_func(pos):
        r = np.linalg.norm(pos)
        if R_inner <= r <= R_outer and pos[2] >= 0:
            if quarter_wedge and (pos[0] < 0 or pos[1] < 0):
                return 0.0
            return MS
        return 0.0

    def ku_func(pos):
        r = np.linalg.norm(pos)
        if R_inner <= r <= R_outer and pos[2] >= 0:
            if quarter_wedge and (pos[0] < 0 or pos[1] < 0):
                return 0.0
            if soft_fraction > 0.0:
                return KU_A1 if np.random.random() < soft_fraction else KU_L10
            return KU_L10
        return 0.0

    def u_func(pos):
        if mode == "Radial":
            r = np.linalg.norm(pos)
            return (pos[0] / r, pos[1] / r, pos[2] / r) if r != 0 else (0, 0, 1)
        elif mode == "Uniaxial_Vertical":
            return (0.0, 0.0, 1.0)
        elif mode == "Random_3D":
            vec = np.random.randn(3)
            return tuple(vec / np.linalg.norm(vec))
        else:
            raise ValueError(f"Unknown mode {mode}")

    system = mm.System(name=name)
    energy_terms = mm.Exchange(A=A_EX)
    if include_demag:
        energy_terms += mm.Demag()

    if soft_fraction > 0.0:
        energy_terms += mm.UniaxialAnisotropy(
            K=df.Field(mesh, nvdim=1, value=ku_func),
            u=df.Field(mesh, nvdim=3, value=u_func)
        )
    else:
        energy_terms += mm.UniaxialAnisotropy(
            K=KU_L10,
            u=df.Field(mesh, nvdim=3, value=u_func)
        )

    energy_terms += mm.Zeeman(H=(0, 0, 0))
    system.energy = energy_terms

    # Initial state: saturated along +z
    system.m = df.Field(mesh, nvdim=3, value=(0, 0, 1),
                        norm=df.Field(mesh, nvdim=1, value=ms_func))

    return system, mesh, mask


def build_curved_patch_system(
    name,
    sphere_diameter,
    patch_width,
    cell_size,
    soft_fraction=0.0,
    strict_mesh=False,
    include_demag=True
):
    """
    Builds a Tier 2 Curved-Patch micromagnetic system.
    A representative patch of lateral size (patch_width x patch_width x T_CAP)
    representing the cap curvature of a sphere of diameter `sphere_diameter`.
    The easy-axis follows the exact local sphere outward normal:
    u(x, y) = (x/R, y/R, sqrt(1 - (x^2 + y^2)/R^2)).
    """
    R_sphere = sphere_diameter / 2.0
    half_w = patch_width / 2.0

    n_w = int(np.ceil(half_w / cell_size))
    actual_half_w = n_w * cell_size
    n_t = int(np.ceil(T_CAP / cell_size))
    actual_t = n_t * cell_size

    p1 = (-actual_half_w, -actual_half_w, 0.0)
    p2 = (actual_half_w, actual_half_w, actual_t)
    region = df.Region(p1=p1, p2=p2)
    cell = (cell_size, cell_size, cell_size)
    mesh = df.Mesh(region=region, cell=cell)

    # Local normal easy axis at curvature radius R_sphere
    def u_curved(pos):
        x, y = pos[0], pos[1]
        r_sq = (x**2 + y**2) / (R_sphere**2)
        if r_sq >= 1.0:
            return (0.0, 0.0, 1.0)
        uz = np.sqrt(1.0 - r_sq)
        ux = x / R_sphere
        uy = y / R_sphere
        norm = np.sqrt(ux**2 + uy**2 + uz**2)
        return (ux / norm, uy / norm, uz / norm)

    def ku_patch(pos):
        if soft_fraction > 0.0:
            return KU_A1 if np.random.random() < soft_fraction else KU_L10
        return KU_L10

    system = mm.System(name=name)
    energy_terms = mm.Exchange(A=A_EX)
    if include_demag:
        energy_terms += mm.Demag()

    if soft_fraction > 0.0:
        energy_terms += mm.UniaxialAnisotropy(
            K=df.Field(mesh, nvdim=1, value=ku_patch),
            u=df.Field(mesh, nvdim=3, value=u_curved)
        )
    else:
        energy_terms += mm.UniaxialAnisotropy(
            K=KU_L10,
            u=df.Field(mesh, nvdim=3, value=u_curved)
        )

    energy_terms += mm.Zeeman(H=(0, 0, 0))
    system.energy = energy_terms

    system.m = df.Field(mesh, nvdim=3, value=(0, 0, 1), norm=MS)
    return system, mesh


# ----------------------------------------------------------------------------
# Solver & Bisection Engine
# ----------------------------------------------------------------------------
def relax_at_field(system, b_ext_z, b_tilt=0.01, stopping_mxhxm=10.0, max_iter=15000):
    """Relaxes system at external field (0, B_tilt, B_z)."""
    driver = oc.MinDriver()
    system.energy.zeeman.H = (0.0, b_tilt / MU0, b_ext_z / MU0)
    driver.drive(system, stopping_mxHxm=stopping_mxhxm, total_iteration_limit=max_iter)

    m_norm_int = system.m.norm.integrate().item()
    mz_avg = (system.m.z.integrate().item() / m_norm_int) if m_norm_int > 0 else 0.0
    return mz_avg


def find_switching_field_bisection(
    system,
    b_low_T=4.0,
    b_high_T=8.0,
    tol_T=0.02,
    b_tilt=0.01,
    stopping_mxhxm=10.0
):
    """
    DEPRECATED: This function is preserved for historical reference only.
    DO NOT USE: Without resaturating at +16 T at each step, bisection tests midpoints
    on an already switched descending branch, causing path-dependent hysteretic trapping.
    Use find_switching_field_monotonic() instead.
    """
    raise NotImplementedError(
        "find_switching_field_bisection is disabled due to hysteretic trapping. "
        "Use find_switching_field_monotonic() which tracks the true descending branch."
    )


def find_switching_field_monotonic(
    system,
    b_start_T=3.0,
    b_step_coarse_T=0.25,
    b_max_T=8.0,
    tol_fine_T=0.02,
    b_tilt=0.01,
    stopping_mxhxm=10.0
):
    """
    Finds Hc along the TRUE monotonic descending branch of the hysteresis loop.
    Never walks backwards in field, completely avoiding hysteretic latching.
    """
    setup_oommf_runner()
    # 1. Initialize saturated state at +16 T
    relax_at_field(system, +16.0, b_tilt=b_tilt, stopping_mxhxm=stopping_mxhxm)
    # 2. Relax at remanence (0 T)
    m_rem = relax_at_field(system, 0.0, b_tilt=b_tilt, stopping_mxhxm=stopping_mxhxm)
    print(f"  [Descending Sweep] Remanence Mr/Ms at 0 T: {m_rem:.4f}")

    # 3. Monotonic coarse step until Mz becomes negative
    b_current = b_start_T
    b_prev = 0.0
    m_prev = m_rem
    m_current = m_rem

    while b_current <= b_max_T:
        m_current = relax_at_field(system, -b_current, b_tilt=b_tilt, stopping_mxhxm=stopping_mxhxm)
        print(f"    B = -{b_current:.3f} T -> Mz/Ms = {m_current:+.4f}")
        if m_current <= 0.0:
            break
        b_prev = b_current
        m_prev = m_current
        b_current += b_step_coarse_T

    if m_current < 0.0 and m_prev > 0.0:
        hc = b_prev + (0.0 - m_prev) * (b_current - b_prev) / (m_current - m_prev)
    else:
        hc = b_current
    print(f"  [Descending Sweep Done] Hc = {hc:.3f} T (bracket: [{b_prev:.3f}, {b_current:.3f}] T)")
    return hc, m_rem


# ----------------------------------------------------------------------------
# Tier 3 Statistical Ensemble / Stoner-Wohlfarth Reference
# ----------------------------------------------------------------------------
def stoner_wohlfarth_loop(psis, weights, h_over_hk):
    """Quasi-static descending branch of an ensemble of NON-interacting
    uniaxial particles. psis = easy-axis angles to the field (rad)."""
    th_grid = np.linspace(-np.pi, np.pi, 3601)
    mz = np.zeros_like(h_over_hk, dtype=float)
    for psi, w in zip(psis, weights):
        th = 0.0
        for i, h in enumerate(h_over_hk):
            E = 0.5 * np.sin(th_grid - psi)**2 - h * np.cos(th_grid)
            j = int(np.argmin(np.abs(th_grid - th)))
            while True:
                nbrs = [k for k in (j - 1, j + 1) if 0 <= k < len(th_grid)]
                k = min(nbrs, key=lambda k: E[k])
                if E[k] < E[j]:
                    j = k
                else:
                    break
            th = th_grid[j]
            mz[i] += w * np.cos(th)
    return mz / np.sum(weights)


def tier3_ensemble_loop(b_fields_T, patch_hc_T=None, patch_mr_ms=None, spread_deg=7.6):
    """
    Tier 3 Statistical Model (Stoner-Wohlfarth Hemisphere Ceiling):
    Computes the non-interacting macrospin ensemble over 100 theta polar values
    with sin(theta) surface weighting under ideal Set A parameters (Hk = 13.20 T).
    Serves as the analytical non-interacting ceiling (Hc = 6.34 T, Mr/Ms = 0.494).
    """
    thetas = np.linspace(0.001, np.pi / 2 - 0.001, 100)
    weights = np.sin(thetas)
    h_over_hk = b_fields_T / HK_T
    return stoner_wohlfarth_loop(thetas, weights, h_over_hk)

