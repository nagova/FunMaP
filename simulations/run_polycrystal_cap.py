"""
run_polycrystal_cap.py  --  FunMaP, Fig. 5 "real cap" model
===========================================================================
Polycrystalline L1_0 FePt cap on the SAME D = 3 um tapered wedge as M1/M6
(geometry and Ms taken unchanged from fept_order_mechanisms.build_wedge_system).

What changes with respect to the ideal cap, and where each number comes from:
  * K_u = S^2 * K_L10, S from XRD (flat reference film): 0.70 +/- 0.09   [measured]
  * 3D Voronoi grains, mean equivalent diameter d_g = 19 nm
    (Scherrer, (111) of flat film, no instrumental correction -> lower bound) [measured]
  * Easy axis of each grain drawn uniformly on the sphere
    (random texture, the same assumption used to extract S)                 [measured/assumed, not fitted]
  * Grain-boundary (GB) layer: the cells within half a cell of each Voronoi bisector plane
    (one cell on each side, i.e. 2 cells across the boundary), with A_gb = x * A; K and the
    easy axis of GB cells are those of their own grain (only the coupling is weakened).
    x = intergrain exchange coupling, the ONLY unmeasured parameter -> swept, never fitted.
    Two-cell layer => the GB-GB link carries exactly x*A whatever averaging
    rule Oxs_ExchangePtwise uses between neighbouring cells.
  * A = 1e-11 J/m, Ms = 1e6 A/m, 0 K, MinDriver (identical to all other runs).

Outputs (one set per run) in OUT_DIR:
  <run_id>.csv            B_ext_T,Mx,My,Mz  (descending branch, sim_intervals.py-compatible)
  <run_id>_params.json    every parameter, RNG seed, grain count, GB volume fraction, interval
  <run_id>_slice_B<..>.npz y-midplane slice (x, z, m, norm, grain_id, gb) at 0 T, last
                          unswitched, first switched field  -> used to draw the "real cap"

Usage:
  python run_polycrystal_cap.py --phase control      # C1, C2 (must pass before anything else)
  python run_polycrystal_cap.py --phase main         # S = 0.70, x in {1, 0.3, 0.1, 0}
  python run_polycrystal_cap.py --phase band         # S = 0.61 and 0.79, same x
  python run_polycrystal_cap.py --phase seed2        # second grain realisation, S = 0.70
  python run_polycrystal_cap.py --phase fine --x X   # best x at cell 1.3 nm (final check)
"""
import os, sys, json, time, argparse
import numpy as np, pandas as pd
from scipy.spatial import cKDTree

for p in [r"C:\Users\admin\miniforge3\envs\ubermag_env\Library\bin", r"C:\Users\admin\miniforge3\Library\bin"]:
    if os.path.exists(p) and p not in os.environ.get("PATH", ""):
        os.environ["PATH"] = p + ";" + os.environ.get("PATH", "")

EXPLORATORY_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "archive", "exploratory_scripts"))
sys.path.insert(0, EXPLORATORY_DIR)
import fept_order_mechanisms as fom
import discretisedfield as df
import micromagneticmodel as mm
import oommfc as oc

OUT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results", "simulations", "loops_polycrystal"))
D_SPHERE, WY0 = 3.0e-6, 80.0e-9          # identical to M1
A_EX = 1.0e-11
B_SAT, B_TILT, STOPPING = 12.0, 0.01, 10.0
M1_SCHEDULE = [-0.5, -1.5, -2.5, -3.5, -4.5, -5.1, -5.25, -5.4, -5.55, -5.7, -5.85,
               -6.0, -6.15, -6.3, -6.45, -6.6, -6.8]
POLY_SCHEDULE = ([round(-0.05 * k, 2) for k in range(1, 81)]          # 0.05 T steps to -4.0 T
                 + [-4.5, -5.0, -5.5, -6.0, -7.0])                     # full branch for SQUID overlay


# ---------------------------------------------------------------- helpers
def cell_centres(mesh):
    p = getattr(mesh.region, "pmin", None)
    p = np.array(p if p is not None else mesh.region.p1, float)
    n, c = np.array(mesh.n), np.array(mesh.cell, float)
    ax = [p[i] + (np.arange(n[i]) + 0.5) * c[i] for i in range(3)]
    return np.stack(np.meshgrid(*ax, indexing="ij"), -1), ax, c


def voronoi_grains(pts_mag, box_lo, box_hi, d_g, rng, w_gb):
    """3D Poisson-Voronoi seeds with mean equivalent-sphere diameter d_g.
    Returns grain id and GB flag (distance to bisector plane < w_gb/2)."""
    vol = np.prod(box_hi - box_lo)
    n_seed = max(1, int(round(vol / (np.pi / 6 * d_g ** 3))))
    seeds = box_lo + rng.random((n_seed, 3)) * (box_hi - box_lo)
    d, i = cKDTree(seeds).query(pts_mag, k=2)
    s1, s2 = seeds[i[:, 0]], seeds[i[:, 1]]
    dist_plane = (d[:, 1] ** 2 - d[:, 0] ** 2) / (2 * np.linalg.norm(s2 - s1, axis=1))
    return i[:, 0], dist_plane < w_gb / 2, n_seed


def random_axes(n, rng):
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def build_system(run_id, S, x, cell_nm, d_g_nm, seed, axes):
    """axes = 'random' (polycrystal) or 'radial' (keep the M1 easy axis in every grain)."""
    base, mesh = fom.build_wedge_system(D_sphere=D_SPHERE, cell_nm=cell_nm, Wy0=WY0,
                                        S_pole=1.0, p=0.0, name=run_id + "_base")
    norm = base.m.norm                                  # Ms(r): geometry exactly as M1
    ms = norm.array[..., 0]
    mag = ms > 0
    # sanity: the constructor must not modify A or Ms with S (the old M2 A->0 issue)
    assert np.allclose(ms[mag], ms[mag].max()), "Ms is not uniform inside the wedge"
    pts, _, c = cell_centres(mesh)
    rng = np.random.default_rng(seed)
    lo, hi = pts[mag].min(0) - c / 2, pts[mag].max(0) + c / 2
    gid, gb, n_grain = voronoi_grains(pts[mag], lo, hi, d_g_nm * 1e-9, rng, c[0])

    u = np.zeros(pts.shape); u[..., 2] = 1.0
    if axes == "random":
        u[mag] = random_axes(n_grain, rng)[gid]
    else:                                               # radial axis field of M1
        u0 = base.energy.uniaxialanisotropy.u
        u = (u0 if isinstance(u0, df.Field) else df.Field(mesh, nvdim=3, value=u0)).array.copy()
    K = np.zeros(ms.shape); K[mag] = S ** 2 * fom.KU_L10
    A = np.full(ms.shape, A_EX)
    a_mag = A[mag]; a_mag[gb] = x * A_EX if x > 0 else 1e-25; A[mag] = a_mag

    gid_full = -np.ones(ms.shape, int); gid_full[mag] = gid
    gb_full = np.zeros(ms.shape, bool); gb_full[mag] = gb

    sysm = mm.System(name=run_id)
    sysm.energy = (mm.Exchange(A=df.Field(mesh, nvdim=1, value=A[..., None]))
                   + mm.Demag()
                   + mm.UniaxialAnisotropy(K=df.Field(mesh, nvdim=1, value=K[..., None]),
                                           u=df.Field(mesh, nvdim=3, value=u))
                   + mm.Zeeman(H=(0, 0, 0)))
    sysm.m = df.Field(mesh, nvdim=3, value=(0, 0, 1), norm=norm)
    info = dict(n_grain=int(n_grain), gb_volume_fraction=float(gb.mean()),
                n_magnetic_cells=int(mag.sum()), Ku_Jm3=float(S ** 2 * fom.KU_L10))
    return sysm, mesh, gid_full, gb_full, info


def save_slice(sysm, gid, gb, path):
    m = sysm.m.array; n = sysm.m.norm.array[..., 0]; j = m.shape[1] // 2
    _, ax, _ = cell_centres(sysm.m.mesh)
    np.savez_compressed(path, x=ax[0] * 1e9, z=ax[2] * 1e9, m=m[:, j, :, :], norm=n[:, j, :],
                        grain_id=gid[:, j, :], gb=gb[:, j, :])


def run(run_id, S, x, cell_nm, d_g_nm, seed, axes, schedule):
    t0 = time.time(); os.makedirs(OUT_DIR, exist_ok=True)
    sysm, mesh, gid, gb, info = build_system(run_id, S, x, cell_nm, d_g_nm, seed, axes)
    drv = oc.MinDriver(); vol = sysm.m.norm.integrate().item()
    mvec = lambda: [getattr(sysm.m, k).integrate().item() / vol for k in "xyz"]
    def at(B):
        sysm.energy.zeeman.H = (0, B_TILT / fom.MU0, B / fom.MU0)
        drv.drive(sysm, stopping_mxHxm=STOPPING); return mvec()
    rec = [dict(B_ext_T=B_SAT, **dict(zip(["Mx", "My", "Mz"], at(B_SAT))))]
    rec.append(dict(B_ext_T=0.0, **dict(zip(["Mx", "My", "Mz"], at(0.0)))))
    save_slice(sysm, gid, gb, os.path.join(OUT_DIR, f"{run_id}_slice_B+0.00.npz"))
    prev_slice, interval = None, None
    for B in schedule:
        mx, my, mz = at(B)
        rec.append(dict(B_ext_T=round(B, 3), Mx=mx, My=my, Mz=mz))
        print(f"{run_id}  B={B:+6.2f} T  <mz>={mz:+.4f}", flush=True)
        if interval is None and mz < 0:
            Bpre = rec[-2]["B_ext_T"]; interval = [abs(Bpre), abs(B)]
            if prev_slice is not None:
                os.replace(prev_slice, os.path.join(OUT_DIR, f"{run_id}_slice_B{Bpre:+.2f}.npz"))
            save_slice(sysm, gid, gb, os.path.join(OUT_DIR, f"{run_id}_slice_B{B:+.2f}.npz"))
            if axes == "radial" or S == 1.0:           # controls: single jump, stop after one more step
                break
        if interval is None:                            # keep the latest unswitched state
            prev_slice = os.path.join(OUT_DIR, f"{run_id}_tmp_slice.npz"); save_slice(sysm, gid, gb, prev_slice)
        if interval is not None and mz < -0.95:
            break
    if prev_slice and os.path.exists(prev_slice):
        os.remove(prev_slice)
    pd.DataFrame(rec)[["B_ext_T", "Mx", "My", "Mz"]].to_csv(os.path.join(OUT_DIR, run_id + ".csv"), index=False)
    out = dict(run_id=run_id, S=S, x_intergrain=x, cell_nm=cell_nm, d_g_nm=d_g_nm, seed=seed, axes=axes,
               A=A_EX, Ms=1e6, B_sat=B_SAT, B_tilt=B_TILT, stopping_mxHxm=STOPPING,
               Mr_over_Ms=rec[1]["Mz"],
               switching_interval_T=interval if interval else f"no switch down to {abs(schedule[-1])} T",
               minutes=(time.time() - t0) / 60, **info)
    json.dump(out, open(os.path.join(OUT_DIR, run_id + "_params.json"), "w"), indent=2)
    print(json.dumps(out, indent=2)); return out


def rid(S, x, cell, dg, seed, axes):
    return (f"Poly_wedge_D3um_cell{str(cell).replace('.', 'p')}nm_S{int(round(S*100)):03d}"
            f"_x{int(round(x*100)):03d}_dg{dg}_seed{seed}_{axes}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", required=True, choices=["control", "main", "band", "seed2", "fine"])
    ap.add_argument("--cell", type=float, default=2.0)
    ap.add_argument("--dg", type=int, default=19)
    ap.add_argument("--x", type=float, default=None)
    a = ap.parse_args(); fom.setup_oommf_runner()
    X = [1.0, 0.3, 0.1, 0.0]
    xs = [a.x] if a.x is not None else X
    if a.phase == "control":
        r1 = run(rid(1.0, 1.0, a.cell, a.dg, 1, "radial") + "_C1", 1.0, 1.0, a.cell, a.dg, 1, "radial", M1_SCHEDULE)
        print("C1 must switch in [6.45, 6.60] T (M1 f_A1 = 0 %, 1.3 nm) within one 0.15 T step, Mr ~ 0.656")
    elif a.phase == "main":
        for x in xs: run(rid(0.70, x, a.cell, a.dg, 1, "random"), 0.70, x, a.cell, a.dg, 1, "random", POLY_SCHEDULE)
    elif a.phase == "band":
        for S in (0.61, 0.79):
            for x in xs: run(rid(S, x, a.cell, a.dg, 1, "random"), S, x, a.cell, a.dg, 1, "random", POLY_SCHEDULE)
    elif a.phase == "seed2":
        for x in xs: run(rid(0.70, x, a.cell, a.dg, 2, "random"), 0.70, x, a.cell, a.dg, 2, "random", POLY_SCHEDULE)
    elif a.phase == "fine":
        assert a.x is not None, "--x required"
        run(rid(0.70, a.x, 1.3, a.dg, 1, "random"), 0.70, a.x, 1.3, a.dg, 1, "random", POLY_SCHEDULE)
