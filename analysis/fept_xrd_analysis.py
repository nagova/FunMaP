"""
fept_xrd_analysis.py
====================
Quantitative XRD Analysis for FePt Thin Films on SiO2 Spherical Substrates.

Re-fits experimental diffractograms to locate the true:
  - (001) superlattice reflection (~23.5° - 24.2°)
  - (110) superlattice reflection (~32.5° - 33.2°)
  - (111) fundamental reflection  (~40.8° - 41.3°)
  - (200) fundamental reflection  (~46.8° - 47.4°)
  - (002) fundamental reflection  (~48.0° - 49.2°)

Computes:
  - Interplanar d-spacings: d_001, d_111, d_200, d_002
  - Tetragonal lattice parameters:
      * c from d_001 and d_002
      * a from d_200 and d_111 (via 1/d_111² = 2/a² + 1/c²)
      * Axial ratio c/a (c/a = 1.00 for disordered A1, c/a ≈ 0.964 for fully ordered L10)
  - Chemical long-range order parameter S:
      * S_axial = sqrt((1 - c/a) / (1 - (c/a)_ordered))  [structural order]
      * S_intensity_002 = sqrt((I_001 / I_002)_meas / (I_001 / I_002)_calc) [diffraction order]
      * S_intensity_111 = sqrt((I_001 / I_111)_meas / (I_001 / I_111)_calc)

Outputs:
  - Summary CSV table in results directory
  - Multi-panel publication-ready fitting plots for each sample (PNG + SVG)
  - Parameter trend plots across sample series

Part of FunMaP (Functional Magnetic Particles Analysis Pipeline).
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
import sys
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")  # Headless non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import leastsq, minimize

# Physical constants
LAMBDA_CU_KA1 = 1.54056  # Angstrom (Cu Ka1)
C_OVER_A_ORDERED = 0.964  # Bulk fully ordered L10 FePt c/a ratio
C_OVER_A_DISORDERED = 1.000  # Disordered A1 cubic FePt c/a ratio


def pseudo_voigt_height(x: np.ndarray, x0: float, height: float, fwhm: float, eta: float) -> np.ndarray:
    """Standard height-parameterized pseudo-Voigt profile."""
    gamma = fwhm / 2.0
    dx = (x - x0) / max(gamma, 1e-6)
    lorentzian = 1.0 / (1.0 + dx**2)
    gaussian = np.exp(-np.log(2.0) * (dx**2))
    return height * (eta * lorentzian + (1.0 - eta) * gaussian)


def voigt_area(height: float, fwhm: float, eta: float) -> float:
    """Analytical integrated area under a height-parameterized pseudo-Voigt profile."""
    c_lorentz = np.pi / 2.0  # ~1.570796
    c_gauss = np.sqrt(np.pi / (4.0 * np.log(2.0)))  # ~1.064467
    return float(height * fwhm * (eta * c_lorentz + (1.0 - eta) * c_gauss))


def fit_peak_in_window(
    x: np.ndarray,
    y: np.ndarray,
    x_range: Tuple[float, float],
    x0_guess: float,
    fwhm_guess: float = 0.5,
    min_fwhm: float = 0.08,
    max_fwhm: float = 2.5,
) -> Optional[Dict[str, object]]:
    """Fits a single pseudo-Voigt peak with linear background in a specified 2theta window."""
    mask = (x >= x_range[0]) & (x <= x_range[1])
    xw, yw = x[mask], y[mask]
    if len(xw) < 8:
        return None

    bg_c = float(np.min(yw))
    bg_s = float((yw[-1] - yw[0]) / (xw[-1] - xw[0])) if xw[-1] != xw[0] else 0.0
    h0 = max(float(np.max(yw) - bg_c), 1.0)

    p0 = [x0_guess, h0, fwhm_guess, 0.5, bg_c, bg_s]

    def residual(p):
        x0, h, fwhm, eta, c, s = p
        if h < 0 or fwhm < min_fwhm or fwhm > max_fwhm or eta < 0.0 or eta > 1.0:
            return 1e6 * np.ones_like(xw)
        if x0 < x_range[0] or x0 > x_range[1]:
            return 1e6 * np.ones_like(xw)
        y_model = pseudo_voigt_height(xw, x0, h, fwhm, eta) + c + s * (xw - x_range[0])
        return yw - y_model

    try:
        popt, ier = leastsq(residual, p0, full_output=False, maxfev=3000)
        if ier not in (1, 2, 3, 4):
            raise RuntimeError("leastsq did not converge")
    except Exception:
        def loss(p):
            return float(np.sum(residual(p) ** 2))
        res = minimize(loss, p0, method="Nelder-Mead")
        popt = res.x

    x0, h, fwhm, eta, c, s = popt
    x0 = float(np.clip(x0, x_range[0], x_range[1]))
    h = max(float(h), 0.0)
    fwhm = float(np.clip(fwhm, min_fwhm, max_fwhm))
    eta = float(np.clip(eta, 0.0, 1.0))

    y_fit = pseudo_voigt_height(xw, x0, h, fwhm, eta) + c + s * (xw - x_range[0])
    ss_res = np.sum((yw - y_fit) ** 2)
    ss_tot = np.sum((yw - np.mean(yw)) ** 2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    area = voigt_area(h, fwhm, eta)

    return {
        "x0": x0,
        "height": h,
        "fwhm": fwhm,
        "eta": eta,
        "bg_c": float(c),
        "bg_s": float(s),
        "integrated_area": area,
        "r2": float(r2),
        "xw": xw,
        "yw": yw,
        "y_fit": y_fit,
    }


def fit_doublet_peaks(
    x: np.ndarray,
    y: np.ndarray,
    x_range: Tuple[float, float],
    x0_1_guess: float,
    x0_2_guess: float,
) -> Optional[Dict[str, object]]:
    """Fits two adjacent pseudo-Voigt peaks (e.g. 200 and 002) in a common window."""
    mask = (x >= x_range[0]) & (x <= x_range[1])
    xw, yw = x[mask], y[mask]
    if len(xw) < 15:
        return None

    bg_c = float(np.min(yw))
    bg_s = float((yw[-1] - yw[0]) / (xw[-1] - xw[0])) if xw[-1] != xw[0] else 0.0
    h0 = max(float(np.max(yw) - bg_c) / 2.0, 1.0)

    p0 = [x0_1_guess, h0, 0.5, 0.5, x0_2_guess, h0, 0.5, 0.5, bg_c, bg_s]

    def residual(p):
        x0_1, h_1, fwhm_1, eta_1, x0_2, h_2, fwhm_2, eta_2, c, s = p
        if h_1 < 0 or h_2 < 0 or fwhm_1 < 0.08 or fwhm_2 < 0.08 or fwhm_1 > 3.0 or fwhm_2 > 3.0:
            return 1e6 * np.ones_like(xw)
        if eta_1 < 0.0 or eta_1 > 1.0 or eta_2 < 0.0 or eta_2 > 1.0:
            return 1e6 * np.ones_like(xw)
        if x0_1 < x_range[0] or x0_1 > x_range[1] or x0_2 < x_range[0] or x0_2 > x_range[1]:
            return 1e6 * np.ones_like(xw)
        y_model = (
            pseudo_voigt_height(xw, x0_1, h_1, fwhm_1, eta_1)
            + pseudo_voigt_height(xw, x0_2, h_2, fwhm_2, eta_2)
            + c + s * (xw - x_range[0])
        )
        return yw - y_model

    try:
        popt, ier = leastsq(residual, p0, full_output=False, maxfev=4000)
        if ier not in (1, 2, 3, 4):
            raise RuntimeError("leastsq doublet did not converge")
    except Exception:
        def loss(p):
            return float(np.sum(residual(p) ** 2))
        res = minimize(loss, p0, method="Nelder-Mead")
        popt = res.x

    x0_1, h_1, fwhm_1, eta_1, x0_2, h_2, fwhm_2, eta_2, c, s = popt
    x0_1 = float(np.clip(x0_1, x_range[0], x_range[1]))
    x0_2 = float(np.clip(x0_2, x_range[0], x_range[1]))
    h_1 = max(float(h_1), 0.0)
    h_2 = max(float(h_2), 0.0)
    fwhm_1 = float(np.clip(fwhm_1, 0.08, 3.0))
    fwhm_2 = float(np.clip(fwhm_2, 0.08, 3.0))
    eta_1 = float(np.clip(eta_1, 0.0, 1.0))
    eta_2 = float(np.clip(eta_2, 0.0, 1.0))

    y_fit = (
        pseudo_voigt_height(xw, x0_1, h_1, fwhm_1, eta_1)
        + pseudo_voigt_height(xw, x0_2, h_2, fwhm_2, eta_2)
        + c + s * (xw - x_range[0])
    )
    ss_res = np.sum((yw - y_fit) ** 2)
    ss_tot = np.sum((yw - np.mean(yw)) ** 2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    return {
        "x0_1": x0_1,
        "h_1": h_1,
        "fwhm_1": fwhm_1,
        "area_1": voigt_area(h_1, fwhm_1, eta_1),
        "x0_2": x0_2,
        "h_2": h_2,
        "fwhm_2": fwhm_2,
        "area_2": voigt_area(h_2, fwhm_2, eta_2),
        "r2": float(r2),
        "xw": xw,
        "yw": yw,
        "y_fit": y_fit,
    }


def calculate_d_spacing(two_theta_deg: float, wavelength: float = LAMBDA_CU_KA1) -> float:
    """Computes Bragg d-spacing in Angstrom: d = lambda / (2 * sin(theta))."""
    theta_rad = np.radians(two_theta_deg / 2.0)
    return float(wavelength / (2.0 * np.sin(theta_rad)))


def load_diffractogram(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Loads two-column .xy XRD data, ignoring comments and header metadata."""
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith(("#", ";", "//", "*")):
                continue
            parts = s.replace(",", " ").split()
            if len(parts) < 2:
                continue
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    arr = np.array(rows, dtype=float)
    idx = np.argsort(arr[:, 0])
    return arr[idx, 0], arr[idx, 1]


def analyze_fept_file(path: str, output_dir: str) -> Dict[str, object]:
    """
    Performs complete profile re-fitting for a single FePt diffractogram:
      - Locates (001), (111), (200), (002)
      - Computes d_001, d_111, d_200, d_002
      - Calculates a, c, c/a, and chemical order parameters S
    """
    x, y = load_diffractogram(path)
    base_name = os.path.splitext(os.path.basename(path))[0]

    # 1. Fit (001) superlattice peak [22.8 - 25.2 deg]
    fit_001 = fit_peak_in_window(x, y, x_range=(22.8, 25.2), x0_guess=23.8, fwhm_guess=0.6)

    # 2. Fit (110) superlattice peak [31.5 - 34.5 deg]
    fit_110 = fit_peak_in_window(x, y, x_range=(31.5, 34.5), x0_guess=32.9, fwhm_guess=0.7)

    # 3. Fit (111) fundamental peak [40.0 - 42.4 deg]
    fit_111 = fit_peak_in_window(x, y, x_range=(40.0, 42.4), x0_guess=41.1, fwhm_guess=0.5)

    # 4. Fit (200) and (002) in the 46.0 - 49.8 deg region
    fit_doublet = fit_doublet_peaks(
        x, y, x_range=(46.0, 49.8), x0_1_guess=47.1, x0_2_guess=48.5
    )

    pos_001 = fit_001["x0"] if fit_001 else np.nan
    pos_110 = fit_110["x0"] if fit_110 else np.nan
    pos_111 = fit_111["x0"] if fit_111 else np.nan

    d_001 = calculate_d_spacing(pos_001) if np.isfinite(pos_001) else np.nan
    d_110 = calculate_d_spacing(pos_110) if np.isfinite(pos_110) else np.nan
    d_111 = calculate_d_spacing(pos_111) if np.isfinite(pos_111) else np.nan

    if fit_doublet and abs(fit_doublet["x0_2"] - fit_doublet["x0_1"]) > 0.3:
        pos_200 = min(fit_doublet["x0_1"], fit_doublet["x0_2"])
        pos_002 = max(fit_doublet["x0_1"], fit_doublet["x0_2"])
        area_200 = fit_doublet["area_1"] if fit_doublet["x0_1"] < fit_doublet["x0_2"] else fit_doublet["area_2"]
        area_002 = fit_doublet["area_2"] if fit_doublet["x0_1"] < fit_doublet["x0_2"] else fit_doublet["area_1"]
    else:
        fit_200 = fit_peak_in_window(x, y, x_range=(46.5, 47.8), x0_guess=47.1)
        fit_002 = fit_peak_in_window(x, y, x_range=(47.8, 49.5), x0_guess=48.5)
        pos_200 = fit_200["x0"] if fit_200 else np.nan
        pos_002 = fit_002["x0"] if fit_002 else np.nan
        area_200 = fit_200["integrated_area"] if fit_200 else np.nan
        area_002 = fit_002["integrated_area"] if fit_002 else np.nan

    d_200 = calculate_d_spacing(pos_200) if np.isfinite(pos_200) else np.nan
    d_002 = calculate_d_spacing(pos_002) if np.isfinite(pos_002) else np.nan

    # --- Compute Lattice Parameters a, c, and c/a ---
    if np.isfinite(d_001) and np.isfinite(d_002):
        c_val = 0.5 * (d_001 + 2.0 * d_002)
    elif np.isfinite(d_001):
        c_val = d_001
    elif np.isfinite(d_002):
        c_val = 2.0 * d_002
    else:
        c_val = np.nan

    a_from_200 = 2.0 * d_200 if np.isfinite(d_200) else np.nan
    if np.isfinite(d_111) and np.isfinite(c_val):
        term = (1.0 / (d_111 ** 2)) - (1.0 / (c_val ** 2))
        a_from_111 = np.sqrt(2.0 / term) if term > 0 else np.nan
    else:
        a_from_111 = np.nan

    if np.isfinite(a_from_200) and np.isfinite(a_from_111):
        a_val = 0.5 * (a_from_200 + a_from_111)
    elif np.isfinite(a_from_200):
        a_val = a_from_200
    elif np.isfinite(a_from_111):
        a_val = a_from_111
    else:
        a_val = np.nan

    c_over_a = (c_val / a_val) if (np.isfinite(c_val) and np.isfinite(a_val) and a_val > 0) else np.nan

    # --- Compute Chemical Order Parameter S ---
    if np.isfinite(c_over_a) and c_over_a <= 1.00:
        s_axial = np.sqrt(max(0.0, (1.0 - c_over_a) / (1.0 - C_OVER_A_ORDERED)))
        s_axial = min(s_axial, 1.0)
    elif np.isfinite(c_over_a) and c_over_a > 1.00:
        s_axial = 0.0
    else:
        s_axial = np.nan

    area_001 = fit_001["integrated_area"] if fit_001 else np.nan
    area_111 = fit_111["integrated_area"] if fit_111 else np.nan

    if np.isfinite(area_001) and np.isfinite(area_002) and area_002 > 0:
        ratio_meas = area_001 / area_002
        s_intensity_002 = np.sqrt(max(0.0, ratio_meas / 0.485))
        s_intensity_002 = min(s_intensity_002, 1.0)
    else:
        s_intensity_002 = np.nan

    if np.isfinite(area_001) and np.isfinite(area_111) and area_111 > 0:
        ratio_111_meas = area_001 / area_111
        s_intensity_111 = np.sqrt(max(0.0, ratio_111_meas / 0.090))
        s_intensity_111 = min(s_intensity_111, 1.0)
    else:
        s_intensity_111 = np.nan

    # --- Generate Multi-Panel Fitting Diagnostic Plot ---
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(f"XRD Peak Fitting & Lattice Analysis -- {base_name}", fontsize=12, fontweight="bold")

    # Panel 1: (001) region
    if fit_001:
        axes[0].plot(fit_001["xw"], fit_001["yw"], "k.", label="Observed", alpha=0.6)
        axes[0].plot(fit_001["xw"], fit_001["y_fit"], "r-", lw=1.8, label="Fit")
        axes[0].axvline(pos_001, color="tab:orange", ls="--", label=f"(001): {pos_001:.3f} deg")
        axes[0].set_title(f"(001) Superlattice\n2theta = {pos_001:.3f} deg (d = {d_001:.3f} A)")
    else:
        m = (x >= 22.5) & (x <= 25.5)
        axes[0].plot(x[m], y[m], "k.-", alpha=0.5)
        axes[0].set_title("(001) Region (no peak resolved)")
    axes[0].set_xlabel("2theta (deg)")
    axes[0].set_ylabel("Intensity (counts)")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=8, loc="upper right")

    # Panel 2: (111) region
    if fit_111:
        axes[1].plot(fit_111["xw"], fit_111["yw"], "k.", label="Observed", alpha=0.6)
        axes[1].plot(fit_111["xw"], fit_111["y_fit"], "b-", lw=1.8, label="Fit")
        axes[1].axvline(pos_111, color="tab:blue", ls="--", label=f"(111): {pos_111:.3f} deg")
        axes[1].set_title(f"(111) Fundamental\n2theta = {pos_111:.3f} deg (d = {d_111:.3f} A)")
    else:
        m = (x >= 40.0) & (x <= 42.5)
        axes[1].plot(x[m], y[m], "k.-", alpha=0.5)
        axes[1].set_title("(111) Region (no peak resolved)")
    axes[1].set_xlabel("2theta (deg)")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=8, loc="upper right")

    # Panel 3: (200) / (002) region
    if fit_doublet:
        axes[2].plot(fit_doublet["xw"], fit_doublet["yw"], "k.", label="Observed", alpha=0.6)
        axes[2].plot(fit_doublet["xw"], fit_doublet["y_fit"], "g-", lw=1.8, label="Fit")
        axes[2].axvline(pos_200, color="tab:purple", ls="--", label=f"(200): {pos_200:.2f} deg")
        axes[2].axvline(pos_002, color="tab:green", ls="--", label=f"(002): {pos_002:.2f} deg")
        axes[2].set_title(f"(200)/(002) Splitting\na = {a_val:.3f} A, c = {c_val:.3f} A\nc/a = {c_over_a:.3f}, S = {s_axial:.2f}")
    else:
        m = (x >= 46.0) & (x <= 49.5)
        axes[2].plot(x[m], y[m], "k.-", alpha=0.5)
        axes[2].set_title(f"a = {a_val:.3f} A, c = {c_val:.3f} A\nc/a = {c_over_a:.3f}, S = {s_axial:.2f}")
    axes[2].set_xlabel("2theta (deg)")
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(fontsize=8, loc="upper right")

    plt.tight_layout()
    fits_dir = os.path.join(output_dir, "fits")
    os.makedirs(fits_dir, exist_ok=True)
    out_png = os.path.join(fits_dir, f"{base_name}_peak_fits.png")
    out_svg = os.path.join(fits_dir, f"{base_name}_peak_fits.svg")
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_svg)
    plt.close(fig)

    return {
        "sample": base_name,
        "two_theta_001_deg": pos_001,
        "two_theta_110_deg": pos_110,
        "two_theta_111_deg": pos_111,
        "two_theta_200_deg": pos_200,
        "two_theta_002_deg": pos_002,
        "d_001_A": d_001,
        "d_110_A": d_110,
        "d_111_A": d_111,
        "d_200_A": d_200,
        "d_002_A": d_002,
        "lattice_a_A": a_val,
        "lattice_c_A": c_val,
        "ratio_c_over_a": c_over_a,
        "order_param_S_axial": s_axial,
        "order_param_S_int_002": s_intensity_002,
        "order_param_S_int_111": s_intensity_111,
        "area_001": area_001,
        "area_111": area_111,
        "area_002": area_002,
        "fit_plot_png": out_png,
    }


def run_batch_analysis(input_dir: str, output_dir: str) -> None:
    """Runs profile fitting and order parameter analysis across all .xy files in input_dir."""
    all_files = sorted(glob.glob(os.path.join(input_dir, "*.xy")))
    if not all_files:
        print(f"No .xy files found in {input_dir}")
        return

    os.makedirs(output_dir, exist_ok=True)
    results = []

    print("=" * 80)
    print("FePt QUANTITATIVE XRD ANALYSIS & CHEMICAL ORDER PARAMETER FITTING")
    print(f"Input Directory : {input_dir}")
    print(f"Output Directory: {output_dir}")
    print("=" * 80)

    selected_files = []
    for fp in all_files:
        fname = os.path.basename(fp)
        if "pre_rockin" in fname.lower() or "rockingcurve" in fname.lower():
            continue
        if "onlydata" in fname.lower():
            continue
        if fname.startswith("SiO_10mic_on_Si"):
            continue
        selected_files.append(fp)

    for fp in selected_files:
        fname = os.path.basename(fp)
        print(f"\nProcessing: {fname[:50]} ...")
        res = analyze_fept_file(fp, output_dir)
        results.append(res)
        pos001 = f"{res['two_theta_001_deg']:.2f} deg" if np.isfinite(res['two_theta_001_deg']) else "N/A"
        pos111 = f"{res['two_theta_111_deg']:.2f} deg" if np.isfinite(res['two_theta_111_deg']) else "N/A"
        la = f"{res['lattice_a_A']:.3f} A" if np.isfinite(res['lattice_a_A']) else "N/A"
        lc = f"{res['lattice_c_A']:.3f} A" if np.isfinite(res['lattice_c_A']) else "N/A"
        ca = f"{res['ratio_c_over_a']:.4f}" if np.isfinite(res['ratio_c_over_a']) else "N/A"
        sa = f"{res['order_param_S_axial']:.2f}" if np.isfinite(res['order_param_S_axial']) else "N/A"
        print(f"  (001) = {pos001} | (111) = {pos111} | a = {la} | c = {lc} | c/a = {ca} | S = {sa}")

    df = pd.DataFrame(results)
    cols_to_export = [
        "sample",
        "two_theta_001_deg",
        "two_theta_110_deg",
        "two_theta_111_deg",
        "two_theta_200_deg",
        "two_theta_002_deg",
        "lattice_a_A",
        "lattice_c_A",
        "ratio_c_over_a",
        "order_param_S_axial",
        "order_param_S_int_111",
        "order_param_S_int_002",
    ]
    summary_csv = os.path.join(output_dir, "fept_lattice_order_parameters.csv")
    df[cols_to_export].to_csv(summary_csv, index=False)
    print("\n" + "=" * 80)
    print(f"SUMMARY SAVED TO: {summary_csv}")
    print(df[cols_to_export].to_string(index=False))
    print("=" * 80)

    # Generate trend figure for sample series
    valid = df[df["ratio_c_over_a"].notna()].copy()
    if len(valid) >= 2:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        x_indices = np.arange(len(valid))
        sample_labels = [
            s.replace("SiO_10mic_+FePt_", "")
            .replace("_on_Si_011_Slit_2mm_s2p3dpmin_20to100_", "")
            .replace("_on_Si_011_Slit_2mm_20to80d_2p3dpmin", "")
            .replace("_on_Si_011_Slit_2mm_20to100d_2p3dpmin", "")
            .replace("_on_Si_011_Slit_2mm_20to90d_2p3dpmin", "")
            .replace("GENBB-Nat_FePt_SiO_Film_20_80_2q_BB_s2p1dpm_2mmSL-20250509-132210", "Continuous_Film")
            for s in valid["sample"]
        ]

        ax1.plot(x_indices, valid["lattice_a_A"], "s-", color="tab:blue", lw=2, label="a (A)")
        ax1.plot(x_indices, valid["lattice_c_A"], "o-", color="tab:orange", lw=2, label="c (A)")
        ax1.set_xticks(x_indices)
        ax1.set_xticklabels(sample_labels, rotation=25, ha="right", fontsize=9)
        ax1.set_ylabel("Lattice Parameter (A)", fontsize=11)
        ax1.set_title("Lattice Parameters a and c across FePt Series", fontsize=11, fontweight="bold")
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        ax2.plot(x_indices, valid["ratio_c_over_a"], "^-", color="tab:purple", lw=2, label="c/a ratio")
        ax2.axhline(0.964, color="gray", ls="--", label="Fully Ordered L10 (0.964)")
        ax2.axhline(1.000, color="gray", ls=":", label="Disordered A1 (1.000)")
        ax2_tw = ax2.twinx()
        ax2_tw.plot(x_indices, valid["order_param_S_axial"], "d--", color="tab:red", lw=2, label="Order Parameter S")
        ax2_tw.set_ylabel("Chemical Order Parameter S", color="tab:red", fontsize=11)
        ax2_tw.tick_params(axis="y", labelcolor="tab:red")
        ax2_tw.set_ylim(-0.05, 1.05)

        ax2.set_xticks(x_indices)
        ax2.set_xticklabels(sample_labels, rotation=25, ha="right", fontsize=9)
        ax2.set_ylabel("Axial Ratio c/a", color="tab:purple", fontsize=11)
        ax2.tick_params(axis="y", labelcolor="tab:purple")
        ax2.set_title("Tetragonality c/a and Chemical Order S", fontsize=11, fontweight="bold")
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        trend_png = os.path.join(output_dir, "fept_ordering_trends.png")
        trend_svg = os.path.join(output_dir, "fept_ordering_trends.svg")
        fig.savefig(trend_png, dpi=300)
        fig.savefig(trend_svg)
        plt.close(fig)
        print(f"Trend plot saved: {trend_png} / .svg")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Quantitative FePt XRD Analysis & Chemical Order Parameter.")
    default_in = os.path.join("MPI-IS", "data", "xrd_films")
    default_out = os.path.join("MPI-IS", "results")
    parser.add_argument("--input", default=default_in, help="Input directory containing .xy files.")
    parser.add_argument("--output", default=default_out, help="Output directory for CSV and figures.")
    args = parser.parse_args()
    run_batch_analysis(args.input, args.output)
