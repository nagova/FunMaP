import os
import glob
import numpy as np
import pandas as pd
from scipy.optimize import leastsq, minimize

LAMBDA = 1.54056  # Cu Ka1 (Angstrom)
C_A_ORDERED = 0.964
C_A_DISORDERED = 1.000

# Atomic scattering factors at sin(theta)/lambda for Cu Ka1
# Fe (Z=26), Pt (Z=78)
def atomic_factors(two_theta_deg):
    s = np.sin(np.radians(two_theta_deg / 2.0)) / LAMBDA
    # Analytical Cromer-Mann approx or effective table
    # At typical angles (2th ~ 24-50 deg, s ~ 0.13 - 0.27 A^-1):
    f_Pt = 78.0 * np.exp(-1.5 * s**2)
    f_Fe = 26.0 * np.exp(-2.0 * s**2)
    return f_Fe, f_Pt

def lorentz_polarization(two_theta_deg):
    th = np.radians(two_theta_deg / 2.0)
    return (1.0 + np.cos(2.0 * th)**2) / (np.sin(th)**2 * np.cos(th))

def theoretical_intensity(hkl, two_theta_deg):
    f_Fe, f_Pt = atomic_factors(two_theta_deg)
    lp = lorentz_polarization(two_theta_deg)
    if hkl == (0, 0, 1):
        F = 2.0 * (f_Pt - f_Fe)  # superlattice
        mult = 2
    elif hkl == (1, 1, 0):
        F = 2.0 * (f_Pt - f_Fe)  # superlattice
        mult = 4
    elif hkl == (1, 1, 1):
        F = 2.0 * (f_Pt + f_Fe)  # fundamental
        mult = 8
    elif hkl == (2, 0, 0):
        F = 2.0 * (f_Pt + f_Fe)  # fundamental
        mult = 4
    elif hkl == (0, 0, 2):
        F = 2.0 * (f_Pt + f_Fe)  # fundamental
        mult = 2
    else:
        F = 1.0
        mult = 1
    return mult * (F**2) * lp

def pseudo_voigt_height(x, x0, height, fwhm, eta):
    gamma = fwhm / 2.0
    dx = (x - x0) / max(gamma, 1e-6)
    lorentzian = 1.0 / (1.0 + dx**2)
    gaussian = np.exp(-np.log(2.0) * (dx**2))
    return height * (eta * lorentzian + (1.0 - eta) * gaussian)

def voigt_area(height, fwhm, eta):
    c_lorentz = np.pi / 2.0
    c_gauss = np.sqrt(np.pi / (4.0 * np.log(2.0)))
    return float(height * fwhm * (eta * c_lorentz + (1.0 - eta) * c_gauss))

def load_diffractogram(path):
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

def fit_single_peak(x, y, x_range, x0_guess, fwhm_guess=0.5):
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
        if h < 0 or fwhm < 0.06 or fwhm > 3.0 or eta < 0 or eta > 1:
            return 1e6 * np.ones_like(xw)
        if x0 < x_range[0] or x0 > x_range[1]:
            return 1e6 * np.ones_like(xw)
        y_m = pseudo_voigt_height(xw, x0, h, fwhm, eta) + c + s * (xw - x_range[0])
        return yw - y_m

    try:
        popt, ier = leastsq(residual, p0, full_output=False, maxfev=3000)
        if ier not in (1, 2, 3, 4):
            raise RuntimeError()
    except Exception:
        def loss(p):
            return float(np.sum(residual(p)**2))
        res = minimize(loss, p0, method="Nelder-Mead")
        popt = res.x

    x0, h, fwhm, eta, c, s = popt
    y_fit = pseudo_voigt_height(xw, x0, h, fwhm, eta) + c + s * (xw - x_range[0])
    ss_res = np.sum((yw - y_fit)**2)
    ss_tot = np.sum((yw - np.mean(yw))**2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    return {
        "x0": float(x0),
        "height": float(h),
        "fwhm": float(fwhm),
        "eta": float(eta),
        "area": voigt_area(h, fwhm, eta),
        "r2": float(r2),
        "xw": xw, "yw": yw, "y_fit": y_fit
    }

def fit_doublet_peak(x, y, x_range, x0_1_guess, x0_2_guess):
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
        if eta_1 < 0 or eta_1 > 1 or eta_2 < 0 or eta_2 > 1:
            return 1e6 * np.ones_like(xw)
        if x0_1 < x_range[0] or x0_1 > x_range[1] or x0_2 < x_range[0] or x0_2 > x_range[1]:
            return 1e6 * np.ones_like(xw)
        y_m = (
            pseudo_voigt_height(xw, x0_1, h_1, fwhm_1, eta_1)
            + pseudo_voigt_height(xw, x0_2, h_2, fwhm_2, eta_2)
            + c + s * (xw - x_range[0])
        )
        return yw - y_m

    try:
        popt, ier = leastsq(residual, p0, full_output=False, maxfev=4000)
        if ier not in (1, 2, 3, 4):
            raise RuntimeError()
    except Exception:
        def loss(p):
            return float(np.sum(residual(p)**2))
        res = minimize(loss, p0, method="Nelder-Mead")
        popt = res.x

    x0_1, h_1, fwhm_1, eta_1, x0_2, h_2, fwhm_2, eta_2, c, s = popt
    y_fit = (
        pseudo_voigt_height(xw, x0_1, h_1, fwhm_1, eta_1)
        + pseudo_voigt_height(xw, x0_2, h_2, fwhm_2, eta_2)
        + c + s * (xw - x_range[0])
    )
    ss_res = np.sum((yw - y_fit)**2)
    ss_tot = np.sum((yw - np.mean(yw))**2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    return {
        "x0_1": float(x0_1), "fwhm_1": float(fwhm_1), "area_1": voigt_area(h_1, fwhm_1, eta_1),
        "x0_2": float(x0_2), "fwhm_2": float(fwhm_2), "area_2": voigt_area(h_2, fwhm_2, eta_2),
        "r2": float(r2), "xw": xw, "yw": yw, "y_fit": y_fit
    }

def calc_d(two_theta):
    if not np.isfinite(two_theta):
        return np.nan
    th = np.radians(two_theta / 2.0)
    return float(LAMBDA / (2.0 * np.sin(th)))

def scherrer_size(two_theta_deg, fwhm_deg, K=0.94, b_inst=0.08):
    th = np.radians(two_theta_deg / 2.0)
    b_obs_rad = np.radians(fwhm_deg)
    b_inst_rad = np.radians(b_inst)
    if b_obs_rad > b_inst_rad:
        b_sample = np.sqrt(b_obs_rad**2 - b_inst_rad**2)
    else:
        b_sample = b_obs_rad
    return float((K * LAMBDA) / (b_sample * np.cos(th)))

data_dir = r"c:\Users\admin\Documents\Coding Projects\FunMaP\MPI-IS\data\xrd_films"
sub_file = os.path.join(data_dir, "SiO_10mic_on_Si_011_Slit_2mm_s2p3dpmin_20to100_.xy")
x_sub, y_sub = load_diffractogram(sub_file)

samples = [
    ("Continuous Film (GENBB)", "GENBB-Nat_FePt_SiO_Film_20_80_2q_BB_s2p1dpm_2mmSL-20250509-132210.xy", 0.0),
    ("FCC Reference (as-deposited)", "SiO_10mic_+FePt_FCC_on_Si_011_Slit_2mm_s2p3dpmin_20to100_.xy", 10.0),
    ("Annealed 500C 1k_3 (3um)", "SiO_10mic_+FePt_annealled_500dC_1k_3_on_Si_011_Slit_2mm_20to100d_2p3dpmin.xy", 3.0),
    ("Annealed 500C 1k_5 (5um)", "SiO_10mic_+FePt_annealled_500dC_1k_5_on_Si_011_Slit_2mm_20to100d_2p3dpmin.xy", 5.0),
    ("Annealed 500C 1k_8 (8um)", "SiO_10mic_+FePt_annealled_500dC_1k_8_on_Si_011_Slit_2mm_20to90d_2p3dpmin.xy", 8.0),
    ("Annealed 500C 1k_10 (10um)", "SiO_10mic_+FePt_annealled_500dC_1k_10_on_Si_011_Slit_2mm_20to80d_2p3dpmin.xy", 10.0),
]

full_results = []

for label, fname, diameter in samples:
    fp = os.path.join(data_dir, fname)
    if not os.path.exists(fp):
        continue
    x, y = load_diffractogram(fp)

    # Substrate normalization & background correction in (001) region
    # Interpolate substrate onto x grid
    y_sub_interp = np.interp(x, x_sub, y_sub)
    # Match scale around silica halo (21-22 deg)
    m_halo = (x >= 21.0) & (x <= 22.0)
    scale = np.median(y[m_halo]) / max(np.median(y_sub_interp[m_halo]), 1.0)
    y_subtracted = np.maximum(y - scale * y_sub_interp, 0.0)

    # Fit peaks on raw and background-subtracted data
    f001 = fit_single_peak(x, y, (22.8, 25.2), 23.8)
    f001_sub = fit_single_peak(x, y_subtracted, (22.8, 25.2), 23.8)
    f110 = fit_single_peak(x, y, (31.5, 34.5), 32.9)
    f111 = fit_single_peak(x, y, (40.0, 42.4), 41.1)
    fd = fit_doublet_peak(x, y, (46.0, 49.8), 47.1, 48.5)

    pos_001 = f001["x0"] if f001 else np.nan
    pos_110 = f110["x0"] if f110 else np.nan
    pos_111 = f111["x0"] if f111 else np.nan

    fwhm_001 = f001["fwhm"] if f001 else np.nan
    fwhm_111 = f111["fwhm"] if f111 else np.nan

    d001 = calc_d(pos_001)
    d110 = calc_d(pos_110)
    d111 = calc_d(pos_111)

    if fd and abs(fd["x0_2"] - fd["x0_1"]) > 0.3:
        pos_200 = min(fd["x0_1"], fd["x0_2"])
        pos_002 = max(fd["x0_1"], fd["x0_2"])
        fwhm_200 = fd["fwhm_1"] if fd["x0_1"] < fd["x0_2"] else fd["fwhm_2"]
        fwhm_002 = fd["fwhm_2"] if fd["x0_1"] < fd["x0_2"] else fd["fwhm_1"]
        area_200 = fd["area_1"] if fd["x0_1"] < fd["x0_2"] else fd["area_2"]
        area_002 = fd["area_2"] if fd["x0_1"] < fd["x0_2"] else fd["area_1"]
    else:
        f200 = fit_single_peak(x, y, (46.5, 47.8), 47.1)
        f002 = fit_single_peak(x, y, (47.8, 49.5), 48.5)
        pos_200 = f200["x0"] if f200 else np.nan
        pos_002 = f002["x0"] if f002 else np.nan
        fwhm_200 = f200["fwhm"] if f200 else np.nan
        fwhm_002 = f002["fwhm"] if f002 else np.nan
        area_200 = f200["area"] if f200 else np.nan
        area_002 = f002["area"] if f002 else np.nan

    d200 = calc_d(pos_200)
    d002 = calc_d(pos_002)

    # Tetragonal lattice parameters
    if np.isfinite(d001) and np.isfinite(d002):
        c_val = 0.5 * (d001 + 2.0 * d002)
    elif np.isfinite(d001):
        c_val = d001
    elif np.isfinite(d002):
        c_val = 2.0 * d002
    else:
        c_val = np.nan

    a_from_200 = 2.0 * d200 if np.isfinite(d200) else np.nan
    if np.isfinite(d111) and np.isfinite(c_val):
        term = (1.0 / (d111**2)) - (1.0 / (c_val**2))
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

    # Order parameters
    if np.isfinite(c_over_a) and c_over_a <= 1.00:
        s_axial = np.sqrt(max(0.0, (1.0 - c_over_a) / (1.0 - C_A_ORDERED)))
        s_axial = min(s_axial, 1.0)
    elif np.isfinite(c_over_a) and c_over_a > 1.00:
        s_axial = 0.0
    else:
        s_axial = np.nan

    area_001 = f001["area"] if f001 else np.nan
    area_111 = f111["area"] if f111 else np.nan

    # Theoretical intensities for S_int
    I_001_calc = theoretical_intensity((0,0,1), pos_001)
    I_111_calc = theoretical_intensity((1,1,1), pos_111)
    I_002_calc = theoretical_intensity((0,0,2), pos_002)

    ratio_calc_111 = I_001_calc / max(I_111_calc, 1e-6)
    ratio_calc_002 = I_001_calc / max(I_002_calc, 1e-6)

    s_int_111 = min(np.sqrt(max(0.0, (area_001 / area_111) / ratio_calc_111)), 1.0) if (area_001 > 0 and area_111 > 0) else np.nan
    s_int_002 = min(np.sqrt(max(0.0, (area_001 / area_002) / ratio_calc_002)), 1.0) if (area_001 > 0 and area_002 > 0) else np.nan

    # Crystallite size via Scherrer (111 and 001)
    D_111_nm = scherrer_size(pos_111, fwhm_111) / 10.0 if np.isfinite(pos_111) and np.isfinite(fwhm_111) else np.nan
    D_001_nm = scherrer_size(pos_001, fwhm_001) / 10.0 if np.isfinite(pos_001) and np.isfinite(fwhm_001) else np.nan
    D_200_nm = scherrer_size(pos_200, fwhm_200) / 10.0 if np.isfinite(pos_200) and np.isfinite(fwhm_200) else np.nan

    # Williamson-Hall microstrain & size estimation using multiple reflections (001, 111, 200, 002)
    peaks_wh = []
    for p_2th, p_fwhm in [(pos_001, fwhm_001), (pos_111, fwhm_111), (pos_200, fwhm_200), (pos_002, fwhm_002)]:
        if np.isfinite(p_2th) and np.isfinite(p_fwhm):
            th = np.radians(p_2th / 2.0)
            b = np.radians(p_fwhm)
            peaks_wh.append((4.0 * np.sin(th), b * np.cos(th)))

    if len(peaks_wh) >= 3:
        arr_wh = np.array(peaks_wh)
        x_w, y_w = arr_wh[:, 0], arr_wh[:, 1]
        x_bar, y_bar = np.mean(x_w), np.mean(y_w)
        slope = np.sum((x_w - x_bar) * (y_w - y_bar)) / max(np.sum((x_w - x_bar)**2), 1e-12)
        y_int = y_bar - slope * x_bar
        microstrain_eps = float(slope)
        wh_size_nm = float((0.94 * LAMBDA / max(y_int, 1e-5)) / 10.0) if y_int > 0 else np.nan
    else:
        microstrain_eps = np.nan
        wh_size_nm = np.nan

    # Texture / orientation index (I_001 / I_111)
    texture_001_111 = (area_001 / area_111) if (area_001 > 0 and area_111 > 0) else np.nan
    doublet_split_deg = (pos_002 - pos_200) if (np.isfinite(pos_002) and np.isfinite(pos_200)) else np.nan

    full_results.append({
        "Sample": label,
        "Sphere_Diam_um": diameter,
        "2th_001": round(pos_001, 3),
        "2th_110": round(pos_110, 3),
        "2th_111": round(pos_111, 3),
        "2th_200": round(pos_200, 3),
        "2th_002": round(pos_002, 3),
        "Doublet_Split_deg": round(doublet_split_deg, 3),
        "a_A": round(a_val, 4),
        "c_A": round(c_val, 4),
        "c_over_a": round(c_over_a, 4),
        "S_axial": round(s_axial, 3),
        "S_int_111": round(s_int_111, 3),
        "S_int_002": round(s_int_002, 3),
        "D_111_nm": round(D_111_nm, 1),
        "D_001_nm": round(D_001_nm, 1),
        "WH_Size_nm": round(wh_size_nm, 1) if np.isfinite(wh_size_nm) else None,
        "Microstrain_pct": round(microstrain_eps * 100, 3) if np.isfinite(microstrain_eps) else None,
        "Texture_001_to_111": round(texture_001_111, 3),
    })

df_res = pd.DataFrame(full_results)
out_csv = r"c:\Users\admin\Documents\Coding Projects\FunMaP\MPI-IS\results\fept_in_depth_xrd_analysis.csv"
df_res.to_csv(out_csv, index=False)
print("=== IN-DEPTH XRD ANALYSIS RESULTS ===")
print(df_res.to_string(index=False))
