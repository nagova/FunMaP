"""
generate_raw_peak_fits.py
=========================
Generates raw peak-fit outputs from FePt XRD scans without downstream derivations.
Outputs:
  - FunMaP/MPI-IS/results/fept_peak_fits.csv
  - FunMaP/MPI-IS/results/fept_instrument_and_fit_metadata.csv
  - FunMaP/MPI-IS/results/fept_sample_log.csv
  - FunMaP/MPI-IS/results/fept_raw_patterns/<sample_id>.csv
"""

import os
import numpy as np
import pandas as pd
from scipy.optimize import leastsq

BASE_DIR = r"c:\Users\admin\Documents\Coding Projects\FunMaP"
DATA_DIR = os.path.join(BASE_DIR, "MPI-IS", "data", "xrd_films")
RESULTS_DIR = os.path.join(BASE_DIR, "MPI-IS", "results")
RAW_PATTERNS_DIR = os.path.join(RESULTS_DIR, "fept_raw_patterns")
os.makedirs(RAW_PATTERNS_DIR, exist_ok=True)

SAMPLES = [
    {
        "id": "film_flat_annealed",
        "file": "GENBB-Nat_FePt_SiO_Film_20_80_2q_BB_s2p1dpm_2mmSL-20250509-132210.xy",
        "sphere_diam_um": 0.0,
        "anneal_state": "annealed_500C",
        "deposited_thickness_nm": 60.0,
        "anneal_temp_C": 500,
        "anneal_time_min": 16.7,
        "substrate": "SiO2/Si(011) wafer",
        "is_same_physical_sample_as": "None",
        "step_size_deg": 0.0050,
        "dwell_time_s": 0.1429,
        "scan_range_two_theta_deg": "20.0-80.0",
        "notes": "Continuous flat FePt film control annealed at 500 C; Bragg-Brentano geometry"
    },
    {
        "id": "spheres_fcc_asdep",
        "file": "SiO_10mic_+FePt_FCC_on_Si_011_Slit_2mm_s2p3dpmin_20to100_.xy",
        "sphere_diam_um": 10.0,
        "anneal_state": "as_deposited",
        "deposited_thickness_nm": 60.0,
        "anneal_temp_C": 25,
        "anneal_time_min": 0,
        "substrate": "SiO2 microspheres (10 um) on Si(011) chip",
        "is_same_physical_sample_as": "spheres_1k_10 (pre-anneal state of same sphere batch)",
        "step_size_deg": 0.0100,
        "dwell_time_s": 0.2609,
        "scan_range_two_theta_deg": "20.0-100.0",
        "notes": "As-deposited FePt cap reference on 10 um spheres; chemically disordered A1 phase (S=0); serves as unannealed control"
    },
    {
        "id": "spheres_1k_3",
        "file": "SiO_10mic_+FePt_annealled_500dC_1k_3_on_Si_011_Slit_2mm_20to100d_2p3dpmin.xy",
        "sphere_diam_um": 3.0,
        "anneal_state": "annealed_500C",
        "deposited_thickness_nm": 60.0,
        "anneal_temp_C": 500,
        "anneal_time_min": 16.7,
        "substrate": "SiO2 microspheres (3 um) on Si(011) chip",
        "is_same_physical_sample_as": "None",
        "step_size_deg": 0.0100,
        "dwell_time_s": 0.2609,
        "scan_range_two_theta_deg": "20.0-100.0",
        "notes": "FePt caps on 3 um SiO2 microspheres annealed at 500 C (1000 s anneal)"
    },
    {
        "id": "spheres_1k_5",
        "file": "SiO_10mic_+FePt_annealled_500dC_1k_5_on_Si_011_Slit_2mm_20to100d_2p3dpmin.xy",
        "sphere_diam_um": 5.0,
        "anneal_state": "annealed_500C",
        "deposited_thickness_nm": 60.0,
        "anneal_temp_C": 500,
        "anneal_time_min": 16.7,
        "substrate": "SiO2 microspheres (5 um) on Si(011) chip",
        "is_same_physical_sample_as": "None",
        "step_size_deg": 0.0100,
        "dwell_time_s": 0.2609,
        "scan_range_two_theta_deg": "20.0-100.0",
        "notes": "FePt caps on 5 um SiO2 microspheres annealed at 500 C (1000 s anneal)"
    },
    {
        "id": "spheres_1k_8",
        "file": "SiO_10mic_+FePt_annealled_500dC_1k_8_on_Si_011_Slit_2mm_20to90d_2p3dpmin.xy",
        "sphere_diam_um": 8.0,
        "anneal_state": "annealed_500C",
        "deposited_thickness_nm": 60.0,
        "anneal_temp_C": 500,
        "anneal_time_min": 16.7,
        "substrate": "SiO2 microspheres (8 um) on Si(011) chip",
        "is_same_physical_sample_as": "None",
        "step_size_deg": 0.0100,
        "dwell_time_s": 0.2609,
        "scan_range_two_theta_deg": "20.0-90.0",
        "notes": "FePt caps on 8 um SiO2 microspheres annealed at 500 C (1000 s anneal)"
    },
    {
        "id": "spheres_1k_10",
        "file": "SiO_10mic_+FePt_annealled_500dC_1k_10_on_Si_011_Slit_2mm_20to80d_2p3dpmin.xy",
        "sphere_diam_um": 10.0,
        "anneal_state": "annealed_500C",
        "deposited_thickness_nm": 60.0,
        "anneal_temp_C": 500,
        "anneal_time_min": 16.7,
        "substrate": "SiO2 microspheres (10 um) on Si(011) chip",
        "is_same_physical_sample_as": "spheres_fcc_asdep (post-anneal counterpart of same sphere batch)",
        "step_size_deg": 0.0100,
        "dwell_time_s": 0.2609,
        "scan_range_two_theta_deg": "20.0-80.0",
        "notes": "FePt caps on 10 um SiO2 microspheres annealed at 500 C (1000 s anneal)"
    },
]

def load_xy(p):
    rows = []
    with open(p, 'r', encoding='utf-8', errors='ignore') as f:
        for l in f:
            if l.startswith(('#', ';', '*', '/')) or not l.strip():
                continue
            parts = l.replace(',', ' ').split()
            if len(parts) >= 2:
                try:
                    rows.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    pass
    arr = np.array(rows, dtype=float)
    idx = np.argsort(arr[:, 0])
    return arr[idx, 0], arr[idx, 1]

def pseudo_voigt(x, x0, h, fwhm, eta):
    gamma = fwhm / 2.0
    dx = (x - x0) / max(gamma, 1e-6)
    lorentz = 1.0 / (1.0 + dx**2)
    gauss = np.exp(-np.log(2.0) * (dx**2))
    return h * (eta * lorentz + (1.0 - eta) * gauss)

def voigt_area(h, fwhm, eta):
    c_lorentz = np.pi / 2.0
    c_gauss = np.sqrt(np.pi / (4.0 * np.log(2.0)))
    return float(h * fwhm * (eta * c_lorentz + (1.0 - eta) * c_gauss))

def pure_inv(A):
    """Pure Python Gauss-Jordan matrix inversion avoiding broken BLAS DLL."""
    n = len(A)
    M = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(A)]
    for i in range(n):
        max_row = max(range(i, n), key=lambda r: abs(M[r][i]))
        if abs(M[max_row][i]) < 1e-14:
            raise np.linalg.LinAlgError("Singular matrix")
        M[i], M[max_row] = M[max_row], M[i]
        pivot = M[i][i]
        for j in range(2 * n):
            M[i][j] /= pivot
        for r in range(n):
            if r != i:
                factor = M[r][i]
                for j in range(2 * n):
                    M[r][j] -= factor * M[i][j]
    return [row[n:] for row in M]

def get_cov_and_errors(res_func, popt, N):
    P = len(popt)
    eps = 1e-5
    J = []
    for j in range(P):
        p_plus = list(popt)
        p_plus[j] += eps
        p_minus = list(popt)
        p_minus[j] -= eps
        col = [(res_func(p_plus)[k] - res_func(p_minus)[k]) / (2.0 * eps) for k in range(N)]
        J.append(col)
    JT_J = [[sum(J[i][k] * J[j][k] for k in range(N)) for j in range(P)] for i in range(P)]
    try:
        cov = pure_inv(JT_J)
        resid = res_func(popt)
        s_sq = sum(r**2 for r in resid) / max(N - P, 1)
        perr = [np.sqrt(max(0.0, cov[i][i] * s_sq)) for i in range(P)]
        cov_scaled = [[cov[i][j] * s_sq for j in range(P)] for i in range(P)]
        return cov_scaled, perr, s_sq
    except Exception:
        return None, [None]*P, None

def fit_single_reflection(x, y, x_range, x0_guess, fwhm_guess=0.5):
    m = (x >= x_range[0]) & (x <= x_range[1])
    xw, yw = x[m], y[m]
    N = len(xw)
    if N < 8:
        return None
    bg_c = float(np.min(yw))
    bg_s = float((yw[-1] - yw[0]) / (xw[-1] - xw[0])) if xw[-1] != xw[0] else 0.0
    h0 = max(float(np.max(yw) - bg_c), 1.0)
    p0 = [x0_guess, h0, fwhm_guess, 0.5, bg_c, bg_s]

    def res(p):
        x0, h, fwhm, eta, c, s = p
        pen = 0.0
        if h < 0: pen += 1e5 * (abs(h) + 1)**2
        if fwhm < 0.04: pen += 1e5 * ((0.04 - fwhm)*100)**2
        if fwhm > 3.0: pen += 1e5 * ((fwhm - 3.0)*100)**2
        if eta < 0.0: pen += 1e5 * (abs(eta)*100)**2
        if eta > 1.0: pen += 1e5 * ((eta - 1.0)*100)**2
        if x0 < x_range[0]: pen += 1e5 * ((x_range[0] - x0)*100)**2
        if x0 > x_range[1]: pen += 1e5 * ((x0 - x_range[1])*100)**2
        
        ym = pseudo_voigt(xw, x0, h, fwhm, eta) + c + s * (xw - x_range[0])
        r = yw - ym
        if pen > 0:
            r = r + np.sqrt(pen / N)
        return r

    popt, ier = leastsq(res, p0, full_output=False, maxfev=6000)
    cov, perr, s_sq = get_cov_and_errors(res, popt, N)
    
    resid = res(popt)
    noise = float(np.std(resid))
    snr = float(popt[1] / max(noise, 1e-6))
    area = voigt_area(popt[1], popt[2], popt[3])
    
    area_err = None
    if cov is not None:
        h_val, f_val, e_val = popt[1], popt[2], popt[3]
        c_l = np.pi / 2.0
        c_g = np.sqrt(np.pi / (4.0 * np.log(2.0)))
        dA_dh = f_val * (e_val * c_l + (1.0 - e_val) * c_g)
        dA_df = h_val * (e_val * c_l + (1.0 - e_val) * c_g)
        dA_de = h_val * f_val * (c_l - c_g)
        grad = [0.0, dA_dh, dA_df, dA_de, 0.0, 0.0]
        var_A = sum(grad[i] * cov[i][j] * grad[j] for i in range(6) for j in range(6))
        area_err = float(np.sqrt(max(0.0, var_A)))
        
    ym_final = pseudo_voigt(xw, popt[0], popt[1], popt[2], popt[3]) + popt[4] + popt[5] * (xw - x_range[0])
    chi2 = float(np.sum((yw - ym_final)**2) / max(N - 6, 1))
    bg_at_peak = float(popt[4] + popt[5] * (popt[0] - x_range[0]))
    
    base_linear = yw[0] + (yw[-1] - yw[0]) * (xw - xw[0]) / (xw[-1] - xw[0])
    meas_area_trapz = float(np.trapezoid(yw - base_linear, xw))

    return {
        'x0': float(popt[0]), 'x0_err': perr[0],
        'h': float(popt[1]), 'h_err': perr[1],
        'fwhm': float(popt[2]), 'fwhm_err': perr[2],
        'eta': float(popt[3]),
        'area': area, 'area_err': area_err,
        'bg_at_peak': bg_at_peak,
        'chi2': chi2,
        'n_points': N,
        'snr': snr,
        'meas_area_trapz': meas_area_trapz,
        'noise': noise
    }

def fit_doublet_reflections(x, y, x_range=(46.0, 49.8), x0_1_guess=47.25, x0_2_guess=48.5):
    m = (x >= x_range[0]) & (x <= x_range[1])
    xw, yw = x[m], y[m]
    N = len(xw)
    bg_c = float(np.min(yw))
    bg_s = float((yw[-1] - yw[0]) / (xw[-1] - xw[0])) if xw[-1] != xw[0] else 0.0
    h_max = float(np.max(yw) - bg_c)
    
    p0 = [x0_1_guess, h_max, 0.1, 0.5, x0_2_guess, h_max * 0.01, 0.5, 0.5, bg_c, bg_s]
    
    def res(p):
        x0_1, h1, fwhm1, eta1, x0_2, h2, fwhm2, eta2, c, s = p
        pen = 0.0
        if h1 < 0: pen += 1e6 * (abs(h1) + 1)**2
        if h2 < 0: pen += 1e6 * (abs(h2) + 1)**2
        if fwhm1 < 0.02: pen += 1e6 * ((0.02 - fwhm1)*100)**2
        if fwhm2 < 0.05: pen += 1e6 * ((0.05 - fwhm2)*100)**2
        if fwhm1 > 3.0: pen += 1e6 * ((fwhm1 - 3.0)*100)**2
        if fwhm2 > 3.0: pen += 1e6 * ((fwhm2 - 3.0)*100)**2
        if x0_1 < 46.5 or x0_1 > 47.8: pen += 1e6 * 10
        if x0_2 < 47.8 or x0_2 > 49.5: pen += 1e6 * 10
        ym = (pseudo_voigt(xw, x0_1, h1, fwhm1, eta1)
              + pseudo_voigt(xw, x0_2, h2, fwhm2, eta2)
              + c + s * (xw - x_range[0]))
        r = yw - ym
        if pen > 0:
            r = r + np.sqrt(pen / N)
        return r

    popt, ier = leastsq(res, p0, full_output=False, maxfev=10000)
    cov, perr, s_sq = get_cov_and_errors(res, popt, N)
    area1 = voigt_area(popt[1], popt[2], popt[3])
    area2 = voigt_area(popt[5], popt[6], popt[7])
    
    def calc_area_err(h_idx, f_idx, e_idx):
        if cov is None: return None
        h_val, f_val, e_val = popt[h_idx], popt[f_idx], popt[e_idx]
        c_l = np.pi / 2.0
        c_g = np.sqrt(np.pi / (4.0 * np.log(2.0)))
        dA_dh = f_val * (e_val * c_l + (1.0 - e_val) * c_g)
        dA_df = h_val * (e_val * c_l + (1.0 - e_val) * c_g)
        dA_de = h_val * f_val * (c_l - c_g)
        grad = [0.0]*10
        grad[h_idx] = dA_dh
        grad[f_idx] = dA_df
        grad[e_idx] = dA_de
        var = sum(grad[i] * cov[i][j] * grad[j] for i in range(10) for j in range(10))
        return float(np.sqrt(max(0.0, var)))

    a1_err = calc_area_err(1, 2, 3)
    a2_err = calc_area_err(5, 6, 7)
    
    ym_final = (pseudo_voigt(xw, popt[0], popt[1], popt[2], popt[3])
                + pseudo_voigt(xw, popt[4], popt[5], popt[6], popt[7])
                + popt[8] + popt[9] * (xw - x_range[0]))
    chi2 = float(np.sum((yw - ym_final)**2) / max(N - 10, 1))

    bg1 = float(popt[8] + popt[9] * (popt[0] - x_range[0]))
    bg2 = float(popt[8] + popt[9] * (popt[4] - x_range[0]))

    return {
        'peak1': {
            'x0': float(popt[0]), 'x0_err': perr[0],
            'h': float(popt[1]), 'h_err': perr[1],
            'fwhm': float(popt[2]), 'fwhm_err': perr[2],
            'eta': float(popt[3]),
            'area': area1, 'area_err': a1_err,
            'bg_at_peak': bg1,
        },
        'peak2': {
            'x0': float(popt[4]), 'x0_err': perr[4],
            'h': float(popt[5]), 'h_err': perr[5],
            'fwhm': float(popt[6]), 'fwhm_err': perr[6],
            'eta': float(popt[7]),
            'area': area2, 'area_err': a2_err,
            'bg_at_peak': bg2,
        },
        'chi2': chi2,
        'n_points': N,
        'fit_range': f"{x_range[0]}-{x_range[1]}",
        'bg_c': float(popt[8]),
        'bg_s': float(popt[9]),
    }

def main():
    peak_fit_rows = []

    for s in SAMPLES:
        sid = s["id"]
        fp = os.path.join(DATA_DIR, s["file"])
        x, y = load_xy(fp)
        
        # 1. Export raw pattern to CSV
        raw_df = pd.DataFrame({"two_theta_deg": x, "intensity_counts": y})
        raw_out = os.path.join(RAW_PATTERNS_DIR, f"{sid}.csv")
        raw_df.to_csv(raw_out, index=False)
        print(f"Exported raw pattern: {raw_out}")
        
        # 2. Reflection 001 [22.8 - 25.2]
        f001 = fit_single_reflection(x, y, (22.8, 25.2), 23.9)
        if sid == "spheres_fcc_asdep":
            peak_fit_rows.append({
                "sample_id": sid,
                "sphere_diam_um": s["sphere_diam_um"],
                "anneal_state": s["anneal_state"],
                "reflection": "001",
                "two_theta_fit_deg": f"{f001['x0']:.4f}",
                "two_theta_err_deg": "NOT_AVAILABLE",
                "fwhm_obs_deg": f"{f001['fwhm']:.4f}",
                "fwhm_err_deg": "NOT_AVAILABLE",
                "integrated_area_counts_deg": f"{f001['meas_area_trapz']:.2f}",
                "integrated_area_err": "NOT_AVAILABLE",
                "peak_height_counts": f"{f001['h']:.1f}",
                "background_at_peak_counts": f"{f001['bg_at_peak']:.1f}",
                "profile_function": "NOT_FITTED: forbidden in disordered fcc A1 (signal is amorphous SiO2 halo)",
                "background_model": "linear",
                "fit_range_two_theta_deg": "22.8-25.2",
                "rwp_or_chi2": f"{f001['chi2']:.2f}",
                "n_points_in_fit": f001["n_points"]
            })
        elif sid in ("spheres_1k_8", "spheres_1k_10") and f001['snr'] < 1.0:
            peak_fit_rows.append({
                "sample_id": sid,
                "sphere_diam_um": s["sphere_diam_um"],
                "anneal_state": s["anneal_state"],
                "reflection": "001",
                "two_theta_fit_deg": f"{f001['x0']:.4f}",
                "two_theta_err_deg": "NOT_AVAILABLE",
                "fwhm_obs_deg": f"{f001['fwhm']:.4f}",
                "fwhm_err_deg": "NOT_AVAILABLE",
                "integrated_area_counts_deg": f"{f001['meas_area_trapz']:.2f}",
                "integrated_area_err": "NOT_AVAILABLE",
                "peak_height_counts": f"{f001['h']:.1f}",
                "background_at_peak_counts": f"{f001['bg_at_peak']:.1f}",
                "profile_function": "NOT_FITTED: submerged in amorphous SiO2 halo and background noise (SNR < 1)",
                "background_model": "linear",
                "fit_range_two_theta_deg": "22.8-25.2",
                "rwp_or_chi2": f"{f001['chi2']:.2f}",
                "n_points_in_fit": f001["n_points"]
            })
        else:
            two_th_err = f"{f001['x0_err']:.4f}" if f001['x0_err'] is not None else "NOT_AVAILABLE"
            fwhm_err = f"{f001['fwhm_err']:.4f}" if f001['fwhm_err'] is not None else "NOT_AVAILABLE"
            area_err = f"{f001['area_err']:.2f}" if f001['area_err'] is not None else "NOT_AVAILABLE"
            peak_fit_rows.append({
                "sample_id": sid,
                "sphere_diam_um": s["sphere_diam_um"],
                "anneal_state": s["anneal_state"],
                "reflection": "001",
                "two_theta_fit_deg": f"{f001['x0']:.4f}",
                "two_theta_err_deg": two_th_err,
                "fwhm_obs_deg": f"{f001['fwhm']:.4f}",
                "fwhm_err_deg": fwhm_err,
                "integrated_area_counts_deg": f"{f001['area']:.2f}",
                "integrated_area_err": area_err,
                "peak_height_counts": f"{f001['h']:.1f}",
                "background_at_peak_counts": f"{f001['bg_at_peak']:.1f}",
                "profile_function": "pseudo_voigt",
                "background_model": "linear",
                "fit_range_two_theta_deg": "22.8-25.2",
                "rwp_or_chi2": f"{f001['chi2']:.2f}",
                "n_points_in_fit": f001["n_points"]
            })

        # 3. Reflection 110 [31.5 - 34.5]
        f110 = fit_single_reflection(x, y, (31.5, 34.5), 32.9)
        if sid == "spheres_fcc_asdep" or (sid in ("spheres_1k_8", "spheres_1k_10") and f110['snr'] < 1.5):
            reason = "forbidden in disordered fcc A1 (signal within background noise)" if sid == "spheres_fcc_asdep" else "submerged in background noise (SNR < 1.5)"
            peak_fit_rows.append({
                "sample_id": sid,
                "sphere_diam_um": s["sphere_diam_um"],
                "anneal_state": s["anneal_state"],
                "reflection": "110",
                "two_theta_fit_deg": f"{f110['x0']:.4f}",
                "two_theta_err_deg": "NOT_AVAILABLE",
                "fwhm_obs_deg": f"{f110['fwhm']:.4f}",
                "fwhm_err_deg": "NOT_AVAILABLE",
                "integrated_area_counts_deg": f"{f110['meas_area_trapz']:.2f}",
                "integrated_area_err": "NOT_AVAILABLE",
                "peak_height_counts": f"{f110['h']:.1f}",
                "background_at_peak_counts": f"{f110['bg_at_peak']:.1f}",
                "profile_function": f"NOT_FITTED: {reason}",
                "background_model": "linear",
                "fit_range_two_theta_deg": "31.5-34.5",
                "rwp_or_chi2": f"{f110['chi2']:.2f}",
                "n_points_in_fit": f110["n_points"]
            })
        else:
            two_th_err = f"{f110['x0_err']:.4f}" if f110['x0_err'] is not None else "NOT_AVAILABLE"
            fwhm_err = f"{f110['fwhm_err']:.4f}" if f110['fwhm_err'] is not None else "NOT_AVAILABLE"
            area_err = f"{f110['area_err']:.2f}" if f110['area_err'] is not None else "NOT_AVAILABLE"
            peak_fit_rows.append({
                "sample_id": sid,
                "sphere_diam_um": s["sphere_diam_um"],
                "anneal_state": s["anneal_state"],
                "reflection": "110",
                "two_theta_fit_deg": f"{f110['x0']:.4f}",
                "two_theta_err_deg": two_th_err,
                "fwhm_obs_deg": f"{f110['fwhm']:.4f}",
                "fwhm_err_deg": fwhm_err,
                "integrated_area_counts_deg": f"{f110['area']:.2f}",
                "integrated_area_err": area_err,
                "peak_height_counts": f"{f110['h']:.1f}",
                "background_at_peak_counts": f"{f110['bg_at_peak']:.1f}",
                "profile_function": "pseudo_voigt",
                "background_model": "linear",
                "fit_range_two_theta_deg": "31.5-34.5",
                "rwp_or_chi2": f"{f110['chi2']:.2f}",
                "n_points_in_fit": f110["n_points"]
            })

        # 4. Reflection 111 [40.0 - 42.0]
        f111 = fit_single_reflection(x, y, (40.0, 42.0), 41.1)
        if sid != "film_flat_annealed" and (f111['h'] <= 1.0 or f111['snr'] < 1.5):
            peak_fit_rows.append({
                "sample_id": sid,
                "sphere_diam_um": s["sphere_diam_um"],
                "anneal_state": s["anneal_state"],
                "reflection": "111",
                "two_theta_fit_deg": f"{f111['x0']:.4f}",
                "two_theta_err_deg": "NOT_AVAILABLE",
                "fwhm_obs_deg": f"{f111['fwhm']:.4f}",
                "fwhm_err_deg": "NOT_AVAILABLE",
                "integrated_area_counts_deg": f"{f111['meas_area_trapz']:.2f}",
                "integrated_area_err": "NOT_AVAILABLE",
                "peak_height_counts": f"{f111['h']:.1f}",
                "background_at_peak_counts": f"{f111['bg_at_peak']:.1f}",
                "profile_function": "NOT_FITTED: peak below detection limit in sphere cap geometry / obscured by rising substrate background",
                "background_model": "linear",
                "fit_range_two_theta_deg": "40.0-42.0",
                "rwp_or_chi2": f"{f111['chi2']:.2f}",
                "n_points_in_fit": f111["n_points"]
            })
        else:
            two_th_err = f"{f111['x0_err']:.4f}" if f111['x0_err'] is not None else "NOT_AVAILABLE"
            fwhm_err = f"{f111['fwhm_err']:.4f}" if f111['fwhm_err'] is not None else "NOT_AVAILABLE"
            area_err = f"{f111['area_err']:.2f}" if f111['area_err'] is not None else "NOT_AVAILABLE"
            peak_fit_rows.append({
                "sample_id": sid,
                "sphere_diam_um": s["sphere_diam_um"],
                "anneal_state": s["anneal_state"],
                "reflection": "111",
                "two_theta_fit_deg": f"{f111['x0']:.4f}",
                "two_theta_err_deg": two_th_err,
                "fwhm_obs_deg": f"{f111['fwhm']:.4f}",
                "fwhm_err_deg": fwhm_err,
                "integrated_area_counts_deg": f"{f111['area']:.2f}",
                "integrated_area_err": area_err,
                "peak_height_counts": f"{f111['h']:.1f}",
                "background_at_peak_counts": f"{f111['bg_at_peak']:.1f}",
                "profile_function": "pseudo_voigt",
                "background_model": "linear",
                "fit_range_two_theta_deg": "40.0-42.0",
                "rwp_or_chi2": f"{f111['chi2']:.2f}",
                "n_points_in_fit": f111["n_points"]
            })

        # 5 & 6. Doublet (200) and (002) in [46.0 - 49.8]
        fd = fit_doublet_reflections(x, y, (46.0, 49.8), 47.25, 48.5)
        p1 = fd['peak1']
        p2 = fd['peak2']
        
        # 200 row
        p1_th_err = f"{p1['x0_err']:.4f}" if p1['x0_err'] is not None and p1['x0_err'] < 10.0 else "NOT_AVAILABLE"
        p1_fwhm_err = f"{p1['fwhm_err']:.4f}" if p1['fwhm_err'] is not None and p1['fwhm_err'] < 10.0 else "NOT_AVAILABLE"
        p1_area_err = f"{p1['area_err']:.2f}" if p1['area_err'] is not None and p1['area_err'] < 1e7 else "NOT_AVAILABLE"
        
        p1_note = "pseudo_voigt (WARNING: feature at 47.30 deg is dominated by single-crystal Si(220) substrate reflection, >2.5M counts)"
        peak_fit_rows.append({
            "sample_id": sid,
            "sphere_diam_um": s["sphere_diam_um"],
            "anneal_state": s["anneal_state"],
            "reflection": "200",
            "two_theta_fit_deg": f"{p1['x0']:.4f}",
            "two_theta_err_deg": p1_th_err,
            "fwhm_obs_deg": f"{p1['fwhm']:.4f}",
            "fwhm_err_deg": p1_fwhm_err,
            "integrated_area_counts_deg": f"{p1['area']:.2f}",
            "integrated_area_err": p1_area_err,
            "peak_height_counts": f"{p1['h']:.1f}",
            "background_at_peak_counts": f"{p1['bg_at_peak']:.1f}",
            "profile_function": p1_note,
            "background_model": "linear (shared doublet window)",
            "fit_range_two_theta_deg": fd["fit_range"],
            "rwp_or_chi2": f"{fd['chi2']:.2f}",
            "n_points_in_fit": fd["n_points"]
        })

        # 002 row
        is_p2_unphysical = (p2['x0_err'] is None or p2['x0_err'] > 5.0 or p2['fwhm'] <= 0.03 or p2['fwhm'] > 2.5 or p2['h'] < 10.0 or p2['area'] < 0)
        if is_p2_unphysical:
            peak_fit_rows.append({
                "sample_id": sid,
                "sphere_diam_um": s["sphere_diam_um"],
                "anneal_state": s["anneal_state"],
                "reflection": "002",
                "two_theta_fit_deg": f"{p2['x0']:.4f}",
                "two_theta_err_deg": "NOT_AVAILABLE",
                "fwhm_obs_deg": f"{abs(p2['fwhm']):.4f}",
                "fwhm_err_deg": "NOT_AVAILABLE",
                "integrated_area_counts_deg": f"{p2['area']:.2f}",
                "integrated_area_err": "NOT_AVAILABLE",
                "peak_height_counts": f"{p2['h']:.1f}",
                "background_at_peak_counts": f"{p2['bg_at_peak']:.1f}",
                "profile_function": "NOT_FITTED: not genuinely resolved; fit unconstrained onto Si(220) substrate tail or unphysical",
                "background_model": "linear (shared doublet window)",
                "fit_range_two_theta_deg": fd["fit_range"],
                "rwp_or_chi2": f"{fd['chi2']:.2f}",
                "n_points_in_fit": fd["n_points"]
            })
        else:
            p2_th_err = f"{p2['x0_err']:.4f}" if p2['x0_err'] is not None and p2['x0_err'] < 10.0 else "NOT_AVAILABLE"
            p2_fwhm_err = f"{p2['fwhm_err']:.4f}" if p2['fwhm_err'] is not None and p2['fwhm_err'] < 10.0 else "NOT_AVAILABLE"
            p2_area_err = f"{p2['area_err']:.2f}" if p2['area_err'] is not None and p2['area_err'] < 1e7 else "NOT_AVAILABLE"
            peak_fit_rows.append({
                "sample_id": sid,
                "sphere_diam_um": s["sphere_diam_um"],
                "anneal_state": s["anneal_state"],
                "reflection": "002",
                "two_theta_fit_deg": f"{p2['x0']:.4f}",
                "two_theta_err_deg": p2_th_err,
                "fwhm_obs_deg": f"{p2['fwhm']:.4f}",
                "fwhm_err_deg": p2_fwhm_err,
                "integrated_area_counts_deg": f"{p2['area']:.2f}",
                "integrated_area_err": p2_area_err,
                "peak_height_counts": f"{p2['h']:.1f}",
                "background_at_peak_counts": f"{p2['bg_at_peak']:.1f}",
                "profile_function": "pseudo_voigt (WARNING: weak feature on descending tail of Si(220) peak)",
                "background_model": "linear (shared doublet window)",
                "fit_range_two_theta_deg": fd["fit_range"],
                "rwp_or_chi2": f"{fd['chi2']:.2f}",
                "n_points_in_fit": fd["n_points"]
            })

    # Save File 1
    df_peak_fits = pd.DataFrame(peak_fit_rows)
    csv_peak_fits = os.path.join(RESULTS_DIR, "fept_peak_fits.csv")
    df_peak_fits.to_csv(csv_peak_fits, index=False)
    print(f"Generated {csv_peak_fits} with {len(df_peak_fits)} rows.")

    # Save File 2: fept_instrument_and_fit_metadata.csv
    metadata_rows = []
    for s in SAMPLES:
        metadata_rows.append({
            "sample_id": s["id"],
            "wavelength_used_A": 1.54059,
            "k_alpha2_stripped": "False",
            "mono_or_filter": "Cu_K-beta filter, no monochromator",
            "instrumental_fwhm_deg_at_20": "NOT_MEASURED",
            "instrumental_fwhm_deg_at_50": "NOT_MEASURED",
            "instrumental_fwhm_source": "NOT_MEASURED",
            "scherrer_k_used_previously": 0.94,
            "step_size_deg": s["step_size_deg"],
            "dwell_time_s": s["dwell_time_s"],
            "scan_range_two_theta_deg": s["scan_range_two_theta_deg"],
            "zero_offset_correction_deg": "NOT_APPLIED",
            "sample_displacement_correction_applied": "NOT_APPLIED",
            "standard_used_for_instrument_function": "NOT_MEASURED"
        })
    df_metadata = pd.DataFrame(metadata_rows)
    csv_metadata = os.path.join(RESULTS_DIR, "fept_instrument_and_fit_metadata.csv")
    df_metadata.to_csv(csv_metadata, index=False)
    print(f"Generated {csv_metadata} with {len(df_metadata)} rows.")

    # Save File 3: fept_sample_log.csv
    sample_log_rows = []
    for s in SAMPLES:
        sample_log_rows.append({
            "sample_id": s["id"],
            "sphere_diam_um": s["sphere_diam_um"],
            "deposited_thickness_nm": s["deposited_thickness_nm"],
            "anneal_temp_C": s["anneal_temp_C"],
            "anneal_time_min": s["anneal_time_min"],
            "substrate": s["substrate"],
            "is_same_physical_sample_as": s["is_same_physical_sample_as"],
            "notes": s["notes"]
        })
    df_sample_log = pd.DataFrame(sample_log_rows)
    csv_sample_log = os.path.join(RESULTS_DIR, "fept_sample_log.csv")
    df_sample_log.to_csv(csv_sample_log, index=False)
    print(f"Generated {csv_sample_log} with {len(df_sample_log)} rows.")

    print("\n=== VERIFYING PROSE & SELF-CHECKS ===")
    # Self-Check 1: Area(001) / Area(110)
    print("\n--- Self-Check 1: Superlattice ratio Area(001)/Area(110) ---")
    for s in SAMPLES:
        sid = s["id"]
        r001 = df_peak_fits[(df_peak_fits['sample_id'] == sid) & (df_peak_fits['reflection'] == '001')].iloc[0]
        r110 = df_peak_fits[(df_peak_fits['sample_id'] == sid) & (df_peak_fits['reflection'] == '110')].iloc[0]
        a001 = float(r001['integrated_area_counts_deg'])
        a110 = float(r110['integrated_area_counts_deg'])
        ratio = a001 / a110 if abs(a110) > 1e-3 else np.nan
        print(f"  {sid}: Area(001)={a001:.2f}, Area(110)={a110:.2f}, Ratio={ratio:.4f}")

    # Self-Check 2: Area(001) / Area(002) vs powder limit 1.87
    print("\n--- Self-Check 2: Order-parameter sanity check Area(001)/Area(002) vs 1.87 ---")
    for s in SAMPLES:
        sid = s["id"]
        r001 = df_peak_fits[(df_peak_fits['sample_id'] == sid) & (df_peak_fits['reflection'] == '001')].iloc[0]
        r002 = df_peak_fits[(df_peak_fits['sample_id'] == sid) & (df_peak_fits['reflection'] == '002')].iloc[0]
        a001 = float(r001['integrated_area_counts_deg'])
        a002 = float(r002['integrated_area_counts_deg'])
        ratio = a001 / a002 if abs(a002) > 1e-3 else np.nan
        print(f"  {sid}: Area(001)={a001:.2f}, Area(002)={a002:.2f}, Ratio={ratio:.4f}")

    # Self-Check 3: Width consistency (FWHM of fundamentals)
    print("\n--- Self-Check 3: Width consistency & Williamson-Hall slope of fundamentals (111, 200, 002) ---")
    for s in SAMPLES:
        sid = s["id"]
        sub = df_peak_fits[df_peak_fits['sample_id'] == sid]
        fwhms = {r['reflection']: float(r['fwhm_obs_deg']) if r['fwhm_obs_deg'] != 'NOT_AVAILABLE' else np.nan for _, r in sub.iterrows()}
        positions = {r['reflection']: float(r['two_theta_fit_deg']) if r['two_theta_fit_deg'] != 'NOT_AVAILABLE' else np.nan for _, r in sub.iterrows()}
        print(f"  {sid} FWHMs: 001={fwhms.get('001')}, 110={fwhms.get('110')}, 111={fwhms.get('111')}, 200={fwhms.get('200')}, 002={fwhms.get('002')}")
        
        # WH slope on fundamentals only: (111), (200), (002)
        wh_pts = []
        for refl in ['111', '200', '002']:
            tth = positions.get(refl)
            fw = fwhms.get(refl)
            if tth is not None and fw is not None and np.isfinite(tth) and np.isfinite(fw) and fw > 0:
                th = np.radians(tth / 2.0)
                x_val = 4.0 * np.sin(th)
                y_val = np.radians(fw) * np.cos(th)
                wh_pts.append((x_val, y_val))
        if len(wh_pts) >= 2:
            xs = [p[0] for p in wh_pts]
            ys = [p[1] for p in wh_pts]
            mx = sum(xs) / len(xs)
            my = sum(ys) / len(ys)
            denom = sum((x_i - mx)**2 for x_i in xs)
            numer = sum((x_i - mx)*(y_i - my) for x_i, y_i in zip(xs, ys))
            slope = numer / denom if denom > 1e-12 else 0.0
            print(f"    WH slope (fundamentals only): {slope:.6f} ({'NEGATIVE' if slope < 0 else 'POSITIVE/FLAT'})")
        else:
            print(f"    WH slope: NOT ENOUGH POINTS")

    # Self-Check 4: Resolvability check of (200) and (002)
    print("\n--- Self-Check 4: Resolvability check of (200) and (002) ---")
    for s in SAMPLES:
        sid = s["id"]
        r200 = df_peak_fits[(df_peak_fits['sample_id'] == sid) & (df_peak_fits['reflection'] == '200')].iloc[0]
        r002 = df_peak_fits[(df_peak_fits['sample_id'] == sid) & (df_peak_fits['reflection'] == '002')].iloc[0]
        fw200 = float(r200['fwhm_obs_deg'])
        fw002 = float(r002['fwhm_obs_deg'])
        tth200 = float(r200['two_theta_fit_deg'])
        tth002 = float(r002['two_theta_fit_deg'])
        split = tth002 - tth200
        resolved = "YES" if (fw200 < split and fw002 < split and fw002 < 1.2) else "NO (NOT RESOLVED / ARTIFACT)"
        print(f"  {sid}: 2th_200={tth200:.3f} (FWHM={fw200:.3f}), 2th_002={tth002:.3f} (FWHM={fw002:.3f}), Split={split:.3f} deg -> Resolved: {resolved}")

    print("\nALL RUNS COMPLETE!")

if __name__ == '__main__':
    main()
