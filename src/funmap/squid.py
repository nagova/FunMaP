"""
SQUID magnetometry processing module.

Implements exact standardized baseline subtraction, branch separation,
coercivity zero-crossing, remanence plateau evaluation, and high-field slope fitting.
"""

import numpy as np


CHI_BG_DEFAULT = 1.33e-4  # emu/T empirical positive linear baseline


def split_branches(field_T, moment):
    """
    Split a hysteresis loop into descending and ascending branches.
    
    Descending: start -> index of minimum field.
    Ascending: index of minimum field -> end.
    
    Parameters
    ----------
    field_T : np.ndarray
        Applied field in Tesla.
    moment : np.ndarray
        Moment values (raw or corrected).
        
    Returns
    -------
    desc_H, desc_M, asc_H, asc_M : np.ndarray
    """
    min_idx = np.argmin(field_T)
    desc_H = field_T[:min_idx + 1]
    desc_M = moment[:min_idx + 1]
    asc_H = field_T[min_idx:]
    asc_M = moment[min_idx:]
    return desc_H, desc_M, asc_H, asc_M


def find_zero_crossing(h, m):
    """
    Find zero-moment field crossing by linear interpolation.
    """
    for i in range(len(m) - 1):
        if (m[i] >= 0 and m[i+1] <= 0) or (m[i] <= 0 and m[i+1] >= 0):
            if m[i+1] == m[i]:
                return h[i]
            h_cross = h[i] + (0.0 - m[i]) * (h[i+1] - h[i]) / (m[i+1] - m[i])
            return h_cross
    return np.nan


def interpolate_moment_at_field(h, m, target_h):
    """
    Linearly interpolate moment at a given field target_h within a branch.
    """
    if len(h) < 2:
        return np.nan
        
    # Check if branch is decreasing or increasing in field
    if h[0] > h[-1]:  # Descending
        idx = np.searchsorted(-h, -target_h)
    else:  # Ascending
        idx = np.searchsorted(h, target_h)
        
    if idx == 0:
        return m[0]
    if idx >= len(h):
        return m[-1]
        
    h0, h1 = h[idx - 1], h[idx]
    m0, m1 = m[idx - 1], m[idx]
    if h1 == h0:
        return m0
    return m0 + (target_h - h0) * (m1 - m0) / (h1 - h0)


def compute_high_field_slope(field_T, moment, h_min=4.0, h_max=7.0):
    """
    Linear fit of raw M vs H on descending branch between h_min and h_max.
    """
    desc_H, desc_M, _, _ = split_branches(field_T, moment)
    mask = (desc_H >= h_min) & (desc_H <= h_max)
    if np.sum(mask) < 2:
        return np.nan
    p = np.polyfit(desc_H[mask], desc_M[mask], 1)
    return p[0]  # slope in emu/T


def process_loop(field_T, moment_raw, chi_bg=CHI_BG_DEFAULT):
    """
    Process a single hysteresis loop with the standard fixed-baseline pipeline.
    
    Parameters
    ----------
    field_T : np.ndarray
        Applied field in Tesla.
    moment_raw : np.ndarray
        Raw moment in emu.
    chi_bg : float
        Baseline slope to subtract (emu/T).
        
    Returns
    -------
    metrics : dict
        Calculated metrics matching Table 1.
    corrected_data : dict
        Corrected field and moments for plotting.
    """
    # Baseline correction: M_c = M - chi_bg * H
    moment_corr = moment_raw - chi_bg * field_T
    
    # Split branches
    desc_H, desc_M_c, asc_H, asc_M_c = split_branches(field_T, moment_corr)
    
    # M_corr(7 T) = mean |M_c| over points with |H| > 6.9 T
    mask_7T = np.abs(field_T) > 6.9
    M_corr_7T = np.mean(np.abs(moment_corr[mask_7T]))
    
    # Coercivity: zero crossing on descending and ascending branches
    hc_desc = abs(find_zero_crossing(desc_H, desc_M_c))
    hc_asc = abs(find_zero_crossing(asc_H, asc_M_c))
    H_c = (hc_desc + hc_asc) / 2.0
    
    # Remanence at zero field: mean of desc M_c(0) and -asc M_c(0), / M_corr
    mr_zero_desc = interpolate_moment_at_field(desc_H, desc_M_c, 0.0)
    mr_zero_asc = interpolate_moment_at_field(asc_H, asc_M_c, 0.0)
    Mr_zero = (mr_zero_desc - mr_zero_asc) / 2.0
    Mr_zero_ratio = Mr_zero / M_corr_7T
    
    # Remanence at +50 mT: mean of desc M_c(+0.05) and -asc M_c(-0.05), / M_corr
    mr_50_desc = interpolate_moment_at_field(desc_H, desc_M_c, 0.05)
    mr_50_asc = interpolate_moment_at_field(asc_H, asc_M_c, -0.05)
    Mr_50 = (mr_50_desc - mr_50_asc) / 2.0
    Mr_50_ratio = Mr_50 / M_corr_7T
    
    # High-field slope on raw moment (4 to 7 T)
    dmdh = compute_high_field_slope(field_T, moment_raw, 4.0, 7.0)
    
    metrics = {
        'H_c_T': H_c,
        'H_c_desc_T': hc_desc,
        'H_c_asc_T': hc_asc,
        'M_corr_7T_emu': M_corr_7T,
        'Mr_zero_emu': Mr_zero,
        'Mr_zero_ratio': Mr_zero_ratio,
        'Mr_50_emu': Mr_50,
        'Mr_50_ratio': Mr_50_ratio,
        'dmdh_emu_per_T': dmdh
    }
    
    corrected_data = {
        'field_T': field_T,
        'moment_raw': moment_raw,
        'moment_corr': moment_corr,
        'moment_norm': moment_corr / M_corr_7T,
        'desc_H': desc_H,
        'desc_M_corr': desc_M_c,
        'asc_H': asc_H,
        'asc_M_corr': asc_M_c
    }
    
    return metrics, corrected_data


def process_loop_tailfit(field_T, moment_raw):
    """
    Process a loop using the tail-fit pipeline (Table S4).
    
    Fits the high-field slope (4 to 7 T) on the loop itself,
    and subtracts that fitted slope.
    """
    chi_tail = compute_high_field_slope(field_T, moment_raw, 4.0, 7.0)
    return process_loop(field_T, moment_raw, chi_bg=chi_tail)
