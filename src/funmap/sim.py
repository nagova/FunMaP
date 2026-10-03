"""
Micromagnetic simulation processing module.

Processes R1 benchmark loops, sub-exchange mesh convergence,
and 3D magnetization state snapshots.
"""

import os
import glob
import numpy as np
from .io import read_r1_loop_csv, read_r1_snapshot_npz


def compute_switching_field(B_ext_T, Mz):
    """
    Compute nucleation threshold H_c by linear interpolation of Mz across zero.
    
    For R1 benchmark caps, switching occurs between -6.20 T and -6.30 T.
    
    Returns
    -------
    H_c_T : float
        Interpolated switching field.
    bracket : tuple of float
        (B_pre, B_post) field steps straddling switching.
    """
    for i in range(len(B_ext_T) - 1):
        if (Mz[i] >= 0 and Mz[i+1] <= 0) or (Mz[i] <= 0 and Mz[i+1] >= 0):
            b0, b1 = B_ext_T[i], B_ext_T[i+1]
            m0, m1 = Mz[i], Mz[i+1]
            if m1 == m0:
                hc = abs(b0)
            else:
                hc = abs(b0 + (0.0 - m0) * (b1 - b0) / (m1 - m0))
            bracket = (min(abs(b0), abs(b1)), max(abs(b0), abs(b1)))
            return float(hc), bracket
            
    return np.nan, (np.nan, np.nan)


def parse_mesh_convergence_table(sim_loops_dir):
    """
    Parse mesh convergence loop files and calculate interpolated H_c.
    
    Expected cell sizes: 1.0, 1.5, 2.0, 3.0, 4.0 nm.
    
    Returns
    -------
    table : list of dict
    """
    pattern = os.path.join(sim_loops_dir, 'R1_hemi_d200nm_cell_*.csv')
    files = glob.glob(pattern)
    table = []
    
    # Material constants
    Ku = 6.6e6       # J/m^3
    Ms = 1.0e6       # A/m
    A = 10.0e-12     # J/m
    mu0 = 4.0 * np.pi * 1e-7
    l_ex = np.sqrt(2.0 * A / (mu0 * Ms**2)) * 1e9   # 3.9894 nm -> ~3.99 nm
    delta_0 = np.sqrt(A / Ku) * 1e9                 # 1.2309 nm -> ~1.23 nm
    
    for f in sorted(files):
        # Extract cell size from filename: e.g. R1_hemi_d200nm_cell_1p0nm.csv
        basename = os.path.basename(f)
        cell_str = basename.replace('R1_hemi_d200nm_cell_', '').replace('.csv', '').replace('nm', '')
        cell_nm = float(cell_str.replace('p', '.'))
        
        b, mx, my, mz = read_r1_loop_csv(f)
        hc, bracket = compute_switching_field(b, mz)
        
        # Bounding box 320 x 320 x 160 nm
        nx = int(np.round(320.0 / cell_nm))
        ny = int(np.round(320.0 / cell_nm))
        nz = int(np.round(160.0 / cell_nm))
        total_cells = nx * ny * nz
        
        # Remanence at zero field
        idx_zero = np.argmin(np.abs(b))
        mz_zero = float(mz[idx_zero])
        
        table.append({
            'filename': basename,
            'cell_size_nm': cell_nm,
            'mesh_size_nm': cell_nm,
            'cell_over_delta0': cell_nm / delta_0,
            'cell_over_lex': cell_nm / l_ex,
            'nx': nx,
            'ny': ny,
            'nz': nz,
            'total_bounding_cells': total_cells,
            'Hc_interpolated_T': hc,
            'switching_field_T': hc,
            'remanence_mz': mz_zero,
            'bracket_lower_T': bracket[0],
            'bracket_upper_T': bracket[1],
            'field_range_min_T': float(np.min(b)),
            'field_range_max_T': float(np.max(b))
        })
        
    # Sort by cell size ascending
    table.sort(key=lambda x: x['cell_size_nm'])
    return table


def extract_midplane_slice_y0(snapshot_npz_path):
    """
    Extract midplane cross-section (y = 0) from 3D snapshot NPZ.
    
    Normalizes m by norm: m_norm = m / norm.
    
    Returns
    -------
    slice_data : dict
        Contains x_coords, z_coords, mx, my, mz, norm, mask
    """
    snap = read_r1_snapshot_npz(snapshot_npz_path)
    m = snap['m']          # (nx, ny, nz, 3)
    norm = snap['norm']    # (nx, ny, nz, 1)
    
    nx, ny, nz, _ = m.shape
    # Midplane index in y: ny // 2
    y_idx = ny // 2
    
    m_y0 = m[:, y_idx, :, :]
    norm_y0 = norm[:, y_idx, :, 0]
    
    # Active cells mask
    active_mask = norm_y0 > 0.0
    
    # Normalized components
    mx = np.zeros_like(norm_y0)
    my = np.zeros_like(norm_y0)
    mz = np.zeros_like(norm_y0)
    
    mx[active_mask] = m_y0[:, :, 0][active_mask] / norm_y0[active_mask]
    my[active_mask] = m_y0[:, :, 1][active_mask] / norm_y0[active_mask]
    mz[active_mask] = m_y0[:, :, 2][active_mask] / norm_y0[active_mask]
    
    return {
        'B_ext_T': snap['B'],
        'mz_mean': snap['mz'],
        'mx': mx,
        'my': my,
        'mz': mz,
        'norm': norm_y0,
        'mask': active_mask,
        'shape': (nx, nz)
    }
