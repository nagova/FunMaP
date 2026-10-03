"""
X-ray diffraction (XRD) processing and theoretical line calculation module.

Rigaku SmartLab Bragg-Brentano symmetric scan analysis.
Computes exact theoretical Bragg angles from physical wavelengths and lattice constants.
"""

import os
import numpy as np
from .io import read_xrd_xy


# Physical emission lines in Angstroms
LAMBDA_CU_KA = 1.540598    # Cu K-alpha weighted average
LAMBDA_CU_KB = 1.392250    # Cu K-beta
LAMBDA_W_LA1 = 1.476390    # Tungsten L-alpha 1
LAMBDA_W_LA2 = 1.487430    # Tungsten L-alpha 2
LAMBDA_W_LB1 = 1.281810    # Tungsten L-beta 1

# Standard lattice parameters in Angstroms
A_SI = 5.4309              # Silicon diamond cubic

# Stoichiometric L1_0 FePt (tetragonal)
A_FEPT_L10 = 3.852
C_FEPT_L10 = 3.713

# Disordered A1 FePt (FCC cubic)
A_FEPT_A1 = 3.820


def bragg_2theta_deg(d_spacing_angstrom, wavelength_angstrom):
    """
    Compute Bragg angle 2theta in degrees from d-spacing and wavelength.
    lambda = 2 d sin(theta) -> 2theta = 2 arcsin(lambda / (2d))
    """
    sin_th = wavelength_angstrom / (2.0 * d_spacing_angstrom)
    if sin_th > 1.0:
        return np.nan
    theta_rad = np.arcsin(sin_th)
    return float(np.degrees(2.0 * theta_rad))


def get_theoretical_reflections():
    """
    Calculate theoretical 2theta reflections for all candidate phases and lines.
    
    Returns
    -------
    lines : dict
        Dictionary of line names, theoretical 2theta (deg), and physical description.
    """
    # Si (220): d = a / sqrt(8)
    d_si220 = A_SI / np.sqrt(8.0)
    
    # FePt L1_0 d-spacings
    d_l10_001 = C_FEPT_L10
    d_l10_110 = A_FEPT_L10 / np.sqrt(2.0)
    d_l10_111 = 1.0 / np.sqrt(2.0 / (A_FEPT_L10**2) + 1.0 / (C_FEPT_L10**2))
    d_l10_200 = A_FEPT_L10 / 2.0
    d_l10_002 = C_FEPT_L10 / 2.0
    
    # FePt A1 d-spacings
    d_a1_111 = A_FEPT_A1 / np.sqrt(3.0)
    d_a1_200 = A_FEPT_A1 / 2.0
    
    lines = {
        # Substrate and instrumental lines diffracted by Si(220)
        'Si(220) [Cu Ka]': {
            '2theta': bragg_2theta_deg(d_si220, LAMBDA_CU_KA),
            'type': 'substrate',
            'desc': 'Primary substrate reflection of Si(110) wafer (47.30 deg)'
        },
        'Si(220) [Cu Kb]': {
            '2theta': bragg_2theta_deg(d_si220, LAMBDA_CU_KB),
            'type': 'instrumental',
            'desc': 'Cu K-beta line diffracted by Si(220) (42.52 deg)'
        },
        'Si(220) [W La1]': {
            '2theta': bragg_2theta_deg(d_si220, LAMBDA_W_LA1),
            'type': 'instrumental',
            'desc': 'Tungsten cathode contamination L-alpha 1 diffracted by Si(220) (45.18 deg)'
        },
        'Si(220) [W La2]': {
            '2theta': bragg_2theta_deg(d_si220, LAMBDA_W_LA2),
            'type': 'instrumental',
            'desc': 'Tungsten cathode contamination L-alpha 2 diffracted by Si(220) (45.55 deg)'
        },
        # L1_0 FePt thin film lines
        'FePt L1_0 (001) [Cu Ka]': {
            '2theta': bragg_2theta_deg(d_l10_001, LAMBDA_CU_KA),
            'type': 'L10_superlattice',
            'desc': 'Diagnostic chemical order superlattice reflection (23.95 deg)'
        },
        'FePt L1_0 (110) [Cu Ka]': {
            '2theta': bragg_2theta_deg(d_l10_110, LAMBDA_CU_KA),
            'type': 'L10_superlattice',
            'desc': 'Chemical order superlattice reflection (32.86 deg)'
        },
        'FePt L1_0 (111) [Cu Ka]': {
            '2theta': bragg_2theta_deg(d_l10_111, LAMBDA_CU_KA),
            'type': 'L10_fundamental',
            'desc': 'Strongest fundamental reflection of L1_0 FePt (41.08 deg)'
        },
        'FePt L1_0 (200) [Cu Ka]': {
            '2theta': bragg_2theta_deg(d_l10_200, LAMBDA_CU_KA),
            'type': 'L10_fundamental',
            'desc': 'Tetragonal fundamental reflection (47.15 deg)'
        },
        'FePt L1_0 (002) [Cu Ka]': {
            '2theta': bragg_2theta_deg(d_l10_002, LAMBDA_CU_KA),
            'type': 'L10_fundamental',
            'desc': 'Tetragonal fundamental reflection (49.03 deg)'
        },
        # A1 FePt thin film lines
        'FePt A1 (111) [Cu Ka]': {
            '2theta': bragg_2theta_deg(d_a1_111, LAMBDA_CU_KA),
            'type': 'A1_fundamental',
            'desc': 'Disordered cubic fundamental reflection (40.85 deg)'
        },
        'FePt A1 (200) [Cu Ka]': {
            '2theta': bragg_2theta_deg(d_a1_200, LAMBDA_CU_KA),
            'type': 'A1_fundamental',
            'desc': 'Disordered cubic fundamental reflection (47.57 deg)'
        }
    }
    return lines


def parse_xrd_metadata_table(raw_xrd_dir):
    """
    Parse metadata headers for all .xy files in raw directory.
    Returns list of metadata records.
    """
    records = []
    for root, dirs, files in os.walk(raw_xrd_dir):
        for f in sorted(files):
            if f.endswith('.xy'):
                path = os.path.join(root, f)
                _, _, meta = read_xrd_xy(path)
                records.append({
                    'filename': f,
                    'start_time': meta.get('MEAS_SCAN_START_TIME', 'NOT AVAILABLE'),
                    'end_time': meta.get('MEAS_SCAN_END_TIME', 'NOT AVAILABLE'),
                    'step_deg': meta.get('MEAS_SCAN_STEP', 'NOT AVAILABLE'),
                    'speed_deg_per_min': meta.get('MEAS_SCAN_SPEED', 'NOT AVAILABLE'),
                    'wave_type': meta.get('MEAS_COND_XG_WAVE_TYPE', 'NOT AVAILABLE'),
                    'axis_x': meta.get('MEAS_SCAN_AXIS_X', 'NOT AVAILABLE'),
                    'data_count': meta.get('MEAS_DATA_COUNT', 'NOT AVAILABLE')
                })
    return records
