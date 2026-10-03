"""
I/O module for reading raw experimental and simulation files.

Follows strict non-negotiable rules:
- Never edits raw data.
- Never guesses or estimates values.
- Parses instrument headers directly from raw files.
"""

import os
import re
import numpy as np
from PIL import Image


def read_squid_dat(filepath):
    """
    Parse a Quantum Design MPMS3 VSM .dat file.
    
    Extracts Magnetic Field (Oe) converted to Tesla (Oe / 1e4)
    and Moment (emu).
    
    Parameters
    ----------
    filepath : str or Path
        Path to .dat file.
        
    Returns
    -------
    field_T : np.ndarray
        Applied field in Tesla.
    moment_emu : np.ndarray
        Measured moment in emu.
    header_info : dict
        Extracted header information.
    """
    with open(filepath, 'r', encoding='latin1') as f:
        lines = f.readlines()
        
    data_idx = -1
    header_info = {}
    for i, line in enumerate(lines):
        line_s = line.strip()
        if 'FILEOPENTIME' in line_s.replace(' ', '').upper():
            header_info['file_open_time'] = line_s
        elif '[Data]' in line_s:
            data_idx = i + 1
            break
            
    if data_idx == -1 or data_idx >= len(lines):
        raise ValueError(f"Could not find [Data] section in {filepath}")
        
    headers = [h.strip() for h in lines[data_idx].split(',')]
    field_col = None
    moment_col = None
    for j, h in enumerate(headers):
        if 'Magnetic Field (Oe)' in h:
            field_col = j
        elif 'Moment (emu)' in h:
            moment_col = j
            
    if field_col is None or moment_col is None:
        raise ValueError(f"Required columns not found in {filepath} headers: {headers}")
        
    data = []
    for line in lines[data_idx+1:]:
        parts = line.strip().split(',')
        if len(parts) > max(field_col, moment_col):
            try:
                f_val = float(parts[field_col])
                m_val = float(parts[moment_col])
                data.append((f_val, m_val))
            except ValueError:
                continue
                
    arr = np.array(data, dtype=float)
    if len(arr) == 0:
        raise ValueError(f"No valid numeric data found in {filepath}")
        
    field_T = arr[:, 0] / 1e4
    moment_emu = arr[:, 1]
    return field_T, moment_emu, header_info


def read_xrd_xy(filepath):
    """
    Parse a Rigaku SmartLab .xy diffractometer file.
    
    Extracts 2theta (deg), intensity (cps), and metadata header tags.
    
    Parameters
    ----------
    filepath : str or Path
        Path to .xy file.
        
    Returns
    -------
    two_theta : np.ndarray
        Diffraction angle 2theta in degrees.
    intensity : np.ndarray
        Measured intensity in counts per second (cps).
    metadata : dict
        Header metadata dictionary.
    """
    metadata = {}
    data = []
    
    with open(filepath, 'r', encoding='latin1') as f:
        for line in f:
            line_s = line.strip()
            if not line_s:
                continue
            if line_s.startswith('*'):
                # Metadata line: *TAG "value" or *TAG value
                match = re.match(r'^\*([A-Za-z0-9_\-]+)\s+"?([^"]*)"?$', line_s)
                if match:
                    tag, val = match.groups()
                    metadata[tag] = val
                else:
                    parts = line_s.split(None, 1)
                    tag = parts[0][1:]
                    val = parts[1].strip('"') if len(parts) > 1 else ''
                    metadata[tag] = val
            else:
                # Numeric data line: 2theta intensity
                parts = line_s.split()
                if len(parts) >= 2:
                    try:
                        th = float(parts[0])
                        inten = float(parts[1])
                        data.append((th, inten))
                    except ValueError:
                        continue
                        
    arr = np.array(data, dtype=float)
    two_theta = arr[:, 0]
    intensity = arr[:, 1]
    return two_theta, intensity, metadata


def read_sem_tiff_metadata(filepath):
    """
    Read Zeiss SmartSEM metadata from TIFF tag 34118.
    
    Parameters
    ----------
    filepath : str or Path
        Path to TIFF file.
        
    Returns
    -------
    meta : dict
        Dictionary of parsed SEM metadata including:
        - pixel_size_nm
        - mag
        - eht_kv
        - date
        - time
        - width_px
        - height_px
    """
    im = Image.open(filepath)
    width_px, height_px = im.size
    meta = {
        'filepath': filepath,
        'filename': os.path.basename(filepath),
        'width_px': width_px,
        'height_px': height_px,
        'pixel_size_nm': None,
        'mag': None,
        'eht_kv': None,
        'date': None,
        'time': None
    }
    
    if 34118 in im.tag_v2:
        val = im.tag_v2[34118]
        if isinstance(val, bytes):
            text = val.decode('latin1', errors='ignore')
        else:
            text = str(val)
            
        lines = text.split('\r\n')
        for i, line in enumerate(lines):
            line_str = line.strip()
            if line_str == 'AP_PIXEL_SIZE' and i + 1 < len(lines):
                # e.g. "Pixel Size = 32.13 nm"
                m = re.search(r'=\s*([0-9\.]+)\s*(nm|um|µm)?', lines[i+1])
                if m:
                    pval = float(m.group(1))
                    unit = m.group(2)
                    if unit in ('um', 'µm'):
                        pval *= 1000.0
                    meta['pixel_size_nm'] = pval
            elif line_str == 'AP_MAG' and i + 1 < len(lines):
                # e.g. "Mag = 3.48 K X"
                m = re.search(r'=\s*([0-9\.]+)\s*(K\s*X|X)?', lines[i+1], re.IGNORECASE)
                if m:
                    mval = float(m.group(1))
                    unit = m.group(2)
                    if unit and 'k' in unit.lower():
                        mval *= 1000.0
                    meta['mag'] = mval
            elif line_str in ('AP_ACTUALKV', 'AP_MANUALKV') and i + 1 < len(lines):
                # e.g. "EHT = 6.00 kV"
                m = re.search(r'=\s*([0-9\.]+)\s*kV', lines[i+1], re.IGNORECASE)
                if m and meta['eht_kv'] is None:
                    meta['eht_kv'] = float(m.group(1))
            elif line_str == 'AP_DATE' and i + 1 < len(lines):
                meta['date'] = lines[i+1].replace('Date:', '').strip()
            elif line_str == 'AP_TIME' and i + 1 < len(lines):
                meta['time'] = lines[i+1].replace('Time:', '').strip()
                
    return meta


def read_r1_loop_csv(filepath):
    """
    Read R1 micromagnetic hysteresis loop CSV file.
    
    Columns: B_ext_T, Mx, My, Mz
    """
    data = np.genfromtxt(filepath, delimiter=',', skip_header=1)
    B_ext_T = data[:, 0]
    Mx = data[:, 1]
    My = data[:, 2]
    Mz = data[:, 3]
    return B_ext_T, Mx, My, Mz


def read_r1_snapshot_npz(filepath):
    """
    Read R1 3D magnetization state snapshot (.npz).
    
    Contains 'm', 'norm', 'B', 'mz'.
    'm' is in A/m, and 'norm' is in A/m.
    Normalized magnetization: m / norm.
    """
    data = np.load(filepath)
    m = data['m']
    norm = data['norm']
    B = float(data['B'])
    mz = float(data['mz'])
    return {'m': m, 'norm': norm, 'B': B, 'mz': mz}
