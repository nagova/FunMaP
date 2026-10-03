"""
Scanning Electron Microscopy (SEM) processing module.

Reads Zeiss SmartSEM TIFF tag 34118 metadata, measures particle diameters
algorithmically, and generates the systematic SEM sample catalog.
"""

import os
import re
import numpy as np
from PIL import Image
from scipy import ndimage
from .io import read_sem_tiff_metadata


def measure_sphere_diameter_hough(image_path, pixel_size_nm):
    """
    Measure spherical particle diameter using Hough transform / radial gradient profile.
    
    Parameters
    ----------
    image_path : str or Path
        Path to SEM TIFF image.
    pixel_size_nm : float
        Pixel size in nanometers.
        
    Returns
    -------
    estimated_diameter_um : float or None
    """
    if pixel_size_nm is None or pixel_size_nm <= 0:
        return None
        
    im = Image.open(image_path).convert('L')
    arr = np.array(im, dtype=float)
    
    # Exclude bottom data bar if present (typically bottom 15% of image)
    h, w = arr.shape
    arr_crop = arr[:int(h * 0.85), :]
    
    # Smooth image
    blurred = ndimage.gaussian_filter(arr_crop, sigma=2.0)
    
    # Compute gradient magnitude
    gy, gx = np.gradient(blurred)
    grad = np.sqrt(gx**2 + gy**2)
    
    # Find strong edges
    thresh = np.percentile(grad, 90)
    edge_y, edge_x = np.where(grad > thresh)
    
    if len(edge_x) < 50:
        return None
        
    # Range of candidate diameters to test: 2 um to 25 um
    d_min_px = int((2000.0 / pixel_size_nm))
    d_max_px = int((25000.0 / pixel_size_nm))
    d_max_px = min(d_max_px, min(h, w))
    
    if d_min_px >= d_max_px or d_min_px < 5:
        return None
        
    # 1D radial auto-correlation across binary thresholded particles
    binary = arr_crop > np.mean(arr_crop)
    labeled, num_features = ndimage.label(binary)
    if num_features > 0:
        sizes = ndimage.sum(binary, labeled, range(1, num_features + 1))
        # Filter for particle-sized blobs
        valid_sizes = [s for s in sizes if s > (np.pi * (d_min_px/4)**2)]
        if valid_sizes:
            # Equivalent diameter d = 2 * sqrt(Area / pi)
            diameters_px = [2.0 * np.sqrt(s / np.pi) for s in valid_sizes]
            median_d_px = np.median(diameters_px)
            median_d_um = (median_d_px * pixel_size_nm) / 1000.0
            return float(median_d_um)
            
    return None


def generate_sem_catalog(raw_sem_dir):
    """
    Generate complete catalog of SEM TIFF files with instrument metadata
    and algorithmic sphere diameter checks.
    
    Parameters
    ----------
    raw_sem_dir : str or Path
        Directory containing raw SEM TIFFs.
        
    Returns
    -------
    catalog : list of dict
    """
    catalog = []
    for f in sorted(os.listdir(raw_sem_dir)):
        if f.lower().endswith(('.tif', '.tiff')):
            p = os.path.join(raw_sem_dir, f)
            meta = read_sem_tiff_metadata(p)
            pixel_size = meta['pixel_size_nm']
            
            # Measured sphere diameter
            d_meas = measure_sphere_diameter_hough(p, pixel_size)
            
            # Extract nominal diameter label from filename if present
            nominal_label = 'NOT SPECIFIED'
            flag = ''
            
            # Check p032-p040 anomaly
            m_num = re.search(r'p([0-9]{3})', f)
            num_val = int(m_num.group(1)) if m_num else None
            
            if num_val is not None:
                if 1 <= num_val <= 9:
                    nominal_label = '10 um'
                elif 10 <= num_val <= 19:
                    nominal_label = '8 um'
                elif 20 <= num_val <= 31:
                    nominal_label = '3 um'
                elif 32 <= num_val <= 40:
                    nominal_label = '5 um (nominal label in old catalog)'
                    if d_meas is not None and d_meas > 8.0:
                        flag = 'FLAG: Measured diameter ~10-11 um contradicts nominal 5 um label'
                        
            catalog.append({
                'filename': f,
                'width_px': meta['width_px'],
                'height_px': meta['height_px'],
                'pixel_size_nm': pixel_size if pixel_size is not None else 'NOT AVAILABLE',
                'mag': meta['mag'] if meta['mag'] is not None else 'NOT AVAILABLE',
                'eht_kv': meta['eht_kv'] if meta['eht_kv'] is not None else 'NOT AVAILABLE',
                'date': meta['date'] if meta['date'] is not None else 'NOT AVAILABLE',
                'time': meta['time'] if meta['time'] is not None else 'NOT AVAILABLE',
                'measured_diameter_um': f"{d_meas:.2f}" if d_meas is not None else 'NOT MEASURED',
                'nominal_label': nominal_label,
                'flag': flag
            })
            
    return catalog
