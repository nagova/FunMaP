"""
FunMaP: Functional Magnetic Particles Analysis Pipeline

A reproducible research pipeline for FePt thin films on spherical SiO2 substrates.
"""

__version__ = "1.0.0"
__author__ = "Natalia Gonzalez-Vazquez, Andrew K. Schulz, Gunther Richter"

from .io import read_squid_dat, read_xrd_xy, read_sem_tiff_metadata, read_r1_loop_csv, read_r1_snapshot_npz
from .squid import process_loop, process_loop_tailfit, split_branches, CHI_BG_DEFAULT
from .stats import run_ols_joint_regression, run_categorical_anova, compute_batch_summary
from .xrd import get_theoretical_reflections, parse_xrd_metadata_table
from .sem import generate_sem_catalog, measure_sphere_diameter_hough
from .sim import compute_switching_field, parse_mesh_convergence_table, extract_midplane_slice_y0
from .style import apply_publication_style, save_figure, COLORS
