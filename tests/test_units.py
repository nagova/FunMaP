"""
Unit tests for FunMaP core pipeline.
Designed for CI: runs on demo_data/ or in-memory fixtures.
Does not require external experimental raw data.
"""

import os
import sys
from pathlib import Path
import pytest
import numpy as np

# Ensure src is on sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from funmap.io import (
    read_squid_dat, read_xrd_xy, read_r1_loop_csv, read_r1_snapshot_npz
)
from funmap.squid import (
    split_branches, find_zero_crossing, interpolate_moment_at_field,
    compute_high_field_slope, process_loop, process_loop_tailfit
)
from funmap.stats import (
    run_ols_joint_regression, run_categorical_anova, compute_batch_summary
)
from funmap.sim import (
    compute_switching_field, parse_mesh_convergence_table
)


def test_split_branches():
    # Loop from 7 to -7 and back to 7
    h = np.array([7.0, 3.0, 0.0, -3.0, -7.0, -3.0, 0.0, 3.0, 7.0])
    m = np.array([1.0, 0.8, 0.6, -0.2, -1.0, -0.8, -0.6, 0.2, 1.0])
    desc_h, desc_m, asc_h, asc_m = split_branches(h, m)
    
    assert len(desc_h) == 5
    assert desc_h[0] == 7.0
    assert desc_h[-1] == -7.0
    assert len(asc_h) == 5
    assert asc_h[0] == -7.0
    assert asc_h[-1] == 7.0


def test_find_zero_crossing():
    # Linear crossing between x=1 and x=3, crossing at x=2
    h = np.array([3.0, 1.0])
    m = np.array([1.0, -1.0])
    h_cross = find_zero_crossing(h, m)
    assert np.isclose(h_cross, 2.0)


def test_interpolate_moment_at_field():
    h = np.array([5.0, 0.0, -5.0])
    m = np.array([1.0, 0.5, -1.0])
    m_zero = interpolate_moment_at_field(h, m, 0.0)
    assert np.isclose(m_zero, 0.5)


def test_process_loop_synthetic():
    h = np.linspace(7.0, -7.0, 141)
    h_full = np.concatenate([h, h[::-1]])
    hc = 1.15
    m_down = np.tanh((h + hc) / 0.5)
    m_up = np.tanh((h[::-1] - hc) / 0.5)
    m_full = np.concatenate([m_down, m_up]) * 1.2e-3
    
    metrics, corr = process_loop(h_full, m_full, chi_bg=0.0)
    assert np.isclose(metrics['H_c_T'], hc, atol=0.05)
    assert 'M_corr_7T_emu' in metrics
    assert 'Mr_50_ratio' in metrics


def test_ols_joint_regression():
    diameters = np.array([3, 5, 8, 10, 3, 5, 8, 10], dtype=float)
    hc = np.array([1.18, 1.13, 1.12, 1.14, 1.10, 1.05, 1.00, 1.08])
    batches = ['Batch 1']*4 + ['Batch 2']*4
    
    res = run_ols_joint_regression(diameters, hc, batches)
    assert 'slope_inv_r' in res
    assert 'p_value_slope' in res
    assert 'delta_hc_10_to_3' in res
    assert res['delta_hc_10_to_3'] > 0


def test_categorical_anova():
    diameters = np.array([3, 5, 8, 10, 3, 5, 8, 10], dtype=float)
    hc = np.array([1.18, 1.13, 1.12, 1.14, 1.10, 1.05, 1.00, 1.08])
    batches = ['Batch 1']*4 + ['Batch 2']*4
    
    res = run_categorical_anova(diameters, hc, batches)
    assert 'F_statistic' in res
    assert 'p_value' in res
    assert res['df_diam'] == 3
    assert res['df_res'] == 3


def test_read_demo_squid():
    demo_file = Path('demo_data/squid_dat/2025.06.Sample3.FePt.SiO2.7T.100Oe.s.dat')
    if demo_file.exists():
        h, m, meta = read_squid_dat(demo_file)
        assert len(h) > 0
        assert len(m) == len(h)
        assert 'file_open_time' in meta
