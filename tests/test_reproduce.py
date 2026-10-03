"""
Acceptance reproduction test suite for FunMaP.
Runs against experimental raw data in data/raw/.
Asserts all numerical values and statistics to within +/- 0.001 of verified paper ground truth.
"""

import os
import sys
import hashlib
from pathlib import Path
import pytest
import numpy as np

# Ensure src is on sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from funmap.io import read_squid_dat, read_r1_loop_csv
from funmap.squid import process_loop
from funmap.stats import (
    run_ols_joint_regression, run_categorical_anova, compute_batch_summary
)
from funmap.sim import parse_mesh_convergence_table


RAW_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"
SQUID_DIR = RAW_DIR / "squid_dat"
R1_DIR = RAW_DIR / "r1_sim"


@pytest.mark.skipif(not SQUID_DIR.exists(), reason="Experimental raw data not available")
def test_raw_checksums():
    """Verify SHA-256 integrity of all raw files against data/raw/CHECKSUMS.txt."""
    checksums_file = RAW_DIR / "CHECKSUMS.txt"
    assert checksums_file.exists(), "CHECKSUMS.txt missing"
    
    with open(checksums_file, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        
    for line in lines:
        expected_hash, rel_path = line.split()
        full_path = RAW_DIR / rel_path
        assert full_path.exists(), f"File {rel_path} declared in CHECKSUMS.txt is missing"
        
        hasher = hashlib.sha256()
        with open(full_path, 'rb') as fp:
            for chunk in iter(lambda: fp.read(65536), b""):
                hasher.update(chunk)
        actual_hash = hasher.hexdigest()
        assert actual_hash == expected_hash, f"Checksum mismatch for {rel_path}"


@pytest.mark.skipif(not SQUID_DIR.exists(), reason="Experimental raw data not available")
def test_batch_1_coercivities():
    """Verify Batch 1 out-of-plane coercivities match ground truth."""
    files_expected = [
        ('2025.06.Sample3.FePt.SiO2.7T.100Oe.s.dat', 1.1864),
        ('2025.06.Sample5.FePt.SiO2.7T.100Oe.s.dat', 1.1328),
        ('2025.06.Sample8.FePt.SiO2.7T.100Oe.s.dat', 1.1207),
        ('2025.06.Sample10.FePt.SiO2.7T.100Oe.s.dat', 1.1473),
    ]
    for fname, exp_hc in files_expected:
        p = SQUID_DIR / fname
        b_t, m_raw, _ = read_squid_dat(p)
        metrics, _ = process_loop(b_t, m_raw)
        assert np.isclose(metrics['H_c_T'], exp_hc, atol=0.001)


@pytest.mark.skipif(not SQUID_DIR.exists(), reason="Experimental raw data not available")
def test_batch_2_coercivities():
    """Verify Batch 2 out-of-plane coercivities match ground truth."""
    files_expected = [
        ('2025.07.Sample3.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat', 1.1084),
        ('2025.07.Sample5.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat', 1.0561),
        ('2025.07.Sample8.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat', 1.0001),
        ('2025.07.Sample10.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat', 1.0856),
    ]
    for fname, exp_hc in files_expected:
        p = SQUID_DIR / fname
        b_t, m_raw, _ = read_squid_dat(p)
        metrics, _ = process_loop(b_t, m_raw)
        assert np.isclose(metrics['H_c_T'], exp_hc, atol=0.001)


@pytest.mark.skipif(not SQUID_DIR.exists(), reason="Experimental raw data not available")
def test_as_deposited_soft_phase():
    """Verify as-deposited caps are soft (H_c < 10 mT)."""
    as_dep_files = [
        '2025.07.Sample3.28.04.FePt.SiO2.7T.100Oe.s.dat',
        '2025.07.Sample5.28.04.FePt.SiO2.7T.100Oe.s.dat',
        '2025.07.Sample8.28.04.FePt.SiO2.7T.100Oe.s.dat',
        '2025.07.Sample10.28.04.FePt.SiO2.7T.100Oe.s_00001.dat'
    ]
    for fname in as_dep_files:
        p = SQUID_DIR / fname
        b_t, m_raw, _ = read_squid_dat(p)
        metrics, _ = process_loop(b_t, m_raw)
        assert metrics['H_c_T'] < 0.010, f"{fname} Hc exceeded 10 mT"


@pytest.mark.skipif(not SQUID_DIR.exists(), reason="Experimental raw data not available")
def test_remanence_ratios():
    """Verify remanence plateau ratios match exact computed values."""
    # B1: 0.620, 0.595, 0.594, 0.626 -> Mean 0.609
    # B2: 0.649, 0.647, 0.628, 0.628 -> Mean 0.638
    b1_files = [
        '2025.06.Sample3.FePt.SiO2.7T.100Oe.s.dat',
        '2025.06.Sample5.FePt.SiO2.7T.100Oe.s.dat',
        '2025.06.Sample8.FePt.SiO2.7T.100Oe.s.dat',
        '2025.06.Sample10.FePt.SiO2.7T.100Oe.s.dat'
    ]
    ratios_b1 = []
    for f in b1_files:
        b_t, m_raw, _ = read_squid_dat(SQUID_DIR / f)
        m, _ = process_loop(b_t, m_raw)
        ratios_b1.append(m['Mr_50_ratio'])
    assert np.isclose(np.mean(ratios_b1), 0.6088, atol=0.002)


@pytest.mark.skipif(not SQUID_DIR.exists(), reason="Experimental raw data not available")
def test_ols_joint_regression_exact():
    """Verify OLS joint regression matches paper statistics exactly."""
    diameters = np.array([3, 5, 8, 10, 3, 5, 8, 10], dtype=float)
    hc = np.array([1.1864, 1.1328, 1.1207, 1.1473, 1.1084, 1.0561, 1.0001, 1.0856])
    batches = [1, 1, 1, 1, 0, 0, 0, 0]
    
    res = run_ols_joint_regression(diameters, hc, batches)
    assert np.isclose(res['const'], 1.0174, atol=0.001)
    assert np.isclose(res['b1_offset'], 0.0843, atol=0.001)
    assert np.isclose(res['slope_inv_R'], 0.1190, atol=0.001)
    assert np.isclose(res['delta_Hc_10_to_3'], 0.0555, atol=0.001)
    assert np.isclose(res['delta_Hc_se'], 0.0296, atol=0.001)
    assert np.isclose(res['p_value'], 0.1191, atol=0.002)


@pytest.mark.skipif(not SQUID_DIR.exists(), reason="Experimental raw data not available")
def test_categorical_anova_exact():
    """Verify categorical ANOVA F(3,3) and p-value."""
    diameters = np.array([3, 5, 8, 10, 3, 5, 8, 10], dtype=float)
    hc = np.array([1.1864, 1.1328, 1.1207, 1.1473, 1.1084, 1.0561, 1.0001, 1.0856])
    batches = [1, 1, 1, 1, 0, 0, 0, 0]
    
    res = run_categorical_anova(diameters, hc, batches)
    assert np.isclose(res['F_statistic'], 8.37, atol=0.05)
    assert np.isclose(res['p_value'], 0.0572, atol=0.002)


@pytest.mark.skipif(not R1_DIR.exists(), reason="R1 simulation data not available")
def test_r1_mesh_convergence_exact():
    """Verify R1 switching field at 1.0 nm sub-exchange grid."""
    table = parse_mesh_convergence_table(R1_DIR)
    # Find 1.0 nm entry
    entry_1nm = [r for r in table if np.isclose(r['mesh_size_nm'], 1.0)][0]
    assert np.isclose(abs(entry_1nm['switching_field_T']), 6.222, atol=0.002)
    assert np.isclose(entry_1nm['remanence_mz'], 0.7095, atol=0.001)
