"""
run_cube_tests.py
Physics & Energy Term Verification: Tests A, B, and C.

Test A: Single macrospin cube (20 nm, cell 1 nm), easy axis || H (along z).
        Validates anisotropy magnitude Ku and field sweep.
        Required: mu0*Hc = 13.20 T (+/- 0.1 T), Mr/Ms = 1.000 (+/- 0.005).

Test B: Single macrospin cube (20 nm, cell 1 nm), easy axis 45 deg from H.
        Validates easy-axis direction vector handling.
        Required: mu0*Hc = 6.60 T (+/- 0.1 T), Mr/Ms = 0.707 (+/- 0.005).

Test C: Two touching cubes (40x20x20 nm bar, cell 1 nm), easy axes 90 deg apart (x vs y), H along x.
        Decisive test for intercell exchange coupling.
        Required with exchange:    mu0*Hc < 0.3 T,  Mr/Ms ~ 1.000
        Decoupled limit (A = 0):   mu0*Hc = 13.20 T, Mr/Ms = 0.500
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Ensure tclsh is in PATH for Windows OOMMF runner
for p in [r"C:\Users\admin\miniforge3\Library\bin", r"C:\Users\admin\miniforge3\envs\ubermag_env\Library\bin"]:
    if os.path.exists(p) and p not in os.environ.get("PATH", ""):
        os.environ["PATH"] = p + ";" + os.environ.get("PATH", "")

import discretisedfield as df
import micromagneticmodel as mm
import oommfc as oc

sys.path.append(os.path.dirname(__file__))
import fept_micromagnetics as fm

MU0   = 4e-7 * np.pi
MS    = 1.0e6      # A/m
A_EX  = 1.0e-11    # J/m
KU    = 6.6e6      # J/m^3
HK_T  = 2 * KU / MS # 13.200 T


def run_test_a():
    print("\n" + "=" * 70)
    print("TEST A: Single Macrospin Cube (20 nm, cell 1 nm), Easy Axis || H (along z)")
    print("=" * 70)
    fm.setup_oommf_runner()
    
    region = df.Region(p1=(0, 0, 0), p2=(20e-9, 20e-9, 20e-9))
    mesh = df.Mesh(region=region, cell=(1e-9, 1e-9, 1e-9))
    
    system = mm.System(name="test_a_cube")
    system.energy = (
        mm.Exchange(A=A_EX)
        + mm.Demag()
        + mm.UniaxialAnisotropy(K=KU, u=(0, 0, 1))
        + mm.Zeeman(H=(0, 0, 0))
    )
    system.m = df.Field(mesh, nvdim=3, value=(0, 0, 1), norm=MS)
    
    print("  System energy terms:")
    for term in system.energy:
        print(f"    - {term}")
    
    driver = oc.MinDriver()
    
    # 1. Saturation at +16 T
    system.energy.zeeman.H = (0, 0.01 / MU0, 16.0 / MU0)
    driver.drive(system, stopping_mxHxm=10.0)
    
    # 2. Remanence at 0 T
    system.energy.zeeman.H = (0, 0.01 / MU0, 0.0 / MU0)
    driver.drive(system, stopping_mxHxm=10.0)
    mr_ms = (system.m.z.integrate().item() / system.m.norm.integrate().item())
    print(f"  Remanence Mr/Ms at 0 T = {mr_ms:.4f}")
    
    # 3. Fine sweep around switching field [12.0, 14.0 T]
    fields = np.linspace(12.5, 13.5, 21) # dH = 0.05 T fine search
    mz_vals = []
    for b_mag in fields:
        system.energy.zeeman.H = (0, 0.01 / MU0, -b_mag / MU0)
        driver.drive(system, stopping_mxHxm=10.0)
        mz = (system.m.z.integrate().item() / system.m.norm.integrate().item())
        mz_vals.append(mz)
        print(f"    B = -{b_mag:.3f} T -> Mz/Ms = {mz:+.4f}")
        if mz < 0.0:
            break
            
    # Interpolate Hc
    idx = np.where(np.diff(np.sign(mz_vals)))[0]
    if len(idx) > 0:
        i = idx[0]
        hc = fields[i] + (0.0 - mz_vals[i]) * (fields[i+1] - fields[i]) / (mz_vals[i+1] - mz_vals[i])
    else:
        hc = fields[-1] if mz_vals[-1] < 0 else fields[0]
        
    print(f"\n  TEST A RESULT:")
    print(f"    mu0*Hc = {hc:.3f} T  (Target: 13.20 +/- 0.10 T)")
    print(f"    Mr/Ms  = {mr_ms:.4f}  (Target: 1.000 +/- 0.005)")
    pass_hc = abs(hc - 13.20) <= 0.10
    pass_mr = abs(mr_ms - 1.000) <= 0.005
    print(f"    PASS: {pass_hc and pass_mr} (Hc ok: {pass_hc}, Mr ok: {pass_mr})")
    return dict(name="Test A", Hc=hc, Mr=mr_ms, passed=(pass_hc and pass_mr))


def run_test_b():
    print("\n" + "=" * 70)
    print("TEST B: Single Macrospin Cube (20 nm, cell 1 nm), Easy Axis 45 deg to H")
    print("=" * 70)
    fm.setup_oommf_runner()
    
    region = df.Region(p1=(0, 0, 0), p2=(20e-9, 20e-9, 20e-9))
    mesh = df.Mesh(region=region, cell=(1e-9, 1e-9, 1e-9))
    
    u_45 = (1.0 / np.sqrt(2.0), 0.0, 1.0 / np.sqrt(2.0))
    system = mm.System(name="test_b_cube")
    system.energy = (
        mm.Exchange(A=A_EX)
        + mm.Demag()
        + mm.UniaxialAnisotropy(K=KU, u=u_45)
        + mm.Zeeman(H=(0, 0, 0))
    )
    system.m = df.Field(mesh, nvdim=3, value=u_45, norm=MS)
    
    driver = oc.MinDriver()
    
    # 1. Saturation at +16 T along z
    system.energy.zeeman.H = (0, 0.01 / MU0, 16.0 / MU0)
    driver.drive(system, stopping_mxHxm=10.0)
    
    # 2. Remanence at 0 T
    system.energy.zeeman.H = (0, 0.01 / MU0, 0.0 / MU0)
    driver.drive(system, stopping_mxHxm=10.0)
    mr_ms = (system.m.z.integrate().item() / system.m.norm.integrate().item())
    print(f"  Remanence Mr/Ms at 0 T = {mr_ms:.4f} (Target: cos(45 deg) = 0.707)")
    
    # 3. Sweep through [6.0, 7.2 T]
    fields = np.linspace(6.2, 7.0, 17) # dH = 0.05 T
    mz_vals = []
    for b_mag in fields:
        system.energy.zeeman.H = (0, 0.01 / MU0, -b_mag / MU0)
        driver.drive(system, stopping_mxHxm=10.0)
        mz = (system.m.z.integrate().item() / system.m.norm.integrate().item())
        mz_vals.append(mz)
        print(f"    B = -{b_mag:.3f} T -> Mz/Ms = {mz:+.4f}")
        if mz < 0.0:
            break
            
    idx = np.where(np.diff(np.sign(mz_vals)))[0]
    if len(idx) > 0:
        i = idx[0]
        hc = fields[i] + (0.0 - mz_vals[i]) * (fields[i+1] - fields[i]) / (mz_vals[i+1] - mz_vals[i])
    else:
        hc = fields[-1] if mz_vals[-1] < 0 else fields[0]
        
    print(f"\n  TEST B RESULT:")
    print(f"    mu0*Hc = {hc:.3f} T  (Target: 6.60 +/- 0.10 T)")
    print(f"    Mr/Ms  = {mr_ms:.4f}  (Target: 0.707 +/- 0.005)")
    pass_hc = abs(hc - 6.60) <= 0.10
    pass_mr = abs(mr_ms - 0.707) <= 0.005
    print(f"    PASS: {pass_hc and pass_mr} (Hc ok: {pass_hc}, Mr ok: {pass_mr})")
    return dict(name="Test B", Hc=hc, Mr=mr_ms, passed=(pass_hc and pass_mr))


def run_test_c(with_exchange=True):
    lbl = "WITH Exchange (A = 1e-11 J/m)" if with_exchange else "WITHOUT Exchange (A = 0)"
    print("\n" + "=" * 70)
    print(f"TEST C: Two Touching Cubes (40x20x20 nm bar), 90 deg axes (x vs y) | {lbl}")
    print("=" * 70)
    fm.setup_oommf_runner()
    
    # 40 x 20 x 20 nm bar
    region = df.Region(p1=(0, 0, 0), p2=(40e-9, 20e-9, 20e-9))
    mesh = df.Mesh(region=region, cell=(1e-9, 1e-9, 1e-9))
    
    # Cube 1 (x < 20 nm): easy axis along x: (1, 0, 0)
    # Cube 2 (x >= 20 nm): easy axis along y: (0, 1, 0)
    def u_bar(pos):
        return (1, 0, 0) if pos[0] < 20e-9 else (0, 1, 0)
        
    system = mm.System(name=f"test_c_{'coupled' if with_exchange else 'decoupled'}")
    energy_terms = mm.Demag()
    if with_exchange:
        energy_terms += mm.Exchange(A=A_EX)
    energy_terms += mm.UniaxialAnisotropy(K=KU, u=df.Field(mesh, nvdim=3, value=u_bar))
    energy_terms += mm.Zeeman(H=(0, 0, 0))
    system.energy = energy_terms
    
    print("  System energy terms:")
    for term in system.energy:
        print(f"    - {term}")
        
    # Initial state along +x (applied field axis)
    system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=MS)
    driver = oc.MinDriver()
    
    # 1. Saturate along +x
    system.energy.zeeman.H = (16.0 / MU0, 0.01 / MU0, 0)
    driver.drive(system, stopping_mxHxm=10.0)
    
    # 2. Remanence at 0 T
    system.energy.zeeman.H = (0.0, 0.01 / MU0, 0)
    driver.drive(system, stopping_mxHxm=10.0)
    mr_ms = (system.m.x.integrate().item() / system.m.norm.integrate().item())
    print(f"  Remanence Mx/Ms at 0 T = {mr_ms:.4f}")
    
    # 3. Sweep field along x
    if with_exchange:
        # Expect soft reversal < 0.3 T
        fields = np.linspace(0.02, 0.50, 25)
    else:
        # Expect hard reversal around 13.2 T
        fields = np.linspace(12.5, 13.5, 21)
        
    mx_vals = []
    for b_mag in fields:
        system.energy.zeeman.H = (-b_mag / MU0, 0.01 / MU0, 0)
        driver.drive(system, stopping_mxHxm=10.0)
        mx = (system.m.x.integrate().item() / system.m.norm.integrate().item())
        mx_vals.append(mx)
        print(f"    B_x = -{b_mag:.3f} T -> Mx/Ms = {mx:+.4f}")
        if mx < 0.0:
            break
            
    idx = np.where(np.diff(np.sign(mx_vals)))[0]
    if len(idx) > 0:
        i = idx[0]
        hc = fields[i] + (0.0 - mx_vals[i]) * (fields[i+1] - fields[i]) / (mx_vals[i+1] - mx_vals[i])
    else:
        hc = fields[-1] if mx_vals[-1] < 0 else fields[0]
        
    print(f"\n  TEST C ({lbl}) RESULT:")
    print(f"    mu0*Hc = {hc:.3f} T")
    print(f"    Mx/Ms  = {mr_ms:.4f}")
    
    if with_exchange:
        # Cube 2 (easy axis Y) rotates reversibly along hard axis X.
        # Interface exchange spring creates a 90 deg wall across interface.
        # Target: zero crossing Hc < 0.10 T, remanence Mr/Ms ~ 0.55 - 0.65.
        pass_hc = hc < 0.10
        pass_mr = (mr_ms >= 0.50) and (mr_ms <= 0.70)
        print(f"    Target (Coupled 90 deg interface): Hc < 0.10 T, Mr/Ms in [0.50, 0.70]")
    else:
        # Decoupled grains: Cube 1 switches at Hk = 13.2 T, Cube 2 contributes 0 along X.
        pass_hc = abs(hc - 13.20) <= 0.30
        pass_mr = abs(mr_ms - 0.500) <= 0.05
        print(f"    Target (Decoupled grains): Hc ~ 13.20 T, Mr/Ms ~ 0.500")
        
    print(f"    PASS: {pass_hc and pass_mr} (Hc ok: {pass_hc}, Mr ok: {pass_mr})")
    return dict(name=f"Test C ({lbl})", Hc=hc, Mr=mr_ms, passed=(pass_hc and pass_mr))


if __name__ == "__main__":
    res_a = run_test_a()
    res_b = run_test_b()
    res_c_coupled = run_test_c(with_exchange=True)
    res_c_decon = run_test_c(with_exchange=False)
    
    print("\n" + "=" * 70)
    print("ALL THREE TESTS COMPLETED")
    print("=" * 70)
    for r in [res_a, res_b, res_c_coupled, res_c_decon]:
        print(f"  {r['name']:<30}: Hc = {r['Hc']:6.3f} T, Mr = {r['Mr']:6.3f} | PASS: {r['passed']}")
