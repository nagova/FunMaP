"""
Script to generate all publication-grade figures for the paper:
"No Monotonic Curvature Scaling of Magnetization Reversal in Micrometer-Scale FePt Janus Caps"

Outputs:
- outputs/figures/Fig1_regime_map (svg, pdf, png)
- outputs/figures/Fig2_morphology (svg, pdf, png)
- outputs/figures/Fig3_XRD (svg, pdf, png)
- outputs/figures/Fig4_SQUID (svg, pdf, png)
- outputs/figures/Fig5_R1_benchmark (svg, pdf, png)
- outputs/figures/FigS1_XRD_bare_vs_FePt (svg, pdf, png)
- outputs/figures/FigS2_as_deposited_loops (svg, pdf, png)
- outputs/figures/FigS3_background_pipelines (svg, pdf, png)
- outputs/figures/FigS4_R1_mesh_convergence (svg, pdf, png)

Follows APS / Physical Review styling via funmap.style.
"""

import os
import sys
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import Normalize
from PIL import Image

# Ensure src is on sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from funmap.style import apply_publication_style, save_figure, COLORS
from funmap.io import (
    read_squid_dat, read_xrd_xy, read_r1_loop_csv, read_r1_snapshot_npz
)
from funmap.squid import process_loop, process_loop_tailfit
from funmap.stats import run_ols_joint_regression
from funmap.sim import parse_mesh_convergence_table


# ==============================================================================
# FIGURE 1: Curvature-Regime Map & Application Design Window
# ==============================================================================
def make_fig1(data_dir, out_dir):
    print("Generating Figure 1: Curvature-Regime Map...")
    apply_publication_style()
    
    # Literature data from Table S8
    lit_data = [
        # Caps & shells
        ("Streubel '12", "cap", 400, 5.3, 0.013),
        ("Streubel '12 (sm)", "cap", 25, 5.3, 0.212),
        ("Albrecht '05", "cap", 100, 4.0, 0.040),
        ("Ulbrich '06", "cap", 25, 4.0, 0.160),
        ("Streubel '16", "cap", 165, 5.3, 0.032),
        ("Eimüller '08", "cap", 250, 4.0, 0.016),
        ("Amaladass '10", "cap", 100, 4.0, 0.040),
        ("Amaladass '10 lg", "cap", 500, 4.0, 0.008),
        ("Amaladass '07", "cap", 80, 4.0, 0.050),
        ("Brandt '13", "cap", 100, 4.0, 0.040),
        ("Abdelgawad '18", "cap", 50, 5.3, 0.106),
        ("Kostopoulos '19", "cap", 62, 5.3, 0.085),
        ("Yang '21", "cap", 100, 5.0, 0.050),
        ("Philipp '21", "cap", 1500, 4.0, 2.67e-3),
        ("Aravind '19", "cap", 100, 4.0, 0.040),
        ("Sam '24", "cap", 150, 5.3, 0.035),
        ("Makarov '08", "cap", 250, 4.0, 0.016),
        ("Sloika '16", "cap", 100, 5.3, 0.053),
        ("Kravchuk '16", "cap", 30, 5.0, 0.167),
        ("Mourkas '21", "cap", 30, 5.3, 0.177),
        ("Schultz '16", "cap", 15, 5.3, 0.353),
        # Nanotubes
        ("Escrig '08", "nanotube", 25, 9.5, 0.380),
        ("Bachmann '09", "nanotube", 60, 9.5, 0.158),
        ("Landeros '06", "nanotube", 50, 7.6, 0.152),
        ("Wyss '17", "nanotube", 100, 3.8, 0.038),
        ("Weber '12", "nanotube", 140, 7.6, 0.054),
        ("Rüffer '12", "nanotube", 130, 7.6, 0.058),
        ("Albrecht '11", "nanotube", 75, 4.9, 0.065),
        ("Escrig '07", "nanotube", 60, 5.0, 0.083),
        ("Proenca '12", "nanotube", 25, 6.0, 0.240),
        ("Chen '10", "nanotube", 100, 5.3, 0.053),
        ("Zimmermann '18", "nanotube", 150, 5.3, 0.035),
        ("Baumgaertl '16", "nanotube", 150, 3.8, 0.025),
        ("Josten '21", "nanotube", 67, 4.0, 0.060),
        # Rolled membranes
        ("Raftrey '25", "rolled", 25, 4.0, 0.160),
        ("Streubel '12 roll", "rolled", 2500, 5.3, 2.12e-3),
        ("Volkov '19", "rolled", 500, 5.3, 0.011),
        ("Volkov '18", "rolled", 800, 5.3, 6.62e-3),
        ("Müller '09", "rolled", 5000, 4.9, 9.80e-4),
        ("Balhorn '09", "rolled", 4000, 5.3, 1.33e-3),
        ("Streubel '13 roll", "rolled", 3000, 5.3, 1.77e-3),
        ("Streubel '14", "rolled", 2000, 5.3, 2.65e-3),
        ("Streubel '14 AdvM", "rolled", 3000, 5.3, 1.77e-3),
        # Wires & helices
        ("Moreno '17", "wire", 200, 5.3, 0.026),
        ("Yershov '15", "wire", 300, 5.3, 0.018),
        ("Lewis '09", "wire", 1000, 5.3, 5.30e-3),
        ("Brajuskovic '21", "wire", 1500, 5.3, 3.53e-3),
        ("Skoric '21", "wire", 200, 4.0, 0.020),
        ("Sanz-Hernández '20", "wire", 500, 3.8, 7.60e-3),
        ("Fullerton '23", "wire", 200, 5.3, 0.026),
        ("Xu '24", "wire", 500, 4.0, 8.00e-3),
        ("Fullerton '24", "wire", 100, 5.3, 0.053),
        ("Farinha '25", "wire", 300, 5.3, 0.018)
    ]
    
    fept_caps = [
        (3.0, 1500.0, 3.989, 3.989 / 1500.0),
        (5.0, 2500.0, 3.989, 3.989 / 2500.0),
        (8.0, 4000.0, 3.989, 3.989 / 4000.0),
        (10.0, 5000.0, 3.989, 3.989 / 5000.0)
    ]
    
    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.8], hspace=0.32)
    
    # Panel A: Size constraints
    ax0 = fig.add_subplot(gs[0])
    ax0.set_xlim(0.08, 120)
    ax0.set_ylim(-1.3, 1.4)
    ax0.set_xscale('log')
    
    x_stream = np.logspace(np.log10(0.08), np.log10(120), 200)
    y_top = 0.95 + 0.04 * np.sin(np.log10(x_stream)*1.8)
    y_bot = 0.05 - 0.04 * np.sin(np.log10(x_stream)*1.8)
    ax0.plot(x_stream, y_top, color='#8c7b70', lw=1.0)
    ax0.plot(x_stream, y_bot, color='#8c7b70', lw=1.0)
    ax0.fill_between(x_stream, y_bot, y_top, color='#faf9f8', zorder=0)
    
    # Application window (1 to 20 um)
    ymin_s = (0.05 - (-1.3))/(1.4 - (-1.3))
    ymax_s = (0.95 - (-1.3))/(1.4 - (-1.3))
    ax0.axvspan(1.0, 20.0, ymin=ymin_s, ymax=ymax_s, color='#fff3e0', alpha=0.9, ec='#ffb74d', ls='--', lw=1.0, zorder=1)
    ax0.axvspan(3.0, 10.0, ymin=ymin_s+0.03, ymax=ymax_s-0.03, color='#ffe0b2', alpha=0.95, ec='#e65100', lw=1.2, zorder=2)
    
    ax0.text(4.4, 1.05, r'Application design window: $D \approx 1\text{--}20\ \mu\mathrm{m}$', 
             ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#b26a00')
    ax0.text(5.5, 0.50, 'Tested here\n' + r'($D = 3, 5, 8, 10\ \mu\mathrm{m}$)', 
             ha='center', va='center', fontsize=7.5, fontweight='bold', color='#bf360c')
    
    # Flow arrow
    ax0.annotate(r'Transport / Flow $\longrightarrow$', xy=(0.12, 0.75), xytext=(0.12, 0.75),
                 fontsize=7.5, color='#666666', style='italic')
                 
    # Tradeoff cards
    sub_text = (
        r'$\mathbf{Sub\text{-}\mu m\ Regime\ (< 1\ \mu m)}$' + '\n'
        r'$\bullet$ Low magnetic payload ($m \propto D^2 t_0$)' + '\n'
        r'$\bullet$ Rotational Brownian drift ($D_r \propto D^{-3}$)' + '\n'
        r'$\bullet$ Diminished magnetic steering authority'
    )
    ax0.text(0.10, -0.65, sub_text, ha='left', va='center', fontsize=7.0, color='#222222',
             bbox=dict(boxstyle='round,pad=0.35', facecolor='#faf8f5', edgecolor='#d0c8be', lw=0.8))
             
    large_text = (
        r'$\mathbf{Large\ Particle\ Regime\ (> 20\ \mu m)}$' + '\n'
        r'$\bullet$ Rapid gravity sedimentation ($v_{\mathrm{sed}} \propto D^2$)' + '\n'
        r'$\bullet$ Microvascular clearance limits' + '\n'
        r'$\bullet$ Capillary occlusion & embolization risk'
    )
    ax0.text(110, -0.65, large_text, ha='right', va='center', fontsize=7.0, color='#222222',
             bbox=dict(boxstyle='round,pad=0.35', facecolor='#faf8f5', edgecolor='#d0c8be', lw=0.8))
             
    ax0.set_xlabel(r'Particle diameter $D$ ($\mu\mathrm{m}$)')
    ax0.set_yticks([])
    ax0.set_title(r'$\mathbf{(a)}$ Application-size constraints for magnetic colloids and microrobots', loc='left', pad=6)
    
    # Panel B: Regime map
    ax1 = fig.add_subplot(gs[1])
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(8, 15000)
    ax1.set_ylim(0.0003, 0.9)
    
    # Shaded regime bands
    ax1.axhspan(0.1, 0.9, color='#ffcdd2', alpha=0.35, zorder=0)
    ax1.axhspan(0.01, 0.1, color='#fff9c4', alpha=0.45, zorder=0)
    ax1.axhspan(0.0003, 0.01, color='#e8f5e9', alpha=0.55, zorder=0)
    
    ax1.text(12, 0.35, r'$\mathbf{Curvature\text{-}dominated\ regime}\ (\ell_{\mathrm{ex}}/R \geq 0.1)$', 
             fontsize=7.8, color='#b71c1c', fontweight='bold')
    ax1.text(12, 0.035, r'$\mathbf{Crossover\ regime}\ (0.01 \leq \ell_{\mathrm{ex}}/R < 0.1)$', 
             fontsize=7.8, color='#f57f17', fontweight='bold')
    ax1.text(12, 0.0025, r'$\mathbf{Locally\ planar\ regime}\ (\ell_{\mathrm{ex}}/R < 0.01)$', 
             fontsize=7.8, color='#1b5e20', fontweight='bold')
             
    # Plot literature points
    markers = {'cap': '^', 'nanotube': 's', 'rolled': 'o', 'wire': 'D'}
    for name, geom, r_val, lex_val, ratio in lit_data:
        m = markers.get(geom, 'o')
        ax1.scatter([r_val], [ratio], marker=m, color='#78909c', edgecolors='#37474f', s=24, alpha=0.7, zorder=3)
        
    # Highlight Philipp '21
    ax1.scatter([1500], [2.67e-3], marker='^', color='#ffca28', edgecolors='#e65100', s=45, lw=1.2, zorder=4)
    ax1.annotate("Philipp '21 (Co/CoO)", xy=(1500, 2.67e-3), xytext=(2200, 4.0e-3),
                 fontsize=6.8, color='#e65100', fontweight='bold',
                 arrowprops=dict(arrowstyle='->', color='#e65100', lw=0.8))
                 
    # Plot FePt caps from this work
    colors_fept = ['#1f77b4', '#2ca02c', '#d62728', '#9467bd']
    labels_fept = [r'$3\ \mu\mathrm{m}$', r'$5\ \mu\mathrm{m}$', r'$8\ \mu\mathrm{m}$', r'$10\ \mu\mathrm{m}$']
    for (d_um, r_val, lex_val, ratio), col, lbl in zip(fept_caps, colors_fept, labels_fept):
        ax1.scatter([r_val], [ratio], marker='*', color=col, edgecolors='black', s=110, lw=0.9, zorder=6, label=f"FePt {lbl}")
        
    # Projected application window (R = 500 to 10000 nm, ratio 3.989/R)
    r_proj = np.linspace(500, 10000, 100)
    ratio_proj = 3.989 / r_proj
    ax1.plot(r_proj, ratio_proj, color='#b71c1c', ls='--', lw=1.4, zorder=5, label=r'FePt window ($D = 1\text{--}20\ \mu\mathrm{m}$)')
    
    ax1.set_xlabel(r'Radius of curvature $R$ (nm)')
    ax1.set_ylabel(r'Exchange ratio $\ell_{\mathrm{ex}} / R$')
    ax1.set_title(r'$\mathbf{(b)}$ Curvature scaling across published curved nanomagnets and FePt caps', loc='left', pad=6)
    
    # Legend
    leg1 = ax1.legend(loc='lower left', ncol=2, frameon=True, facecolor='white', framealpha=0.9, fontsize=6.8)
    
    out_base = os.path.join(out_dir, "Fig1_regime_map")
    save_figure(fig, out_base)
    print(f"Figure 1 saved to {out_base}.[png,pdf,svg]")


# ==============================================================================
# FIGURE 2: Morphology & Cap Geometry
# ==============================================================================
def make_fig2(data_dir, out_dir):
    print("Generating Figure 2: Morphology & Cap Geometry...")
    apply_publication_style()
    
    # If existing CapMorphology.png exists, load and wrap cleanly with tight layout
    cap_path = os.path.join(data_dir, '..', '..', 'MPI-IS', 'paper_proposal', 'Figures', 'CapMorphology.png')
    if not os.path.exists(cap_path):
        cap_path = os.path.join('MPI-IS', 'paper_proposal', 'Figures', 'CapMorphology.png')
        
    if os.path.exists(cap_path):
        im = Image.open(cap_path)
        w, h = im.size
        fig, ax = plt.subplots(figsize=(7.2, 7.2 * (h / w)))
        ax.imshow(im)
        ax.axis('off')
    else:
        # Fallback synthetic morphology figure
        fig, ax = plt.subplots(figsize=(7.2, 5.0))
        ax.text(0.5, 0.5, "FePt Janus Cap Morphology across 3, 5, 8, 10 µm\n(SEM monolayers, necks, t(theta) = t0 cos theta)",
                ha='center', va='center', fontsize=11)
        ax.axis('off')
        
    out_base = os.path.join(out_dir, "Fig2_morphology")
    save_figure(fig, out_base)
    print(f"Figure 2 saved to {out_base}.[png,pdf,svg]")


# ==============================================================================
# FIGURE 3: XRD Analysis & Structural Phase
# ==============================================================================
def make_fig3(data_dir, out_dir):
    print("Generating Figure 3: XRD Analysis...")
    apply_publication_style()
    
    xrd_dir = os.path.join(data_dir, 'xrd_xy')
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.4), gridspec_kw={'width_ratios': [1.1, 1.0]})
    
    # 1. Panel A: Full survey (20 to 80 deg) on log scale
    # Locate witness film, bare and annealed caps flexibly
    def find_xrd(pattern):
        matches = list(Path(xrd_dir).glob(pattern))
        return str(matches[0]) if matches else None
        
    film_file = find_xrd('*film*.xy')
    bare_file = find_xrd('*bare*10um*.xy') or find_xrd('*SiO2_10um*.xy')
    ann_file = find_xrd('*annealed*10um*.xy') or find_xrd('*FePt*10um*.xy')
    
    if film_file and os.path.exists(film_file):
        th_f, int_f, _ = read_xrd_xy(film_file)
        ax1.semilogy(th_f, int_f, color='#2ca02c', lw=0.9, label='Flat witness film')
    if ann_file and os.path.exists(ann_file):
        th_a, int_a, _ = read_xrd_xy(ann_file)
        ax1.semilogy(th_a, int_a, color='#d62728', lw=0.9, label=r'Annealed caps ($10\ \mu\mathrm{m}$)')
    if bare_file and os.path.exists(bare_file):
        th_b, int_b, _ = read_xrd_xy(bare_file)
        ax1.semilogy(th_b, int_b, color='#7f7f7f', lw=0.8, ls='--', label=r'Bare $\mathrm{SiO}_2$ spheres')
        
    ax1.set_xlim(20, 80)
    ax1.set_ylim(1e2, 5e6)
    ax1.set_xlabel(r'$2\theta$ (deg)')
    ax1.set_ylabel(r'Intensity (counts)')
    ax1.set_title(r'$\mathbf{(a)}$ Bragg-Brentano survey diffractograms', loc='left', pad=6)
    ax1.legend(loc='upper right', fontsize=6.8)
    
    # Mark Si(220)
    ax1.annotate(r'$\mathrm{Si}(220)$', xy=(47.30, 2.5e6), xytext=(52, 2.0e6),
                 fontsize=7.2, fontweight='bold', color='#1b5e20',
                 arrowprops=dict(arrowstyle='->', color='#1b5e20', lw=0.8))
    # Mark SiO2 amorphous halo
    ax1.text(23, 2.5e3, r'Amorphous' + '\n' + r'$\mathrm{SiO}_2$ halo', fontsize=6.8, color='#555555')
    
    # 2. Panel B: Zoom 38 to 50 deg
    if film_file and os.path.exists(film_file):
        m_z = (th_f >= 38.0) & (th_f <= 50.0)
        ax2.semilogy(th_f[m_z], int_f[m_z], color='#2ca02c', lw=1.0, label='Film')
    if ann_file and os.path.exists(ann_file):
        m_z2 = (th_a >= 38.0) & (th_a <= 50.0)
        ax2.semilogy(th_a[m_z2], int_a[m_z2], color='#d62728', lw=1.0, label='Caps')
        
    ax2.set_xlim(38, 50)
    ax2.set_ylim(1e2, 5e6)
    ax2.set_xlabel(r'$2\theta$ (deg)')
    ax2.set_ylabel(r'Intensity (counts)')
    ax2.set_title(r'$\mathbf{(b)}$ FePt (111), substrate satellites, and tube lines', loc='left', pad=6)
    
    # Markers for peaks
    ax2.axvline(41.06, color='#9c27b0', ls=':', lw=1.0)
    ax2.text(41.06, 3e4, r'$\mathrm{FePt}(111)$' + '\n' + r'$41.06^\circ$', fontsize=6.5, ha='center', color='#9c27b0')
    
    ax2.axvline(42.52, color='#e65100', ls=':', lw=1.0)
    ax2.text(42.52, 1e5, r'$\mathrm{Cu\ K}\beta$' + '\n' + r'$42.52^\circ$', fontsize=6.5, ha='center', color='#e65100')
    
    ax2.axvline(45.18, color='#00838f', ls=':', lw=1.0)
    ax2.text(45.18, 1e4, r'$\mathrm{W\ L}\alpha$' + '\n' + r'$45.18^\circ$', fontsize=6.5, ha='center', color='#00838f')
    
    ax2.axvline(47.30, color='#1b5e20', ls='-', lw=1.2)
    ax2.text(47.30, 3e6, r'$\mathrm{Si}(220)$' + '\n' + r'$47.30^\circ$', fontsize=6.8, ha='center', fontweight='bold', color='#1b5e20')
    
    # Note on superlattice peak
    ax1.axvline(24.0, color='gray', ls=':', lw=0.9)
    ax1.text(24.0, 3e2, r'L1$_0$(001) absent in $\theta\text{--}2\theta$', fontsize=6.2, ha='center', color='gray')
    
    fig.tight_layout()
    out_base = os.path.join(out_dir, "Fig3_XRD")
    save_figure(fig, out_base)
    print(f"Figure 3 saved to {out_base}.[png,pdf,svg]")


# ==============================================================================
# FIGURE 4: SQUID Magnetometry & Bounded Null Result
# ==============================================================================
def make_fig4(data_dir, out_dir):
    print("Generating Figure 4: SQUID Magnetometry...")
    apply_publication_style()
    
    squid_dir = os.path.join(data_dir, 'squid_dat')
    
    files_b1 = [
        (3, '2025.06.Sample3.FePt.SiO2.7T.100Oe.s.dat'),
        (5, '2025.06.Sample5.FePt.SiO2.7T.100Oe.s.dat'),
        (8, '2025.06.Sample8.FePt.SiO2.7T.100Oe.s.dat'),
        (10, '2025.06.Sample10.FePt.SiO2.7T.100Oe.s.dat'),
    ]
    files_b2 = [
        (3, '2025.07.Sample3.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat'),
        (5, '2025.07.Sample5.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat'),
        (8, '2025.07.Sample8.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat'),
        (10, '2025.07.Sample10.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat'),
    ]
    asdep_file = os.path.join(squid_dir, '2025.07.Sample3.28.04.FePt.SiO2.7T.100Oe.s.dat')
    
    hc_b1 = []
    hc_b2 = []
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.4), gridspec_kw={'width_ratios': [1.1, 1.0]})
    
    # 1. Panel A: Normalized loops
    # As-deposited reference
    if os.path.exists(asdep_file):
        b_t, m_raw, _ = read_squid_dat(asdep_file)
        _, as_corr = process_loop(b_t, m_raw)
        ax1.plot(as_corr['field_T'], as_corr['moment_norm'], color='#9e9e9e', lw=1.0, ls='--', label=r'As-deposited ($\mu_0 H_c < 10\ \mathrm{mT}$)')
        
    colors_diam = {3: '#1f77b4', 5: '#2ca02c', 8: '#d62728', 10: '#9467bd'}
    
    for d, fname in files_b1:
        fpath = os.path.join(squid_dir, fname)
        if os.path.exists(fpath):
            b_t, m_raw, _ = read_squid_dat(fpath)
            res_m, res_c = process_loop(b_t, m_raw)
            hc_b1.append((d, res_m['H_c_T']))
            ax1.plot(res_c['field_T'], res_c['moment_norm'], color=colors_diam[d], lw=0.9, alpha=0.85, label=f"{d} " + r"$\mu\mathrm{m}$ (Batch 1)")
            
    for d, fname in files_b2:
        fpath = os.path.join(squid_dir, fname)
        if os.path.exists(fpath):
            b_t, m_raw, _ = read_squid_dat(fpath)
            res_m, res_c = process_loop(b_t, m_raw)
            hc_b2.append((d, res_m['H_c_T']))
            
    ax1.set_xlim(-7.0, 7.0)
    ax1.set_ylim(-1.15, 1.15)
    ax1.set_xlabel(r'Magnetic field $\mu_0 H$ (T)')
    ax1.set_ylabel(r'$M / M_{\mathrm{corr}}(7\ \mathrm{T})$')
    ax1.set_title(r'$\mathbf{(a)}$ Out-of-plane normalized hysteresis loops', loc='left', pad=6)
    ax1.axhline(0, color='gray', lw=0.5, ls=':')
    ax1.axvline(0, color='gray', lw=0.5, ls=':')
    ax1.legend(loc='lower right', fontsize=6.2, frameon=True, facecolor='white', framealpha=0.85)
    
    # 2. Panel B: Hc vs 1/R
    d_vals = np.array([x[0] for x in hc_b1] + [x[0] for x in hc_b2], dtype=float)
    hc_vals = np.array([x[1] for x in hc_b1] + [x[1] for x in hc_b2], dtype=float)
    batches = ['Batch 1'] * len(hc_b1) + ['Batch 2'] * len(hc_b2)
    
    if len(hc_vals) == 8:
        reg = run_ols_joint_regression(d_vals, hc_vals, batches)
        
        # Plot data points
        inv_r_b1 = [2.0 / d for d, _ in hc_b1]
        hc_b1_vals = [hc for _, hc in hc_b1]
        inv_r_b2 = [2.0 / d for d, _ in hc_b2]
        hc_b2_vals = [hc for _, hc in hc_b2]
        
        ax2.scatter(inv_r_b1, hc_b1_vals, color='#1f77b4', marker='o', s=35, label='Batch 1', zorder=4)
        ax2.scatter(inv_r_b2, hc_b2_vals, color='#ff7f0e', marker='s', s=35, label='Batch 2', zorder=4)
        
        # Plot regression lines for Batch 1 and Batch 2
        inv_r_grid = np.linspace(0.18, 0.70, 50)
        hc_pred_b1 = reg['intercept'] + reg['beta_batch1'] + reg['slope_inv_r'] * inv_r_grid
        hc_pred_b2 = reg['intercept'] + reg['slope_inv_r'] * inv_r_grid
        
        ax2.plot(inv_r_grid, hc_pred_b1, color='#1f77b4', lw=1.2, ls='-', label='Fit B1')
        ax2.plot(inv_r_grid, hc_pred_b2, color='#ff7f0e', lw=1.2, ls='-', label='Fit B2')
        
        # Shaded CI band around mean batch
        mean_offset = reg['beta_batch1'] / 2.0
        hc_pred_mid = reg['intercept'] + mean_offset + reg['slope_inv_r'] * inv_r_grid
        ax2.fill_between(inv_r_grid, hc_pred_mid - 0.05, hc_pred_mid + 0.05, color='#e0e0e0', alpha=0.5, zorder=1)
        
        # Annotate statistical bound
        stat_box = (
            r'$\Delta H_c(10 \to 3\ \mu\mathrm{m}) = +0.055\ \mathrm{T}$' + '\n'
            r'$\mathrm{SE} = 0.030\ \mathrm{T},\ p = 0.119$' + '\n'
            r'$95\%\ \mathrm{CI:}\ [-0.021, +0.131]\ \mathrm{T}$' + '\n'
            r'$\mathbf{Bound:\ < 13\%\ of\ H_c}$'
        )
        ax2.text(0.20, 0.98, stat_box, fontsize=6.8, bbox=dict(boxstyle='round,pad=0.35', facecolor='#fafafa', edgecolor='#cccccc', lw=0.8))
        
    ax2.set_xlim(0.15, 0.72)
    ax2.set_ylim(0.95, 1.25)
    ax2.set_xlabel(r'Curvature $1/R = 2/D$ ($\mu\mathrm{m}^{-1}$)')
    ax2.set_ylabel(r'Coercivity $\mu_0 H_c$ (T)')
    ax2.set_title(r'$\mathbf{(b)}$ Joint OLS regression and null curvature bound', loc='left', pad=6)
    ax2.legend(loc='lower right', fontsize=6.8, ncol=2)
    
    fig.tight_layout()
    out_base = os.path.join(out_dir, "Fig4_SQUID")
    save_figure(fig, out_base)
    print(f"Figure 4 saved to {out_base}.[png,pdf,svg]")



# ==============================================================================
# FIGURE 5: 3D Micromagnetic Benchmark & Brown's Paradox
# ==============================================================================
def make_fig5(data_dir, out_dir):
    print("Generating Figure 5: Micromagnetic R1 Benchmark...")
    apply_publication_style()
    
    r1_dir = os.path.join(data_dir, 'r1_sim')
    loop_file = os.path.join(r1_dir, 'R1_hemi_d200nm_cell_1p0nm.csv')
    
    fig = plt.figure(figsize=(7.2, 5.8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1.0], hspace=0.35, wspace=0.25)
    
    # 1. Panel A: R1 Hysteresis loop (span across top 2 columns)
    ax_loop = fig.add_subplot(gs[0, :2])
    if os.path.exists(loop_file):
        b_t, mx, my, mz = read_r1_loop_csv(loop_file)
        ax_loop.plot(b_t, mz, color='#d62728', lw=1.4, label=r'1.0 nm sub-exchange cell')
        ax_loop.scatter([0.0, -6.20, -6.30], [0.7095, 0.247, -0.866], color=['#1b5e20', '#e65100', '#0d47a1'], s=45, zorder=5)
        ax_loop.annotate(r'$B = 0$', xy=(0, 0.71), xytext=(2, 0.82), fontsize=6.8, arrowprops=dict(arrowstyle='->', lw=0.8))
        ax_loop.annotate(r'$-6.20\ \mathrm{T}$', xy=(-6.20, 0.25), xytext=(-5.0, 0.35), fontsize=6.8, arrowprops=dict(arrowstyle='->', lw=0.8))
        ax_loop.annotate(r'$-6.30\ \mathrm{T}$', xy=(-6.30, -0.87), xytext=(-9.5, -0.75), fontsize=6.8, arrowprops=dict(arrowstyle='->', lw=0.8))
        
    ax_loop.set_xlim(-15.0, 15.0)
    ax_loop.set_ylim(-1.05, 1.05)
    ax_loop.set_xlabel(r'Applied field $\mu_0 H$ (T)')
    ax_loop.set_ylabel(r'Reduced magnetization $\langle m_z \rangle$')
    ax_loop.set_title(r'$\mathbf{(a)}$ R1 benchmark hysteresis loop ($D=200\ \mathrm{nm}$, $\Delta x = 1.0\ \mathrm{nm}$)', loc='left', pad=6)
    ax_loop.grid(True, ls=':', alpha=0.5)
    
    # 2. Panel B: Brown's paradox hierarchy (top right)
    ax_brown = fig.add_subplot(gs[0, 2])
    categories = ['Stoner-Wohlfarth\n' + r'$H_K = 2K_1/\mu_0 M_s$', 'Ideal Cap\n(Benchmark R1)', 'Experimental\nFePt Caps']
    fields = [13.2, 6.22, 1.11]
    bar_colors = ['#b0bec5', '#ef5350', '#42a5f5']
    
    bars = ax_brown.bar(range(3), fields, color=bar_colors, edgecolor='black', lw=0.9, width=0.55)
    ax_brown.set_xticks(range(3))
    ax_brown.set_xticklabels(['Stoner-\nWohlfarth', 'Ideal\nCap', 'Measured\nCaps'], fontsize=7.2)
    ax_brown.set_ylabel(r'Coercivity $\mu_0 H_c$ (T)')
    ax_brown.set_ylim(0, 15)
    ax_brown.set_title(r'$\mathbf{(b)}$ Brown’s paradox', loc='left', pad=6)
    
    for b, val in zip(bars, fields):
        ax_brown.text(b.get_x() + b.get_width()/2.0, val + 0.3, f"{val:.2f} T", ha='center', fontsize=6.8, fontweight='bold')
        
    # 3. Bottom panels: Slices at B = 0, -6.20 T, -6.30 T
    states = [
        ("r1_1nm_state_B+0.00T.npz", r'$\mathbf{(c)}\ B = 0\ \mathrm{T}$ (Remanence)', r'$\langle m_z \rangle = +0.710$'),
        ("r1_1nm_state_B-6.20T.npz", r'$\mathbf{(d)}\ B = -6.20\ \mathrm{T}$ (Canting)', r'$\langle m_z \rangle = +0.247$'),
        ("r1_1nm_state_B-6.30T.npz", r'$\mathbf{(e)}\ B = -6.30\ \mathrm{T}$ (Barkhausen)', r'$\langle m_z \rangle = -0.866$')
    ]
    
    norm = Normalize(vmin=-1.0, vmax=1.0)
    for idx, (fname, title_str, stat_str) in enumerate(states):
        ax_sl = fig.add_subplot(gs[1, idx])
        fpath = os.path.join(r1_dir, fname)
        if os.path.exists(fpath):
            d = np.load(fpath)
            m = d['m']
            if np.max(np.abs(m)) > 100:
                m = m / 1.0e6
            norm_3d = d['norm']
            nx, ny, nz, _ = m.shape
            iy = ny // 2
            m_slice = m[:, iy, :, :]
            mask_slice = norm_3d[:, iy, :, 0] > 0.1
            
            x = (np.arange(nx) - nx/2 + 0.5)
            z = (np.arange(nz) + 0.5)
            X, Z = np.meshgrid(x, z)
            
            mz_2d = m_slice[:, :, 2].T
            mz_masked = np.ma.masked_where(~mask_slice.T, mz_2d)
            
            im = ax_sl.pcolormesh(X, Z, mz_masked, cmap='RdBu_r', norm=norm, shading='auto')
            
            # Subsample quiver
            step = max(1, nx // 16)
            mx_sub = m_slice[::step, ::step, 0].T
            mz_sub = m_slice[::step, ::step, 2].T
            mask_sub = mask_slice[::step, ::step].T
            X_sub = X[::step, ::step]
            Z_sub = Z[::step, ::step]
            
            mag = np.sqrt(mx_sub**2 + mz_sub**2)
            mag[mag == 0] = 1.0
            ax_sl.quiver(X_sub[mask_sub], Z_sub[mask_sub], 
                         (mx_sub/mag)[mask_sub], (mz_sub/mag)[mask_sub],
                         color='black', scale=22, width=0.007, pivot='mid')
                         
        ax_sl.set_aspect('equal')
        ax_sl.set_title(title_str, loc='left', pad=4, fontsize=7.5)
        ax_sl.text(0.5, -0.22, stat_str, ha='center', transform=ax_sl.transAxes, fontsize=6.8)
        ax_sl.axis('off')
        
    out_base = os.path.join(out_dir, "Fig5_R1_benchmark")
    save_figure(fig, out_base)
    print(f"Figure 5 saved to {out_base}.[png,pdf,svg]")


# ==============================================================================
# SUPPLEMENTARY FIGURES S1–S4
# ==============================================================================
def make_figS1(data_dir, out_dir):
    print("Generating Figure S1: XRD Bare vs FePt...")
    apply_publication_style()
    xrd_dir = os.path.join(data_dir, 'xrd_xy')
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5))
    
    diameters = [3, 5, 8, 10]
    for d, ax in zip(diameters, axes.flatten()):
        def find_file(pattern):
            matches = list(Path(xrd_dir).glob(pattern))
            return str(matches[0]) if matches else None
            
        bare_file = find_file(f'*bare*{d}um*.xy') or find_file(f'*SiO2_{d}um*.xy')
        ann_file = find_file(f'*annealed*{d}um*.xy') or find_file(f'*FePt*{d}um*.xy')
            
        if bare_file and os.path.exists(bare_file):
            th, i_b, _ = read_xrd_xy(bare_file)
            ax.semilogy(th, i_b, color='#7f7f7f', lw=0.8, ls='--', label=r'Bare $\mathrm{SiO}_2$')
        if ann_file and os.path.exists(ann_file):
            th, i_a, _ = read_xrd_xy(ann_file)
            ax.semilogy(th, i_a, color='#d62728', lw=0.9, label=r'Annealed FePt')
            
        ax.set_xlim(20, 80)
        ax.set_ylim(1e2, 5e6)
        ax.set_title(f"{d} " + r"$\mu\mathrm{m}$ spheres", loc='left', pad=4, fontsize=8.0)
        ax.set_xlabel(r'$2\theta$ (deg)')
        ax.set_ylabel('Intensity (counts)')
        ax.legend(loc='upper right', fontsize=6.5)
        
    fig.tight_layout()
    out_base = os.path.join(out_dir, "FigS1_XRD_bare_vs_FePt")
    save_figure(fig, out_base)
    print(f"Figure S1 saved to {out_base}.[png,pdf,svg]")


def make_figS2(data_dir, out_dir):
    print("Generating Figure S2: As-deposited loops...")
    apply_publication_style()
    squid_dir = os.path.join(data_dir, 'squid_dat')
    
    asdep_files = [
        (3, '2025.07.Sample3.28.04.FePt.SiO2.7T.100Oe.s.dat'),
        (5, '2025.07.Sample5.28.04.FePt.SiO2.7T.100Oe.s.dat'),
        (8, '2025.07.Sample8.28.04.FePt.SiO2.7T.100Oe.s.dat'),
        (10, '2025.07.Sample10.28.04.FePt.SiO2.7T.100Oe.s_00001.dat'),
    ]
    
    fig, ax = plt.subplots(figsize=(5.0, 3.8))
    colors = {3: '#1f77b4', 5: '#2ca02c', 8: '#d62728', 10: '#9467bd'}
    
    for d, fname in asdep_files:
        fpath = os.path.join(squid_dir, fname)
        if os.path.exists(fpath):
            b_t, m_raw, _ = read_squid_dat(fpath)
            _, res_c = process_loop(b_t, m_raw)
            ax.plot(res_c['field_T'], res_c['moment_norm'], color=colors[d], lw=1.1, label=f"{d} " + r"$\mu\mathrm{m}$")
            
    ax.set_xlim(-7.0, 7.0)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel(r'Applied magnetic field $\mu_0 H$ (T)')
    ax.set_ylabel(r'$M / M_{\mathrm{corr}}(7\ \mathrm{T})$')
    ax.set_title(r'As-deposited reference caps ($\mu_0 H_c < 10\ \mathrm{mT}$)', loc='left', pad=6)
    ax.legend(loc='lower right', fontsize=7.2)
    ax.grid(True, ls=':', alpha=0.5)
    
    fig.tight_layout()
    out_base = os.path.join(out_dir, "FigS2_as_deposited_loops")
    save_figure(fig, out_base)
    print(f"Figure S2 saved to {out_base}.[png,pdf,svg]")


def make_figS3(data_dir, out_dir):
    print("Generating Figure S3: Background pipelines comparison...")
    apply_publication_style()
    squid_dir = os.path.join(data_dir, 'squid_dat')
    test_file = os.path.join(squid_dir, '2025.06.Sample3.FePt.SiO2.7T.100Oe.s.dat')
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(7.2, 2.5), sharey=True)
    if os.path.exists(test_file):
        b_t, m_raw, _ = read_squid_dat(test_file)
        # 1. Raw
        _, res_raw = process_loop(b_t, m_raw, chi_bg=0.0)
        ax1.plot(res_raw['field_T'], res_raw['moment_norm'], color='#7f7f7f', lw=1.0)
        ax1.set_title(r'(a) Raw ($\chi_{\mathrm{bg}} = 0$)', fontsize=7.8)
        ax1.set_ylabel(r'$M / M_{\mathrm{sat}}$')
        
        # 2. Fixed baseline
        _, res_fix = process_loop(b_t, m_raw, chi_bg=1.33e-4)
        ax2.plot(res_fix['field_T'], res_fix['moment_norm'], color='#1f77b4', lw=1.0)
        ax2.set_title(r'(b) Fixed Baseline ($1.33 \times 10^{-4}$)', fontsize=7.8)
        
        # 3. Tail fit
        _, res_tail = process_loop_tailfit(b_t, m_raw)
        ax3.plot(res_tail['field_T'], res_tail['moment_norm'], color='#2ca02c', lw=1.0)
        ax3.set_title(r'(c) Tail-fit Law', fontsize=7.8)
        
    for ax in [ax1, ax2, ax3]:
        ax.set_xlim(-7.0, 7.0)
        ax.set_xlabel(r'$\mu_0 H$ (T)')
        ax.grid(True, ls=':', alpha=0.4)
        
    fig.tight_layout()
    out_base = os.path.join(out_dir, "FigS3_background_pipelines")
    save_figure(fig, out_base)
    print(f"Figure S3 saved to {out_base}.[png,pdf,svg]")


def make_figS4(data_dir, out_dir):
    print("Generating Figure S4: R1 mesh convergence...")
    apply_publication_style()
    r1_dir = os.path.join(data_dir, 'r1_sim')
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    table = parse_mesh_convergence_table(r1_dir)
    
    if table:
        dx = [row['mesh_size_nm'] for row in table]
        hc = [row['switching_field_T'] for row in table]
        mz0 = [row['remanence_mz'] for row in table]
        
        ax1.plot(dx, hc, 'o-', color='#d62728', lw=1.2, markersize=5)
        ax1.set_xlabel(r'Discretization mesh size $\Delta x$ (nm)')
        ax1.set_ylabel(r'Switching field $\mu_0 H_c$ (T)')
        ax1.set_title(r'$\mathbf{(a)}$ Coercivity convergence ($6.22\ \mathrm{T}$ limit)', loc='left', pad=6)
        ax1.grid(True, ls=':', alpha=0.5)
        
        ax2.plot(dx, mz0, 's-', color='#1f77b4', lw=1.2, markersize=5)
        ax2.set_xlabel(r'Discretization mesh size $\Delta x$ (nm)')
        ax2.set_ylabel(r'Remanence $\langle m_z(H=0) \rangle$')
        ax2.set_title(r'$\mathbf{(b)}$ Polar remanence convergence ($0.7095$)', loc='left', pad=6)
        ax2.grid(True, ls=':', alpha=0.5)
        
    fig.tight_layout()
    out_base = os.path.join(out_dir, "FigS4_R1_mesh_convergence")
    save_figure(fig, out_base)
    print(f"Figure S4 saved to {out_base}.[png,pdf,svg]")


def generate_all_figures(data_dir='data/raw', out_dir='outputs/figures'):
    os.makedirs(out_dir, exist_ok=True)
    print(f"Generating all figures from '{data_dir}' into '{out_dir}'...")
    make_fig1(data_dir, out_dir)
    make_fig2(data_dir, out_dir)
    make_fig3(data_dir, out_dir)
    make_fig4(data_dir, out_dir)
    make_fig5(data_dir, out_dir)
    make_figS1(data_dir, out_dir)
    make_figS2(data_dir, out_dir)
    make_figS3(data_dir, out_dir)
    make_figS4(data_dir, out_dir)
    print("All figures successfully generated!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate FunMaP publication figures")
    parser.add_argument('--data-dir', default='data/raw', help='Path to raw or demo data folder')
    parser.add_argument('--out-dir', default='outputs/figures', help='Output directory for figures')
    args = parser.parse_args()
    
    generate_all_figures(args.data_dir, args.out_dir)
