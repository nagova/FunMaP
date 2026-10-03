"""
Script to generate all tables from raw data files:
- outputs/tables/sample_map.csv
- outputs/tables/Table_1_SQUID_summary.csv
- outputs/tables/Table_S4_tailfit_summary.csv
- outputs/tables/Table_S6_mesh_convergence.csv
- outputs/tables/xrd_metadata.csv
- outputs/tables/sem_catalog.csv
"""

import os
import sys
import csv
from pathlib import Path
import numpy as np

# Ensure src is on sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from funmap.io import read_squid_dat
from funmap.squid import process_loop, process_loop_tailfit
from funmap.stats import compute_batch_summary
from funmap.sim import parse_mesh_convergence_table
from funmap.xrd import parse_xrd_metadata_table
from funmap.sem import generate_sem_catalog


def make_sample_map(squid_dir, out_path):
    """
    Generate outputs/tables/sample_map.csv mapping SQUID files to batch/diameter.
    """
    # 12 expected files
    mapping = [
        ('2025.06.Sample3.FePt.SiO2.7T.100Oe.s.dat', 'Batch 1', 3, 'Annealed (500 C, 1 h)'),
        ('2025.06.Sample5.FePt.SiO2.7T.100Oe.s.dat', 'Batch 1', 5, 'Annealed (500 C, 1 h)'),
        ('2025.06.Sample8.FePt.SiO2.7T.100Oe.s.dat', 'Batch 1', 8, 'Annealed (500 C, 1 h)'),
        ('2025.06.Sample10.FePt.SiO2.7T.100Oe.s.dat', 'Batch 1', 10, 'Annealed (500 C, 1 h)'),
        ('2025.07.Sample3.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat', 'Batch 2', 3, 'Annealed (500 C, 1 h)'),
        ('2025.07.Sample5.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat', 'Batch 2', 5, 'Annealed (500 C, 1 h)'),
        ('2025.07.Sample8.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat', 'Batch 2', 8, 'Annealed (500 C, 1 h)'),
        ('2025.07.Sample10.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat', 'Batch 2', 10, 'Annealed (500 C, 1 h)'),
        ('2025.07.Sample3.28.04.FePt.SiO2.7T.100Oe.s.dat', 'As-deposited', 3, 'As-deposited A1'),
        ('2025.07.Sample5.28.04.FePt.SiO2.7T.100Oe.s.dat', 'As-deposited', 5, 'As-deposited A1'),
        ('2025.07.Sample8.28.04.FePt.SiO2.7T.100Oe.s.dat', 'As-deposited', 8, 'As-deposited A1'),
        ('2025.07.Sample10.28.04.FePt.SiO2.7T.100Oe.s_00001.dat', 'As-deposited', 10, 'As-deposited A1')
    ]
    
    rows = []
    for fname, batch, d, state in mapping:
        fpath = os.path.join(squid_dir, fname)
        if os.path.exists(fpath):
            _, _, meta = read_squid_dat(fpath)
            open_time = meta.get('file_open_time', 'NOT AVAILABLE')
        else:
            open_time = 'FILE NOT FOUND'
        rows.append({
            'filename': fname,
            'batch': batch,
            'diameter_um': d,
            'state': state,
            'file_open_time': open_time
        })
        
    with open(out_path, 'w', newline='', encoding='utf-8') as fp:
        writer = csv.DictWriter(fp, fieldnames=['filename', 'batch', 'diameter_um', 'state', 'file_open_time'])
        writer.writeheader()
        writer.writerows(rows)
    return rows


def make_table_1_and_s4(squid_dir, out_tab1_path, out_tabs4_path):
    """
    Process SQUID files and output Table 1 (fixed baseline) and Table S4 (tail-fit).
    """
    batches = {
        'Annealed Batch 1 (June 2025)': [
            ('Sample 3', 3, '2025.06.Sample3.FePt.SiO2.7T.100Oe.s.dat'),
            ('Sample 5', 5, '2025.06.Sample5.FePt.SiO2.7T.100Oe.s.dat'),
            ('Sample 8', 8, '2025.06.Sample8.FePt.SiO2.7T.100Oe.s.dat'),
            ('Sample 10', 10, '2025.06.Sample10.FePt.SiO2.7T.100Oe.s.dat'),
        ],
        'Annealed Batch 2 (July 2025)': [
            ('Sample 3', 3, '2025.07.Sample3.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat'),
            ('Sample 5', 5, '2025.07.Sample5.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat'),
            ('Sample 8', 8, '2025.07.Sample8.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat'),
            ('Sample 10', 10, '2025.07.Sample10.28.04.500.1h.FePt.SiO2.7T.100Oe.s.dat'),
        ],
        'As-deposited Reference (July 2025)': [
            ('Sample 3 (A1)', 3, '2025.07.Sample3.28.04.FePt.SiO2.7T.100Oe.s.dat'),
            ('Sample 5 (A1)', 5, '2025.07.Sample5.28.04.FePt.SiO2.7T.100Oe.s.dat'),
            ('Sample 8 (A1)', 8, '2025.07.Sample8.28.04.FePt.SiO2.7T.100Oe.s.dat'),
            ('Sample 10 (A1)', 10, '2025.07.Sample10.28.04.FePt.SiO2.7T.100Oe.s_00001.dat'),
        ]
    }
    
    rows_tab1 = []
    rows_tabs4 = []
    
    for bname, sample_list in batches.items():
        b_hc = []
        b_mr0 = []
        b_mr50 = []
        b_mcorr = []
        b_slope = []
        
        b_hc_tf = []
        b_mr50_tf = []
        
        for sname, d, fname in sample_list:
            fpath = os.path.join(squid_dir, fname)
            f_T, m_emu, _ = read_squid_dat(fpath)
            
            # Fixed baseline (Table 1)
            m_fixed, _ = process_loop(f_T, m_emu)
            # Tail-fit (Table S4)
            m_tf, _ = process_loop_tailfit(f_T, m_emu)
            
            b_hc.append(m_fixed['H_c_T'])
            b_mr0.append(m_fixed['Mr_zero_ratio'])
            b_mr50.append(m_fixed['Mr_50_ratio'])
            b_mcorr.append(m_fixed['M_corr_7T_emu'] * 1e3)  # 10^-3 emu
            b_slope.append(m_fixed['dmdh_emu_per_T'] * 1e4) # 10^-4 emu/T
            
            b_hc_tf.append(m_tf['H_c_T'])
            b_mr50_tf.append(m_tf['Mr_50_ratio'])
            
            rows_tab1.append({
                'Group': bname,
                'Sample': sname,
                'Diameter_um': d,
                'Hc_T': f"{m_fixed['H_c_T']:.4f}",
                'Mr_H0_ratio': f"{m_fixed['Mr_zero_ratio']:.3f}",
                'Mr_50mT_ratio': f"{m_fixed['Mr_50_ratio']:.3f}",
                'Mcorr_7T_memu': f"{m_fixed['M_corr_7T_emu']*1e3:.3f}",
                'dmdh_1e-4_emu_T': f"{m_fixed['dmdh_emu_per_T']*1e4:.2f}"
            })
            
            rows_tabs4.append({
                'Group': bname,
                'Sample': sname,
                'Diameter_um': d,
                'Hc_fixed_T': f"{m_fixed['H_c_T']:.4f}",
                'Hc_tailfit_T': f"{m_tf['H_c_T']:.4f}",
                'Mr50_fixed': f"{m_fixed['Mr_50_ratio']:.3f}",
                'Mr50_tailfit': f"{m_tf['Mr_50_ratio']:.3f}"
            })
            
        # Summary row
        rows_tab1.append({
            'Group': bname,
            'Sample': f"{bname} Mean +/- SD",
            'Diameter_um': '-',
            'Hc_T': f"{np.mean(b_hc):.4f} +/- {np.std(b_hc, ddof=1):.4f}",
            'Mr_H0_ratio': f"{np.mean(b_mr0):.3f} +/- {np.std(b_mr0, ddof=1):.3f}",
            'Mr_50mT_ratio': f"{np.mean(b_mr50):.3f} +/- {np.std(b_mr50, ddof=1):.3f}",
            'Mcorr_7T_memu': f"{np.mean(b_mcorr):.3f} +/- {np.std(b_mcorr, ddof=1):.3f}",
            'dmdh_1e-4_emu_T': f"{np.mean(b_slope):.2f} +/- {np.std(b_slope, ddof=1):.2f}"
        })
        
    # Write Table 1
    with open(out_tab1_path, 'w', newline='', encoding='utf-8') as fp:
        fields = ['Group', 'Sample', 'Diameter_um', 'Hc_T', 'Mr_H0_ratio', 'Mr_50mT_ratio', 'Mcorr_7T_memu', 'dmdh_1e-4_emu_T']
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows_tab1)
        
    # Write Table S4
    with open(out_tabs4_path, 'w', newline='', encoding='utf-8') as fp:
        fields_s4 = ['Group', 'Sample', 'Diameter_um', 'Hc_fixed_T', 'Hc_tailfit_T', 'Mr50_fixed', 'Mr50_tailfit']
        writer = csv.DictWriter(fp, fieldnames=fields_s4)
        writer.writeheader()
        writer.writerows(rows_tabs4)


def make_all_tables(data_raw_dir='data/raw', out_dir='outputs/tables', out_tables_dir=None):
    if out_tables_dir is not None:
        out_dir = out_tables_dir
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. Sample Map
    squid_dir = os.path.join(data_raw_dir, 'squid_dat')
    make_sample_map(squid_dir, os.path.join(out_dir, 'sample_map.csv'))
    
    # 2. Table 1 & Table S4
    make_table_1_and_s4(
        squid_dir,
        os.path.join(out_dir, 'Table_1_SQUID_summary.csv'),
        os.path.join(out_dir, 'Table_S4_tailfit_summary.csv')
    )
    
    # 3. Table S6 Mesh Convergence
    sim_dir = os.path.join(data_raw_dir, 'r1_sim')
    mesh_table = parse_mesh_convergence_table(sim_dir)
    with open(os.path.join(out_dir, 'Table_S6_mesh_convergence.csv'), 'w', newline='', encoding='utf-8') as fp:
        if mesh_table:
            writer = csv.DictWriter(fp, fieldnames=list(mesh_table[0].keys()))
            writer.writeheader()
            writer.writerows(mesh_table)
            
    # 4. XRD Metadata
    xrd_dir = os.path.join(data_raw_dir, 'xrd_xy')
    xrd_meta = parse_xrd_metadata_table(xrd_dir)
    with open(os.path.join(out_dir, 'xrd_metadata.csv'), 'w', newline='', encoding='utf-8') as fp:
        if xrd_meta:
            writer = csv.DictWriter(fp, fieldnames=list(xrd_meta[0].keys()))
            writer.writeheader()
            writer.writerows(xrd_meta)
            
    # 5. SEM Catalog
    sem_dir = os.path.join(data_raw_dir, 'sem_tiff')
    if os.path.exists(sem_dir):
        sem_cat = generate_sem_catalog(sem_dir)
        with open(os.path.join(out_dir, 'sem_catalog.csv'), 'w', newline='', encoding='utf-8') as fp:
            if sem_cat:
                writer = csv.DictWriter(fp, fieldnames=list(sem_cat[0].keys()))
                writer.writeheader()
                writer.writerows(sem_cat)
                
    print("All tables successfully generated in", out_dir)


generate_all_tables = make_all_tables


if __name__ == '__main__':
    make_all_tables()
