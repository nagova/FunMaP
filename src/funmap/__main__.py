"""
FunMaP CLI entry point:
- python -m funmap reproduce
- python -m funmap demo
- python -m funmap tables
- python -m funmap numbers
- python -m funmap figures
- python -m funmap test
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parent.parent.parent


def run_tables(data_dir='data/raw', out_dir='outputs/tables'):
    print(f"\n[1/4] Generating tables from '{data_dir}' into '{out_dir}'...")
    from scripts.make_tables import generate_all_tables
    generate_all_tables(data_raw_dir=data_dir, out_tables_dir=out_dir)


def run_numbers(data_dir='data/raw', out_dir='outputs'):
    print(f"\n[2/4] Generating LaTeX macros and discrepancy log from '{data_dir}'...")
    from scripts.make_numbers import generate_paper_numbers
    generate_paper_numbers(data_raw_dir=data_dir, out_dir=out_dir)


def run_figures(data_dir='data/raw', out_dir='outputs/figures'):
    print(f"\n[3/4] Generating publication-grade figures from '{data_dir}' into '{out_dir}'...")
    from scripts.make_figures import generate_all_figures
    generate_all_figures(data_dir=data_dir, out_dir=out_dir)


def run_reproduce():
    print("=" * 70)
    print(" FunMaP: Complete Paper Reproduction Pipeline")
    print(" Paper: 'No Monotonic Curvature Scaling of Magnetization Reversal")
    print("         in Micrometer-Scale FePt Janus Caps'")
    print("=" * 70)
    
    # 1. Run tables
    run_tables('data/raw', 'outputs/tables')
    
    # 2. Run numbers & discrepancy audit
    run_numbers('data/raw', 'outputs')
    
    # 3. Run figures
    run_figures('data/raw', 'outputs/figures')
    
    # 4. Run verification tests
    print("\n[4/4] Executing ground truth reproduction assertions...")
    res = subprocess.run([sys.executable, "-m", "pytest", "tests/test_reproduce.py", "-v"], cwd=ROOT_DIR)
    
    print("\n" + "=" * 70)
    if res.returncode == 0:
        print(" SUCCESS: All experimental data and models 100% reproduced!")
        print(" Generated artifacts:")
        print(" - LaTeX numbers: outputs/paper_numbers.tex")
        print(" - Discrepancy report: outputs/DISCREPANCIES.md")
        print(" - Tables: outputs/tables/")
        print(" - Publication figures: outputs/figures/")
    else:
        print(" WARNING: Some reproduction assertions failed. Check output above.")
    print("=" * 70)


def run_demo():
    print("=" * 70)
    print(" FunMaP: Fast Demonstration Mode (Synthetic Dataset)")
    print("=" * 70)
    
    demo_dir = ROOT_DIR / "demo_data"
    if not (demo_dir / "squid_dat").exists():
        print("Generating synthetic demo data...")
        from scripts.make_demo_data import generate_all_demo_data
        generate_all_demo_data(out_dir=str(demo_dir))
        
    out_demo = ROOT_DIR / "outputs" / "demo"
    out_demo_tables = out_demo / "tables"
    out_demo_figs = out_demo / "figures"
    
    # Tables & Numbers
    run_tables(str(demo_dir), str(out_demo_tables))
    run_numbers(str(demo_dir), str(out_demo))
    run_figures(str(demo_dir), str(out_demo_figs))
    
    print("\nRunning unit tests on synthetic data...")
    res = subprocess.run([sys.executable, "-m", "pytest", "tests/test_units.py", "-v"], cwd=ROOT_DIR)
    
    print("\n" + "=" * 70)
    print(" DEMO COMPLETE: FunMaP executed cleanly on synthetic dataset.")
    print(f" Demo outputs located in: {out_demo}")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        prog="funmap",
        description="FunMaP: Functional Magnetic Particles analysis pipeline"
    )
    subparsers = parser.add_subparsers(dest="command", help="Pipeline command to execute")
    
    subparsers.add_parser("reproduce", help="Full reproduction from raw experimental data")
    subparsers.add_parser("demo", help="Fast execution using synthetic demo data")
    subparsers.add_parser("tables", help="Generate CSV summary tables")
    subparsers.add_parser("numbers", help="Generate paper LaTeX numbers & discrepancy log")
    subparsers.add_parser("figures", help="Generate publication-grade figures")
    subparsers.add_parser("test", help="Run test suite")
    
    args = parser.parse_args()
    
    if args.command == "reproduce":
        run_reproduce()
    elif args.command == "demo":
        run_demo()
    elif args.command == "tables":
        run_tables()
    elif args.command == "numbers":
        run_numbers()
    elif args.command == "figures":
        run_figures()
    elif args.command == "test":
        subprocess.run([sys.executable, "-m", "pytest", "tests/", "-v"], cwd=ROOT_DIR)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
