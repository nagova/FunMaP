"""Switching intervals for all simulated loops -> SI table (CSV + LaTeX).
H_c is reported as [|H_pre|, |H_post|]: last field with <m_z> > 0 and first with <m_z> < 0
on the descending branch. Step = interval width. Runs that never switch are flagged.
Usage: python sim_intervals.py [--data DIR] [--out DIR]"""
import os, csv, glob, argparse
HERE=os.path.dirname(os.path.abspath(__file__))
ap=argparse.ArgumentParser(); ap.add_argument('--data',default=os.path.join(HERE,'..','data')); ap.add_argument('--out',default=os.path.join(HERE,'..','output'))
A=ap.parse_args(); os.makedirs(A.out,exist_ok=True)
rows=[]
for f in sorted(glob.glob(os.path.join(A.data,'sim_loops','*.csv'))):
    r=list(csv.reader(open(f))); h=r[0]; d=[list(map(float,x)) for x in r[1:] if x]
    iz=h.index('Mz'); mr=[x[iz] for x in d if x[0]==0.0][0]; iv=None
    for a,b in zip(d,d[1:]):
        if a[iz]>0 and b[iz]<0: iv=(abs(a[0]),abs(b[0])); break
    rows.append(dict(run=os.path.basename(f)[:-4],Mr_over_Ms=round(mr,3),
        H_pre_T=iv[0] if iv else '',H_post_T=iv[1] if iv else '',step_T=round(iv[1]-iv[0],3) if iv else '',
        note='' if iv else f'no switch down to {abs(d[-1][0])} T'))
with open(os.path.join(A.out,'sim_switching_intervals.csv'),'w',newline='') as fo:
    w=csv.DictWriter(fo,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
with open(os.path.join(A.out,'sim_switching_intervals.tex'),'w') as fo:
    fo.write('\\begin{tabular}{lccc}\n\\hline\nRun & $M_r/M_s$ & $\\mu_0H_c$ interval (T) & step (T)\\\\\n\\hline\n')
    for x in rows:
        iv=f"[{x['H_pre_T']:.2f}, {x['H_post_T']:.2f}]" if x['H_pre_T']!='' else x['note']
        fo.write(f"{x['run'].replace('_','\\_')} & {x['Mr_over_Ms']:.3f} & {iv} & {x['step_T']}\\\\\n")
    fo.write('\\hline\n\\end{tabular}\n')
for x in rows: print(x)
