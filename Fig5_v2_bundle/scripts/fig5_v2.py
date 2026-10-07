"""Fig. 5 (v2): R1 benchmark + curvature test + chemical-order test.
All simulated coercive fields are switching INTERVALS [|H_pre|, |H_post|]
(last unswitched / first switched field on the descending branch); no interpolation.
Usage: python fig5_v2.py [--data DIR] [--out DIR]   (defaults: ../data, ../output)"""
import sys, os, argparse, matplotlib.ticker
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
ap=argparse.ArgumentParser(); ap.add_argument('--data',default=os.path.join(HERE,'..','data')); ap.add_argument('--out',default=os.path.join(HERE,'..','output'))
A=ap.parse_args(); DATA=os.path.abspath(A.data); OUT=os.path.abspath(A.out); os.makedirs(OUT,exist_ok=True)
from style import *; from icons import *
import numpy as np, pandas as pd
from matplotlib.colors import TwoSlopeNorm
R=os.path.join(DATA,'R1_Benchmark')+'/'
br=pd.read_csv(R+'descending_branch_B_vs_Mz_cell1p0nm.csv')
fig=plt.figure(figsize=(DOUBLE,170*MM))
gs=fig.add_gridspec(3,5,height_ratios=[1,1.0,0.66],width_ratios=[1.25,1,1,1,0.05],hspace=0.62,wspace=0.32)
ax=fig.add_subplot(gs[0,0]); plabel(ax,'a',-0.3)
b=br[br.B_ext_T<=0]
ax.plot(b.B_ext_T,b.Mz,'-',color=INK,lw=0.8); ax.plot(b.B_ext_T,b.Mz,'o',ms=1.8,color=INK)
for B,l in [(0,'b'),(-6.2,'c'),(-6.3,'d')]:
    y=br.loc[np.isclose(br.B_ext_T,B),'Mz'].values[0]; ax.plot(B,y,'o',ms=4,mfc='none',mec=ACC,mew=1)
    ax.annotate(l,(B,y),(B+(0.5 if B<0 else -0.9),y+(0.12 if l!='d' else 0.15)),fontsize=7,fontweight='bold',color=ACC)
ax.set_xlim(-7.2,0.4); ax.set_ylim(-1.05,1.0); ax.axhline(0,color='0.8',lw=0.4)
ax.set_xlabel('$\\mu_0H$ (T)'); ax.set_ylabel('$\\langle m_z\\rangle$')
ax.text(0.47,0.56,'$D$ = 200 nm' + '\n' + r'$\Delta x$ = 1.0 nm',transform=ax.transAxes,fontsize=5.8,color='0.35',ha='left')
norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1)
for k,(B,l,lab) in enumerate([('+0.00','b','0 T'),('-6.20','c','−6.20 T'),('-6.30','d','−6.30 T')]):
    d=np.load(R+f'midplane_slices_y0/r1_midplane_slice_y0_B{B}T.npz'); x=d['x']; z=d['z']; msk=d['norm']>0; m=d['m']/np.where(msk,d['norm'],1)[...,None]
    mz=np.where(msk,m[...,2],np.nan)
    a=fig.add_subplot(gs[0,k+1]); plabel(a,l,-0.08)
    im=a.pcolormesh(x,z,mz.T,cmap='RdBu_r',norm=norm,shading='nearest',rasterized=True)
    st=10; X,Z=np.meshgrid(x[::st],z[::st],indexing='ij'); mm=msk[::st,::st]
    a.quiver(X[mm],Z[mm],m[::st,::st,0][mm],m[::st,::st,2][mm],scale=1/8.5,scale_units='xy',angles='xy',width=0.008,headwidth=3.2,headlength=3.5,pivot='mid',color=INK)
    a.set_aspect('equal'); a.set_xlim(-120,120); a.set_ylim(0,170); a.set_title(lab,fontsize=7,pad=2)
    a.set_xlabel('$x$ (nm)'); 
    if k==0: a.set_ylabel('$z$ (nm)')
    else: a.set_yticklabels([])
    mzv=br.loc[np.isclose(br.B_ext_T,float(B)),'Mz'].values[0]
    a.text(0.03,0.04,f'$\\langle m_z\\rangle$ = {mzv:+.3f}',transform=a.transAxes,fontsize=5.8)
    if k==1: a.annotate('',(108,108),(108,158),arrowprops=dict(arrowstyle='->',lw=0.8,color=ACC)); a.text(101,135,'$B$',color=ACC,fontsize=6.5,ha='right')
cax=fig.add_subplot(gs[0,4]); cb=fig.colorbar(im,cax=cax,ticks=[-1,0,1]); cb.set_label('$m_z$',labelpad=1); cb.outline.set_linewidth(0.5)
# ---- (e) Hc vs D : simulation (wedges, R1 3D) and experiment
import glob as _g, csv as _csv
L=os.path.join(DATA,'sim_loops')+'/'
def bracket(f):
    r=list(_csv.reader(open(f))); d=[list(map(float,x)) for x in r[1:] if x]
    for a_,b_ in zip(d,d[1:]):
        if a_[3]>0 and b_[3]<0: return abs(a_[0]),abs(b_[0])
    return None
WD={0.2:'200nm',0.5:'500nm',1:'1um',3:'3um',5:'5um',8:'8um',10:'10um'}
ev=fig.add_subplot(gs[1,0:2]); plabel(ev,'e',-0.16)
for Dm,tag in WD.items():
    lo,hi=bracket(L+f'L10_pure_wedge_D{tag}_cell1p3nm.csv')
    ev.errorbar(Dm,(lo+hi)/2,yerr=(hi-lo)/2,fmt='s',ms=3.6,color=ACC,elinewidth=0.8,capsize=1.5)
lo,hi=bracket(L+'R1_hemi_d200nm_cell_1p0nm.csv'); ev.errorbar(0.2*0.85,(lo+hi)/2,yerr=(hi-lo)/2,fmt='D',ms=3.6,color=INK,elinewidth=0.8,capsize=1.5,zorder=4)
for f,Dm in [('M1_wedge_D3um_cell1p3nm_fA1_10pct',3),('M1_wedge_D3um_cell1p3nm_fA1_20pct',3),('M6_snapshot_rerun_D10um_cell1p3nm_15pctA1',10)]:
    lo,hi=bracket(L+f+'.csv'); ev.errorbar(Dm*1.13,(lo+hi)/2,yerr=(hi-lo)/2,fmt='^',ms=3.6,mfc='white',mec='#c0605a',color='#c0605a',elinewidth=0.8,capsize=1.5)
X=pd.read_csv(os.path.join(DATA,'squid','fig4_extracted_values.csv'))
for _,r in X[X.b!='AD'].iterrows():
    ev.plot(r.d*(0.95 if r.b=='B1' else 1.05),r.Hc,'o',ms=4.2,mfc=D_COL[r.d] if r.b=='B1' else 'white',mec=D_COL[r.d],mew=1.0,zorder=3)
ev.axhspan(1.00,1.19,color=D_COL[10],alpha=0.13,lw=0)
ev.set_xscale('log'); ev.set_yscale('log'); ev.set_xlim(0.13,14); ev.set_ylim(0.6,10)
ev.set_xticks([0.2,0.5,1,2,5,10]); ev.set_xticklabels(['0.2','0.5','1','2','5','10'])
ev.set_yticks([1,2,5,10]); ev.set_yticklabels(['1','2','5','10']); ev.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ev.set_xlabel('sphere diameter $D$ (µm)'); ev.set_ylabel('$\\mu_0H_c$ (T)')
ev.text(0.62,7.6,'pure L1$_0$, simulated',color=ACC,fontsize=6.3,ha='center')
ev.text(0.15,4.9,'3D hemisphere (R1)',color=INK,fontsize=5.8,ha='left',va='top')
ev.text(4.2,3.75,'10–20 % soft\nphase (sim.)',color='#c0605a',fontsize=5.8,ha='center',va='top',linespacing=1.0)
ev.text(0.62,0.8,'measured caps (Batch 1 ●, Batch 2 ○)',color=D_COL[10],fontsize=6.3,ha='center')
ev.annotate('',(0.3,1.22),(0.3,5.9),arrowprops=dict(arrowstyle='->',lw=0.9,color='0.35'))
ev.text(0.32,2.6,'×5.6',fontsize=6.5,color='0.3',fontweight='bold')
# ---- (f) Hc vs mean anisotropy K/K0 = 1 - f_A1
fv=fig.add_subplot(gs[1,2:4]); plabel(fv,'f',-0.14)
fr=[0,10,20,30,40,50,66]; mid=[]
for p in fr:
    lo,hi=bracket(L+f'M1_wedge_D3um_cell1p3nm_fA1_{p:02d}pct.csv'); k=1-p/100; mid.append((k,(lo+hi)/2))
    fv.errorbar(k,(lo+hi)/2,yerr=(hi-lo)/2,fmt='s',ms=3.6,color=ACC,elinewidth=0.8,capsize=1.5,zorder=3)
mid=np.array(mid); cf=np.polyfit(mid[:,0],mid[:,1],1)
kk=np.linspace(0,1,50); fv.plot(kk,np.polyval(cf,kk),'--',color=ACC,lw=0.7)
S,dS=0.70,0.09; k0=S**2; klo,khi=(S-dS)**2,(S+dS)**2
fv.axvspan(klo,khi,color=CONCEPT['struct'],alpha=0.16,lw=0); fv.axvline(k0,color=CONCEPT['struct'],lw=0.8)
fv.axhspan(1.00,1.19,color=D_COL[10],alpha=0.18,lw=0)
hS=np.polyval(cf,k0); hlo,hhi=np.polyval(cf,klo),np.polyval(cf,khi)
kb=[(1.00-cf[1])/cf[0],(1.19-cf[1])/cf[0]]
fv.plot(k0,hS,'o',ms=4.5,mfc='white',mec=INK,mew=1,zorder=4)
fv.set_xlim(0,1.04); fv.set_ylim(0,7.2)
fv.set_xlabel('mean anisotropy $\\langle K_u\\rangle/K_{\\mathrm{L1_0}}$  (= 1 − soft fraction)'); fv.set_ylabel('$\\mu_0H_c$ (T)')
fv.text(k0+0.015,6.55,'XRD: $S$ = 0.70 ± 0.09\n$K_u\\propto S^2$ → 0.49',fontsize=5.9,color=CONCEPT['struct'],va='top',linespacing=1.05)
fv.text(k0-0.02,hS+0.1,f'{hlo:.1f}–{hhi:.1f} T',fontsize=6.2,va='bottom',ha='right',color=INK)
fv.text(1.0,0.55,'measured 1.00–1.19 T',fontsize=6.2,color=D_COL[10],ha='right')
fv.text(0.98,2.3,'simulated, $D$ = 3 µm,\n0 K, soft cells\nbelow exchange length',fontsize=5.8,color=ACC,ha='right',linespacing=1.0)
fv.annotate('',(k0+0.045,1.25),(k0+0.045,hS-0.35),arrowprops=dict(arrowstyle='->',lw=0.9,color='0.35'))
fv.text(k0+0.06,1.75,f'×{hS/1.10:.1f} not explained\nby measured order',fontsize=5.8,color='0.3',ha='left',va='center',linespacing=1.0)
print('fit',cf,'H(S^2)',hS,hlo,hhi,'band crossing K/K0',kb)
# (g) ladder
e=fig.add_subplot(gs[2,1:4]); plabel(e,'g',-0.62,1.0); 
rows=[('anisotropy field $H_K = 2K_u/\\mu_0M_s$',13.2,INK),('Stoner–Wohlfarth minimum $H_K/2$',6.6,INK),('R1 ideal cap (simulation)',6.22,ACC),
      ('annealed caps, all $D$ (Batch 1 ●, Batch 2 ○)',None,D_COL[10]),('as-deposited A1',None,GREY)]
ys=np.arange(len(rows))[::-1]
for (lab,v,c),y in zip(rows,ys):
    e.text(-0.02,y,lab,va='center',ha='right',fontsize=6.5,transform=e.get_yaxis_transform())
    if v: e.plot(v,y,'o',ms=4.5,color=c); e.text(v*1.2,y,f'{v:g} T',va='center',fontsize=6.3,color=c)
ya=ys[3]; e.fill_betweenx([ya-0.22,ya+0.22],1.00,1.19,color=D_COL[10],alpha=0.25,lw=0)
e.plot([1.147],[ya],'o',ms=4.5,color=INK); e.plot([1.063],[ya],'o',ms=4.5,mfc='white',mec=INK); e.text(1.19*1.2,ya,'1.00–1.19 T',va='center',fontsize=6.3,color=D_COL[10])
yd=ys[4]; e.annotate('',(0.010,yd),(0.002,yd),arrowprops=dict(arrowstyle='<-',lw=0.9,color=GREY,ls='--')); e.text(0.011,yd,'< 10 mT',va='center',fontsize=6.3,color='0.4')
e.set_xscale('log'); e.set_xlim(0.001,20); e.set_ylim(-0.7,len(rows)-0.4); e.set_yticks([])
for s in ['left','right','top']: e.spines[s].set_visible(False)
e.tick_params(which='both',top=False,right=False); e.set_xlabel('$\\mu_0H_c$ (T, log scale)')
e.set_xticks([0.001,0.01,0.1,1,10]); e.set_xticklabels(['1 mT','10 mT','0.1 T','1 T','10 T'])
def arr(x1,x2,y,c,t):
    e.annotate('',(x2,y),(x1,y),arrowprops=dict(arrowstyle='->',lw=1.0,color=c)); e.text(np.sqrt(x1*x2),y+0.12,t,ha='center',va='bottom',fontsize=6.3,color=c,fontweight='bold')
arr(0.006,1.0,(ys[3]+ys[4])/2,'#1f6f4a','annealing: > 100×')
yg=(ys[2]+ys[3])/2; e.annotate('',(6.0,yg),(1.12,yg),arrowprops=dict(arrowstyle='->',lw=1.0,color=ACC)); e.text(0.95,yg,'gap to ideal cap ~5.6×',ha='right',va='center',fontsize=6.3,color=ACC,fontweight='bold')
c3=ax.inset_axes([0.42,0.04,0.56,0.46]); c3.set_xlim(-2.0,2.7); c3.set_ylim(-0.6,2.2); c3.set_aspect('equal'); c3.axis('off')
cap3d(c3,0,0,1.0,t=0.6)
for a_ in np.deg2rad([-60,-30,0,30,60]):
    r0=1.0+0.6*np.cos(a_)*0.55; arrow(c3,(r0*np.sin(a_)*0.85,r0*np.cos(a_)*0.85),((r0+0.32)*np.sin(a_)*0.85,(r0+0.32)*np.cos(a_)*0.85),color='white',lw=0.7,ms=4,z=6)

c3.text(0.25,-0.45,'radial easy axis $\\hat{u}=\\hat{r}$, $t_0$ = 60 nm',ha='center',va='top',fontsize=5.6,color='0.3')
import numpy as _np
lx,hx=_np.log10(0.001),_np.log10(20); ylo,yhi=e.get_ylim()
def fx(v): return (_np.log10(v)-lx)/(hx-lx)
def fy(y): return (y-ylo)/(yhi-ylo)
for (v,y,kind) in [(13.2,ys[0],'cell'),(6.6,ys[1],'grain'),(6.22,ys[2],'cap'),(1.19,ys[3],'loop'),(0.010,ys[4],'loop0')]:
    xo=fx(v)-0.075 if kind!='loop0' else fx(v)+0.115
    if kind=='loop': xo=fx(1.0)-0.075
    ii=e.inset_axes([xo,fy(y)-0.07,0.045,0.14]); ii.set_xlim(-1,1); ii.set_ylim(-1,1); ii.set_aspect('equal'); ii.axis('off')
    if kind=='cell': unit_cell(ii,-0.75,-0.85,1.2,'L10',r=0.12)
    elif kind=='grain': grain_icon(ii,0,0,1.5,angle=45)
    elif kind=='cap': glyph_mesh_cap(ii,0,-0.6,0.9,CONCEPT['sim'])
    elif kind=='loop': glyph_loop(ii,0,0,1.8,1.4,D_COL[10],hc=0.25)
    else: glyph_loop(ii,0,0,1.8,1.4,'0.45',hc=0.02)
save(fig,os.path.join(OUT,'Fig5_simulation_v2.svg'))
