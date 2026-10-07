"""Shared vector icons for FunMaP figures (matplotlib patches -> editable SVG)."""
import numpy as np
from matplotlib.patches import Circle, Polygon, Rectangle, FancyBboxPatch, FancyArrowPatch, Ellipse
from style import D_COL, INK, GREY
CAP="#C07A3A"; SPHERE="#E4E6EA"; SPHERE_EDGE="#9AA0A6"
CONCEPT={"fab":"#8C6D4F","struct":"#1B7F79","mag":"#B3261E","sim":"#4A6FA5","stats":"#555555"}

def cap_icon(ax,x,y,r,angle=0,t=0.28,ring=None,z=3,lw=0.5):
    """Janus sphere: SiO2 core + FePt cap with t(θ)=t0 cosθ, pointing along `angle` (deg, 0=up)."""
    ax.add_patch(Circle((x,y),r,fc=SPHERE,ec=SPHERE_EDGE,lw=lw,zorder=z))
    ax.add_patch(Ellipse((x-0.35*r,y+0.35*r),0.5*r,0.32*r,angle=35,fc='white',ec='none',alpha=0.8,zorder=z+0.1))
    th=np.linspace(-np.pi/2,np.pi/2,60); a=np.deg2rad(angle)
    outer=np.c_[(r+t*r*np.cos(th))*np.sin(th),(r+t*r*np.cos(th))*np.cos(th)]
    inner=np.c_[r*np.sin(th),r*np.cos(th)][::-1]
    P=np.r_[outer,inner]; R=np.array([[np.cos(a),np.sin(a)],[-np.sin(a),np.cos(a)]])
    P=P@R.T+[x,y]; ax.add_patch(Polygon(P,closed=True,fc=CAP,ec='#8A5226',lw=lw*0.6,zorder=z+0.2))
    if ring: ax.add_patch(Circle((x,y),r*(1+t)*1.12,fc='none',ec=ring,lw=1.0,zorder=z-0.1))

def monolayer(ax,x0,y0,w,n=6,r=None,caps=True,z=3):
    r=r or w/(2*n); dy=r*np.sqrt(3)
    for j in range(3):
        for i in range(n-(j%2)):
            xx=x0+r+i*2*r+(r if j%2 else 0); yy=y0+r+j*dy
            if caps: cap_icon(ax,xx,yy,r*0.92,t=0.25,z=z,lw=0.3)
            else:
                ax.add_patch(Circle((xx,yy),r*0.92,fc=SPHERE,ec=SPHERE_EDGE,lw=0.3,zorder=z))

def box(ax,x,y,w,h,color,title,z=1):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.0,rounding_size=1.2",fc=color,ec=color,alpha=0.10,lw=0,zorder=z))
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.0,rounding_size=1.2",fc='none',ec=color,lw=0.7,zorder=z))
    ax.text(x+w/2,y+h-1.2,title,ha='center',va='top',fontsize=6.6,fontweight='bold',color=color,zorder=z+1)

def arrow(ax,p0,p1,color=INK,lw=0.8,ms=6,z=5,style='-|>'):
    ax.add_patch(FancyArrowPatch(p0,p1,arrowstyle=style,mutation_scale=ms,lw=lw,color=color,zorder=z))

def glyph_loop(ax,x,y,w,h,color,hc=0.3):
    H=np.linspace(-1,1,60)
    ax.plot(x+H*w/2,y+np.tanh((H+hc)/0.15)*h/2,color=color,lw=0.8)
    ax.plot(x+H*w/2,y+np.tanh((H-hc)/0.15)*h/2,color=color,lw=0.8)
    ax.plot([x-w/2,x+w/2],[y,y],color='0.6',lw=0.3); ax.plot([x,x],[y-h/2,y+h/2],color='0.6',lw=0.3)

def glyph_xrd(ax,x,y,w,h,color):
    t=np.linspace(0,1,120); s=0.15*np.exp(-t*2)+0.9*np.exp(-((t-0.62)/0.02)**2)+0.25*np.exp(-((t-0.45)/0.015)**2)
    ax.plot(x-w/2+t*w,y-h/2+s*h,color=color,lw=0.8)

def glyph_sem(ax,x,y,w,h,color):
    ax.add_patch(Rectangle((x-w/2,y-h/2),w,h,fc='#6E7378',ec=color,lw=0.7,zorder=3))
    for i in range(3):
        for j in range(2):
            ax.add_patch(Circle((x-w/2+w*(0.2+0.3*i)+(0.15*w if j else 0),y-h/2+h*(0.3+0.42*j)),h*0.19,fc='#D9DCE0',ec='none',zorder=3.1))

def glyph_regression(ax,x,y,w,h,color):
    xs=np.linspace(-0.4,0.4,4)
    for k,dy in enumerate([0.12,-0.12]):
        ys=0.15*xs+dy+np.array([0.03,-0.04,0.02,-0.01])
        ax.plot(x+xs*w,y+ys*h*2,'o',ms=2.2,mfc=color if k==0 else 'white',mec=color,mew=0.6)
        ax.plot(x+np.array([-0.45,0.45])*w,y+(0.15*np.array([-0.45,0.45])+dy)*h*2,'-',color=color,lw=0.6)

def glyph_mesh_cap(ax,x,y,r,color):
    th=np.linspace(-np.pi/2,np.pi/2,40)
    for f in [1.0,1.12,1.24]:
        ax.plot(x+r*f*np.sin(th)*np.sqrt(np.clip(np.cos(th),0,1))**0,y+r*f*np.cos(th),color=color,lw=0.4)
    for a in np.linspace(-np.pi/2,np.pi/2,9):
        ax.plot([x+r*np.sin(a),x+1.24*r*np.sin(a)],[y+r*np.cos(a),y+1.24*r*np.cos(a)],color=color,lw=0.4)

def furnace(ax,x,y,w,h,color):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0,rounding_size=0.8",fc='#F3E3D3',ec=color,lw=0.7,zorder=3))
    zz=np.linspace(-w/2+1,w/2-1,40); ax.plot(x+zz,y-h/2+1+0.6*np.abs(((zz*2)%2)-1),color='#D9531E',lw=0.6,zorder=3.2)
    ax.plot(x+zz,y+h/2-1-0.6*np.abs(((zz*2)%2)-1),color='#D9531E',lw=0.6,zorder=3.2)

def icon_tube(ax,x,y,s,color='0.45'):
    ax.add_patch(Ellipse((x-s*0.9,y),s*0.5,s*1.0,fc='white',ec=color,lw=0.6,zorder=4))
    ax.plot([x-s*0.9,x+s*0.9],[y+s*0.5,y+s*0.5],color=color,lw=0.6,zorder=4); ax.plot([x-s*0.9,x+s*0.9],[y-s*0.5,y-s*0.5],color=color,lw=0.6,zorder=4)
    ax.add_patch(Ellipse((x+s*0.9,y),s*0.5,s*1.0,fc='0.85',ec=color,lw=0.6,zorder=4))
    ax.add_patch(Ellipse((x+s*0.9,y),s*0.25,s*0.55,fc='white',ec=color,lw=0.5,zorder=4.1))

def icon_helix(ax,x,y,s,color='0.45'):
    t=np.linspace(0,4*np.pi,120); ax.plot(x+(t/(4*np.pi)-0.5)*2*s,y+0.45*s*np.sin(t),color=color,lw=0.8,zorder=4)

def icon_roll(ax,x,y,s,color='0.45'):
    t=np.linspace(0,4.2*np.pi,200); rr=s*0.12*(1+t/(2*np.pi)); ax.plot(x+rr*np.cos(t),y+rr*np.sin(t),color=color,lw=0.7,zorder=4)

# ---------------- crystal structures (oblique projection) ----------------
FE="#B5541C"; PT="#8E949B"
def _proj(p,ox,oy,s,c_over_a=1.0):
    x,y,z=p; return ox+s*(x+0.42*y), oy+s*(c_over_a*z+0.30*y)
def unit_cell(ax,ox,oy,s,kind='L10',c_over_a=0.964,r=0.13,z0=5):
    from matplotlib.patches import Wedge
    P=lambda p:_proj(p,ox,oy,s,c_over_a)
    edges=[((0,0,0),(1,0,0)),((0,0,0),(0,1,0)),((0,0,0),(0,0,1)),((1,0,0),(1,1,0)),((1,0,0),(1,0,1)),((0,1,0),(1,1,0)),((0,1,0),(0,1,1)),((0,0,1),(1,0,1)),((0,0,1),(0,1,1)),((1,1,0),(1,1,1)),((1,0,1),(1,1,1)),((0,1,1),(1,1,1))]
    for a,b in edges:
        (x1,y1),(x2,y2)=P(a),P(b); hidden=(a==(0,1,0) or b==(0,1,0)) and not (a==(0,1,1) or b==(0,1,1))
        ax.plot([x1,x2],[y1,y2],color='0.55',lw=0.5,ls=':' if hidden else '-',zorder=z0)
    corners=[(i,j,k) for i in (0,1) for j in (0,1) for k in (0,1)]
    faces_z=[(0.5,0.5,0),(0.5,0.5,1)]; faces_side=[(0.5,0,0.5),(0.5,1,0.5),(0,0.5,0.5),(1,0.5,0.5)]
    atoms=[(p,'Pt') for p in corners+faces_z]+[(p,'Fe') for p in faces_side]
    atoms.sort(key=lambda t:-t[0][1])   # draw back to front
    for p,el in atoms:
        x,y=P(p); rr=r*s
        if kind=='L10':
            col=PT if el=='Pt' else FE
            ax.add_patch(Circle((x,y),rr,fc=col,ec='white',lw=0.4,zorder=z0+1+(1-p[1])))
        else:  # A1: random occupancy -> half/half
            ax.add_patch(Wedge((x,y),rr,90,270,fc=FE,ec='none',zorder=z0+1+(1-p[1])))
            ax.add_patch(Wedge((x,y),rr,270,90,fc=PT,ec='none',zorder=z0+1+(1-p[1])))
            ax.add_patch(Circle((x,y),rr,fc='none',ec='white',lw=0.4,zorder=z0+1.01+(1-p[1])))

def bragg_brentano(ax,x,y,s,color='0.3'):
    """Tube, sample (wafer + caps), detector at symmetric angle θ."""
    th=np.deg2rad(28); L=s
    ax.add_patch(Rectangle((x-0.45*s,y-0.04*s),0.9*s,0.04*s,fc='#BFC4CA',ec='none',zorder=4))
    for k in range(7): cap_icon(ax,x-0.36*s+k*0.12*s,y+0.045*s,0.05*s,t=0.3,z=5,lw=0.2)
    tx,ty=x-L*np.cos(th),y+L*np.sin(th); dx,dy=x+L*np.cos(th),y+L*np.sin(th)
    arrow(ax,(tx,ty),(x-0.02*s,y+0.06*s),color=color,lw=0.7,ms=5); arrow(ax,(x+0.02*s,y+0.06*s),(dx,dy),color=color,lw=0.7,ms=5)
    ax.add_patch(Rectangle((tx-0.1*s,ty-0.04*s),0.12*s,0.1*s,angle=-28,fc='0.35',ec='none',zorder=4))
    ax.add_patch(Rectangle((dx-0.02*s,dy-0.04*s),0.12*s,0.1*s,angle=28,fc='0.35',ec='none',zorder=4))
    ax.text(tx-0.05*s,ty+0.12*s,'Cu Kα',ha='center',fontsize=5.8); ax.text(dx+0.05*s,dy+0.12*s,'detector',ha='center',fontsize=5.8)
    ax.text(x-0.28*s,y+0.2*s,'θ',fontsize=6.5,style='italic'); ax.text(x+0.22*s,y+0.2*s,'θ',fontsize=6.5,style='italic')

def grain_icon(ax,x,y,s,angle=45,color='0.35'):
    ax.add_patch(Ellipse((x,y),s,0.7*s,fc='#EDE7DF',ec=color,lw=0.6,zorder=4))
    a=np.deg2rad(angle); arrow(ax,(x-0.35*s*np.sin(a),y-0.35*s*np.cos(a)),(x+0.35*s*np.sin(a),y+0.35*s*np.cos(a)),color=color,lw=0.7,ms=4,style='<|-|>')

def cap3d(ax,x,y,R,t=0.3,slice_plane=True,field=True):
    """Pseudo-3D hemispherical cap (shell) with optional y=0 slice plane and field arrow."""
    ang=np.linspace(0,np.pi,100)
    # outer shell silhouette
    ax.fill(np.r_[x+R*(1+t)*np.cos(ang)],np.r_[y+R*(1+t)*np.sin(ang)],color=CAP,alpha=0.95,zorder=3,lw=0)
    ax.add_patch(Ellipse((x,y),2*R*(1+t),0.5*R*(1+t),fc='#D99A62',ec='#8A5226',lw=0.5,zorder=3.5))
    ax.add_patch(Ellipse((x,y),2*R,0.5*R,fc=SPHERE,ec=SPHERE_EDGE,lw=0.5,zorder=3.6))
    for f in [0.35,0.65]:
        ax.plot(x+R*(1+t)*np.cos(ang)*np.sqrt(1-f**2),y+R*(1+t)*f+0.12*R*np.sin(ang)*0,color='#9E5F2C',lw=0.4,zorder=3.4)
    ax.add_patch(Ellipse((x-0.35*R,y+0.75*R),0.5*R,0.28*R,angle=30,fc='white',alpha=0.5,ec='none',zorder=3.7))
    if slice_plane:
        ax.add_patch(Polygon([(x-1.55*R,y-0.25*R),(x+1.55*R,y-0.25*R),(x+1.55*R,y+1.55*R),(x-1.55*R,y+1.55*R)],closed=True,fc='#4A6FA5',alpha=0.10,ec='#4A6FA5',lw=0.6,ls='--',zorder=4))
        ax.text(x+1.5*R,y+1.45*R,'$y$ = 0 slice',ha='right',va='top',fontsize=5.8,color='#4A6FA5')
    if field:
        arrow(ax,(x+1.85*R,y+1.4*R),(x+1.85*R,y+0.3*R),color=CONCEPT['mag'],lw=1.0,ms=6); ax.text(x+1.95*R,y+0.85*R,'$\\mathbf{B}$',fontsize=7,color=CONCEPT['mag'])

def shell_patch_pts(cx,cy,R,t,th1,th2,f1=0.0,f2=1.0,n=12):
    """Polygon for shell sector between polar angles th1..th2 (rad, 0=up) and radial fractions f1..f2 of local thickness t*cos(th)."""
    th=np.linspace(th1,th2,n)
    rin=R+f1*t*R*np.cos(th); rout=R+f2*t*R*np.cos(th)
    P=np.r_[np.c_[cx+rout*np.sin(th),cy+rout*np.cos(th)],np.c_[cx+rin[::-1]*np.sin(th[::-1]),cy+rin[::-1]*np.cos(th[::-1])]]
    return P

def ideal_cap(ax,cx,cy,R,t=0.6):
    ax.add_patch(Circle((cx,cy),R,fc=SPHERE,ec=SPHERE_EDGE,lw=0.5,zorder=2))
    ax.add_patch(Polygon(shell_patch_pts(cx,cy,R,t,-np.deg2rad(85),np.deg2rad(85),n=80),closed=True,fc=CAP,ec='#8A5226',lw=0.5,zorder=3))
    for a in np.deg2rad(np.linspace(-70,70,9)):
        rm=R+0.5*t*R*np.cos(a); x0,y0=cx+rm*np.sin(a),cy+rm*np.cos(a); L=0.17*R
        arrow(ax,(x0-L*np.sin(a)/2,y0-L*np.cos(a)/2),(x0+L*np.sin(a)/2,y0+L*np.cos(a)/2),color='white',lw=0.7,ms=4,z=4)

def poly_cap(ax,cx,cy,R,t=0.6,seed=4,neck=True):
    rng=np.random.default_rng(seed)
    ax.add_patch(Circle((cx,cy),R,fc=SPHERE,ec=SPHERE_EDGE,lw=0.5,zorder=2))
    edges=np.deg2rad(np.sort(np.r_[-85,85,rng.uniform(-80,80,13)]))
    for a1,a2 in zip(edges[:-1],edges[1:]):
        split=rng.uniform(0.35,0.65) if (a2-a1)>np.deg2rad(8) else None
        parts=[(0,split),(split,1)] if split else [(0,1)]
        for f1,f2 in parts:
            shade=rng.uniform(-0.12,0.12); base=np.array([0xC0,0x7A,0x3A])/255; col=np.clip(base*(1+shade),0,1)
            ax.add_patch(Polygon(shell_patch_pts(cx,cy,R,t,a1,a2,f1,f2,n=8),closed=True,fc=col,ec='#6B3F1A',lw=0.35,zorder=3))
            am=(a1+a2)/2; rm=R+(f1+f2)/2*t*R*np.cos(am)
            if t*R*np.cos(am)*(f2-f1)>0.08*R:
                ang=am+np.deg2rad(rng.uniform(-60,60)); L=0.12*R; x0,y0=cx+rm*np.sin(am),cy+rm*np.cos(am)
                arrow(ax,(x0-L*np.sin(ang)/2,y0-L*np.cos(ang)/2),(x0+L*np.sin(ang)/2,y0+L*np.cos(ang)/2),color='white',lw=0.6,ms=3.5,z=4)
    if neck:
        nx=cx+1.02*R; ax.add_patch(Polygon([(nx-0.02*R,cy+0.02*R),(nx+0.16*R,cy+0.06*R),(nx+0.16*R,cy-0.02*R),(nx-0.02*R,cy-0.04*R)],closed=True,fc='#8A5226',ec='none',zorder=3))
        ang=np.linspace(np.pi*0.6,np.pi*1.25,40); ax.plot(cx+2.2*R+R*np.cos(ang),cy+R*np.sin(ang)*0.95,color=SPHERE_EDGE,lw=0.5,zorder=2)
