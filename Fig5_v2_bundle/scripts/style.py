import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, re
MM=1/25.4
DOUBLE=180*MM; SINGLE=86*MM
D_COL={3:"#E0A030",5:"#C8553D",8:"#7E3F8F",10:"#2E3A87"}
GREY="#9AA0A6"; INK="#1F2328"; ACC="#B3261E"
plt.rcParams.update({
 "font.family":"sans-serif","font.sans-serif":["Liberation Sans","Arial","Helvetica"],
 "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.labelsize":8,"axes.titlesize":8,
 "xtick.labelsize":7,"ytick.labelsize":7,"legend.fontsize":7,
 "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
 "xtick.minor.width":0.4,"ytick.minor.width":0.4,"xtick.major.size":2.5,"ytick.major.size":2.5,
 "xtick.direction":"in","ytick.direction":"in","xtick.top":True,"ytick.right":True,
 "lines.linewidth":1.0,"legend.frameon":False,"savefig.dpi":600,"image.interpolation":"none",
 "mathtext.fontset":"custom","mathtext.rm":"Liberation Sans","mathtext.it":"Liberation Sans:italic","mathtext.bf":"Liberation Sans:bold",
})
def plabel(ax,s,x=-0.16,y=1.04):
    s=s if s.startswith("(") else f"({s})"
    ax.text(x,y,s,transform=ax.transAxes,fontsize=9,fontweight="bold",va="bottom",ha="left")
def save(fig,path):
    fig.savefig(path,format="svg",bbox_inches="tight",pad_inches=0.02)
    s=open(path).read()
    s=re.sub(r"font-family:\s*'?Liberation Sans'?","font-family:Helvetica, Arial, sans-serif",s)
    s=s.replace("'Liberation Sans'","Helvetica, Arial, sans-serif").replace("Liberation Sans","Helvetica")
    open(path,"w").write(s)
    import svgfix; svgfix.fix(path)
    fig.savefig(path.replace(".svg",".png"),dpi=300,bbox_inches="tight",pad_inches=0.02)
    fig.savefig(path.replace(".svg",".pdf"),bbox_inches="tight",pad_inches=0.02)
