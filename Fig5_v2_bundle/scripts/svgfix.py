"""Replace TeX-encoded mathtext glyphs (cmsy10/cmex10/STIX) in matplotlib SVGs with Unicode in a common font."""
import re,sys
MAP={'cmsy10':{'/':'∝','¢':'·','Ã':'×','Ä':'±','°':'→','Ê':'≤','Ë':'≥'},'cmex10':{'\xad':'⟨','\xae':'⟩','p':'√'}}
FONT="'Cambria Math', 'Segoe UI Symbol', 'STIX Two Math', 'Apple Symbols', 'DejaVu Math TeX Gyre', 'DejaVu Sans', sans-serif"
def fix(path):
    s=open(path,encoding='utf-8').read(); n=0
    def rep(m):
        nonlocal n
        fam,body=m.group(1),m.group(3)
        mp=MAP.get(fam,{}); new=''.join(mp.get(c,c) for c in body); n+=1
        return m.group(0).replace(f"'{fam}'",FONT).replace('>'+body+'<','>'+new+'<')
    s=re.sub(r"font-family: '(cmsy10|cmex10|cmmi10|cmr10|STIXGeneral)'([^>]*)>([^<]*)<",rep,s)
    open(path,'w',encoding='utf-8').write(s); return n
if __name__=='__main__':
    for p in sys.argv[1:]: print(p,fix(p))
