# Zamkniecie kanalu E1 wg 296b w postaci OGOLNEJ: J_phi + n <= 1
# (warunek konieczny dla najlepszego kierunku).  Mierzony moment to czas
# symulowany pierwszego checkpointu, na ktorym warunek zachodzi.
import re,glob,os
D='/tmp/claude-1000/-home-teddy-Projekty-Positronium/5981d9ee-dfbe-485a-8f93-5c6d87fd7653/scratchpad/close/'
def load():
    res={}
    for f in sorted(glob.glob(D+'*.err')):
        cfg,ph,seed=os.path.basename(f)[:-4].split('-')
        if not os.path.exists(f[:-4]+'.out'): continue
        out=open(f[:-4]+'.out').read()
        m=re.search(r'fotony=(\d+) odmowy ceiling=(\d+)',out)
        if not m: continue              # przebieg jeszcze trwa
        rows=[]
        for ln in open(f):
            if not ln.startswith('ACTION'): continue
            d=dict(re.findall(r'(\w+)=([-\d.e+]+)',ln))
            rows.append((float(d['t']),float(d['n']),float(d['Jphi_over_h'])))
        close=next(((t,n,J) for t,n,J in rows if J+n<=1.0),None)
        res[(cfg,int(ph),int(seed))]=dict(close=close,photons=int(m.group(1)),
            refused=int(m.group(2)),last=rows[-1])
    return res
if __name__=='__main__':
    for k,v in sorted(load().items()):
        c=v['close']
        s=('ZAMKN t=%9.3f ps n=%.4f J=%.4f'%(c[0]*1e12,c[1],c[2]) if c
           else 'otwarty  (n=%.4f J=%.4f)'%(v['last'][1],v['last'][2]))
        print('%-4s %d %2d  %s  fot=%d odm=%d'%(k[0],k[1],k[2],s,v['photons'],v['refused']))
