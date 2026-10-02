import sys,re
sys.path.insert(0,'/tmp/claude-1000/-home-teddy-Projekty-Positronium/5981d9ee-dfbe-485a-8f93-5c6d87fd7653/scratchpad/post')
import an
eps1=3.2625e-06
for f in sys.argv[1:]:
    r=an.rows(f)
    post=[x for x in r if x['n']<0.5 and abs(x['dtot'])<0.01]
    lsall=[l for l in open(f) if 'CREM_LS ' in l]
    dl=[abs(float(re.search(r'L_after=([-\d.e+]+)',l).group(1))-float(re.search(r'L_before=([-\d.e+]+)',l).group(1))) for l in lsall]
    de=[abs(float(re.search(r'dE_demanded_eV=([-\d.e+]+)',l).group(1))) for l in lsall]
    print('==',f.split('/')[-1])
    print('  Q1 max|dL_transport| na polkrok = %.3e hbar (%d polkrokow)'%(max(dl),len(dl)))
    print('  Q5 max|dE_demanded|            = %.3e eV'%max(de))
    if post:
        dev=max(abs(x['J']/x['n']-1) for x in post)
        print('  Q2 po zamknieciu: %d chk, max|J/n - 1| = %.3e, koncowe J/n = %.6f'%(len(post),dev,post[-1]['J']/post[-1]['n']))
        last=post[-100:]
        ratio=[(-x['dtot']-0) for x in last]
        # dn/n na checkpoint z kolejnych n
        rs=[]
        for a,b in zip(last,last[1:]):
            nm=(a['n']+b['n'])/2; rs.append(((a['n']-b['n'])/nm)/(eps1/2*nm**-3))
        print('  Q3 dn/n / ((e1/2) n^-3) w ostatnich %d chk: srednio %.4f'%(len(rs),sum(rs)/len(rs)))
        print('     dL_kredyt/chk = %.4e   dL_transport/chk = %.4e'%(sum(x['dcr'] for x in post)/len(post),sum(x['dtr'] for x in post)/len(post)))
