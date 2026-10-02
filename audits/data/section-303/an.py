import re,sys,statistics as st
def load(f):
    chk=[]; cur=None
    for ln in open(f):
        if ' ACTION ' in ln:
            d=dict(re.findall(r'(\w+)=([-\d.e+]+)',ln))
            cur={'k':int(d['checkpoint']),'J':float(d['Jphi_over_h']),'n':float(d['n']),'t':float(d['t']),'ls':[]}
            chk.append(cur)
        elif 'CREM_LS ' in ln and cur is not None:
            d=dict(re.findall(r'(\w+)=([-\d.e+]+)',ln))
            cur['ls'].append(d)
    return chk
def rows(f):
    c=load(f); out=[]
    for a,b in zip(c,c[1:]):
        if len(a['ls'])!=2: continue
        dtr=sum(float(x['L_after'])-float(x['L_before']) for x in a['ls'])
        dtot=b['J']-a['J']
        dem=sum(float(x['dE_demanded_eV']) for x in a['ls'])
        U=float(a['ls'][0]['U_dipole_eV'])
        out.append(dict(k=a['k'],n=a['n'],J=a['J'],dtot=dtot,dtr=dtr,dcr=dtot-dtr,dem=dem,U=U,dt=b['t']-a['t'],
                        Lc=float(a['ls'][0]['L_circ'])))
    return out
if __name__=='__main__':
    for f in sys.argv[1:]:
        r=rows(f); post=[x for x in r if x['n']<0.5 and abs(x['dtot'])<0.01]
        pre=[x for x in r if x['n']>0.5 and abs(x['dtot'])<0.01]
        print('==',f.split('/')[-1],' checkpointow po zamknieciu:',len(post))
        for name,g in (('przed zamknieciem (n>1/2, bez fotonu)',pre),('po zamknieciu',post)):
            if not g: continue
            m=lambda key: st.mean(x[key] for x in g)
            print('  %s: n=%d'%(name,len(g)))
            print('    dL_transport/chk = %+.4e hbar   dL_kredyt/chk = %+.4e   dL_calk = %+.4e'%(m('dtr'),m('dcr'),m('dtot')))
            print('    dE_demanded/chk  = %+.4e eV     U_dipole = %+.4e eV   iloraz %.3g'%(m('dem'),m('U'),m('dem')/m('U') if m('U') else float('nan')))
            print('    dt/chk = %.4e s   L/L_circ (J/n) = %.5f -> %.5f'%(m('dt'),g[0]['J']/g[0]['n'],g[-1]['J']/g[-1]['n']))
