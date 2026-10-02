import re,glob,os,math,json
B='/tmp/claude-1000/-home-teddy-Projekty-Positronium/5981d9ee-dfbe-485a-8f93-5c6d87fd7653/scratchpad/grp/'
out=[]
for src in ('s298','s300','s301','s302'):
    for f in sorted(glob.glob(B+src+'/*.err')):
        rows=[dict(re.findall(r'(\w+)=([-\d.e+]+)',l)) for l in open(f) if l.lstrip('0123456789. ').startswith('ACTION')]
        k=0
        for a,c in zip(rows,rows[1:]):
            na,nc=float(a['n']),float(c['n']); Ja,Jc=float(a['Jphi_over_h']),float(c['Jphi_over_h'])
            if nc<na-1e-3 and abs(Jc-Ja)>1e-3:
                k+=1
                out.append(dict(src=src,run=os.path.basename(f)[:-4],k=k,nb=na,Jb=Ja,eb=float(a['e']),na=nc,Ja=Jc,ea=float(c['e']),t=float(c['t'])))
json.dump(out,open(B+'photons.json','w'))
print('fotonow:',len(out),' z',len(set((p['src'],p['run']) for p in out)),'przebiegow')
