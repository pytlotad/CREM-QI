import re,glob,os,math,sys
D=sys.argv[1]
print('%5s %4s %7s %5s %6s %11s | %5s %8s %8s %8s %5s | %s'%('J0','kan','outcome','fot','odm','t_sym ps','chk1','n przed','n po','granica','H2','koniec n/J/e  chk'))
for f in sorted(glob.glob(D+'*.out')):
    b=os.path.basename(f)[:-4]; J,ph=b.split('-'); J=float(J)
    out=open(f).read()
    m=re.search(r'fotony=(\d+) odmowy ceiling=(\d+)',out); o=re.search(r'outcome=(\d) .*calib=([\d.e+-]+)',out)
    rows=[dict(re.findall(r'(\w+)=([-\d.e+]+)',l)) for l in open(D+b+'.err') if l.startswith('ACTION')]
    f1=None
    for a,c in zip(rows,rows[1:]):
        if abs(float(c['Jphi_over_h'])-float(a['Jphi_over_h']))>0.02: f1=(a,c); break
    bound=max(1-J,1/math.sqrt(3)); last=rows[-1] if rows else None
    h2='-'
    if f1: h2='TAK' if float(f1[1]['n'])>=bound-1e-3 else 'NIE'
    print('%5.2f %4s %7s %5s %6s %11s | %5s %8s %8s %8.4f %5s | %s %d'%(J,'para' if ph=='1' else 'orto',
      o.group(1) if o else 'zabity', m.group(1) if m else '?', m.group(2) if m else '?',
      '%.4f'%(float(o.group(2))*1e12) if o else '-',
      f1[1]['checkpoint'] if f1 else '-', '%.4f'%float(f1[0]['n']) if f1 else '-','%.4f'%float(f1[1]['n']) if f1 else '-', bound, h2,
      '%.4f/%.4f/%.4f'%(float(last['n']),float(last['Jphi_over_h']),float(last['e'])) if last else '-', len(rows)))
