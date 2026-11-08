#!/usr/bin/env python3
"""Audit 741: F' of the 1s state after the last 2p -> 1s photon grouped by the HELD F of the 2p state (CREM_DEBUG_BORN line
"BORN photon F held") and the line of each group; compares with the 740 run (same cascades).  usage: ana741.py new741.err [new740.err]"""
import sys, math, importlib.util, os
here=os.path.dirname(os.path.abspath(__file__))
argv=sys.argv; sys.argv=[argv[0]]
spec=importlib.util.spec_from_file_location('a',os.path.join(here,'../section-739/ana739d.py')); a=importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
sys.argv=argv
e=1.602176634e-19; h=6.62607015e-34; hc=1239.8419843320026
QM={(0.5,0):0.0,(0.5,1):1/3,(1.5,1):2/3,(1.5,2):0.0}
def parse(path):
    out={}; idx=None; held=None
    for line in open(path):
        if line.startswith('START idx'): idx=int(line.split()[2]); held=None
        elif line.startswith('BORN photon F held'):
            t=line.split(); held=(int(t[4].split('=')[1]),int(t[5].split('=')[1]),int(t[6].split('=')[1]))
        elif line.startswith('CREM_LEDGER t') and idx is not None:
            parts=line.split('|'); t=parts[0].split(); B=a.side(parts[1][3:]); A=a.side(parts[2][3:])
            if 1.0<=B['L']/a.hbar<2.0 and A['L']/a.hbar<1.0:
                out[idx]=(float(t[t.index('Eg')+1])/e,B,A,held,line)
    return out
r=parse(argv[1]); r0=parse(argv[2]) if len(argv)>2 else {}
per={}; lines={}; nflip=0; same=0
for i,(E,B,A,held,ln) in sorted(r.items()):
    j,Fnow,c=a.basins(B,A); Fp=0 if c<-0.5 else 1
    Fh=held[2] if held else Fnow
    if held and held[1]!=held[2]: nflip+=1; print(f'  run {i}: 2j={held[0]} F_now={held[1]} F_held={held[2]} F\'={Fp} line {E:.10f} eV')
    elif i in r0 and r0[i][4]==ln: same+=1
    per.setdefault((j,Fh),[0,0]); per[(j,Fh)][0]+=1; per[(j,Fh)][1]+=(Fp==0)
    lines.setdefault((j,Fh,Fp),[]).append(E)
print(f'{len(r)} runs; F_held != F_now: {nflip}; other runs with the photon ledger line identical to the 740 run: {same}/{len(r)-nflip}')
for k in sorted(per):
    m,k0=per[k]; q=QM[k]; s=math.sqrt(q*(1-q)/m) if 0<q<1 else float('nan')
    print(f'  2P{"1/2" if k[0]==0.5 else "3/2"} F_held={k[1]}: {m:3d} runs, F\'=0 {k0:3d} ({k0/m:.3f}; 6j {q:.3f}; {(k0/m-q)/s if s==s else 0:+.2f} sigma)')
print(f'{"group":26s} {"n":>3s} {"mean [eV]":>16s} {"spread [eV]":>12s} {"nm":>14s} {"MHz":>16s}')
for k in sorted(lines):
    v=lines[k]; m=sum(v)/len(v)
    print(f'2P{"1/2" if k[0]==0.5 else "3/2"} F={k[1]} -> 1S F\'={k[2]}    {len(v):3d} {m:16.10f} {max(v)-min(v):12.3e} {hc/m:14.7f} {m*e/h/1e6:16.1f}')
