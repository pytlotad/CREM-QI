# Audyt 800 POST HOC: srednia kanalow wodor / antywodor 2P na ziarno (z out/ana800.txt).
import re, statistics as st
t=open('out/ana800.txt').read()
F1,F0=t.split('=== F = 0 (PHEN=1)')[0],t.split('=== F = 0 (PHEN=1)')[1].split('H800p')[0]
for blk,name in [(F1,'F = 1'),(F0,'F = 0')]:
    d={(m.group(1),m.group(2)):float(m.group(3)) for m in re.finditer(r'(proton|antiproton) ziarno (\d+).*?2P ([\d.]+) MHz',blk)}
    ms=[]
    for s in sorted({k[1] for k in d},key=int):
        a=(d[('proton',s)]+d[('antiproton',s)])/2; ms.append(a)
        print('%s ziarno %s: (H + antyH)/2 = %.2f MHz = %.5f pomiaru; H - antyH = %+.2f MHz'%(name,s,a,a/10969.045,d[('proton',s)]-d[('antiproton',s)]))
    print('%s: rozstep srednich par %.3f MHz, odchylenie std %.3f MHz'%(name,max(ms)-min(ms),st.pstdev(ms)))
