# Audit 298: test znakow z cenzura prawostronna.  Dla przebiegu otwartego
# czas zamkniecia jest >= ostatniego czasu symulowanego (dolne ograniczenie).
import sys, math
sys.path.insert(0,'/tmp/claude-1000/-home-teddy-Projekty-Positronium/5981d9ee-dfbe-485a-8f93-5c6d87fd7653/scratchpad/close')
from parse2 import load
R=load()
def obs(k):
    v=R.get(k)
    if v is None: return None
    if v['close']: return (v['close'][0],True)
    return (v['last'][0],False)          # cenzura: >= last t
def order(a,b):
    # +1 gdy a pozniej niz b, -1 gdy wczesniej, 0 nierozstrzygalne
    if a is None or b is None: return None
    (ta,ca),(tb,cb)=a,b
    if ca and cb: return (ta>tb)-(ta<tb)
    if ca and not cb: return -1 if ta<tb else 0
    if cb and not ca: return 1 if tb<ta else 0
    return 0
def binom_two_sided(k,n):
    if n==0: return 1.0
    p=sum(math.comb(n,i) for i in range(0,n+1) if abs(i-n/2)>=abs(k-n/2))/2**n
    return min(1.0,p)
def table(label,keyA,keyB):
    plus=minus=undec=0; rows=[]
    for s in range(40,50):
        a,b=obs(keyA(s)),obs(keyB(s)); o=order(a,b)
        fa='%8.1f%s'%(a[0]*1e12,'' if a[1] else '+') if a else '   -    '
        fb='%8.1f%s'%(b[0]*1e12,'' if b[1] else '+') if b else '   -    '
        rows.append('   %d  %s  %s  %s'%(s,fa,fb,{1:'A pozniej',-1:'A wczesniej',0:'nierozstrz.',None:'brak'}[o]))
        if o==1: plus+=1
        elif o==-1: minus+=1
        else: undec+=1
    n=plus+minus
    print(label); print('\n'.join(rows))
    print('   rozstrzygalnych %d: A pozniej %d, A wczesniej %d, nierozstrz. %d;  p (znaki, dwustronnie) = %.3f\n'
          %(n,plus,minus,undec,binom_two_sided(plus,n)))
print('czasy w ps; "+" = cenzura (zamkniecie POZNIEJ niz podany czas)\n')
for cfg in ('zach','ret'):
    table('%s: A = orto, B = para'%cfg, lambda s:(cfg,2,s), lambda s:(cfg,1,s))
for ph,name in ((1,'para'),(2,'orto')):
    table('%s: A = retardowana, B = zachowawcza'%name, lambda s:('ret',ph,s), lambda s:('zach',ph,s))
tot=len(R); cens=sum(1 for v in R.values() if not v['close'])
print('przebiegow %d, ocenzurowanych %d'%(tot,cens))

print('\nWZGLEDNE ROZNICE SPAROWANE (tylko pary, w ktorych oba przebiegi zamkniete)')
import statistics as st
for label,kA,kB in [('zach orto-para',lambda s:('zach',2,s),lambda s:('zach',1,s)),
                    ('ret  orto-para',lambda s:('ret',2,s),lambda s:('ret',1,s)),
                    ('para ret-zach ',lambda s:('ret',1,s),lambda s:('zach',1,s)),
                    ('orto ret-zach ',lambda s:('ret',2,s),lambda s:('zach',2,s))]:
    d=[]
    for s in range(40,50):
        a,b=obs(kA(s)),obs(kB(s))
        if a and b and a[1] and b[1]: d.append((s,(a[0]-b[0])/b[0]))
    if not d: print('  %s: brak par'%label); continue
    vals=[x for _,x in d]
    print('  %s: n=%d  %s'%(label,len(d),'  '.join('%d:%+.4f%%'%(s,100*x) for s,x in d)))
    if len(vals)>1:
        print('      srednia %+.4f%%  mediana %+.4f%%  max|.| %.4f%%'%(100*st.mean(vals),100*st.median(vals),100*max(abs(x) for x in vals)))
