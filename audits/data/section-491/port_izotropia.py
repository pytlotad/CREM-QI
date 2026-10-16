# literal pure-Python port of izotropia_r3.py
import math
for e in [0,0.5,0.66,0.9]:
    N=200000; num=0; den=0
    for i in range(N):
        f=2*math.pi*i/N; w=1+e*math.cos(f); num+=math.cos(f)**2*w; den+=w
    print(e, num/den)
