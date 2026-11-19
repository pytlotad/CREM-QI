# Domyslna miara precesji na obieg radialny (od audytu 804): srednia kata delta miedzy kolejnymi apocentrami tylko z odcinkow |dt/mediana - 1| <= 0,1,
# Omega_obr = srednia / (2 pi T_n), T_n = period z linii sum.  Detektor apsyd w silniku ma od 804 histereze (0,1 period); filtr jest zabezpieczeniem.
# Zwraca (Omega_obr, liczba odrzuconych odcinkow, liczba odcinkow, liczba apocentrow).  Uzycie: from per_orbit import om
import math
import numpy as np
def _v3(w, tag): i = w.index(tag); return np.array([float(x) for x in w[i + 1:i + 4]])
def om(path, tag='mechB0'):
    T = []; L = []; S = []; per = None
    for ln in open(path):
        if not ln.startswith(f'PAIR TORQUE799 {tag} '): continue
        w = ln.split()
        if w[3] == 'apo': T.append(float(w[w.index('t') + 1])); L.append(_v3(w, 'L')); S.append(_v3(w, 's2'))
        elif w[3] == 'sum': per = float(w[w.index('period') + 1])
    d = np.array([np.linalg.norm(S[j] - S[j - 1])/np.linalg.norm(np.cross(L[j - 1], S[j - 1])) for j in range(1, len(S))])
    dt = np.diff(T); ok = np.abs(dt/np.median(dt) - 1) <= 0.1
    return float(np.mean(d[ok]))/(2*math.pi*per), int((~ok).sum()), len(d), len(T)
