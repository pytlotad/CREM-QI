# Audit 792 POST HOC (after the scored ana792 output): does the orbit energy explain the channel gap on n = 3, L = 1, seed 5?
# Omega ~ <r^-3> L ~ |E|^{3/2} / L^2 at fixed L; Kepler energy eps at the end of each one-period piece (T traces), mean over 200 pieces.
import io, contextlib, numpy as np
with contextlib.redirect_stdout(io.StringIO()): exec(open(__file__.replace('post792.py', 'ana792.py')).read())
mu = 9.1093837015e-31*1.67262192369e-27/(9.1093837015e-31 + 1.67262192369e-27)
EN = -13.605693122994/9*1.602176634e-19/mu
pr(f'eps nominalne (n = 3, Bohr, masa zredukowana): {EN:.6e} J/kg')
res = {}
for t in ('proton', 'antiproton'):
    rows = chunks(f'{O}/T/{t}_n3L1_5.err'); ep = np.array([r['eps'] for r in rows]); om = omega(f'{O}/T/{t}_n3L1_5.err')[0]
    f = (ep.mean()/EN)**1.5; res[t] = (om, om/f)
    pr(f'  {t}: eps/nom - 1 obieg 1 {ep[0]/EN - 1:+.4f}, obieg 200 {ep[-1]/EN - 1:+.4f}, srednio {ep.mean()/EN - 1:+.4f} | Omega {om/1e6:.2f} MHz -> przeliczone na E_nom {om/f/1e6:.2f} MHz'
       f' | Omega n^3 L^2 {om/f*27/1e6:.3f} MHz')
pr(f"  luka proton / antyproton: surowa {res['proton'][0]/res['antiproton'][0] - 1:+.4f}, po przeliczeniu na E_nom {res['proton'][1]/res['antiproton'][1] - 1:+.4f}")
pr('  (eps w chwili konca kawalka, tuz przy peryapsie: zawiera chwilowe czlony nie-keplerowskie; przeliczenie tylko orientacyjne)')
