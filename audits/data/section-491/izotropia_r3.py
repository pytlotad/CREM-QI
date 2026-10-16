# Izotropia kierunku r z waga r^-3 na orbicie Keplera (<cos^2 f> = 1/2 przy kazdym e).
import numpy as np
# r^-3 weighted in-plane direction tensor on Kepler orbit
for e in [0,0.5,0.66,0.9]:
    f=np.linspace(0,2*np.pi,200001)[:-1]; w=(1+e*np.cos(f))  # (1/r^3)*r^2 df ∝ (1+e cos f)
    print(e, (np.cos(f)**2*w).sum()/w.sum())
