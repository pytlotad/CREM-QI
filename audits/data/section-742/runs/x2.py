# audit 742, dead-path test: run from a copy of the tree (not the repository); multiplies the out-of-window midpoint rates by 2
p='modules/crem_trajectory.hpp'; s=open(p).read()
s=s.replace("rates.push_back(2.0*(2.0/3.0)*pairCoulombStrength","rates.push_back(2.0*2.0*(2.0/3.0)*pairCoulombStrength",1)
s=s.replace("rates.push_back(weight*(2.0/3.0)*pairCoulombStrength","rates.push_back(2.0*weight*(2.0/3.0)*pairCoulombStrength",1)
s=s.replace("rates.push_back(4.0*power/(hbar*omega));","rates.push_back(2.0*4.0*power/(hbar*omega));",1)
open(p,'w').write(s)
# then: p585 PAIR=proton,electron NLEV=3 LFRAC=0.16666666666666666 CREM_EMISSION_LEDGER=1, cmp with section-741/out/probes/new/H3s.*
