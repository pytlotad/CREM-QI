#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(int argc,char**argv){ applyPairFromOption(argv[1]);
  const double aB=pairBohrRadius(activePair), mu=firstMass*secondMass/(firstMass+secondMass), a=4*aB, e=std::sqrt(1-0.75*0.75);
  const double n=std::sqrt(pairCoulombStrength/(mu*a*a*a)); const double P=larmorOrbitAveragedPower(a,e), S=eccentricOrbitHazardSuppression(e);
  const double Eph=pairCoulombStrength/(2*aB)*(1-0.25);
  std::printf("%s: P %.4e W, S(e) %.5f, hbar omega %.4e J, E_photon %.4e J; tau = hbar omega/(P S) %.5e s; E_photon/P %.5e s; orbits %.5e\n",argv[1],P,S,hbar*n,Eph,hbar*n/(P*S),Eph/P,hbar*n/(P*S)*n/(2*pi)); }
