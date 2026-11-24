// Audit 851 (static probe): N photon increments applyStochasticDipolePhoton on a frozen Kepler state at n = 4 (no dynamics);
// per-increment mismatch P = d(kin - kin_COM) + share, the DRIFT850 measure, in J and in n (n = 3 conversion n^3/(2B), as ana850).
// CREM_PHOTON_INCREMENT_STABLE in the environment selects the 851 path.  argv: pair N
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char**argv){
  if(argc<3){ std::fprintf(stderr,"usage\n"); return 1; }
  applyPairFromOption(argv[1]); const int N=std::atoi(argv[2]);
  const double M=firstMass+secondMass, mu=firstMass*secondMass/M, k=pairCoulombStrength;
  const double B=mu*k*k/(2*hbar*hbar), a=16*hbar*hbar/(mu*k), ecc=0.5, Ean=0.3;
  const double nmot=std::sqrt(k/(mu*a*a*a)), ec=std::sqrt(1-ecc*ecc), w=1-ecc*std::cos(Ean);
  const Vec3 r{a*(std::cos(Ean)-ecc),a*ec*std::sin(Ean),0}, v=Vec3{-a*nmot*std::sin(Ean),a*nmot*ec*std::cos(Ean),0}*(1/w);
  State s{}; s.firstPosition=r*(secondMass/M); s.secondPosition=r*(-firstMass/M); s.firstVelocity=v*(secondMass/M); s.secondVelocity=v*(-firstMass/M);
  const auto kr=[&](const State& x){ const Vec3 P=x.firstVelocity*firstMass+x.secondVelocity*secondMass;
    return (long double)kineticEnergy(x.firstVelocity,firstMass)+(long double)kineticEnergy(x.secondVelocity,secondMass)-(long double)(P.squaredNorm()/(2*M)); };
  const double share=(B/9-B/16)/N, conv=27.0/(2*B);
  long double sum=0, sum2=0; int neg=0, emitted=0;
  for(int j=0;j<N;++j){
    const double th=std::acos(1-2*(j+0.5)/N), ph=2.399963229728653*j;
    const Vec3 dir{std::sin(th)*std::cos(ph),std::sin(th)*std::sin(ph),std::cos(th)};
    const long double before=kr(s);
    const StochasticPhotonRecoil inc=applyStochasticDipolePhoton(s,share,dir);
    if(!inc.emitted) continue; ++emitted;
    const long double P=kr(s)-before+share; sum+=P; sum2+=P*P; if(P<0) ++neg;
  }
  const double mean=(double)(sum/emitted), rms=(double)std::sqrt(sum2/emitted);
  std::printf("S851 pair %s stable %d N %d emitted %d share %.6e J meanP %.6e J = %.4e n  rmsP %.4e n  sumP %.4e n  neg %d\n",
    argv[1],std::getenv("CREM_PHOTON_INCREMENT_STABLE")?1:0,N,emitted,share,mean,mean*conv,rms*conv,(double)sum*conv,neg);
}
