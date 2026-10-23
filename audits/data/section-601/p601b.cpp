// Audit 601: shares of the 1S field eigenstates after the photon from p-Ps 2^1P, engine functions (ps_field_eigen.hpp), isotropic L^.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
int main(){
  const double fields[2]={0.03,0.1};
  const double semi=4.0*pairBohrRadius(activePair);       // n = 2
  const double m1=std::abs(firstMagneticMoment), m2=std::abs(secondMagneticMoment);
  double E[3]; positroniumTripletLevels(semi,1,m1,m2,E);
  const double h=2*pi*hbar;
  // model 1S splitting: the engine's level difference (spinCouplingEnergy at n = 1, L = hbar/2), from the HFS 203.381 GHz of audit 575
  const double dE=203.3811e9*h;
  std::printf("E_J [GHz] %.6f %.6f %.6f  dE1S [GHz] %.4f  mu1 %.6e mu2 %.6e\n",E[0]/h/1e9,E[1]/h/1e9,E[2]/h/1e9,dE/h/1e9,m1,m2);
  std::mt19937_64 rng(6011); std::uniform_real_distribution<double> U(0,1);
  const double sgn1=firstCharge<0?-1.0:1.0, sgn2=secondCharge<0?-1.0:1.0;
  for(double B: fields){
    const Vec3 field{0,0,B};
    long cnt[4]={0,0,0,0}; const int N=2000000;
    for(int i=0;i<N;++i){
      const double cz=2*U(rng)-1, ph=2*pi*U(rng), sz=std::sqrt(1-cz*cz);
      const Vec3 Lh{sz*std::cos(ph),sz*std::sin(ph),cz};
      psFieldEigenAtFormation(E,-sgn1*m1*B,-sgn2*m2*B,field,Lh,Vec3{0,0,1},true,U(rng),rng());
      const PsFieldEigenImage im=psFieldEigenAtPhoton(dE,U(rng),2*pi*U(rng));
      const int c=im.m==1?0:im.m==-1?1:(im.singletWeight<0.5?2:3); ++cnt[c];
    }
    std::printf("B %.2f T  T+1 %.5f  T-1 %.5f  T0' %.5f  S' %.5f  (N %d)\n",B,cnt[0]/double(N),cnt[1]/double(N),cnt[2]/double(N),cnt[3]/double(N),N);
  }
}
