// Audit 739, H739a: spin energy of hydrogen 1s on a grid of cos(s_p, s_e); the two basin values, A = F1 - F0 in MHz, the isotropic F = 0 share.
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  applyPairFromOption("proton,electron");
  gSingletChannelEligible=true;
  const double a=pairBohrRadius(activePair), L=0.5*hbar, MHz=1.0e-6/(2.0*pi*hbar), ueV=1.0e6/elementaryCharge;
  const double ecc=std::sqrt(std::max(0.0,1.0-0.25));
  const Vec3 z{0,0,1}, x{1,0,0};
  const Vec3 mp=z*firstMagneticMoment;           // proton: s_p || m_p
  double E0=0,E1=0,maxDev0=0,maxDev1=0; int n0=0,n1=0;
  const int N=2001; double ref0=0,ref1=0;
  for(int i=0;i<N;++i){
    const double c=-1.0+2.0*i/(N-1), s=std::sqrt(std::max(0.0,1.0-c*c));
    const Vec3 se{s,0,c};                          // electron spin direction
    const Vec3 me=se*(-secondMagneticMoment);      // electron moment antiparallel to spin
    const double U=spinCouplingEnergy(a,ecc,z,x,mp,me,L);
    if(c<-0.5){ if(!n0) ref0=U; maxDev0=std::max(maxDev0,std::abs(U/ref0-1)); E0+=U; ++n0; }
    else      { if(!n1) ref1=U; maxDev1=std::max(maxDev1,std::abs(U/ref1-1)); E1+=U; ++n1; }
    if(i%250==0||i==N-1) std::printf("cos %+.4f  U %.9f ueV\n",c,U*ueV);
  }
  E0/=n0; E1/=n1;
  std::printf("F=0 basin: %d points, U %.9f ueV (max rel dev %.2e); F=1 basin: %d points, U %.9f ueV (max rel dev %.2e)\n",n0,E0*ueV,maxDev0,n1,E1*ueV,maxDev1);
  std::printf("U1/U0 = %.12f (Born: -1/3)\n",E1/E0);
  std::printf("A = U1 - U0 = %.6f MHz = %.6f ueV (measured 1420.405751768 MHz, Hellwig 1970): %+.6f MHz\n",(E1-E0)*MHz,(E1-E0)*ueV,(E1-E0)*MHz-1420.405751768);
  const int M=20000; int f0=0;
  for(int k=0;k<M;++k){ const double cz=1.0-2.0*(k+0.5)/M; if(cz<-0.5) ++f0; }   // cos(s_p, s_e) of a Fibonacci sphere about s_p = z
  std::printf("isotropic F=0 share %.5f (QM 1/4)\n",double(f0)/M);
}
