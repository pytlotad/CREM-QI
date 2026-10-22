#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <complex>
// Audit 559: the engine's S' image (psSingletImageInField) for Ps 1S in a field -- w against Breit-Rabi P_S(S'), and the survival over 0.5 ns from it in the transport (rules 541) against QM.
int main(){ std::setvbuf(stdout,nullptr,_IOLBF,0);
  applyPairFromOption("electron,positron");
  const double a=pairBohrRadius(activePair), L=0.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass), h=2*pi*hbar;
  const double n=quiggaRosnerContactDensity(a,L), Gp=4*pi*classicalElectronRadius*classicalElectronRadius*c*n, eps=4*(pi*pi-9)*fineStructureConstant/(9*pi), Go=Gp*eps, Gd=Gp-Go;
  const double dE=204.6577e9*h, Kc=dE/4, m1=firstMagneticMoment, m2=secondMagneticMoment; const Vec3 z{0,0,1}, x{1,0,0};
  for(double B:{0.1,1.0,3.0}){
    gExternalMagneticField=z*B; const double Z=m1*B, X=Z/Kc, PSbr=0.5*(1+1/std::sqrt(1+X*X));
    const PsFieldEigenImage im=psSingletImageInField(dE,m1*B,-m2*B,gExternalMagneticField,0.3);
    SecularSpinOrbitState s; s.orbitalAngularMomentum=z*L; s.periapsisDirection=x; gSingletChannelEligible=true; gSingletSpinTransport=false;
    s.firstDipole=im.n1*(-m1); s.secondDipole=im.n2*m2;
    const double w0=contactTwoPhotonWeight(s.firstDipole,s.secondDipole); double lnS=0;
    for(int i=1;i<=2000;++i){ s=advanceCoupledSecularSpinOrbit(s,a,mu,2.5e-13).state; lnS-=(Go+contactTwoPhotonWeight(s.firstDipole,s.secondDipole)*Gd)*2.5e-13; }
    using C=std::complex<double>; const double th=0.5*std::atan2(4*Z,4*Kc), vS=std::cos(th), vT=-std::sin(th);
    C H[2][2]={{C(-3*Kc,-0.5*hbar*Gp),C(2*Z,0)},{C(2*Z,0),C(Kc,-0.5*hbar*Go)}}; const double t=0.5e-9;
    const C tr=H[0][0]+H[1][1], det=H[0][0]*H[1][1]-H[0][1]*H[1][0], disc=std::sqrt(tr*tr-4.0*det), l1=(tr+disc)/2.0, l2=(tr-disc)/2.0;
    const C e1=std::exp(C(0,-1)*l1*t/hbar), e2=std::exp(C(0,-1)*l2*t/hbar); C U[2][2];
    for(int i=0;i<2;++i) for(int j=0;j<2;++j){ const C I=(i==j)?1.0:0.0; U[i][j]=(e1*(H[i][j]-l2*I)-e2*(H[i][j]-l1*I))/(l1-l2); }
    const double qm=std::norm(U[0][0]*vS+U[0][1]*vT)+std::norm(U[1][0]*vS+U[1][1]*vT);
    std::printf("B %.1f T: image m %d, P_S %.7f, w %.7f, Breit-Rabi P_S(S') %.7f; S(0.5 ns) engine %.6e, QM %.6e (ratio %.5f)\n",B,im.m,im.singletWeight,w0,PSbr,std::exp(lnS),qm,std::exp(lnS)/qm);
  }
}
