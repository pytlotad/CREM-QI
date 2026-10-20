#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <complex>
// Audit 541: Ps 1S triplet "T0" start (parallel spins _|_ B) in the engine's transport with the field rules; survival from the model's hazard
// Gamma_o + w (Gamma_p - Gamma_o) vs QM {S, T0} with decay (same Gamma_p, Gamma_o and the model's HFS; also with the measured Delta E).
int main(){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true; gSingletSpinTransport=false;
  const double aB=pairBohrRadius(activePair), a=aB, L=0.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass), h=2*pi*hbar;
  const double n=quiggaRosnerContactDensity(a,L), Gp=4*pi*classicalElectronRadius*classicalElectronRadius*c*n, eps=4*(pi*pi-9)*fineStructureConstant/(9*pi), Go=Gp*eps, Gd=Gp-Go;
  const double muB=9.2740100783e-24, g=std::abs(firstGFactor);
  std::printf("model: Gamma_p %.5e/s (tau %.2f ps), Gamma_o %.5e/s (tau %.2f ns), field rules %s\n",Gp,1e12/Gp,Go,1e9/Go,std::getenv("CREM_NO_PS_FIELD_RULES")?"OFF":"on");
  const Vec3 z{0,0,1}, x{1,0,0};
  for(double B:{0.1,1.0,3.0}){
    gExternalMagneticField=z*B; const double Z=(g/2)*muB*B;
    auto qm=[&](double dE,double t){ using C=std::complex<double>; const double J=dE/4;
      C H[2][2]={{C(-3*J,-0.5*hbar*Gp),C(2*Z,0)},{C(2*Z,0),C(J,-0.5*hbar*Go)}};   // exact 2x2 exponential via eigen
      const C tr=H[0][0]+H[1][1], det=H[0][0]*H[1][1]-H[0][1]*H[1][0], disc=std::sqrt(tr*tr-4.0*det), l1=(tr+disc)/2.0, l2=(tr-disc)/2.0;
      const C e1=std::exp(C(0,-1)*l1*t/hbar), e2=std::exp(C(0,-1)*l2*t/hbar);
      // U = (e1 (H - l2) - e2 (H - l1))/(l1 - l2); psi0 = (0,1)
      const C p0=(e1*H[0][1]-e2*H[0][1])/(l1-l2), p1=(e1*(H[1][1]-l2)-e2*(H[1][1]-l1))/(l1-l2);
      return std::norm(p0)+std::norm(p1); };
    SecularSpinOrbitState s; s.orbitalAngularMomentum=z*L; s.periapsisDirection=x; s.firstDipole=x*(-firstMagneticMoment); s.secondDipole=x*secondMagneticMoment;
    double lnS=0, t=0; const double dt=2.5e-13; double S05=0,S2=0, wtail=0; long ntail=0;
    for(int i=1;i<=8000;++i){ s=advanceCoupledSecularSpinOrbit(s,a,mu,dt).state; t+=dt; const double w=contactTwoPhotonWeight(s.firstDipole,s.secondDipole); lnS-=(Go+w*Gd)*dt;
      if(i==2000) S05=std::exp(lnS); if(i==8000) S2=std::exp(lnS); if(i>7200){ wtail+=w; ++ntail; } }
    const double late=Go+wtail/ntail*Gd;
    const double xm=4*Z/(204.6577e9*h), PSm=(1-1/std::sqrt(1+xm*xm))/2;
    std::printf("B %.1f T: engine S(0.5, 2 ns) %.6f %.6f, late rate %.5e; QM (model HFS) %.6f %.6f, eigen late %.5e (ratio %.4f); QM (measured dE) %.6f %.6f\n",
      B,S05,S2,late,qm(204.6577e9*h,0.5e-9),qm(204.6577e9*h,2e-9),Go+PSm*Gd,late/(Go+PSm*Gd),qm(203.38910e9*h,0.5e-9),qm(203.38910e9*h,2e-9));
  }
}
