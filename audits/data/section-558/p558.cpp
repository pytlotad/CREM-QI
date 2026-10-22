#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <complex>
// Audit 558: Ps 1S singlet-like eigenstate S' in a field in the engine's transport (field rules 541, no singlet transport) vs QM {S, T0} with decay started from S'; and the current
// singlet transport (gSingletSpinTransport = true, antiparallel along B) for comparison.  Survival from the hazard Gamma_o + w (Gamma_p - Gamma_o).
int main(){ std::setvbuf(stdout,nullptr,_IOLBF,0);
  applyPairFromOption("electron,positron");
  const double aB=pairBohrRadius(activePair), a=aB, L=0.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass), h=2*pi*hbar;
  const double n=quiggaRosnerContactDensity(a,L), Gp=4*pi*classicalElectronRadius*classicalElectronRadius*c*n, eps=4*(pi*pi-9)*fineStructureConstant/(9*pi), Go=Gp*eps, Gd=Gp-Go;
  const double muB=9.2740100783e-24, g=std::abs(firstGFactor), dE=204.6577e9*h, Kc=dE/4;
  const Vec3 z{0,0,1}, x{1,0,0};
  for(double B:{0.1,1.0,3.0}){
    gExternalMagneticField=z*B; const double Z=(g/2)*muB*B, X=Z/Kc;
    // QM: H = [[-3K - i hbar Gp/2, 2Z],[2Z, K - i hbar Go/2]] (S, T0); start: S' of the Hermitian part
    using C=std::complex<double>;
    const double th=0.5*std::atan2(4*Z,4*Kc);   // Hermitian [[-3K,2Z],[2Z,K]]: lower eigenvector (cos th, -sin th)
    const double vS=std::cos(th), vT=-std::sin(th), PSqm=vS*vS;
    auto qmS=[&](double t){ C H[2][2]={{C(-3*Kc,-0.5*hbar*Gp),C(2*Z,0)},{C(2*Z,0),C(Kc,-0.5*hbar*Go)}};
      const C tr=H[0][0]+H[1][1], det=H[0][0]*H[1][1]-H[0][1]*H[1][0], disc=std::sqrt(tr*tr-4.0*det), l1=(tr+disc)/2.0, l2=(tr-disc)/2.0;
      const C e1=std::exp(C(0,-1)*l1*t/hbar), e2=std::exp(C(0,-1)*l2*t/hbar);
      C U[2][2]; for(int i=0;i<2;++i) for(int j=0;j<2;++j){ const C I=(i==j)?1.0:0.0; U[i][j]=(e1*(H[i][j]-l2*I)-e2*(H[i][j]-l1*I))/(l1-l2); }
      const C p0=U[0][0]*vS+U[0][1]*vT, p1=U[1][0]*vS+U[1][1]*vT; return std::norm(p0)+std::norm(p1); };
    // late QM rate: the slower-decaying root's imaginary part (2|Im l|/hbar) -- report the faster one (singlet-like)
    C H[2][2]={{C(-3*Kc,-0.5*hbar*Gp),C(2*Z,0)},{C(2*Z,0),C(Kc,-0.5*hbar*Go)}}; const C tr=H[0][0]+H[1][1], det=H[0][0]*H[1][1]-H[0][1]*H[1][0], disc=std::sqrt(tr*tr-4.0*det);
    const double rA=-2*std::imag((tr+disc)/2.0)/hbar, rB=-2*std::imag((tr-disc)/2.0)/hbar, rSing=std::max(rA,rB);
    for(int mode=0;mode<2;++mode){
      SecularSpinOrbitState s; s.orbitalAngularMomentum=z*L; s.periapsisDirection=x;
      gSingletChannelEligible=true; gSingletSpinTransport=(mode==1);
      Vec3 n1,n2;
      if(mode==0){ const double beta=0.5*(std::atan(X)-pi); n1=Vec3{std::cos(beta),0,std::sin(beta)}; n2=Vec3{std::cos(beta),0,-std::sin(beta)}; }
      else { n1=z*-1.0; n2=z; }
      s.firstDipole=n1*(-firstMagneticMoment); s.secondDipole=n2*secondMagneticMoment;
      const double w0=contactTwoPhotonWeight(s.firstDipole,s.secondDipole);
      double lnS=0, S05=0, S2=0, wmin=1, wmax=0, wtail=0; long nt=0; const double dt=2.5e-13; const int steps=8000;
      for(int i=1;i<=steps;++i){ s=advanceCoupledSecularSpinOrbit(s,a,mu,dt).state; const double w=contactTwoPhotonWeight(s.firstDipole,s.secondDipole); lnS-=(Go+w*Gd)*dt;
        wmin=std::min(wmin,w); wmax=std::max(wmax,w); if(i==2000) S05=std::exp(lnS); if(i==8000) S2=std::exp(lnS); if(i>7200){ wtail+=w; ++nt; } }
      const double late=Go+wtail/nt*Gd;
      std::printf("B %.1f T %s: w0 %.6f (P_S(S') %.6f), w range %.6f-%.6f; S(0.5, 2 ns) %.6e %.6e, late rate %.5e/s | QM from S': %.6e %.6e, singlet-like rate %.5e (ratio %.4f)\n",
        B,mode?"singlet transport (now)":"S' + rules 541",w0,PSqm,wmin,wmax,S05,S2,late,qmS(0.5e-9),qmS(2e-9),rSing,late/rSing);
    }
  }
}
