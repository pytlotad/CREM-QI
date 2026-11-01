// Audit 679: Ps 1S m = 0 images in a field (S' and T0') -- the engine's 559 image, the 676 rule and the hybrid (w = P_S with n_z = <sigma_z> where
// possible) -- singlet weight w, survival in the engine's transport (rules 541) against QM, and the Zeeman part of the image energy.
// build (repository root): g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-679/p679.cpp $(root-config --libs) -o p679
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <complex>
int main(){ std::setvbuf(stdout,nullptr,_IOLBF,0);
  applyPairFromOption("electron,positron");
  const double a=pairBohrRadius(activePair), L=0.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass), h=2*pi*hbar, u=1.602176634e-25;
  const double n=quiggaRosnerContactDensity(a,L), Gp=4*pi*classicalElectronRadius*classicalElectronRadius*c*n, eps=4*(pi*pi-9)*fineStructureConstant/(9*pi), Go=Gp*eps, Gd=Gp-Go;
  const double dE=204.6577e9*h, Kc=dE/4, m1=firstMagneticMoment, m2=secondMagneticMoment; const Vec3 z{0,0,1}, x{1,0,0}, y{0,1,0};
  using C=std::complex<double>;
  std::printf("B[T] state P_S s=<s1z> image | w | S(0.5ns) eng/QM | S(5ns) eng/QM | Zeeman part image-QM [ueV]\n");
  for(double B:{0.1,1.0,3.0}){
    gExternalMagneticField=z*B; const double Z=m1*B;
    const double th=0.5*std::atan2(4*Z,4*Kc), vS=std::cos(th), vT=-std::sin(th);
    C H[2][2]={{C(-3*Kc,-0.5*hbar*Gp),C(2*Z,0)},{C(2*Z,0),C(Kc,-0.5*hbar*Go)}};
    const C tr=H[0][0]+H[1][1], det=H[0][0]*H[1][1]-H[0][1]*H[1][0], disc=std::sqrt(tr*tr-4.0*det), l1=(tr+disc)/2.0, l2=(tr-disc)/2.0;
    const auto qm=[&](double cS,double cT,double t){ const C e1=std::exp(C(0,-1)*l1*t/hbar), e2=std::exp(C(0,-1)*l2*t/hbar); C U[2][2];
      for(int i=0;i<2;++i) for(int j=0;j<2;++j){ const C I=(i==j)?1.0:0.0; U[i][j]=(e1*(H[i][j]-l2*I)-e2*(H[i][j]-l1*I))/(l1-l2); }
      return std::norm(U[0][0]*cS+U[0][1]*cT)+std::norm(U[1][0]*cS+U[1][1]*cT); };
    for(int st=0;st<2;++st){
      const double cS=st==0?vS:-vT, cT=st==0?vT:vS;          // S' = vS|S> + vT|T0>, T0' orthogonal
      const double PS=cS*cS, s1=2*cS*cT, s2=-s1, sgn=s1>=0?1:-1, t=std::sqrt(std::max(0.0,1-s1*s1));
      for(int im=0;im<3;++im){
        Vec3 n1,n2; const char* name="";
        const double cosD=(1+2*PS-4*PS*PS)/((1-2*PS)*(1-2*PS));
        if(im==0||(im==2&&cosD>1)){ const double sp=std::sqrt(PS), cp=std::sqrt(1-PS); n1=x*cp+z*(sgn*sp); n2=x*cp-z*(sgn*sp); name=im==0?"I559":"Ihyb(=I559)"; }
        else if(im==1){ const double side=PS<0.5?1:-1; n1=x*t+z*s1; n2=x*(side*t)+z*s2; name="Isz"; }
        else { const double sD=std::sqrt(std::max(0.0,1-cosD*cosD)); n1=x*t+z*s1; n2=(x*cosD+y*sD)*t+z*s2; name="Ihyb"; }
        SecularSpinOrbitState s; s.orbitalAngularMomentum=z*L; s.periapsisDirection=x; gSingletChannelEligible=true; gSingletSpinTransport=false;
        s.firstDipole=n1*(-m1); s.secondDipole=n2*m2;
        const double w0=contactTwoPhotonWeight(s.firstDipole,s.secondDipole); double lnS=0, S05=0;
        for(int i=1;i<=20000;++i){ s=advanceCoupledSecularSpinOrbit(s,a,mu,2.5e-13).state; lnS-=(Go+contactTwoPhotonWeight(s.firstDipole,s.secondDipole)*Gd)*2.5e-13; if(i==2000) S05=std::exp(lnS); }
        const double zim=m1*B*(n1.z-n2.z), zqm=m1*B*(s1-s2);
        std::printf("%.1f %s P_S %.7f s %+.6f %-11s | w %.7f | %.5f | %.5f | %+9.4f\n",B,st==0?"S' ":"T0'",PS,s1,name,w0,S05/qm(cS,cT,0.5e-9),std::exp(lnS)/qm(cS,cT,5e-9),(zim-zqm)/u);
      }
    }
  }
}
