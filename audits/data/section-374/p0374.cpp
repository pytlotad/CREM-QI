#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 374: 2^3P_J spin-orbit + tensor energies under three classical vector-model assignments (orbit: n = 2, L = 3/2 hbar).
int main(){
  const double hP=2*pi*hbar, mu=firstMagneticMoment, mred=firstMass*secondMass/(firstMass+secondMass), GHz=1e-9/hP;
  const double aPs=pairBohrRadius(activePair); const int n=2; const double a=n*n*aPs, Lorb=1.5, eP=std::sqrt(1-std::pow(Lorb/n,2));
  const Vec3 Lh{0,0,1}, P{1,0,0}, Q{0,1,0};
  const auto orb=orbitAveragedBmtAngularVelocities(a,Lh*(Lorb*hbar),Vec3{1,0,0}*(mu*1e-12),Vec3{0,1,0}*(mu*1e-12),mred,0.0,P);
  const double wL=dot(orb.first+orb.second,Lh);
  auto U=[&](const Vec3& s){ return orbitAveragedDipoleEnergy(a,eP,Lh,P,s*mu,s*(-mu)); };
  const double Uiso=(U({1,0,0})+U({0,1,0})+U({0,0,1}))/3.0;
  auto tensor=[&](double ct){ const double st=std::sqrt(std::max(0.0,1-ct*ct)); double t=0; const int M=256;
    for(int j=0;j<M;++j){ const double ph=2*pi*(j+0.5)/M; t+=U(P*(st*std::cos(ph))+Q*(st*std::sin(ph))+Lh*ct)-Uiso; } return t/M; };
  struct V{const char* name; double L,S; bool qJ;} vs[]={{"V1 |L|=3/2 |S|=1 |J|=J+1/2",1.5,1.0,false},{"V2 |L|=|S|=sqrt2 |J|=sqrt(J(J+1))",std::sqrt(2.0),std::sqrt(2.0),true},{"V3 |L|=3/2 |S|=sqrt2 |J|=J+1/2",1.5,std::sqrt(2.0),false}};
  const double m10=5.48723, m21=4.38804;
  for(auto& v: vs){ double E[3],T[3];
    for(int J=0;J<=2;++J){ const double Jm=v.qJ?std::sqrt(J*(J+1.0)):J+0.5; double ct=(Jm*Jm-v.L*v.L-v.S*v.S)/(2*v.L*v.S); ct=std::clamp(ct,-1.0,1.0);
      T[J]=tensor(ct); E[J]=0.5*hbar*wL*v.S*ct+T[J]; }
    std::printf("%s:\n  tensor J0 %.6f J1 %.6f J2 %.6f GHz (J0-J2 %.2e)\n  P1-P0 %.3f (meas %.3f, %+.0f%%)  P2-P1 %.3f (meas %.3f, %+.0f%%)\n",
      v.name,T[0]*GHz,T[1]*GHz,T[2]*GHz,(T[0]-T[2])*GHz,(E[1]-E[0])*GHz,m10,100*((E[1]-E[0])*GHz/m10-1),(E[2]-E[1])*GHz,m21,100*((E[2]-E[1])*GHz/m21-1)); }
}
