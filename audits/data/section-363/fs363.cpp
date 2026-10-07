#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 363: n = 2 fine structure of the classical secular Hamiltonian, Langer L, classical triplet |S| = hbar, |J| = (J+1/2) hbar.
int main(){
  const double hP=2*pi*hbar, m=firstMass, mred=firstMass*secondMass/(firstMass+secondMass), k=pairCoulombStrength, mu=firstMagneticMoment;
  const double alpha=k/(hbar*c), a4=std::pow(alpha,4)*m*c*c, aPs=pairBohrRadius(activePair), psi0=1.0/(pi*aPs*aPs*aPs);
  const int n=2; const double a=n*n*aPs; const double GHz=1e-9/hP;
  auto Erel=[&](double Lh){ // numeric orbit average of the e+e- Darwin Hamiltonian
    const double eta=Lh/n, e=std::sqrt(std::max(0.0,1-eta*eta)), L=Lh*hbar; const int N=200000; double s=0,w=0;
    for(int i=0;i<N;++i){ const double E=2*pi*(i+0.5)/N, tw=1-e*std::cos(E), r=a*tw;
      const double p2=mred*k*(2/r-1/a), pr2=std::max(0.0,p2-L*L/(r*r));
      const double H=-p2*p2/(4*m*m*m*c*c)-k/(2*m*m*c*c*r)*(p2+pr2); s+=H*tw; w+=tw; }
    return s/w; };
  auto ErelAnalytic=[&](double Lh){ return a4*(11.0/(64.0*n*n*n*n)-1.0/(4.0*n*n*n*Lh)); };
  const double ES=Erel(0.5), EP=Erel(1.5);
  std::printf("alpha^4 m c^2 = %.4f GHz\n",a4*GHz);
  std::printf("E_rel 2S: numeric %.6f GHz analytic %.6f ; 2P: numeric %.6f analytic %.6f ; 2P-2S = %.6f GHz (alpha^4mc^2/24 = %.6f)\n",
    ES*GHz,ErelAnalytic(0.5)*GHz,EP*GHz,ErelAnalytic(1.5)*GHz,(EP-ES)*GHz,a4/24*GHz);
  // contact, 2^3S1: moments antiparallel (triplet), n_QR = |psi0|^2 / n^3 at L = hbar/2
  const double Econt=(2*mu0/3)*mu*mu*psi0/(n*n*n);
  std::printf("contact 2^3S1 = %.6f GHz\n",Econt*GHz);
  // spin-orbit and tensor for 2P (L = 1.5 hbar)
  const Vec3 Lh{0,0,1}, P{1,0,0}, Q{0,1,0}; const double Lp=1.5, eP=std::sqrt(1-std::pow(Lp/n,2));
  const auto orb=orbitAveragedBmtAngularVelocities(a,Lh*(Lp*hbar),Vec3{1,0,0}*(mu*1e-12),Vec3{0,1,0}*(mu*1e-12),mred,0.0,P);
  const double wL=dot(orb.first+orb.second,Lh);   // (omega1+omega2).Lhat
  std::printf("orbital precession: omega1.L = %.4e, omega2.L = %.4e rad/s; transverse |.| = %.2e\n",dot(orb.first,Lh),dot(orb.second,Lh),(orb.first-Lh*dot(orb.first,Lh)).norm());
  auto Uax=[&](const Vec3& s){ return orbitAveragedDipoleEnergy(a,eP,Lh,P,s*mu,s*(-mu)); };
  const double Uiso=(Uax({1,0,0})+Uax({0,1,0})+Uax({0,0,1}))/3.0;
  const double Sm=1.0;  // classical triplet |S| = hbar
  double E[3];
  for(int J=0;J<=2;++J){
    const double Jm=J+0.5, ct=(Jm*Jm-Lp*Lp-Sm*Sm)/(2*Lp*Sm), st=std::sqrt(std::max(0.0,1-ct*ct));
    const double Eso=0.5*hbar*wL*ct*Sm*2.0/2.0*2.0/2.0; // S_i = (hbar/2) Shat each: E = (hbar/2)(omega1+omega2).Shat
    double Ut=0; const int M=256; for(int j=0;j<M;++j){ const double ph=2*pi*(j+0.5)/M; Ut+=Uax(P*(st*std::cos(ph))+Q*(st*std::sin(ph))+Lh*ct)-Uiso; } Ut/=M;
    E[J]=EP+Eso+Ut;
    std::printf("2^3P_%d: cos(S,L) = %+.4f  E_SO = %+.4f GHz  E_tensor = %+.4f GHz  E_rel = %.4f GHz  total %.4f GHz\n",J,ct,Eso*GHz,Ut*GHz,EP*GHz,E[J]*GHz);
  }
  const double ES1=ES+Econt;
  const double meas[3]={18499.65,13012.42,8624.38};
  std::printf("2^3S1: E_rel %.4f + contact %.4f = %.4f GHz\n",ES*GHz,Econt*GHz,ES1*GHz);
  // Measured nu_J = E(2^3S1) - E(2^3P_J) > 0: in positronium 2^3S1 lies ABOVE every 2^3P_J (P0 lowest).
  for(int J=0;J<=2;++J) std::printf("nu_%d = E(S1) - E(P%d): model %+.3f GHz  measured %+.5f GHz  model - measured %+.3f GHz\n",J,J,(ES1-E[J])*GHz,meas[J]/1000,(ES1-E[J])*GHz-meas[J]/1000);
  std::printf("P multiplet: P1-P0 model %.3f measured %.3f ; P2-P1 model %.3f measured %.3f ; P2-P0 model %.3f measured %.3f GHz\n",
    (E[1]-E[0])*GHz,(meas[0]-meas[1])/1000,(E[2]-E[1])*GHz,(meas[1]-meas[2])/1000,(E[2]-E[0])*GHz,(meas[0]-meas[2])/1000);
  const double Dar=a4/(8.0*n*n*n);
  std::printf("terms without a classical counterpart, S states only: quantum Darwin contact alpha^4mc^2/(8n^3) = %.3f GHz; virtual annihilation (3S1) alpha^4mc^2/(4n^3) = %.3f GHz\n",Dar*GHz,a4/(4.0*n*n*n)*GHz);
  for(int J=0;J<=2;++J) std::printf("nu_%d with both added to 2^3S1: %+.3f GHz (measured %.3f, remaining %+.3f)\n",J,(ES1+Dar+a4/(4.0*n*n*n)-E[J])*GHz,meas[J]/1000,(ES1+Dar+a4/(4.0*n*n*n)-E[J])*GHz-meas[J]/1000);
  std::printf("missing virtual annihilation of 2^3S1 (alpha^4 m c^2 / (4 n^3)) = %.3f GHz\n",a4/(4.0*n*n*n)*GHz);
}
