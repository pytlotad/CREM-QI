#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 377: spin-spin coupling in the secular transport vs Quigga-Rosner.
// For each (n, L): trace of the partner part of the transport rate (as the s-state isotropy rule extracts it), at default nodes and
// x4 nodes; the Plummer contact term alone averaged over the same Kepler orbit (fine quadrature); the QR value; energies.
static double mred,a1,g1;
static Vec3 ax(int k){ return k==0?Vec3{1,0,0}:(k==1?Vec3{0,1,0}:Vec3{0,0,1}); }

// partner trace / 3 for particle 1 (rad/s per A m^2), engine call with given node override (0 = default)
static double partnerTrace(double a,const Vec3& L,int nodes,int& usedNodes){
  const double tiny=1e-12;
  const auto o=orbitAveragedBmtAngularVelocities(a,L,ax(0)*firstMagneticMoment*tiny,ax(0)*secondMagneticMoment*tiny,mred,0.0,Vec3{1,0,0},nodes);
  usedNodes=o.phaseNodes; double t=0;
  for(int k=0;k<3;++k){ const auto p=orbitAveragedBmtAngularVelocities(a,L,ax(k)*firstMagneticMoment,ax(k)*secondMagneticMoment,mred,0.0,Vec3{1,0,0},nodes);
    t+=dot(ax(k),p.first-o.first)/secondMagneticMoment; }
  return t/3;
}
// <n_Plummer> over the Kepler orbit, uniform in time via E nodes (fine)
static double plummerDensity(double a,double e,int nodes){
  const double eps=magneticDipoleRadius(); double s=0,w=0;
  for(int i=0;i<nodes;++i){ const double E=2*pi*(i+0.5)/nodes, tw=1-e*std::cos(E), r=a*tw;
    s+=tw*3*eps*eps/(4*pi*std::pow(r*r+eps*eps,2.5)); w+=tw; }
  return s/w;
}
int main(){
  mred=firstMass*secondMass/(firstMass+secondMass); a1=pairBohrRadius(activePair); g1=std::abs(firstGyromagneticRatioOf());
  const double mu=firstMagneticMoment, eps=magneticDipoleRadius();
  std::printf("a_pair %.6e m, eps %.6e m (= %.5f r*), |gamma| %.6e, mu %.6e\n",a1,eps,eps/comptonBarrierRadius,g1,mu);
  const double states[][2]={{1,0.5},{2,0.5},{2,1.5},{3,0.5},{3,1.5},{3,2.5}};
  std::printf("\n n  L/hb   e        nodes  trace/3 default   x4 nodes    (2/3)mu0<n>  ratio  <n>/a^-3     n_QR/a^-3   <n>/n_QR    rate transport  rate QR\n");
  for(auto& s: states){
    const double n=s[0], Lh=s[1], a=n*n*a1;
    const double e=std::sqrt(std::max(0.0,1-Lh*Lh/(n*n)));
    const Vec3 L{0,0,Lh*hbar};
    int used=0,used4=0; const double t=partnerTrace(a,L,0,used); const double t4=partnerTrace(a,L,4*used,used4);
    const double nP=plummerDensity(a,e,1<<20), nP2=orbitAveragedContactDensity(a,e);
    const double contact=(2.0/3.0)*mu0*nP;   // trace/3 of mu0 n I + poles (-mu0 n) -> (2/3) mu0 n
    const double nQR=Lh<1.0?1.0/(2*pi*a1*a1*a1*n*n*n*Lh):0.0;
    std::printf("%2.0f  %4.1f  %.5f  %5d  %.6e  %.6e  %.6e  %.4f  %.4e  %.4e  %.4e  %.4e  %.4e\n",n,Lh,e,used,t,t4,contact,t/contact,
      nP*a1*a1*a1,nQR*a1*a1*a1,nQR>0?nP/nQR:0.0,std::abs(t)*mu,(2.0/3.0)*mu0*nQR*g1*mu);
    std::printf("      engine orbitAveragedContactDensity / fine: %.6f\n",nP2/nP);
  }
  // energies at n = 1, L = hbar/2: orbit-averaged U for aligned moments (para, cos +1) and anti-aligned (ortho, cos -1), several axes
  std::printf("\nEnergy (orbitAveragedDipoleEnergy, n = 1, L = hbar/2, e = %.5f):\n",std::sqrt(0.75));
  const double e1=std::sqrt(0.75);
  double meanSplit=0;
  for(int k=0;k<3;++k){
    const Vec3 m1=ax(k)*mu, m2p=ax(k)*secondMagneticMoment;
    const double Up=orbitAveragedDipoleEnergy(a1,e1,Vec3{0,0,1},Vec3{1,0,0},m1,m2p);
    const double Uo=orbitAveragedDipoleEnergy(a1,e1,Vec3{0,0,1},Vec3{1,0,0},m1,m2p*-1.0);
    std::printf("  moments along %c: U(aligned) %.6e J, U(anti) %.6e J, split %.6e J = %.6f GHz\n","xyz"[k],Up,Uo,Uo-Up,(Uo-Up)/(2*pi*hbar)/1e9);
    meanSplit+=(Uo-Up)/3;
  }
  const double nQR1=1.0/(2*pi*a1*a1*a1*0.5);
  const double contactSplitQR=2*(2.0/3.0)*mu0*mu*mu*nQR1;   // U = -(2/3) mu0 n m1.m2 ; aligned - anti = 2 (2/3) mu0 n mu^2
  std::printf("  orientation-averaged split %.6e J = %.6f GHz; Fermi with n_QR %.6e J = %.4f GHz; ratio %.4e\n",meanSplit,meanSplit/(2*pi*hbar)/1e9,
    contactSplitQR,contactSplitQR/(2*pi*hbar)/1e9,meanSplit/contactSplitQR);
}
