#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
// Audit 498: bound-electron g factor in hydrogen 1s from the engine's spin response to a weak external field.
// A: the secular path itself (orbitAveragedBmtAngularVelocities, external part) -- orbit does not respond to B_ext.
// B: the same Thomas-BMT field with the orbit's Larmor velocity omega_L x r_e added (proton at rest), as the full path's Lorentz force
//    would produce; Coulomb field of the proton at the electron.  Ratio to the free-electron rate, linear part from +B and -B.
int main(){
  applyPairFromOption("proton,electron");
  const double aB=pairBohrRadius(activePair), mu=firstMass*secondMass/(firstMass+secondMass), Bz=1.0e-3, ge=secondGFactor;
  const double qe=secondCharge, me=secondMass, ke=pairCoulombStrength;
  const Vec3 zhat{0,0,1};
  const double freeRate=dot(thomasBmtEffectiveField({},{{},zhat*Bz},ge)*(-qe/me),zhat);
  std::mt19937_64 rng(498); std::uniform_real_distribution<double> U(0,2*pi);
  const int nOr=std::getenv("NOR")?std::atoi(std::getenv("NOR")):400;
  for(double Lh: {0.5,1.0}){
    double sumA=0, sumB=0, sumB0=0; int ok=0;
    for(int k=0;k<nOr;++k){
      const double cz=1-2*(k+0.5)/nOr, sz=std::sqrt(1-cz*cz), ph=pi*(3-std::sqrt(5.0))*k;
      const Vec3 n{sz*std::cos(ph),sz*std::sin(ph),cz};
      const Vec3 t0=std::abs(n.z)<0.9?Vec3{0,0,1}:Vec3{1,0,0};
      Vec3 e1=cross(n,t0); e1=e1*(1/e1.norm()); const Vec3 e2=cross(n,e1); const double ps=U(rng);
      const Vec3 peri=e1*std::cos(ps)+e2*std::sin(ps), quad=cross(n,peri);
      const Vec3 Lvec=n*(Lh*hbar);
      const Vec3 m1=zhat*firstMagneticMoment, m2=zhat*(-secondMagneticMoment);   // moments irrelevant for the external part
      double A[2];
      for(int s=0;s<2;++s){ gExternalMagneticField=zhat*(s==0?Bz:-Bz);
        const OrbitAveragedBmtAngularVelocities r=orbitAveragedBmtAngularVelocities(aB,Lvec,m1,m2,mu,0.0,peri);
        A[s]=r.valid?dot(r.secondExternal,zhat):NAN; }
      gExternalMagneticField={};
      if(!std::isfinite(A[0])||!std::isfinite(A[1])) continue;
      sumA+=0.5*(A[0]-A[1])/freeRate;
      // B: own Kepler nodes, same geometry
      const double ecc=std::sqrt(std::max(0.0,1-Lh*Lh*hbar*hbar/(mu*ke*aB))), ec=std::sqrt(1-ecc*ecc), nmot=std::sqrt(ke/(mu*aB*aB*aB));
      const int N=512; double tot[2]={0,0}, tot0=0;
      for(int j=0;j<N;++j){ const double E=2*pi*(j+0.5)/N, w=1-ecc*std::cos(E);
        const Vec3 rr=peri*(aB*(std::cos(E)-ecc))+quad*(aB*ec*std::sin(E));                    // electron - proton
        const Vec3 vv=(peri*(-aB*nmot*std::sin(E))+quad*(aB*nmot*ec*std::cos(E)))*(1/w);
        const Vec3 re=rr*(firstMass/(firstMass+secondMass)), ve=vv*(firstMass/(firstMass+secondMass));
        const double r=rr.norm(); const Vec3 Ec=rr*(firstCharge/(4*pi*8.8541878128e-12*r*r*r));   // proton's field at the electron
        for(int s=0;s<2;++s){ const Vec3 Bx=zhat*(s==0?Bz:-Bz); const Vec3 wL=Bx*(-qe/(2*me));
          const Vec3 v=ve+cross(wL,re);
          tot[s]+=dot(thomasBmtEffectiveField(v,{Ec,Bx},ge)*(-qe/me),zhat)*w; }
        tot0+=dot(thomasBmtEffectiveField(ve,{Ec,{}},ge)*(-qe/me),zhat)*w; }
      sumB+=0.5*(tot[0]-tot[1])/N/freeRate; ++ok;
    }
    const double al=fineStructureConstant;
    std::printf("L = %.1f hbar (e = %.4f): orientations %d\n",Lh,std::sqrt(std::max(0.0,1-Lh*Lh*hbar*hbar/(mu*ke*aB))),ok);
    std::printf("  A  engine secular external rate / free - 1 = %+.4f e-6   (analytic -alpha^2/2 (m_p/M)^2 = %+.4f e-6)\n",(sumA/ok-1)*1e6,
                -0.5*al*al*std::pow(firstMass/(firstMass+secondMass),2)*1e6);
    std::printf("  B  with the orbit's Larmor velocity        = %+.4f e-6   (Breit -alpha^2/3 = %+.4f e-6; measured -17.709(13) e-6)\n",(sumB/ok-1)*1e6,-al*al/3*1e6);
  }
}
