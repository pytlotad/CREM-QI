#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 389: classical tensor dipole-dipole energy of Ps n = 1, L = hbar/2 without the s-state isotropy rule,
// as a function of the spin direction relative to the orbit (axis z, periapsis x); o-Ps (moments antiparallel), p-Ps (parallel).
// Independent quadrature (uniform eccentric anomaly, time weight 1 - e cos E, no softening) vs the engine's
// orbitAveragedDipoleEnergy with its orbit-contact part removed.  Energies in GHz (E/h).
int main(){
  const double mred=firstMass*secondMass/(firstMass+secondMass), a=pairBohrRadius(activePair), L=0.5*hbar;
  const double e=std::sqrt(1.0-L*L/(mred*pairCoulombStrength*a)), h=2*pi*hbar, GHz=1e9*h;
  const Vec3 axis{0,0,1}, peri{1,0,0};
  const double m1=firstMagneticMoment, m2=secondMagneticMoment, nqr=quiggaRosnerContactDensity(a,L);
  std::printf("# a = %.6e m, e = %.6f, n_QR = %.6e m^-3, |psi(0)|^2 = 1/(pi a^3) = %.6e\n",a,e,nqr,1.0/(pi*a*a*a));
  auto quad=[&](const Vec3& u1,const Vec3& u2){ const int N=1<<16; long double s=0;
    for(int i=0;i<N;++i){ const double E=2*pi*(i+0.5)/N, w=1-e*std::cos(E);
      const Vec3 r=peri*(a*(std::cos(E)-e))+Vec3{0,1,0}*(a*std::sqrt(1-e*e)*std::sin(E)); const double rr=r.norm(); const Vec3 n=r/rr;
      s+=w*(mu0/(4*pi))*(m1*m2)*(dot(u1,u2)-3*dot(u1,n)*dot(u2,n))/(rr*rr*rr); }
    return double(s/N); };
  // Same quadrature with the engine's Plummer field (softening = magneticDipoleRadius), contact (magnetization) term excluded.
  auto quadPlummer=[&](const Vec3& u1,const Vec3& u2){ const int N=1<<16; long double s=0; const double soft=magneticDipoleRadius();
    for(int i=0;i<N;++i){ const double E=2*pi*(i+0.5)/N, w=1-e*std::cos(E);
      const Vec3 r=peri*(a*(std::cos(E)-e))+Vec3{0,1,0}*(a*std::sqrt(1-e*e)*std::sin(E));
      s+=w*(-dot(u1*m1,plummerDipoleField(r,u2*m2,soft))); }
    return double(s/N); };
  auto engine=[&](const Vec3& u1,const Vec3& u2){ const Vec3 M1=u1*m1, M2=u2*m2;
    return orbitAveragedDipoleEnergy(a,e,axis,peri,M1,M2)+(2.0/3.0)*mu0*dot(M1,M2)*orbitAveragedContactDensity(a,e); };
  const double dUF=(4.0/3.0)*mu0*m1*m2*nqr;
  std::printf("# contact o-p difference (4/3) mu0 m^2 n_QR = %.4f GHz\n",dUF/GHz);
  for(int ch=0;ch<2;++ch){ const double sgn=ch==0?-1.0:1.0; // o-Ps: antiparallel moments
    double mn=1e300,mx=-1e300,worst=0,worstP=0,worstQP=0,absG=0,absP=0,absS=0; std::printf("# %s: theta(deg from L) phi(deg from periapsis)  E_tensor quad [GHz]  engine [GHz]\n",ch==0?"o-Ps":"p-Ps");
    for(int it=0;it<=6;++it) for(int ip=0;ip<=6;++ip){ const double th=pi/2*it/6, ph=pi/2*ip/6;
      const Vec3 u{std::sin(th)*std::cos(ph),std::sin(th)*std::sin(ph),std::cos(th)};
      const double q=quad(u,u*sgn), g=engine(u,u*sgn); mn=std::min(mn,q); mx=std::max(mx,q); worst=std::max(worst,std::abs(g-q)/std::abs(q));
      const double qp=quadPlummer(u,u*sgn); worstP=std::max(worstP,std::abs(g-qp)/std::abs(q)); worstQP=std::max(worstQP,std::abs(qp-q)/std::abs(q)); absG=std::max(absG,std::abs(g-q)); absP=std::max(absP,std::abs(g-qp)); absS=std::max(absS,std::abs(qp-q));
      if(ip==0||ip==6||it==0) std::printf("  %5.1f %5.1f   %+.6f   %+.6f\n",th*180/pi,ph*180/pi,q/GHz,g/GHz); }
    std::printf("# %s orientation splitting (max-min) = %.4f GHz = %.4f x contact difference; worst |engine/quad - 1| = %.2e\n",ch==0?"o-Ps":"p-Ps",(mx-mn)/GHz,(mx-mn)/dUF,worst);
    std::printf("#   softening (quad Plummer vs quad bare): %.2e;  engine vs quad Plummer (engine's own quadrature): %.2e  [softening radius %.3e m, periapsis %.3e m]\n",worstQP,worstP,magneticDipoleRadius(),a*(1-e));
    std::printf("#   absolute, as a fraction of the splitting: engine-quad %.2e, softening %.2e, engine-quadPlummer %.2e\n",absG/(mx-mn),absS/(mx-mn),absP/(mx-mn)); }
}
