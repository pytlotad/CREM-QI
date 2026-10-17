#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 499: secular Larmor precession -- tests.  argv[1] = pair.  T1: bound-electron spin rate in a weak external field (H 1s);
// T2: L and periapsis turn about B at omega_L; T3 (Ps): omega_L = 0.
int main(int argc,char**argv){
  applyPairFromOption(argc>1?argv[1]:"proton,electron");
  const double aB=pairBohrRadius(activePair), mu=firstMass*secondMass/(firstMass+secondMass);
  const Vec3 z{0,0,1};
  std::printf("pair %s: omega_L/B = %.6e rad/s/T (e/2m_e = %.6e)\n",argc>1?argv[1]:"proton,electron",
    (gExternalMagneticField=z, secularLarmorAngularVelocity().z), eCharge/(2*9.1093837015e-31));
  // T1 (as 498 A, Fibonacci orientations, Langer and circle)
  const double Bz=1e-3; const double freeRate=dot(thomasBmtEffectiveField({},{{},z*Bz},secondGFactor)*(-secondCharge/secondMass),z);
  for(double Lh: {0.5,1.0}){ double sum=0; int ok=0; const int nOr=4000;
    for(int k=0;k<nOr;++k){ const double cz=1-2*(k+0.5)/nOr, sz=std::sqrt(1-cz*cz), ph=pi*(3-std::sqrt(5.0))*k;
      const Vec3 n{sz*std::cos(ph),sz*std::sin(ph),cz}; const Vec3 t0=std::abs(n.z)<0.9?Vec3{0,0,1}:Vec3{1,0,0};
      Vec3 peri=cross(n,t0); peri=peri*(1/peri.norm()); const double ps=2.399963*k; peri=peri*std::cos(ps)+cross(n,peri)*std::sin(ps);
      double A[2]; for(int s=0;s<2;++s){ gExternalMagneticField=z*(s?-Bz:Bz);
        const auto r=orbitAveragedBmtAngularVelocities(aB,n*(Lh*hbar),z*firstMagneticMoment,z*(-secondMagneticMoment),mu,0.0,peri);
        A[s]=r.valid?dot(r.secondExternal,z):NAN; }
      if(std::isfinite(A[0])&&std::isfinite(A[1])){ sum+=0.5*(A[0]-A[1])/freeRate; ++ok; } }
    std::printf("T1 L = %.1f hbar: second particle external rate / free - 1 = %+.4f e-6 (%d orientations)\n",Lh,(sum/ok-1)*1e6,ok); }
  // T2: transport over a short time, B = 1 T, L tilted 60 deg, spins tiny-coupled as prepared
  gExternalMagneticField=z*1.0;
  const Vec3 wL=secularLarmorAngularVelocity();
  SecularSpinOrbitState st{}; const double th=pi/3; const Vec3 Lhat{std::sin(th),0,std::cos(th)};
  st.orbitalAngularMomentum=Lhat*(1.0*hbar); const double sc=std::getenv("TINY")?1e-12:1.0; st.firstDipole=Vec3{0,1,0}*(firstMagneticMoment*sc); st.secondDipole=Vec3{0,1,0}*(-secondMagneticMoment*sc);
  st.periapsisDirection=Vec3{0,1,0};
  const double T=std::abs(wL.z)>0?0.05/std::abs(wL.z):1e-12;
  const SecularSpinOrbitAdvance adv=advanceCoupledSecularSpinOrbit(st,aB,mu,T);
  const auto az=[](Vec3 v){ return std::atan2(v.y,v.x); };
  const double dL=az(adv.state.orbitalAngularMomentum)-az(st.orbitalAngularMomentum), dP=az(adv.state.periapsisDirection)-az(st.periapsisDirection);
  std::printf("T2 B 1 T, t = %.3e s: completed %d; L azimuth rate / omega_L.z = %.6f, periapsis azimuth rate / omega_L.z = %.6f, |L| ratio %.9f, L.z ratio %.9f\n",
    T,(int)adv.completed,wL.z!=0?dL/T/wL.z:0.0,wL.z!=0?dP/T/wL.z:0.0,adv.state.orbitalAngularMomentum.norm()/st.orbitalAngularMomentum.norm(),
    adv.state.orbitalAngularMomentum.z/st.orbitalAngularMomentum.z);
}
