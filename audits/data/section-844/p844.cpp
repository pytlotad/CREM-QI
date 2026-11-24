// Audit 844 (static probe, two-body, from p843b): spin part of the in-plane LRL rotation rate on the Kepler orbit of an 838 node,
// orbit scaled rp -> s rp (e fixed).  Engine: R = a_form(spins) - a_form(spins = 0), forms Z (L_T), ZD (+841), ZDC (+841 + 843).
// Reference BO'C (c^-2, written out in prereg844.txt): mu R = v x rot Acal - grad H_SS, Acal = (C x d)/r^3, C = mu k q1 q2 W/c^2,
// W = sum sigma_p S_p, sigma_p = (g_p - 1)/(2 m_p^2) + g_p/(2 m1 m2), S_p = m_p/gamma_p, gamma_p = g_p q_p/(2 m_p).
// argv: pair rp e m1x m1y m1z m2x m2y m2z N s
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char**argv){
  if(argc<12){ std::fprintf(stderr,"usage\n"); return 1; }
  applyPairFromOption(argv[1]); gMeasurementConservativeForce=true;
  const double sc=std::atof(argv[11]), rp=std::atof(argv[2])*sc, ecc=std::atof(argv[3]);
  const Vec3 m1{std::atof(argv[4]),std::atof(argv[5]),std::atof(argv[6])}, m2{std::atof(argv[7]),std::atof(argv[8]),std::atof(argv[9])};
  const int N=std::atoi(argv[10]);
  const double M=firstMass+secondMass, mu=firstMass*secondMass/M, k=pairCoulombStrength, a=rp/(1-ecc);
  const double nmot=std::sqrt(k/(mu*a*a*a)), ec=std::sqrt(1-ecc*ecc), cc=c*c;
  const Vec3 P{1,0,0}, Q{0,1,0};
  // reference constants (independent of the engine's force code: only charges, masses, g factors)
  const double q1=firstCharge, q2=secondCharge, g1=firstGFactor, g2=secondGFactor;
  const double ga1=g1*q1/(2*firstMass), ga2=g2*q2/(2*secondMass);
  const double s1=(g1-1)/(2*firstMass*firstMass)+g1/(2*firstMass*secondMass), s2=(g2-1)/(2*secondMass*secondMass)+g2/(2*firstMass*secondMass);
  const Vec3 W=m1*(s1/ga1)+m2*(s2/ga2), Cv=W*(mu*coulomb*q1*q2/cc);
  const double km=mu0/(4*pi);
  const auto refAcc=[&](const Vec3& d,const Vec3& v){
    const double r=d.norm(); const Vec3 n=d*(1/r); const double r3=r*r*r, r4=r3*r, r5=r4*r, r7=r5*r*r;
    const Vec3 rotA=(n*(3*dot(Cv,n))-Cv)*(1/r3);
    const double A=dot(m1,m2), a1=dot(m1,d), a2=dot(m2,d);
    const Vec3 gH=(d*(-3*A/r5)-(m1*a2+m2*a1)*(3/r5)+d*(15*a1*a2/r7))*km;
    return (cross(v,rotA)-gH)*(1/mu); };
  const auto acc=[&](const State& s,int form){   // 0 Z, 1 ZD, 2 ZDC
    gSpinOrbitDirac=(form>=1); gSpinOrbitCanonical=(form==2); MutualForces F=allExternalForces(s); gSpinOrbitDirac=false; gSpinOrbitCanonical=false;
    return relativisticAcceleration(s.firstVelocity,F.first,firstMass)-relativisticAcceleration(s.secondVelocity,F.second,secondMass); };
  double om[4]={0,0,0,0}, sw=0;
  for(int j=0;j<N;++j){
    const double E=2*pi*(j+0.5)/N, w=1-ecc*std::cos(E);
    const Vec3 r=P*(a*(std::cos(E)-ecc))+Q*(a*ec*std::sin(E)), v=(P*(-a*nmot*std::sin(E))+Q*(a*nmot*ec*std::cos(E)))*(1/w);
    State s{}; s.firstPosition=r*(secondMass/M); s.secondPosition=r*(-firstMass/M); s.firstVelocity=v*(secondMass/M); s.secondVelocity=v*(-firstMass/M);
    State z=s;
    s.firstProperDipole=m1; s.secondProperDipole=m2; s.firstDipole=m1; s.secondDipole=m2; synchronizeCovariantDipoles(s); synchronizeCovariantDipoles(z);
    const Vec3 hv=cross(r,v);
    for(int f=0;f<4;++f){
      const Vec3 R=(f<3?acc(s,f)-acc(z,f):refAcc(r,v)), dA=cross(R,hv)+cross(v,cross(r,R));
      om[f]+=dot(Q,dA)/(k*ecc/mu)*w; }
    sw+=w;
  }
  std::printf("SPIN844 pair %s s %g e %.12f Z %.15e ZD %.15e ZDC %.15e BOC %.15e\n",argv[1],sc,ecc,om[0]/sw,om[1]/sw,om[2]/sw,om[3]/sw);
}
