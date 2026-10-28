// Audit 643: is L_kan - L_mech of the engine's own spin-orbit Hamiltonian (Breit-Pauli, linear in p) equal to the hidden angular momentum Delta J_hid?
// (A) spin-orbit precession rates of the engine (orbitAveragedBmtAngularVelocities, moments ~0, no field): omega_i . L^ against the Breit-Pauli c_i;
// (B) Delta J_hid (hiddenAngularMomentumOrbitAveraged) against -mu <r^2 W_perp> with the engine's c_i, random spins and orbits.  argv[1] = pair.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
int main(int argc,char**argv){
  applyPairFromOption(argc>1?argv[1]:"electron,positron");
  gExternalMagneticField=Vec3{};
  const double m1=firstMass, m2=secondMass, M=m1+m2, mu=m1*m2/M, g1=std::abs(firstGFactor), g2=std::abs(secondGFactor);
  const double K=pairCoulombStrength;                                // k |q1 q2|
  const double c1=(g1-1)/(2*m1*m1)+g1/(2*m1*m2), c2=(g2-1)/(2*m2*m2)+g2/(2*m1*m2);
  const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1;
  std::printf("pair %s: g %.9f %.9f; Breit-Pauli c1 %.9e c2 %.9e (c1/c2 %.12f); predicted f1 %.9f f2 %.9f\n",argc>1?argv[1]:"electron,positron",
    g1,g2,c1,c2,c1/c2,1-m2/(g1*M),1-m1/(g2*M));
  const Vec3 z{0,0,1}, x{1,0,0}; const double aB=pairBohrRadius(activePair);
  // (A) rates along L, tiny moments (spin-spin off), Kepler orbit n = 2, l = 1
  for(double Lh: {1.5,0.5}){
    const double n=2, a=n*n*aB, L=Lh*hbar, e=std::sqrt(std::max(0.0,1-std::pow(Lh/n,2)));
    const double eps=1e-12;
    const auto r=orbitAveragedBmtAngularVelocities(a,z*L,x*(firstMagneticMoment*eps),x*(secondMagneticMoment*eps),mu,0.0,x);
    const double r3K=1.0/(a*a*a*std::pow(1-e*e,1.5));
    const double w1=dot(r.first,z), w2=dot(r.second,z);
    // H = c_i K/(c^2 r^3) L . S_i  ->  omega_i = c_i K <r^-3> L / c^2
    std::printf("(A) L = %.1f hbar: omega1.L^ %.6e omega2.L^ %.6e rad/s; ratio %.12f (c1/c2 %.12f); omega_i / (c_i K <r^-3>_Kepler L / c^2) = %.9f %.9f\n",
      Lh,w1,w2,w1/w2,c1/c2,w1/(c1*K*r3K*L/(c*c)),w2/(c2*K*r3K*L/(c*c)));
  }
  // (B) Delta J_hid vs canonical shift, same orbit average (Kepler nodes, as hiddenAngularMomentumOrbitAveraged)
  std::mt19937_64 rng(643); std::uniform_real_distribution<double> U(0,1);
  const auto iso=[&]{ const double cz=2*U(rng)-1, ph=2*pi*U(rng), sz=std::sqrt(1-cz*cz); return Vec3{sz*std::cos(ph),sz*std::sin(ph),cz}; };
  double worstSign=1, worstRel=0;
  for(int k=0;k<2000;++k){
    const double n=1+3*U(rng), a=n*n*aB, e=0.95*U(rng);
    Vec3 ax=iso(), pe=cross(ax,iso()); pe=pe*(1/pe.norm());
    const Vec3 sp1=iso(), sp2=iso();                                   // spin unit vectors
    const Vec3 mm1=sp1*(s1*firstMagneticMoment), mm2=sp2*(s2*secondMagneticMoment);
    const Vec3 dJ=hiddenAngularMomentumOrbitAveraged(a,e,ax,pe,mm1,mm2);
    // canonical: -mu K/c^2 sum_i c_i <S_i_perp / r>, S_i = (hbar/2) s_i; the hidden function averages <(m - (m.r^)r^)/r> -- reuse it with
    // "moments" chosen so that k_m q_j m_i' = -mu K c_i S_i / c^2, i.e. m_i' = -mu K c_i (hbar/2) s_i / (c^2 k_m q_j)
    const double km=1e-7;
    const Vec3 mc1=sp1*(-mu*K*c1*0.5*hbar/(c*c*km*secondCharge)), mc2=sp2*(-mu*K*c2*0.5*hbar/(c*c*km*firstCharge));
    const Vec3 dL=hiddenAngularMomentumOrbitAveraged(a,e,ax,pe,mc1,mc2);
    // per-spin comparison: split by setting one spin's moment to zero
    const Vec3 dJ1=hiddenAngularMomentumOrbitAveraged(a,e,ax,pe,mm1,Vec3{}), dL1=hiddenAngularMomentumOrbitAveraged(a,e,ax,pe,mc1,Vec3{});
    const Vec3 dJ2=dJ-dJ1, dL2=dL-dL1;
    worstSign=std::min(worstSign,dot(dL,dJ)/(dL.norm()*dJ.norm()));
    const double f1=1-m2/(g1*M), f2=1-m1/(g2*M);
    worstRel=std::max({worstRel,(dL1-dJ1*f1).norm()/dJ1.norm(),(dL2-dJ2*f2).norm()/dJ2.norm()});
    if(k<3) std::printf("(B) sample %d: |dJ_hid| %.6e hbar, |dL_can| %.6e hbar, cos %.9f, per spin dL_i/dJ_i %.9f %.9f\n",k,dJ.norm()/hbar,dL.norm()/hbar,
      dot(dL,dJ)/(dL.norm()*dJ.norm()),dot(dL1,dJ1)/dJ1.squaredNorm(),dot(dL2,dJ2)/dJ2.squaredNorm());
  }
  std::printf("(B) 2000 configurations: min cos(dL_can, dJ_hid) %.9f, max |dL_i - f_i dJ_i|/|dJ_i| %.3e\n",worstSign,worstRel);
}
