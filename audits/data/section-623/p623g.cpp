// Audit 623, T3: symmetry of the engine's random laws.  (i) E1 photon about the orbital axis: the engine's
// rotatingDipoleCosineFromUniform and its helicity probability (1 + c)^2/(2(1 + c^2)) (copied from crem_collapse.hpp, the
// non-Zeeman branch).  (ii) 3 gamma: drawContactAnnihilationPhotons with the 1 - c^2/3 law about S and with the xi law (588)
// for m = +1, T0 and three random elliptical states.  argv[1] = draws per case.
#include "modules/crem_collapse.hpp"
#include "modules/contact_annihilation.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
struct Acc { double s=0, s2=0; long n=0; void add(double x){ s+=x; s2+=x*x; ++n; }
  double mean() const { return s/n; } double err() const { const double m=mean(); return std::sqrt(std::max(0.0,s2/n-m*m)/n); } };
static void line(const char* name,const Acc& a,double expect){
  std::printf("   %-34s %+.6f +- %.6f  (oczekiwane %+.6f, odchylka %+.2f sigma)\n",name,a.mean(),a.err(),expect,(a.mean()-expect)/a.err()); }
int main(int argc,char**argv){
  const long N=argc>1?std::atol(argv[1]):2000000;
  std::uint64_t st=splitMix64(std::getenv("SEED")?std::strtoull(std::getenv("SEED"),nullptr,10):623623ULL);
  // (i) E1
  { Acc c,h,hc,c2; double worstMirror=0.0;
    for(long i=0;i<N;++i){ const double u1=drawUniformUnit(st), u3=drawUniformUnit(st);
      const double cs=rotatingDipoleCosineFromUniform(u1);
      worstMirror=std::max(worstMirror,std::abs(cs+rotatingDipoleCosineFromUniform(1.0-u1)));
      const double pplus=std::clamp((1.0+cs)*(1.0+cs)/(2.0*(1.0+cs*cs)),0.0,1.0);
      const double hel=u3<pplus?1.0:-1.0;
      c.add(cs); h.add(hel); hc.add(hel*cs); c2.add(cs*cs); }
    std::printf("(i) foton E1 wokol osi L, %ld losowan\n",N);
    line("<c> (P-nieparzysta)",c,0.0); line("<h> (P-nieparzysta)",h,0.0);
    line("<h c> (P-parzysta)",hc,0.5); line("<c^2> (P-parzysta)",c2,0.4);
    std::printf("   max |c(u) + c(1 - u)| = %.3e (odbicie prawa: u -> 1 - u daje c -> -c)\n",worstMirror);
    // Zeeman photon, Delta M = 0 (sin^2 theta about B); Delta M = +-1 uses the same 1 + c^2 law with the sign of Delta M in
    // the helicity probability, (1 + dM c)^2/(2(1 + c^2)), symmetric under (c, dM) -> (-c, -dM) term by term.
    Acc z; double worstZ=0.0;
    for(long i=0;i<N;++i){ const double u=drawUniformUnit(st), cz=sinSquaredCosineFromUniform(u);
      worstZ=std::max(worstZ,std::abs(cz+sinSquaredCosineFromUniform(1.0-u))); z.add(cz); }
    line("Zeeman dM = 0: <c> wokol B",z,0.0);
    std::printf("   Zeeman dM = 0: max |c(u) + c(1 - u)| = %.3e\n\n",worstZ); }
  // (ii) 3 gamma
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  struct Case { const char* name; Vec3 s1, s2; bool xiLaw; };
  std::uint64_t rs=splitMix64(6230001ULL);
  const auto rnd=[&](){ return contactIsotropicDirection(rs); };
  Case cases[]={{"prawo 1 - c^2/3 wokol S = z",{0,0,1},{0,0,1},false},
                {"xi: m = +1 (s1 = s2 = z)",{0,0,1},{0,0,1},true},
                {"xi: T0 (s1 = z, s2 = -z)",{0,0,1},{0,0,-1},true},
                {"xi: eliptyczny 1",rnd(),rnd(),true},{"xi: eliptyczny 2",rnd(),rnd(),true},{"xi: eliptyczny 3",rnd(),rnd(),true}};
  std::printf("(ii) generator 3 gamma, %ld losowan na przypadek; a = os stanu (Im(xi* x xi) albo Re xi przy T0)\n",N);
  for(const Case& cs: cases){
    const ContactTripletPolarization pol=contactTripletPolarization(cs.s1*g1,cs.s2*g2);
    Vec3 axis{0,0,1};
    if(cs.xiLaw){ const auto& x=pol.xi; using C=std::complex<double>;
      const C cx=std::conj(x[1])*x[2]-std::conj(x[2])*x[1], cy=std::conj(x[2])*x[0]-std::conj(x[0])*x[2], cz=std::conj(x[0])*x[1]-std::conj(x[1])*x[0];
      Vec3 v{std::imag(cx),std::imag(cy),std::imag(cz)};
      if(v.norm()<1e-9) v=Vec3{std::real(x[0]),std::real(x[1]),std::real(x[2])};
      if(v.norm()<1e-9) v=Vec3{std::imag(x[0]),std::imag(x[1]),std::imag(x[2])};
      axis=v*(1.0/v.norm()); }
    const Vec3 spin=(cs.s1+cs.s2)*0.5;
    Acc aP,aCPT,aCP,aN2; long three=0;
    for(long i=0;i<N;++i){
      const ContactAnnihilationPhotons p=drawContactAnnihilationPhotons(0.0,spin.norm()>1e-12?spin:axis,st,cs.xiLaw?&pol:nullptr);
      if(p.count!=3) continue; ++three;
      int o[3]={0,1,2}; std::sort(o,o+3,[&](int a,int b){ return p.fraction[a]>p.fraction[b]; });
      const Vec3 k1=p.direction[o[0]], k2=p.direction[o[1]]; Vec3 n=cross(k1,k2); n=n*(1.0/n.norm());
      const double ak=dot(axis,k1), an=dot(axis,n);
      aP.add(ak); aCPT.add(an); aCP.add(ak*an); aN2.add(an*an); }
    std::printf("  %s (os %+.4f %+.4f %+.4f; 3 gamma %ld)\n",cs.name,axis.x,axis.y,axis.z,three);
    line("<a.k1> (P-nieparzysta)",aP,0.0); line("O_CPT <a.n12> (T-nieparzysta)",aCPT,0.0);
    line("O_CP <(a.k1)(a.n12)> (CP-nieparzysta)",aCP,0.0);
    std::printf("   <(a.n12)^2> (P-parzysta, kontrola)   %.6f +- %.6f%s\n",aN2.mean(),aN2.err(),cs.xiLaw?"":"  (prawo 1 - c^2/3: 0,3)");
  }
}
