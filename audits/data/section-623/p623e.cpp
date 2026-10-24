// Audit 623, T2b: are the field-eigenstate draws C-symmetric in distribution?  Under C the field reverses (B -> -B), the
// moments reverse at fixed spins, so a_i -> -a_i and the frame axis b -> -b; the engine's Hamiltonian becomes Pi H Pi
// (up <-> down in both slots) and its eigenvectors come out in another order, so one uniform u picks another index.
// Here the engine's own psFieldEigenAtPhoton (1S, 4 x 4) and psFieldEigenAtFormation (2P, 12 x 12) are scanned over a grid
// of u for random configurations, original and C image, and the outcome distributions are compared in physical terms:
// weight of the drawn state, spin projections on the ORIGINAL field direction, singlet weight.  argv: grid1S grid2P cases.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <map>
#include <tuple>
using Key=std::tuple<long long,long long,long long,long long>;
static long long rq(double x,double s){ return std::llround(x*s); }
int main(int argc,char**argv){
  const int grid1=argc>1?std::atoi(argv[1]):100000, grid2=argc>2?std::atoi(argv[2]):20000, cases=argc>3?std::atoi(argv[3]):60;
  std::uint64_t rs=splitMix64(6230623ULL);
  const auto dir=[&](){ const double z=2*drawUniformUnit(rs)-1, ph=2*pi*drawUniformUnit(rs), s=std::sqrt(1-z*z); return Vec3{s*std::cos(ph),s*std::sin(ph),z}; };
  const double m=firstMagneticMoment;          // electron = positron moment
  const double fields[3]={0.01,0.1,1.0};
  const double hfs=203.39e9*2*pi*hbar;         // 1S splitting scale (J); only sets the mixing, not compared to data
  double worst1=0.0, worst2=0.0, unp1=0.0, unp2=0.0; long classes1=0, classes2=0; int mism1=0, mism2=0;
  for(int c=0;c<cases;++c){
    const double Bm=fields[c%3]; const Vec3 bHat=dir(), s1=dir(), s2=dir(), L=dir();
    const double psi=2*pi*drawUniformUnit(rs);
    // ---- 1S at the photon: rho = product of the two spin coherent states (in the frame of the run)
    std::map<Key,double> f[2];
    for(int img=0;img<2;++img){
      PsFieldEigenState st; const Vec3 B=bHat*(img?-Bm:Bm); psFieldFrame(B,st.e1,st.e2,st.b);
      // original: slot 1 electron (sgn -1), slot 2 positron (+1); image: slot 1 positron, slot 2 electron, same spins
      const double sg1=img?1.0:-1.0, sg2=-sg1; st.a1=-sg1*m*Bm; st.a2=-sg2*m*Bm;
      const auto inF=[&](const Vec3& v){ return Vec3{dot(v,st.e1),dot(v,st.e2),dot(v,st.b)}; };
      std::complex<double> c1[2], c2[2]; psCoherentHalf(inF(s1),c1); psCoherentHalf(inF(s2),c2);
      std::complex<double> p4[4]; for(int a=0;a<2;++a) for(int b=0;b<2;++b) p4[2*a+b]=c1[a]*c2[b];
      for(int a=0;a<4;++a) for(int cc=0;cc<4;++cc) st.rho[a][cc]=std::real(p4[a]*std::conj(p4[cc]));
      gPsFieldEigen=st;
      for(int j=0;j<grid1;++j){ const PsFieldEigenImage im=psFieldEigenAtPhoton(hfs,(j+0.5)/grid1,psi);
        f[img][Key{rq(im.q,1e8),rq(dot(im.n1,bHat),1e6),rq(dot(im.n2,bHat),1e6),rq(im.singletWeight,1e8)}]+=1.0/grid1; }
    }
    for(int img=0;img<2;++img) for(const auto& kv: f[img]){ const auto it=f[1-img].find(kv.first);
      const double d=std::abs(kv.second-(it==f[1-img].end()?0.0:it->second)); worst1=std::max(worst1,d); if(it==f[1-img].end()){ ++mism1; unp1=std::max(unp1,std::get<0>(kv.first)*1e-8); } }
    classes1+=f[0].size();
    // ---- 2P at formation (triplet along s1), the engine's own levels at n = 2
    if(c<cases){ std::map<Key,double> g[2];
      double E[3]; positroniumTripletLevels(4.0*pairBohrRadius(activePair),1,m,m,E);
      for(int img=0;img<2;++img){ const Vec3 B=bHat*(img?-Bm:Bm); const double sg1=img?1.0:-1.0, sg2=-sg1;
        for(int j=0;j<grid2;++j){ psFieldEigenAtFormation(E,-sg1*m*Bm,-sg2*m*Bm,B,L,s1,false,(j+0.5)/grid2,1ULL);
          const PsFieldEigenState& st=gPsFieldEigen;
          // physical <sigma_i . bHat>: diagonal of rho in the run's frame, sign of its b relative to bHat
          const double sb=dot(st.b,bHat);
          const double z1=sb*(st.rho[0][0]+st.rho[1][1]-st.rho[2][2]-st.rho[3][3]);
          const double z2=sb*(st.rho[0][0]-st.rho[1][1]+st.rho[2][2]-st.rho[3][3]);
          g[img][Key{rq(st.probability,1e8),rq(z1,1e6),rq(z2,1e6),rq(st.singletWeight,1e8)}]+=1.0/grid2; } }
      for(int img=0;img<2;++img) for(const auto& kv: g[img]){ const auto it=g[1-img].find(kv.first);
        const double d=std::abs(kv.second-(it==g[1-img].end()?0.0:it->second)); worst2=std::max(worst2,d);
        if(it==g[1-img].end()){ ++mism2; unp2=std::max(unp2,std::get<0>(kv.first)*1e-8); if(std::getenv("VERBOSE")) std::printf("  2P bez pary: konf %d B %.2f obraz %d p %lld z1 %lld z2 %lld s %lld udzial %.2e\n",c,Bm,img,
          std::get<0>(kv.first),std::get<1>(kv.first),std::get<2>(kv.first),std::get<3>(kv.first),kv.second); } }
      classes2+=g[0].size(); }
  }
  std::printf("1S przy fotonie: %d konfiguracji (B = 0,01 / 0,1 / 1 T), siatka u %d: klas %ld, klas bez pary %d (najwieksza waga takiej klasy %.2e), max |p - p_C| = %.3e (rozdzielczosc 1/siatka = %.1e)\n",
              cases,grid1,classes1,mism1,unp1,worst1,1.0/grid1);
  std::printf("2P przy powstaniu: %d konfiguracji, siatka u %d: klas %ld, klas bez pary %d (najwieksza waga takiej klasy %.2e), max |p - p_C| = %.3e (rozdzielczosc %.1e)\n",
              cases,grid2,classes2,mism2,unp2,worst2,1.0/grid2);
}
