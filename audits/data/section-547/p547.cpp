#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
// Audit 547: the engine's import-23 functions -- (1) B = 0 spectrum of the 2P Hamiltonian (Jacobi) = E_J with 2J + 1, singlet 3 x 0; (2) isotropic coherent
// configurations: sum_k p_k s_k vs QM 542 (3P); (3) 1S: from a T0 spin state the draw and the images (w of the m = 0 image vs Breit-Rabi P_S).
int main(){
  applyPairFromOption("electron,positron");
  const double a=4*pairBohrRadius(activePair), m1=firstMagneticMoment, m2=secondMagneticMoment, h=2*pi*hbar;
  double E[3]; positroniumTripletLevels(a,1,m1,m2,E);
  { std::array<double,12> w; std::array<double,144> V; psJacobiEigen<12>(psTwoPHamiltonian(E,0,0),w,V); std::sort(w.begin(),w.end());
    std::printf("E_J [GHz] %.6f %.6f %.6f; B = 0 spectrum:",E[0]/h*1e-9,E[1]/h*1e-9,E[2]/h*1e-9); for(double x:w) std::printf(" %.6f",x/h*1e-9); std::printf("\n"); }
  const double qm[6]={9.46040e-03,6.74008e-02,1.98277e-01,1.95910e-01,1.70929e-01,1.67155e-01}; const double Bs[6]={0.01,0.03,0.1,0.3,1.0,3.0};
  std::mt19937_64 rng(547); std::normal_distribution<double> nd;
  auto iso=[&]{ Vec3 v{nd(rng),nd(rng),nd(rng)}; return v*(1.0/v.norm()); };
  for(int ib=0;ib<6;++ib){ const double B=Bs[ib]; const Vec3 field=Vec3{0.3,-0.4,0.866}*(B/Vec3{0.3,-0.4,0.866}.norm());
    const double a1=m1*B, a2=-m2*B;   // electron first (sgn -1): -sgn|m|B
    std::array<double,12> w; std::array<double,144> V; psJacobiEigen<12>(psTwoPHamiltonian(E,a1,a2),w,V);
    double sk[12]; for(int k=0;k<12;++k){ double s=0; for(int i=0;i<3;++i){ const double v1=V[(4*i+1)*12+k], v2=V[(4*i+2)*12+k]; s+=0.5*(v1-v2)*(v1-v2);} sk[k]=s; }
    Vec3 e1,e2,b; psFieldFrame(field,e1,e2,b); auto fr=[&](const Vec3& v){ return Vec3{dot(v,e1),dot(v,e2),dot(v,b)}; };
    const int N=20000; double sum=0, sum2=0;
    for(int n=0;n<N;++n){ std::complex<double> lc[3], sc[2]; psCoherentL1(fr(iso()),lc); psCoherentHalf(fr(iso()),sc); double tot=0, ps=0;
      for(int k=0;k<12;++k){ std::complex<double> o=0; for(int i=0;i<3;++i) for(int p=0;p<2;++p) for(int q=0;q<2;++q) o+=V[(4*i+2*p+q)*12+k]*lc[i]*sc[p]*sc[q];
        tot+=std::norm(o); ps+=std::norm(o)*sk[k]; }
      sum+=ps/tot; sum2+=(ps/tot)*(ps/tot); }
    const double mean=sum/N, se=std::sqrt((sum2/N-mean*mean)/N);
    // 1S from a pure T0 spin state: rho = |T0><T0|
    gPsFieldEigen=PsFieldEigenState{}; gPsFieldEigen.e1=e1; gPsFieldEigen.e2=e2; gPsFieldEigen.b=b; gPsFieldEigen.a1=a1; gPsFieldEigen.a2=a2;
    for(int x=1;x<=2;++x) for(int y=1;y<=2;++y) gPsFieldEigen.rho[x][y]=0.5;
    const double dE=204.6577e9*h, xB=2*(a1-a2)/dE, PSbr=0.5*(1-1/std::sqrt(1+xB*xB));
    int nS=0, nT=0; double wT=0; std::mt19937_64 r2(7); std::uniform_real_distribution<double> U(0,1);
    for(int n=0;n<4000;++n){ const PsFieldEigenImage im=psFieldEigenAtPhoton(dE,U(r2),2*pi*U(r2)); if(im.singlet) ++nS; else { ++nT; wT=(1-dot(im.n1,im.n2))/2; } }
    std::printf("B %4.2f T: isotropic sum p_k s_k %.5e +- %.1e, QM %.5e, ratio %.4f | 1S from T0: singlet-image fraction %.4f (QM |<S'|T0>|^2 = P_S %.4f), T0' image w %.6f vs Breit-Rabi P_S %.6f\n",
      B,mean,se,qm[ib],mean/qm[ib],double(nS)/4000,PSbr,wT,PSbr);
  }
}
