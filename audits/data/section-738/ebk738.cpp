// Audit 738: C++ order-hbar^2 EBK elements (ebkRadialWave1) vs exact QM and ana737 (H738a).
#include "modules/crem_collapse.hpp"
#include <cstdio>
static std::vector<double> grid(int nmax){ const double rmax=3.0*nmax*nmax+60.0; const int M=60000+4000*nmax;
  std::vector<double> r(M); const double xmin=1.0e-4, xmax=std::sqrt(rmax);
  for(int i=0;i<M;++i){ const double x=xmin+(xmax-xmin)*i/(M-1); r[i]=x*x; } return r; }
static std::vector<double> uqm(int n,int l,const std::vector<double>& r){
  std::vector<double> u(r.size()); const int k=n-l-1; const double a=2*l+1;
  for(size_t i=0;i<r.size();++i){ const double rho=2*r[i]/n; long double L0=1, L1=1+a-rho, Lk=(k==0)?L0:L1;
    for(int j=1;j<k;++j){ const long double L2=((2*j+1+a-rho)*L1-(j+a)*L0)/(j+1); L0=L1; L1=L2; Lk=L2; }
    u[i]=static_cast<double>(r[i]*std::pow(rho,l)*std::exp(-rho/2)*Lk); }
  double nn=0; for(size_t i=1;i<r.size();++i) nn+=0.5*(u[i]*u[i]+u[i-1]*u[i-1])*(r[i]-r[i-1]);
  for(double& x: u) x/=std::sqrt(nn); return u; }
static double el(const std::vector<double>& a,const std::vector<double>& b,const std::vector<double>& r,int p){
  double R=0; for(size_t i=1;i<r.size();++i) R+=0.5*(a[i]*b[i]*std::pow(r[i],p)+a[i-1]*b[i-1]*std::pow(r[i-1],p))*(r[i]-r[i-1]); return std::abs(R); }
int main(){
  applyPairFromOption("proton,electron");
  const int T[7][4]={{2,1,1,0},{3,1,1,0},{3,1,2,0},{3,0,2,1},{3,2,2,1},{4,1,1,0},{4,3,3,2}};
  const double py[7]={1.000145,1.000134,1.000034,1.000106,1.000009,1.000129,1.000002};   // ana737 k1-in (6 digits)
  bool okA=true;
  for(int i=0;i<7;++i){ const auto r=grid(std::max(T[i][0],T[i][2]));
    const double q=el(uqm(T[i][0],T[i][1],r),uqm(T[i][2],T[i][3],r),r,1);
    const double R1=ebkRadialElement(T[i][0],T[i][1],T[i][2],T[i][3],1);
    const double ratio=R1/q; okA&=std::abs(ratio-py[i])<2e-6;
    std::printf("%d%d -> %d%d: R_QM %.9f  R_order1 %.9f  ratio %.7f (ana737 %.6f, diff %+.1e)\n",T[i][0],T[i][1],T[i][2],T[i][3],q,R1,ratio,py[i],ratio-py[i]); }
  const int U[3][4]={{3,2,1,0},{3,2,2,0},{4,0,3,2}};
  bool okB=true;
  for(int i=0;i<3;++i){ const auto r=grid(std::max(U[i][0],U[i][2]));
    const double q=el(uqm(U[i][0],U[i][1],r),uqm(U[i][2],U[i][3],r),r,2);
    const double R1=ebkRadialElement(U[i][0],U[i][1],U[i][2],U[i][3],2);
    okB&=std::abs(R1/q-1)<5e-4;
    std::printf("E2 %d%d -> %d%d: R2_QM %.7f  R2_order1 %.7f  ratio-1 %+.2e\n",U[i][0],U[i][1],U[i][2],U[i][3],q,R1,R1/q-1); }
  const int S[6][2]={{1,0},{2,0},{2,1},{3,0},{3,1},{8,1}}; const double pph[6]={1.0002144,1.0000641,1.0000087,1.0000363,1.0000027,1.0000005};
  bool okP=true;
  for(int i=0;i<6;++i){ double ph=0; ebkRadialWave1(S[i][0],S[i][1],grid(S[i][0]),&ph); okP&=std::abs(ph-pph[i])<1e-6;
    std::printf("phase Q_1 %d%d: %.7f (ana737 %.7f)\n",S[i][0],S[i][1],ph,pph[i]); }
  const double A2p=ebkE1Rate(2,1,1,0); std::printf("A_E1(2p->1s) %.6e 1/s -> tau %.6f ns\n",A2p,1e9/A2p);
  std::printf("H738a: elements %d, E2 %d, phase %d -> %s\n",okA,okB,okP,(okA&&okB&&okP)?"TAK":"NIE");
}
