// Audit 386: classical Kepler orbit averages <r^k> at L = (l+1/2) hbar (also sqrt(l(l+1)) and l) vs exact hydrogenic QM.
// Units a_B = 1, hbar = 1; a = n^2, e = sqrt(1 - (L/n)^2); <f> = (1/2pi) int f(r)(1 - e cos E) dE (trapezoid, periodic).
#include <cmath>
#include <cstdio>
static double avg(int n,double L,int k){
  const double a=double(n)*n, e=std::sqrt(std::max(0.0,1.0-(L/n)*(L/n)));
  const int N=1<<20; long double s=0;
  for(int i=0;i<N;++i){ const double E=2*M_PI*(i+0.5)/N, w=1-e*std::cos(E); s+=std::pow(a*w,double(k))*w; }
  return double(s/N); }
static double closedForm(int n,double L,int k){
  switch(k){ case 1: return 0.5*(3.0*n*n-L*L); case 2: return 0.5*n*n*(5.0*n*n-3.0*L*L);
    case -1: return 1.0/(n*n); case -2: return 1.0/(n*n*n*L); case -3: return 1.0/(n*n*n*L*L*L); } return NAN; }
static double qm(int n,int l,int k){
  switch(k){ case 1: return 0.5*(3.0*n*n-l*(l+1.0)); case 2: return 0.5*n*n*(5.0*n*n+1.0-3.0*l*(l+1.0));
    case -1: return 1.0/(n*n); case -2: return 1.0/(n*n*n*(l+0.5)); case -3: return l>0?1.0/(n*n*n*l*(l+0.5)*(l+1.0)):INFINITY; } return NAN; }
// L that makes the classical closed form equal to QM (NAN if none in (0, n]).
static double matchL(int n,int l,int k){
  const double q=qm(n,l,k); double L2;
  switch(k){ case 1: L2=3.0*n*n-2*q; return L2>0?std::sqrt(L2):NAN; case 2: L2=(5.0*n*n-2*q/(n*n))/3; return L2>0?std::sqrt(L2):NAN;
    case -1: return NAN; case -2: return 1.0/(n*n*n*q); case -3: return std::isfinite(q)?std::cbrt(1.0/(n*n*n*q)):NAN; } return NAN; }
int main(){
  const int ks[5]={1,2,-1,-2,-3}; double worstQuad=0;
  std::printf("# n l  k   class(Langer)       QM                  rel.err(Langer)  rel.err(sqrt(l(l+1)))  rel.err(L=l)   L_match\n");
  for(int n=1;n<=5;++n) for(int l=0;l<n;++l) for(int k:ks){
    const double Lg=l+0.5, c=avg(n,Lg,k), cf=closedForm(n,Lg,k), q=qm(n,l,k);
    worstQuad=std::max(worstQuad,std::abs(c/cf-1));
    const double Ls=std::sqrt(l*(l+1.0)), Ll=l;
    const double es=Ls>0?closedForm(n,Ls,k)/q-1:NAN, el=Ll>0?closedForm(n,Ll,k)/q-1:NAN;
    std::printf("%d %d %+d  %.12e  %.12e  %+.6e    %+.6e           %+.6e  %.6f\n",n,l,k,c,q,std::isfinite(q)?c/q-1:NAN,es,el,matchL(n,l,k)); }
  std::printf("# H386a worst |quadrature/closed - 1| = %.3e\n",worstQuad);
  // H386d check: Ps 2P spin-orbit constant.  Breit, equal masses, g = 2(1+a_e): H_SO = (3/2 + 2 a_e) alpha hbar^3/(m^2 c r^3) L.S/hbar^2,
  // <r^-3> QM 2p = 1/(24 a_Ps^3), a_Ps = 2 a0 -> A_QM = (3/2 + 2 a_e)/192 alpha^4 m c^2; classical Langer: x 8/9.
  { const double alpha=7.2973525643e-3, mc2=510998.95, h=4.135667696e-15, ae=1.15965218e-3;
    const double unit=std::pow(alpha,4)*mc2/h/1e9, Aqm=(1.5+2*ae)/192*unit, Acl=Aqm*8.0/9.0, Amodel=2.437;
    // Measured 2^3P gaps (GHz, as used in audit 374b): P1-P0 = A + 6B, P2-P1 = 2A - (12/5)B.
    const double g10=5.487, g21=4.388, Afit=(g21+0.4*g10)/2.4, Bfit=(g10-Afit)/6;
    std::printf("# H386d Ps 2P: A_QM(Breit) = %.4f GHz  A_class(Langer) = %.4f GHz  A_model(363/374) = %.3f GHz  model/class = %.5f\n",Aqm,Acl,Amodel,Amodel/Acl);
    std::printf("# H386d from measured gaps: A = %.4f GHz, B = %.4f GHz (A/A_QM = %.4f); model/A = %.4f\n",Afit,Bfit,Afit/Aqm,Amodel/Afit);
    std::printf("# H386d P2-P1 = 2A - 2.4B: QM-fit %.3f + %.3f; model 4.873 + (-0.54) [374a]\n",2*Afit,-2.4*Bfit);
    std::printf("# H386d P1-P0 = A + 6B:   QM-fit %.3f + %.3f; model 2.437 + (+0.54) [374a]\n",Afit,6*Bfit); }
  for(int n=1;n<=5;++n){ const double maxCl=avg(n,1e-9,2); std::printf("# H386c n=%d s-state: max_L classical <r^2> (L->0) = %.6f  QM = %.6f  ratio %.6f\n",n,maxCl,qm(n,0,2),maxCl/qm(n,0,2)); }
}
