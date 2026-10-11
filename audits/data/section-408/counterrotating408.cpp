// Audit 408: E1 spectrum of a Kepler orbit split into modes (s, m): s = time harmonic of the mean motion,
// m = +1 co-rotating / -1 counter-rotating about L.  zeta = x + i y; P_{s,+} ~ s^4 |c_s|^2, P_{s,-} ~ s^4 |c_{-s}|^2.
// Checks sum_s (P_{s,+} - P_{s,-})/s = f(e) P, f = (1-e^2)^{3/2}/(1+e^2/2), then the "one photon = one mode" scheme:
// a photon of mode (s, m) carries s hbar of J (energy s hbar omega) and m hbar of L; count rate per mode
// P_{s,m}/(s hbar omega); modes that land outside the level grid (n' >= 1, 0 <= l' <= n'-1) are absent.
#include <cmath>
#include <complex>
#include <cstdio>
#include <map>
#include <vector>
struct Spectrum { std::vector<double> Pp, Pm; double P=0; };
static Spectrum modes(double e,int N=8192,int smax=80){
  std::vector<std::complex<double>> z(N);
  for(int j=0;j<N;++j){ const double M=2*M_PI*j/N; double E=M; for(int it=0;it<100;++it) E-=(E-e*std::sin(E)-M)/(1-e*std::cos(E));
    z[j]={std::cos(E)-e,std::sqrt(1-e*e)*std::sin(E)}; }
  Spectrum S; S.Pp.resize(smax+1); S.Pm.resize(smax+1);
  for(int s=1;s<=smax;++s){ std::complex<double> cp=0,cm=0;
    for(int j=0;j<N;++j){ const double M=2*M_PI*j/N; cp+=z[j]*std::polar(1.0,-s*M); cm+=z[j]*std::polar(1.0,s*M); }
    cp/=N; cm/=N; S.Pp[s]=std::pow(s,4)*std::norm(cp); S.Pm[s]=std::pow(s,4)*std::norm(cm); S.P+=S.Pp[s]+S.Pm[s]; }
  return S; }
static double f(double e){ return std::pow(1-e*e,1.5)/(1+e*e/2); }
static int part2();
int mainOld(){ part2();

  for(double e: {0.0,0.3,0.6614,0.8660}){ const Spectrum S=modes(e); double torque=0;
    for(int s=1;s<(int)S.Pp.size();++s) torque+=(S.Pp[s]-S.Pm[s])/s;
    std::printf("e=%.4f: sum(P+ - P-)/s / P = %.6f  f(e) = %.6f   P(1,+)=%.4f P(1,-)=%.4f P(2,+)=%.4f P(2,-)=%.4f P(3,+)=%.4f (fractions of P)\n",
      e,torque/S.P,f(e),S.Pp[1]/S.P,S.Pm[1]/S.P,S.Pp[2]/S.P,S.Pm[2]/S.P,S.Pp[3]/S.P); }
  struct Level { int n,l; double tau; const char* name; };
  for(const Level L: {Level{2,1,1.600,"H 2p"},Level{3,1,5.58,"H 3p"},Level{3,2,15.6,"H 3d"},Level{4,1,12.3,"H 4p"},Level{4,3,72.7,"H 4f"}}){
    const double e=std::sqrt(1-std::pow((L.l+0.5)/L.n,2)); const Spectrum S=modes(e);
    std::map<std::pair<int,int>,double> ch; double tot=0;
    for(int s=1;s<(int)S.Pp.size();++s) for(int m: {+1,-1}){ const int n2=L.n-s, l2=L.l-m;
      if(n2>=1&&l2>=0&&l2<=n2-1){ const double w=(m>0?S.Pp[s]:S.Pm[s])/s/S.P; ch[{n2,l2}]+=w; tot+=w; } }
    std::printf("%s (n=%d l=%d, e=%.4f): current rate f = %.4f, one-mode rate = %.4f (P/hbar w); tau_one-mode / tau_current = %.4f; branches:",
      L.name,L.n,L.l,e,f(e),tot,f(e)/tot);
    for(auto& [k,v]: ch) std::printf(" ->n%d l%d %.3f",k.first,k.second,v/tot);
    std::printf("\n"); }
}
// ---- part 2: absolute hydrogen lifetimes, three schemes ----
// current (model, audits 352/355/361/370): rate f(e) P/(hbar w); for |L| < hbar (s states) Delta l = +1 (import 3).
// modes-blocked: photon = mode (s,m), Delta n = s, Delta l = -m (L' = L - m hbar), modes off the grid absent.
// modes-floor:   same, but n' = max(n - s, l' + 1): a harmonic too large for the grid lands on the lowest rung that
//                still holds L' (the action rule's floor, audit 355, with the classical L <= J); l' < 0 absent.
// P = e^2 <|a|^2>/(6 pi eps0 c^3), <|a|^2> = (k/mu)^2 (1+e^2/2)/(a^4 (1-e^2)^{5/2}), hbar w = hbar sqrt(k/(mu a^3)), a = n^2 a_mu.

static int part2(){
  const double hbar=1.054571817e-34, eps0=8.8541878128e-12, cl=299792458.0, qe=1.602176634e-19, me=9.1093837015e-31, mp=1.67262192369e-27;
  const double mu=me*mp/(me+mp), k=qe*qe/(4*M_PI*eps0), a0=hbar*hbar/(mu*k);
  struct Level { int n,l; double tau; const char* note; };
  const Level levels[]={{2,0,0,"2s"},{2,1,1.600,"2p meas"},{3,0,158.4,"3s QM"},{3,1,5.58,"3p meas"},{3,2,15.6,"3d QM"},
                        {4,0,226.,"4s QM"},{4,1,12.3,"4p QM"},{4,2,36.5,"4d QM"},{4,3,72.7,"4f QM"}};
  std::printf("\nlevel   tau_ref[ns]  current     modes-blocked  modes-floor   branches (modes-floor)\n");
  for(const Level& L: levels){
    const double J=L.n, Lh=L.l+0.5, e=std::sqrt(1-(Lh/J)*(Lh/J)), a=L.n*L.n*a0;
    const double acc2=(k/mu)*(k/mu)*(1+e*e/2)/(std::pow(a,4)*std::pow(1-e*e,2.5));
    const double P=qe*qe*acc2/(6*M_PI*eps0*cl*cl*cl), hw=hbar*std::sqrt(k/(mu*a*a*a));
    const Spectrum S=modes(e);
    // current: f P/hw, s-states keep the same f (import 3 only changes the sign of Delta l)
    const double rCur=f(e)*P/hw;
    double rBlk=0,rFlr=0; std::map<std::pair<int,int>,double> br;
    for(int s=1;s<(int)S.Pp.size();++s) for(int m: {+1,-1}){
      const double w=(m>0?S.Pp[s]:S.Pm[s])/S.P*P/(s*hw); const int l2=L.l-m;
      { const int n2=L.n-s; if(n2>=1&&l2>=0&&l2<=n2-1) rBlk+=w; }
      { const int n2=std::max(L.n-s,l2+1); if(l2>=0&&n2<L.n){ rFlr+=w; br[{n2,l2}]+=w; } } }
    const auto tau=[](double r){ return r>0?1e9/r:INFINITY; };
    std::printf("%-6s %8.3f   %10.4g   %10.4g   %10.4g   ",L.note,L.tau,tau(rCur),tau(rBlk),tau(rFlr));
    for(auto& [kk,v]: br) std::printf(" ->%d%c %.3f",kk.first,"spdfg"[kk.second],v/rFlr);
    std::printf("\n"); }
  return 0; }
int main(){
  std::printf("e      share of photons in m=-1 modes  [sum P-/s / sum (P+ + P-)/s]   (sum P+/s)/(f P)\n");
  for(double e: {0.1,0.3,0.5,0.6614,0.75,0.8660,0.9270,0.968,0.986}){ const Spectrum S=modes(e,8192,200); double sp=0,sm=0;
    for(int s=1;s<(int)S.Pp.size();++s){ sp+=S.Pp[s]/s; sm+=S.Pm[s]/s; }
    std::printf("%.4f   %.5f                                              %.5f\n",e,sm/(sp+sm),sp/(f(e)*S.P)); }
}
