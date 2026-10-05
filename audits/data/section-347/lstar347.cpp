#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <vector>
// Audit 347: contact <n>_orb at n = 1 versus L, the single L* that reproduces each MEASURED lifetime,
// its tolerance band, and L distributions (microcanonical; uniform in L) as rate / mean-lifetime averages.
static double a1, psi0, sv, epsOP;
static double contact(double L){ // L in hbar, orbit n = 1
  const double e=std::sqrt(std::max(0.0,1.0-L*L));
  return orbitAveragedContactDensity(a1,std::min(e,1.0-1e-12))/psi0; }
static double tauOf(double nrel,int ch){ return 1.0/(sv*nrel*psi0*(ch==1?1.0:epsOP)); }
static double solve(double target,int ch){ // tau increasing in L
  double lo=1e-6,hi=0.999; for(int k=0;k<200;++k){ double m=0.5*(lo+hi); (tauOf(contact(m),ch)<target?lo:hi)=m; } return 0.5*(lo+hi); }
int main(){
  a1=pairBohrRadius(activePair); psi0=1.0/(pi*a1*a1*a1);
  sv=4.0*pi*classicalElectronRadius*classicalElectronRadius*c; epsOP=orePowellSuppression();
  const double tp=125.14e-12, to=142.04e-9;
  std::printf("a1=%.6e m eps=%.6e m  eps/a1=%.4e  sigma v |psi0|^2 -> %.4f ps   1/eps_OP=%.4f\n",
    a1,magneticDipoleRadius(),magneticDipoleRadius()/a1,1e12/(sv*psi0),1.0/epsOP);
  std::printf("\nTABLE n(L) at n = 1\n");
  for(double L: {0.0,1e-3,0.01,0.03,0.05,0.07,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.4,0.5,0.7,0.9,0.99})
    std::printf("  L=%.3f e=%.10f rp/eps=%.4e n/psi0=%.5e tau_p=%.4e s tau_o=%.4e s\n",L,std::sqrt(1-L*L),
      a1*(1-std::sqrt(1-L*L))/magneticDipoleRadius(),contact(L),tauOf(contact(L),1),tauOf(contact(L),2));
  for(int ch=1;ch<=2;++ch){ const double T=ch==1?tp:to; const char*nm=ch==1?"p-Ps":"o-Ps";
    const double Ls=solve(T,ch), lo=solve(0.9*T,ch), hi=solve(1.1*T,ch), lo3=solve(0.999*T,ch), hi3=solve(1.001*T,ch);
    std::printf("\n%s: target tau %.5e s -> n/psi0 = %.5f ; L* = %.6f hbar (e=%.8f, rp/eps=%.4f)\n",nm,T,
      1.0/(sv*psi0*(ch==1?1.0:epsOP)*T),Ls,std::sqrt(1-Ls*Ls),a1*(1-std::sqrt(1-Ls*Ls))/magneticDipoleRadius());
    std::printf("   band tau +-10%%: L in [%.6f, %.6f] width %.2e ; +-0.1%%: [%.7f, %.7f] width %.2e\n",hi,lo,lo-hi,hi3,lo3,lo3-hi3); }
  // distributions: average over L^2 uniform (microcanonical) and L uniform, on [0, Lmax]
  const int N=200000;
  std::vector<double> nn(N+1), LL(N+1);
  for(int k=0;k<=N;++k){ const double u=(k+0.5)/(N+1); LL[k]=u; nn[k]=contact(u); }
  auto avg=[&](double Lmax,bool squared,double& rate,double& invRate){
    double sw=0,sn=0,si=0; for(int k=0;k<=N;++k){ if(LL[k]>Lmax) break; const double w=squared?LL[k]:1.0;
      sw+=w; sn+=w*nn[k]; si+=w/nn[k]; } rate=sn/sw; invRate=si/sw; };
  std::printf("\nDISTRIBUTIONS on [0, Lmax]:  <n>/psi0 (rate average = initial slope)  and  <1/n>*psi0 (mean lifetime)\n");
  for(bool sq: {true,false}) for(double Lm: {1.0,0.5,0.3,0.2,0.15,0.1}){
    double r,ir; avg(Lm,sq,r,ir);
    std::printf("  %s Lmax=%.2f: <n>=%.4e  tau_p(rate)=%.4e s  tau_o(rate)=%.4e s | <1/n>^-1=%.4e  E[T]_p=%.4e s E[T]_o=%.4e s  E[T]/tau_rate=%.3e\n",
      sq?"L^2 unif (microcan.)":"L uniform           ",Lm,r,tauOf(r,1),tauOf(r,2),1.0/ir,tauOf(1.0/ir,1),tauOf(1.0/ir,2),ir*r); }
  std::printf("\nCLASSICAL MOMENTS at n = 1: <r> = a(1+e^2/2) = 1.5 a (QM 1s) only at e = 1, L = 0; <r^2> = a^2(1+3e^2/2) <= 2.5 a^2 (QM 3 a^2)\n");
}
