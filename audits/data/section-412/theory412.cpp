// Audit 412: theory of the rung positions of audit 411's path-memory wave.
// Circular orbit, units hbar = mu = k = 1.  Loop k back returns with phase -2 pi kappa k J and weight q^k
// (q = exp(-1/tauOrbits)); the samples of one loop near the particle add up to N_eff = sqrt(2pi) sigma/Ds with
// sigma = lambda/4 = pi/(2 kappa p), Ds = 2 pi r/256; the tangential force times r gives
//     dJ/dt = -gamma J + A F(theta),  theta = 2 pi kappa J,  A = 64 sqrt(2 pi) eps,
//     F(theta) = sum_k q^k sin(k theta) = q sin(theta)/(1 - 2 q cos(theta) + q^2)
// (kappa, p and r cancel).  Stable fixed points: A F(theta) = gamma J with d/dJ[-gamma J + A F] < 0.
// argv: eps kappa [gamma=2e-4] [tauOrbits=3] [Jmax=6]  -> stable J per rung m (kappa J = m + delta).
#include <cmath>
#include <cstdio>
#include <cstdlib>
int main(int argc,char**argv){
  const double eps=std::atof(argv[1]), kappa=std::atof(argv[2]);
  const double gamma=argc>3?std::atof(argv[3]):2e-4, tauO=argc>4?std::atof(argv[4]):3, Jmax=argc>5?std::atof(argv[5]):6;
  const double q=std::exp(-1.0/tauO), A=64.0*std::sqrt(2*M_PI)*eps;
  const auto F=[&](double th){ return q*std::sin(th)/(1-2*q*std::cos(th)+q*q); };
  const auto g=[&](double J){ return -gamma*J+A*F(2*M_PI*kappa*J); };      // dJ/dt
  std::printf("eps=%g kappa=%g gamma=%g tau=%g orbits: A=%.4e, q=%.4f, F_max=%.4f at theta=%.4f\n",eps,kappa,gamma,tauO,A,q,
    q*std::sin(std::acos(2*q/(1+q*q)))/(1-2*q*2*q/(1+q*q)+q*q),std::acos(2*q/(1+q*q)));
  const int N=600000; double prev=g(0.2);
  for(int i=1;i<=N;++i){ const double J=0.2+(Jmax-0.2)*i/N; const double cur=g(J);
    if(prev>0&&cur<=0){ // + to - crossing: stable fixed point
      double lo=J-(Jmax-0.2)/N, hi=J; for(int it=0;it<100;++it){ const double mid=0.5*(lo+hi); if(g(mid)>0) lo=mid; else hi=mid; }
      const double Js=0.5*(lo+hi), kJ=kappa*Js; std::printf("  stable J = %.5f  (kappa J = %.5f, rung m = %ld, delta = %+.5f)\n",Js,kJ,std::lround(kJ-0.0),kJ-std::floor(kJ+0.5)); }
    prev=cur; }
}
