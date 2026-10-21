#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
// Audit 542: engine transport of Ps 2^3P_J (a = 4 a_B, |L| = 1.5 hbar, parallel spins, |J_e| = (J + 1/2) hbar) in a field, isotropic ensembles (no basin);
// singlet fraction at the 2P decay P_S = Gamma int e^{-Gamma t} w dt, w = contactTwoPhotonWeight, tau = 3.348 ns; compare with m542 variant C.
int main(int argc,char** argv){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true; gSingletSpinTransport=false;
  const double a=4*pairBohrRadius(activePair), mu=firstMass*secondMass/(firstMass+secondMass), tau=3.348e-9, G=1/tau;
  const int N=argc>1?std::atoi(argv[1]):40; const double B=argc>2?std::atof(argv[2]):0.01, dt=argc>3?std::atof(argv[3]):2e-12;
  gExternalMagneticField=Vec3{0,0,1}*B; const long steps=(long)std::ceil(5*tau/dt); double res[3],se[3];
  for(int J=0;J<3;++J){ std::mt19937_64 rng(542+J); std::normal_distribution<double> nd; double sum=0,sum2=0;
    const double ct=((J+0.5)*(J+0.5)-2.25-1)/3, st=std::sqrt(std::max(0.0,1-ct*ct));
    for(int k=0;k<N;++k){ double q[4],qn=0; for(double& v:q){v=nd(rng); qn+=v*v;} qn=std::sqrt(qn); for(double& v:q) v/=qn;
      const Vec3 u{q[1],q[2],q[3]}; auto R=[&](Vec3 v){ const Vec3 t=cross(u,v)*2.0; return v+t*q[0]+cross(u,t); };
      const Vec3 n=R(Vec3{st,0,ct}), lh=R(Vec3{0,0,1}), ph=R(Vec3{1,0,0});
      SecularSpinOrbitState s; s.orbitalAngularMomentum=lh*(1.5*hbar); s.periapsisDirection=ph; s.firstDipole=n*(-firstMagneticMoment); s.secondDipole=n*secondMagneticMoment;
      double acc=0, w=0;
      for(long i=0;i<steps;++i){ s=advanceCoupledSecularSpinOrbit(s,a,mu,dt).state; w=contactTwoPhotonWeight(s.firstDipole,s.secondDipole); acc+=G*std::exp(-G*(i+1)*dt)*dt*w; }
      acc+=std::exp(-5.0)*w; sum+=acc; sum2+=acc*acc; }
    res[J]=sum/N; se[J]=std::sqrt(std::max(0.0,sum2/N-res[J]*res[J])/N); }
  std::printf("B %5.2f T engine (dt %.1e, N %d): P_S J0 %.5e(%.1e) J1 %.5e(%.1e) J2 %.5e(%.1e) 3P %.5e\n",B,dt,N,res[0],se[0],res[1],se[1],res[2],se[2],(res[0]+3*res[1]+5*res[2])/9);
}
