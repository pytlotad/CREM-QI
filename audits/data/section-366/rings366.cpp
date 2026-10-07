#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
// Audit 366: overlapping current loops of the model (charge ring, R = 2 mu/(e c), charge moving at c).
// Volume integrals over the relative position d inside a ball of radius Dmax:
//   magnetic  Int U_mag d^3d,  U_mag = -(mu0/4pi) I1 I2 sum dl1.dl2 / |d + x2 - x1|   (Neumann, constant currents)
//   electric  Int [U_C - U_C,point] d^3d,  U_C = k q1 q2/e^2 <1/|d + x2 - x1|>.
// Both rings share the spin axis z (coaxial when d || z).  Triplet: e- and e+ circulate the same way -> currents opposite.
int main(int argc,char**argv){
  const double k=pairCoulombStrength, mu=firstMagneticMoment, e=std::sqrt(k*4*pi*epsilon0), R=2*mu/(e*c), hP=2*pi*hbar;
  const int N=96; const long M=argc>1?std::atol(argv[1]):400000; const double Dmax=6*R;
  std::vector<Vec3> x(N), t(N); for(int i=0;i<N;++i){ const double p=2*pi*(i+0.5)/N; x[i]={R*std::cos(p),R*std::sin(p),0}; t[i]={-std::sin(p),std::cos(p),0}; }
  const double dl=2*pi*R/N, I=e*c/(2*pi*R);  // |current|; m = I pi R^2 = e c R / 2 = mu
  std::mt19937_64 rng(366); std::uniform_real_distribution<double> U(-1,1);
  double sumM=0,sumM2=0,sumE=0,sumE2=0; long cnt=0;
  for(long s=0;s<M;++s){ Vec3 d{U(rng)*Dmax,U(rng)*Dmax,U(rng)*Dmax}; if(d.norm()>Dmax) continue; ++cnt;
    double neu=0, inv=0;
    for(int i=0;i<N;++i) for(int j=0;j<N;++j){ const double r=(d+x[j]-x[i]).norm(); const double ir=1.0/std::max(r,1e-3*R); neu+=dot(t[i],t[j])*ir; inv+=ir; }
    neu*=dl*dl; inv/=double(N)*N;
    // triplet: currents opposite -> I1 I2 = -I^2 ; U_mag = -(mu0/4pi) I1 I2 neu
    const double um=-(mu0/(4*pi))*(-I*I)*neu;
    const double ue=-k*(inv-1.0/d.norm());   // opposite charges: -k <1/r> + k/|d|
    sumM+=um; sumM2+=um*um; sumE+=ue; sumE2+=ue*ue; }
  const double vol=4.0/3.0*pi*Dmax*Dmax*Dmax;
  const double IM=sumM/cnt*vol, IE=sumE/cnt*vol;
  const double eM=std::sqrt(sumM2/cnt-std::pow(sumM/cnt,2))/std::sqrt(double(cnt))*vol, eE=std::sqrt(sumE2/cnt-std::pow(sumE/cnt,2))/std::sqrt(double(cnt))*vol;
  const double fermiTriplet=(2*mu0/3)*mu*mu;            // -(2 mu0/3) m1.m2 with m1.m2 = -mu^2
  const double finite=(2*pi/3)*k*2*R*R;
  const double aPs=pairBohrRadius(activePair), nc2=1.0/(pi*aPs*aPs*aPs)/8.0;
  std::printf("R = %.3f fm, samples in ball %ld, ring nodes %d\n",R*1e15,cnt,N);
  std::printf("magnetic (triplet): MC %.4e +- %.1e J m^3   Fermi -(2mu0/3)m1.m2 = %.4e   ratio %.4f\n",IM,eM,fermiTriplet,IM/fermiTriplet);
  std::printf("electric:           MC %.4e +- %.1e J m^3   (2pi/3)k 2R^2   = %.4e   ratio %.4f\n",IE,eE,finite,IE/finite);
  std::printf("at n = 2 (x |psi_1s(0)|^2/8): magnetic %.3f +- %.3f GHz (Fermi %.3f), electric %.3f +- %.3f GHz (365: %.3f); singlet magnetic = -triplet by construction (I1 I2 sign)\n",
    IM*nc2/hP*1e-9,eM*nc2/hP*1e-9,fermiTriplet*nc2/hP*1e-9,IE*nc2/hP*1e-9,eE*nc2/hP*1e-9,finite*nc2/hP*1e-9);
}
