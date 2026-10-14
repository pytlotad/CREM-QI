// Audit 438 calibration: on a FROZEN Kepler orbit (J = 2, e = 0.6, hbar = mu = k = 1, theta' = v^2, theta(peri) = psi)
//  orbit drive  D(psi) = <cos theta dg/dt> J^3   (dJ/dt per eps from U = -eps (g - gbar) cos theta)
//  clock pull   P(psi) = <(g - gbar) sin theta>  (pull term: theta' = v^2 - kappa eps (g - gbar) sin theta)
#include <cmath>
#include <cstdio>
int main(){ const double J=2.0, e=0.6, r0=2.0, a=J*J, w=1/(J*J*J), T=2*M_PI/w; const int N=8000, orbits=4;
  auto rOf=[&](double t){ const double M=w*t; double E=M; for(int k=0;k<20;++k) E-=(E-e*std::sin(E)-M)/(1-e*std::cos(E)); return a*(1-e*std::cos(E)); };
  double gbar=0; for(int i=0;i<N;++i) gbar+=std::exp(-rOf((i+0.5)*T/N)/r0); gbar/=N;
  std::printf("# gbar = %.6e\n# psi  D(psi)  P(psi)\n",gbar);
  double Dc=0,Ds=0,Pc=0,Ps=0; const int K=16;
  for(int j=0;j<K;++j){ const double psi=2*M_PI*j/K; double D=0, P=0, gPrev=std::exp(-rOf(0)/r0); const double dt=T/N;
    for(int i=1;i<=orbits*N;++i){ const double t=i*dt, r=rOf(t), g=std::exp(-r/r0); const double Ecc=std::acos(std::max(-1.0,std::min(1.0,(1-r/a)/e))); const double Mm=std::fmod(w*t,2*M_PI); const double Ea=(Mm<=M_PI)?Ecc:2*M_PI-Ecc; const double th=psi+J*(Ea+e*std::sin(Ea))+2*M_PI*J*std::floor(w*t/(2*M_PI)); D+=std::cos(th)*(g-gPrev); P+=(g-gbar)*std::sin(th)*dt; gPrev=g; }
    D=D/(orbits*T)*J*J*J; P/=orbits*T; std::printf("%.4f %+.6e %+.6e\n",psi,D,P);
    Dc+=D*std::cos(psi); Ds+=D*std::sin(psi); Pc+=P*std::cos(psi); Ps+=P*std::sin(psi); }
  const double DA=2*std::hypot(Dc,Ds)/K, PA=2*std::hypot(Pc,Ps)/K;
  std::printf("# D = %.4e sin(psi - %.4f)   P = %.4e sin(psi - %.4f)\n",DA,std::atan2(-Dc,Ds),PA,std::atan2(-Pc,Ps)); }
