#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 380b: the e-dependence of the per-orbit dilation slip is the next order of 1 - sqrt(1 - v^2/c^2).
// slip/(pi J/2hbar) - 1 vs  delta2 = <v_i^4>/(4 c^2 <v_i^2>) = (v_c^2/(16 c^2)) (4/sqrt(1-e^2) - 3),  v_c^2 = K/(mu a),
// and delta3 = 5<v_i^6>/(64 c^4 <v_i^2>) with <v^6> from the same orbit quadrature.
int main(){
  const double mu=firstMass*secondMass/(firstMass+secondMass), K=pairCoulombStrength, a1=pairBohrRadius(activePair);
  const double m1=firstMass, Om=m1*c*c/hbar;
  std::printf(" n  e      measured-1      delta2 (v^4)    residual        delta3 (v^6, numeric)  residual-delta3  4J/L-3\n");
  for(int n=1;n<=4;++n) for(double e: {0.0,0.3,0.6,0.866,0.95}){
    const double a=n*n*a1, w=std::sqrt(K/(mu*a*a*a)), T=2*pi/w, J=std::sqrt(mu*K*a), L=J*std::sqrt(1-e*e);
    const int N=200000; long double slip=0, s2=0, s4=0, s6=0;
    for(int i=0;i<N;++i){
      const double E=2*pi*(i+0.5)/N, tw=1-e*std::cos(E), dt=tw*(T/N), r=a*tw;
      const long double v2=(long double)(K/mu*(2/r-1/a)), vi2=v2/4.0L, x=vi2/((long double)c*c);
      slip+=Om*(-std::expm1(0.5L*std::log1p(-x)))*dt;            // 1 - sqrt(1-x), cancellation-free
      s2+=vi2*dt; s4+=vi2*vi2*dt; s6+=vi2*vi2*vi2*dt;
    }
    const double measured=(double)(slip/(pi*J/(2*hbar)))-1.0;
    const double vc2=K/(mu*a), d2=vc2/(16*c*c)*(4/std::sqrt(1-e*e)-3);
    const double d2num=(double)(s4/(4.0L*c*c*s2)), d3=(double)(5.0L*s6/(64.0L*(long double)c*c*c*c*s2));
    std::printf("%2d %.3f  %.6e  %.6e  %+.3e      %.3e              %+.2e       %.4f   (d2 quadrature/analytic %.10f)\n",
      n,e,measured,d2,measured-d2,d3,measured-d2-d3,4*J/L-3,d2num/d2);
  }
}
