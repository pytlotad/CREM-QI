// Audit 433: doubly retarded echo e- -> e+ -> e- on a circular Ps orbit; smoothness in n_eff; ratio to radiation reaction.
#include <cmath>
#include <cstdio>
#include <vector>
struct V{double x,y,z;};
V operator+(V a,V b){return {a.x+b.x,a.y+b.y,a.z+b.z};} V operator-(V a,V b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
V operator*(V a,double s){return {a.x*s,a.y*s,a.z*s};} double dot(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;}
V cross(V a,V b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};} double norm(V a){return std::sqrt(dot(a,a));}
const double c=299792458.0, e=1.602176634e-19, m=9.1093837015e-31, eps0=8.8541878128e-12, hbar=1.054571817e-34;
const double K=e*e/(4*M_PI*eps0), mu=m/2, alpha=K/(hbar*c), aPs=hbar*hbar/(mu*K);
// Lienard-Wiechert field of charge q at observation point x0, source state (xs, v, a) at its retarded time
void LW(double q,V x0,V xs,V v,V a,V& E,V& B,bool accelOnly){
  V R=x0-xs; double Rn=norm(R); V n=R*(1/Rn); V b=v*(1/c), bd=a*(1/c); double g2=1/(1-dot(b,b)), k=1-dot(n,b);
  V Ev=(n-b)*(1/(g2*k*k*k*Rn*Rn)); V Ea=cross(n,cross(n-b,bd))*(1/(c*k*k*k*Rn));
  E=(accelOnly?Ea:Ev+Ea)*(q/(4*M_PI*eps0)); B=cross(n,E)*(1/c); }
int main(){
  std::printf("# n_eff  delay*c/(2r)  phase_lag/(2alpha/n)  torque_echo[N m]  torque_RR[N m]  ratio\n");
  std::vector<double> ns, tq;
  for(double nf=0.8; nf<=6.0001; nf+=0.005){
    const double r=aPs*nf*nf, w=std::sqrt(K/(mu*r*r*r)), R=r/2;
    auto xe=[&](double t){ return V{R*std::cos(w*t),R*std::sin(w*t),0}; };
    auto ve=[&](double t){ return V{-R*w*std::sin(w*t),R*w*std::cos(w*t),0}; };
    auto ae=[&](double t){ return xe(t)*(-w*w); };
    auto xp=[&](double t){ return xe(t)*(-1.0); }; auto vp=[&](double t){ return ve(t)*(-1.0); };
    const double t3=0; double t2=-r/c; for(int i=0;i<200;++i) t2=t3-norm(xe(t3)-xp(t2))/c;
    double t1=t2-r/c; for(int i=0;i<200;++i) t1=t2-norm(xp(t2)-xe(t1))/c;
    V E1,B1; LW(-e,xp(t2),xe(t1),ve(t1),ae(t1),E1,B1,false);            // electron's field at the positron
    const V ap=(E1+cross(vp(t2),B1))*(e/m);                              // positron acceleration driven by it ("reflection")
    V E2,B2; LW(+e,xe(t3),xp(t2),vp(t2),ap,E2,B2,true);                  // echo: acceleration field of e+ back at e-
    const V F=(E2+cross(ve(t3),B2))*(-e);
    const double torqueEcho=2*cross(xe(t3),F).z;                         // both particles, by symmetry
    // radiation reaction (Abraham-Lorentz, each particle): F = (2/3) K/c^3 * da/dt, da/dt = -w^2 v
    const V jerk=ve(t3)*(-w*w); const V Frr=jerk*((2.0/3.0)*K/(c*c*c));
    const double torqueRR=2*cross(xe(t3),Frr).z;
    const double delay=t3-t1;
    std::printf("%.3f %.9f %.9f %+.6e %+.6e %+.6f\n",nf,delay*c/(2*r),w*delay/(2*alpha/nf),torqueEcho,torqueRR,torqueEcho/torqueRR);
    ns.push_back(nf); tq.push_back(torqueEcho); }
  // smoothness: fit log|tq| = a + b log n + c (log n)^2, then period-1 residual amplitude
  const int N=ns.size(); double S[5][5]={{0}}, Y[3]={0};
  for(int i=0;i<N;++i){ const double x=std::log(ns[i]), y=std::log(std::abs(tq[i])); double p[3]={1,x,x*x};
    for(int a=0;a<3;++a){ Y[a]+=p[a]*y; for(int b=0;b<3;++b) S[a][b]+=p[a]*p[b]; } }
  // solve 3x3
  double A[3][4]; for(int a=0;a<3;++a){ for(int b=0;b<3;++b) A[a][b]=S[a][b]; A[a][3]=Y[a]; }
  for(int i=0;i<3;++i){ for(int j=i+1;j<3;++j){ const double f=A[j][i]/A[i][i]; for(int k=i;k<4;++k) A[j][k]-=f*A[i][k]; } }
  double co[3]; for(int i=2;i>=0;--i){ double s=A[i][3]; for(int k=i+1;k<3;++k) s-=A[i][k]*co[k]; co[i]=s/A[i][i]; }
  double cs=0,sn=0,res2=0; for(int i=0;i<N;++i){ const double x=std::log(ns[i]); const double res=std::log(std::abs(tq[i]))-(co[0]+co[1]*x+co[2]*x*x);
    cs+=res*std::cos(2*M_PI*ns[i]); sn+=res*std::sin(2*M_PI*ns[i]); res2+=res*res; }
  std::printf("# power-law fit: log|torque| = %.6f + %.6f log n + %.6f (log n)^2; rms residual %.3e; period-1 amplitude %.3e\n",
    co[0],co[1],co[2],std::sqrt(res2/N),2*std::sqrt(cs*cs+sn*sn)/N);
}
