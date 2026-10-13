// Audit 435: multiple reflections of the Compton-clock echo between e- and e+ over one orbit.  c = lambdabar_C = 1.
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include <utility>
const double alpha=7.2973525693e-3;
void orbit(double a,double e,double t,double& x,double& y,double& vx,double& vy){
  const double w=std::sqrt(2*alpha/(a*a*a)); double M=w*t, E=M; for(int i=0;i<12;++i) E-=(E-e*std::sin(E)-M)/(1-e*std::cos(E));
  const double b=a*std::sqrt(1-e*e), Ed=w/(1-e*std::cos(E)); x=a*(std::cos(E)-e); y=b*std::sin(E); vx=-a*std::sin(E)*Ed; vy=b*std::cos(E)*Ed; }
void particle(int s,double a,double e,double t,double& x,double& y,double& vx,double& vy){ orbit(a,e,t,x,y,vx,vy); const double f=0.5*s; x*=f; y*=f; vx*=f; vy*=f; }
double tau(int s,double a,double e,double t0,double t1){ const int S=64; double sum=0; for(int i=0;i<=S;++i){ const double t=t0+(t1-t0)*i/S; double x,y,u,w; particle(s,a,e,t,x,y,u,w); sum+=((i==0||i==S)?1:(i%2?4:2))*std::sqrt(1-u*u-w*w); } return sum*(t1-t0)/(3*S); }
int main(){ setvbuf(stdout,nullptr,_IONBF,0);
  for(auto [n,e]: {std::pair{1.0,0.0},{2.0,0.0},{3.0,0.0},{2.0,0.5}}){
    const double a=n*n*2/alpha, T=2*M_PI*std::sqrt(a*a*a/(2*alpha));
    double t=0; int s=+1; int legs=0; double dB[2]={0,0}, propTau[2]={0,0}, lenSum=0;
    double xs,ys,us,ws; particle(s,a,e,t,xs,ys,us,ws);
    while(true){ // leg from particle s at time t to particle -s at arrival time ta
      double ta=t+a; for(int i=0;i<30;++i){ double x,y,u,w; particle(-s,a,e,ta,x,y,u,w); ta=t+std::hypot(x-xs,y-ys); }
      double xr,yr,ur,wr; particle(-s,a,e,ta,xr,yr,ur,wr);
      // de Broglie phase p.Dx of EACH particle over the leg duration [t, ta] (m = 1, trapezoid): index 0 = e- (s = +1), 1 = e+
      for(int q: {+1,-1}){ double xa,ya,ua,wa,xb,yb,ub,wb; particle(q,a,e,t,xa,ya,ua,wa); particle(q,a,e,ta,xb,yb,ub,wb);
        dB[q>0?0:1]+=0.5*((ua+ub)*(xb-xa)+(wa+wb)*(yb-ya)); }
      lenSum+=ta-t; ++legs; t=ta; s=-s; xs=xr; ys=yr; us=ur; ws=wr;
      if(t>=T && s==+1) break; }
    // receivers alternate; sum over all legs of p.Dx for each particle during the legs it receives = its p.Dx over the time it is "waiting"; compare with full-orbit p.Dx
    const double phaseClock=tau(+1,a,e,0,t);           // electron proper time from first emission to final return
    const double Ntrips=legs/2.0;
    std::printf("n %.0f e %.1f: legs %d (round trips %.2f, pi n/alpha %.2f, x%.5f); elapsed/T %.6f\n",n,e,legs,Ntrips,M_PI*n/alpha,Ntrips/(M_PI*n/alpha),t/T);
    std::printf("   clock - wave phase = tau_e(end) - tau_e(0) + legs*pi = %.4f rad (omega_C T = %.4f); rung spacing in n ~ 2pi/(3 omega_C T/n) = %.3e\n",phaseClock+legs*M_PI,T,2*M_PI/(3*T/n));
    std::printf("   sum p.Dx/hbar over legs: e+ %.6f, e- %.6f, total %.6f (2 pi n = %.6f, x%.5f); each / (pi n): %.5f %.5f\n",dB[1],dB[0],dB[0]+dB[1],2*M_PI*n,(dB[0]+dB[1])/(2*M_PI*n),dB[1]/(M_PI*n),dB[0]/(M_PI*n));
    std::printf("   amplitude after the chain: log10 Thomson (r_e/r)^legs = %.1f, resonant (lambdabar/r)^legs = %.1f\n",legs*std::log10(alpha/a),legs*std::log10(1.0/a));
  }
}
