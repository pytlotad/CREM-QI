// Audit 434: Compton clock + reflection.  Units c = lambdabar_C = 1, time 1/omega_C.
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>
const double alpha=7.2973525693e-3;
typedef std::complex<double> C;
int main(){ setvbuf(stdout,nullptr,_IONBF,0);
  // ---- Part A: 1D kinematics.  Source x_s(t) = beta t, proper phase tau = t/gamma.  Mirror at rest at x = D (ahead).
  std::printf("# Part A: envelope period behind the moving clock (backward emission + reflection from a mirror ahead)\n");
  for(double beta: {0.1,0.3,0.5}){
    const double g=1/std::sqrt(1-beta*beta), D=4000, t0=2000;   // observe at lab time t0, region x < x_s(t0)
    // backward wave at x < source: emitted at te with t0 - te = x_s(te) - x  ->  te = (t0 + x)/(1 + beta)... (x_s = beta te)
    // reflected wave: emitted at te forward, reaches mirror at D, returns to x: (D - beta te) + (D - x) = t0 - te -> te = (t0 + x - 2D)/(1 - beta)
    const int N=200000; const double x0=-6000, x1=0; std::vector<double> I(N);
    for(int i=0;i<N;++i){ const double x=x0+(x1-x0)*i/(N-1);
      const double te1=(t0+x)/(1+beta), te2=(t0+x-2*D)/(1-beta);
      const C psi=std::exp(C(0,te1/g))+std::exp(C(0,te2/g+M_PI)); I[i]=std::norm(psi); }
    // period of |psi|^2 by counting maxima
    int maxima=0; double first=0,last=0; for(int i=1;i<N-1;++i) if(I[i]>I[i-1]&&I[i]>=I[i+1]&&I[i]>3.0){ const double x=x0+(x1-x0)*i/(N-1); if(!maxima) first=x; last=x; ++maxima; }
    const double period=(last-first)/(maxima-1);
    std::printf("beta %.2f: envelope period %.5f, 2pi/(gamma beta) %.5f, ratio %.6f\n",beta,period,2*M_PI/(g*beta),period/(2*M_PI/(g*beta)));
  }
  // ---- Part B: Ps, electron clock, echo off the positron, exact double retardation.
  auto orbit=[&](double a,double e,double t,double& x,double& y,double& vx,double& vy){ // relative orbit, mu k = 2 alpha
    const double w=std::sqrt(2*alpha/(a*a*a)); double M=w*t, E=M; for(int i=0;i<12;++i) E-= (E-e*std::sin(E)-M)/(1-e*std::cos(E));
    const double b=a*std::sqrt(1-e*e), Ed=w/(1-e*std::cos(E));
    x=a*(std::cos(E)-e); y=b*std::sin(E); vx=-a*std::sin(E)*Ed; vy=b*std::cos(E)*Ed; };
  auto electron=[&](double a,double e,double t,double& x,double& y,double& vx,double& vy){ orbit(a,e,t,x,y,vx,vy); x*=0.5; y*=0.5; vx*=0.5; vy*=0.5; };
  auto echo=[&](double a,double e,double t3,double& PhiC,double& PhiK,double& dB){
    double x3,y3,u,w; electron(a,e,t3,x3,y3,u,w);
    double t2=t3-a, t1=t2-a;
    for(int i=0;i<20;++i){ double xp,yp,vx,vy; electron(a,e,t2,xp,yp,vx,vy); xp=-xp; yp=-yp; t2=t3-std::hypot(x3-xp,y3-yp); }
    double xp,yp,vx,vy; electron(a,e,t2,xp,yp,vx,vy); xp=-xp; yp=-yp;
    for(int i=0;i<20;++i){ double x1,y1,u1,w1; electron(a,e,t1,x1,y1,u1,w1); t1=t2-std::hypot(xp-x1,yp-y1); }
    PhiC=t3-t1; const int S=32; double s=0, px=0; double xa,ya,ua,wa; electron(a,e,t1,xa,ya,ua,wa);
    for(int i=0;i<=S;++i){ const double t=t1+(t3-t1)*i/S; double xx,yy,uu,ww; electron(a,e,t,xx,yy,uu,ww); const double wt=(i==0||i==S)?1:(i%2?4:2); s+=wt*(-0.5*(uu*uu+ww*ww)); }
    PhiK=s*(t3-t1)/(3*S);
    double xb,yb,ub,wb; electron(a,e,t3,xb,yb,ub,wb); dB=0.5*((ua+ub)*(xb-xa)+(wa+wb)*(yb-ya)); }; // p.Dx/hbar with m = 1 (lambdabar units)
  std::printf("# Part B1: circular orbit, echo phases\n");
  for(double n: {1.0,2.0,3.0}){ const double a=n*n*2/alpha; double PC,PK,dB; echo(a,0.0,0.0,PC,PK,dB);
    std::printf("n %.0f: Phi_C %.6f (4n^2/alpha %.6f, x%.6f)  Phi_K %.6e (-alpha/2 = %.6e, x%.5f)  p.Dx/hbar %.6e (alpha x%.5f)\n",n,PC,4*n*n/alpha,PC/(4*n*n/alpha),PK,-alpha/2,PK/(-alpha/2),dB,dB/alpha); }
  std::printf("# Part B2: orbit-averaged <cos Delta phi> vs n_eff\n");
  for(double e: {0.0,0.5}){
    std::vector<double> ns, cs; const int P=e>0?16:1;
    for(double n=0.9;n<=3.1+1e-9;n+=0.0001){ const double a=n*n*2/alpha, T=2*M_PI*std::sqrt(a*a*a/(2*alpha)); double sum=0;
      for(int j=0;j<P;++j){ double PC,PK,dB; echo(a,e,T*j/P,PC,PK,dB); sum+=std::cos(PC+PK+M_PI); } ns.push_back(n); cs.push_back(sum/P); }
    for(double nc: {1.0,2.0,3.0}){ int z=0; for(size_t i=1;i<ns.size();++i) if(std::abs(ns[i]-nc)<0.05&&std::abs(ns[i-1]-nc)<0.05&&cs[i]*cs[i-1]<0) ++z;
      std::printf("e %.1f, n %.0f: zero crossings per unit n %.1f (prediction 2*4n/(pi alpha)... crossings = 2 per period: %.1f)\n",e,nc,z/0.1,2*4*nc/(M_PI*alpha)); }
    double c1=0,s1=0,rms=0; for(size_t i=0;i<ns.size();++i){ c1+=cs[i]*std::cos(2*M_PI*ns[i]); s1+=cs[i]*std::sin(2*M_PI*ns[i]); rms+=cs[i]*cs[i]; }
    std::printf("e %.1f: period-1 amplitude %.3e, rms %.3e, ratio %.3e\n",e,2*std::hypot(c1,s1)/ns.size(),std::sqrt(rms/ns.size()),2*std::hypot(c1,s1)/ns.size()/std::sqrt(2*rms/ns.size()));
  }
  std::printf("# H434d: Thomson reflection ~ r_e/r: r_e = alpha (lambdabar units); n=1: %.3e, n=2: %.3e\n",alpha/(2/alpha),alpha/(8/alpha));
}
