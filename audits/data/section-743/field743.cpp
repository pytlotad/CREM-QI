// Audit 743 (from field446.cpp, audit 446): matter-wave field of the pair in the relative coordinate, source on the particle.
// Without env JR0 the program is field446 unchanged (main446; control H743f).  With JR0, L0 (audit 743 mode, run743):
//   V = -1/r + beta/r^2 (grid: r clipped at 0.25; particle: analytic), start at periapsis of torus (J_r, L); dt = T_r0/1000;
//   source phase PHASE=L: phi_s = int (v^2/2 - V) dt = S0 - E t (gauge, 691/603); PHASE=K: int v^2/2 dt (446);
//   Gamma = GT/T_r0; ONEP=1: source only during the first radial period; particle: SUB Verlet substeps per field step, field force
//   eps Re[exp(-i (phi_s + chi)) grad psi(x_p)] held over the step.  Output per radial period (periapsis to periapsis):
//   k t J_r L Re<Z> Im<Z> <||psi||^2> |psi|max, Z = exp(-i phi_s) psi(x_p).
// argv (743 mode): periods eps chi A [N=128 Lbox=24]; env BETA (0.2), GT (0.1), GABS (absolute Gamma, overrides GT), PHASE (L),
//   SUB (20), ONEP, EVERY (1).
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
typedef std::complex<double> C;
static void fft1(std::vector<C>& a,bool inv){ const int n=a.size();
  for(int i=1,j=0;i<n;++i){ int bit=n>>1; for(;j&bit;bit>>=1) j^=bit; j^=bit; if(i<j) std::swap(a[i],a[j]); }
  for(int len=2;len<=n;len<<=1){ const double ang=2*M_PI/len*(inv?1:-1); const C wl(std::cos(ang),std::sin(ang));
    for(int i=0;i<n;i+=len){ C w(1); for(int j=0;j<len/2;++j){ const C u=a[i+j], v=a[i+j+len/2]*w; a[i+j]=u+v; a[i+j+len/2]=u-v; w*=wl; } } }
  if(inv) for(auto& x:a) x/=n; }
static void fft2(std::vector<C>& f,int N,bool inv){ std::vector<C> row(N);
  for(int y=0;y<N;++y){ for(int x=0;x<N;++x) row[x]=f[y*N+x]; fft1(row,inv); for(int x=0;x<N;++x) f[y*N+x]=row[x]; }
  for(int x=0;x<N;++x){ for(int y=0;y<N;++y) row[y]=f[y*N+x]; fft1(row,inv); for(int y=0;y<N;++y) f[y*N+x]=row[y]; } }
static int main446(int argc,char** argv){ setvbuf(stdout,nullptr,_IOLBF,0);
  const double J0=std::atof(argv[1]), e0=std::atof(argv[2]); const int orbits=std::atoi(argv[3]);
  const double eps=std::atof(argv[4]), chi=std::atof(argv[5]), gam=std::atof(argv[6]), A=std::atof(argv[7]), Gorb=std::atof(argv[8]);
  const int N=argc>9?std::atoi(argv[9]):128; const double Lb=argc>10?std::atof(argv[10]):24.0, dx=Lb/N, soft=0.1, sigma=0.3;
  const double T0=2*M_PI*J0*J0*J0, Gam=1.0/(Gorb*T0), dt=T0/1000;
  std::vector<C> psi(N*N,0.0); std::vector<double> V(N*N), absb(N*N), k2(N*N);
  for(int iy=0;iy<N;++iy) for(int ix=0;ix<N;++ix){ const double x=-Lb/2+ix*dx, y=-Lb/2+iy*dx, r=std::sqrt(x*x+y*y+soft*soft);
    V[iy*N+ix]=-1/r; const double d=std::min(std::min(x+Lb/2,Lb/2-x),std::min(y+Lb/2,Lb/2-y)), w=2.0;
    absb[iy*N+ix]= d<w? 0.5*std::pow((w-d)/w,2):0.0;
    const int kx=ix<N/2?ix:ix-N, ky=iy<N/2?iy:iy-N; const double qx=2*M_PI*kx/Lb, qy=2*M_PI*ky/Lb; k2[iy*N+ix]=qx*qx+qy*qy; }
  std::vector<C> potStep(N*N), kinStep(N*N);
  for(int i=0;i<N*N;++i){ potStep[i]=std::exp(C(-(Gam+absb[i])*dt*0.5,-V[i]*dt*0.5)); kinStep[i]=std::exp(C(0,-0.5*k2[i]*dt)); }
  const double a=J0*J0, rp=a*(1-e0), vp=std::sqrt((1+e0)/(a*(1-e0)));
  double x=rp,y=0,vx=0,vy=vp, phis=0, t=0, angAcc=0; int count=0; double sJ=0,sL=0; int nb=0;
  auto field=[&](double px,double py,C& f,C& gx,C& gy){ const double fx=(px+Lb/2)/dx, fy=(py+Lb/2)/dx; const int ix=(int)std::floor(fx), iy=(int)std::floor(fy);
    const double ux=fx-ix, uy=fy-iy; auto P=[&](int i,int j){ i=(i%N+N)%N; j=(j%N+N)%N; return psi[j*N+i]; };
    auto bil=[&](auto g){ return (1-ux)*(1-uy)*g(ix,iy)+ux*(1-uy)*g(ix+1,iy)+(1-ux)*uy*g(ix,iy+1)+ux*uy*g(ix+1,iy+1); };
    f=bil([&](int i,int j){return P(i,j);}); gx=bil([&](int i,int j){return (P(i+1,j)-P(i-1,j))/(2*dx);}); gy=bil([&](int i,int j){return (P(i,j+1)-P(i,j-1))/(2*dx);}); };
  auto acc=[&](double px,double py,double qx,double qy,double ph,double& ax,double& ay){ const double r=std::sqrt(px*px+py*py+soft*soft);
    ax=-px/(r*r*r)-gam*qx; ay=-py/(r*r*r)-gam*qy; if(eps!=0){ C f,gx,gy; field(px,py,f,gx,gy); const C ph1=std::exp(C(0,-(ph+chi)));
      ax+=eps*std::real(ph1*gx); ay+=eps*std::real(ph1*gy); } };
  std::printf("# orbit t <J> <L> e |psi|max\n");
  double ax,ay; acc(x,y,vx,vy,phis,ax,ay);
  while(count<orbits){
    // field: half potential, kinetic, half potential (Strang) + source
    for(int i=0;i<N*N;++i) psi[i]*=potStep[i];
    fft2(psi,N,false); for(int i=0;i<N*N;++i) psi[i]*=kinStep[i]; fft2(psi,N,true);
    for(int i=0;i<N*N;++i) psi[i]*=potStep[i];
    if(A!=0){ const C s=C(0,-1)*A*std::exp(C(0,phis))*dt; const int cx=(int)((x+Lb/2)/dx), cy=(int)((y+Lb/2)/dx), R=(int)(3*sigma/dx)+1;
      for(int j=cy-R;j<=cy+R;++j) for(int i=cx-R;i<=cx+R;++i){ if(i<0||j<0||i>=N||j>=N) continue; const double X=-Lb/2+i*dx-x, Y=-Lb/2+j*dx-y;
        psi[j*N+i]+=s*std::exp(-(X*X+Y*Y)/(2*sigma*sigma))/(2*M_PI*sigma*sigma); } }
    // particle: velocity Verlet
    const double vxh=vx+0.5*dt*ax, vyh=vy+0.5*dt*ay; const double x1=x+dt*vxh, y1=y+dt*vyh;
    phis+=0.5*(vx*vx+vy*vy)*dt;
    double ax1,ay1; acc(x1,y1,vxh,vyh,phis,ax1,ay1); const double vx1=vxh+0.5*dt*ax1, vy1=vyh+0.5*dt*ay1;
    double d=std::atan2(y1,x1)-std::atan2(y,x); if(d>M_PI) d-=2*M_PI; if(d<-M_PI) d+=2*M_PI; angAcc+=d;
    x=x1; y=y1; vx=vx1; vy=vy1; ax=ax1; ay=ay1; t+=dt;
    if(x*x+y*y>9.5*9.5){ std::printf("# left grid region at t=%.1f\n",t); break; }
    if(angAcc>=2*M_PI){ angAcc-=2*M_PI; ++count; const double r=std::sqrt(x*x+y*y+soft*soft), E=0.5*(vx*vx+vy*vy)-1/r;
      if(!(E<0)){ std::printf("# unbound\n"); break; } const double J=1/std::sqrt(-2*E), L=x*vy-y*vx; sJ+=J; sL+=L; ++nb;
      if(count%(std::getenv("EVERY")?std::atoi(std::getenv("EVERY")):10)==0){ double m=0; for(auto& c:psi) m=std::max(m,std::abs(c)); std::printf("%d %.2f %.5f %.5f %.4f %.3e\n",count,t,sJ/nb,sL/nb,std::sqrt(std::max(0.0,1-(sL/nb)*(sL/nb)/((sJ/nb)*(sJ/nb)))),m); sJ=sL=0; nb=0; } } } return 0; }
static int run743(int argc,char** argv){ setvbuf(stdout,nullptr,_IOLBF,0);
  auto envd=[](const char* k,double d){ const char* v=std::getenv(k); return v?std::atof(v):d; };
  const int periods=std::atoi(argv[1]); const double eps=std::atof(argv[2]), chi=std::atof(argv[3]), A=std::atof(argv[4]);
  const int N=argc>5?std::atoi(argv[5]):128; const double Lb=argc>6?std::atof(argv[6]):24.0, dx=Lb/N, sigma=0.3, rc=0.25;
  const double beta=envd("BETA",0.2), GT=envd("GT",0.1), Jr0=envd("JR0",0.5), L0=envd("L0",1.0); const int SUB=(int)envd("SUB",20), every=(int)envd("EVERY",1);
  const bool gauge=!(std::getenv("PHASE")&&std::string(std::getenv("PHASE"))=="K"), onep=std::getenv("ONEP")!=nullptr;
  const double Lam0=std::sqrt(L0*L0+2*beta), N0=Jr0+Lam0, Tr0=2*M_PI*N0*N0*N0, Gam=std::getenv("GABS")?envd("GABS",0):GT/Tr0, dt=Tr0/1000;
  std::vector<C> psi(N*N,0.0); std::vector<double> V(N*N), absb(N*N), k2(N*N);
  for(int iy=0;iy<N;++iy) for(int ix=0;ix<N;++ix){ const double x=-Lb/2+ix*dx, y=-Lb/2+iy*dx, r=std::max(std::sqrt(x*x+y*y),rc);
    V[iy*N+ix]=-1/r+beta/(r*r); const double d=std::min(std::min(x+Lb/2,Lb/2-x),std::min(y+Lb/2,Lb/2-y)), w=2.0;
    absb[iy*N+ix]= d<w? 0.5*std::pow((w-d)/w,2):0.0;
    const int kx=ix<N/2?ix:ix-N, ky=iy<N/2?iy:iy-N; const double qx=2*M_PI*kx/Lb, qy=2*M_PI*ky/Lb; k2[iy*N+ix]=qx*qx+qy*qy; }
  std::vector<C> potStep(N*N), kinStep(N*N);
  for(int i=0;i<N*N;++i){ potStep[i]=std::exp(C(-(Gam+absb[i])*dt*0.5,-V[i]*dt*0.5)); kinStep[i]=std::exp(C(0,-0.5*k2[i]*dt)); }
  const double a0=N0*N0, e0=std::sqrt(std::max(0.0,1-Lam0*Lam0/(N0*N0))), rp=a0*(1-e0);
  double x=rp,y=0,vx=0,vy=L0/rp, phis=0, t=0;
  auto Vp=[&](double px,double py){ const double r=std::sqrt(px*px+py*py); return -1/r+beta/(r*r); };
  auto ac=[&](double px,double py,double& ax,double& ay){ const double r2=px*px+py*py, r=std::sqrt(r2), g=1/(r2*r)-2*beta/(r2*r2); ax=-px*g; ay=-py*g; };
  auto field=[&](double px,double py,C& f,C& gx,C& gy){ const double fx=(px+Lb/2)/dx, fy=(py+Lb/2)/dx; const int ix=(int)std::floor(fx), iy=(int)std::floor(fy);
    const double ux=fx-ix, uy=fy-iy; auto P=[&](int i,int j){ i=(i%N+N)%N; j=(j%N+N)%N; return psi[j*N+i]; };
    auto bil=[&](auto g){ return (1-ux)*(1-uy)*g(ix,iy)+ux*(1-uy)*g(ix+1,iy)+(1-ux)*uy*g(ix,iy+1)+ux*uy*g(ix+1,iy+1); };
    f=bil([&](int i,int j){return P(i,j);}); gx=bil([&](int i,int j){return (P(i+1,j)-P(i-1,j))/(2*dx);}); gy=bil([&](int i,int j){return (P(i,j+1)-P(i,j-1))/(2*dx);}); };
  std::printf("# 743 beta %.4f GT %.4f Jr0 %.4f L0 %.4f Tr0 %.6f phase %s onep %d eps %.3e chi %.6f A %.3f N %d Lb %.2f SUB %d\n",beta,GT,Jr0,L0,Tr0,gauge?"L":"K",onep?1:0,eps,chi,A,N,Lb,SUB);
  std::printf("# k t J_r L ReZ ImZ <norm> |psi|max\n");
  int k=0; double vrPrev=0; C Zs=0; double Ns=0; long ns=0; bool srcOn=true;
  while(k<periods){
    for(int i=0;i<N*N;++i) psi[i]*=potStep[i];
    fft2(psi,N,false); for(int i=0;i<N*N;++i) psi[i]*=kinStep[i]; fft2(psi,N,true);
    for(int i=0;i<N*N;++i) psi[i]*=potStep[i];
    if(A!=0&&srcOn){ const C s=C(0,-1)*A*std::exp(C(0,phis))*dt; const int cx=(int)((x+Lb/2)/dx), cy=(int)((y+Lb/2)/dx), R=(int)(3*sigma/dx)+1;
      for(int j=cy-R;j<=cy+R;++j) for(int i=cx-R;i<=cx+R;++i){ if(i<0||j<0||i>=N||j>=N) continue; const double X=-Lb/2+i*dx-x, Y=-Lb/2+j*dx-y;
        psi[j*N+i]+=s*std::exp(-(X*X+Y*Y)/(2*sigma*sigma))/(2*M_PI*sigma*sigma); } }
    C f,gx,gy; field(x,y,f,gx,gy); Zs+=std::exp(C(0,-phis))*f; double nn=0; for(auto& c:psi) nn+=std::norm(c); Ns+=nn*dx*dx; ++ns;
    double fax=0,fay=0; if(eps!=0){ const C ph1=std::exp(C(0,-(phis+chi))); fax=eps*std::real(ph1*gx); fay=eps*std::real(ph1*gy); }
    const double h=dt/SUB; double ax,ay; ac(x,y,ax,ay);
    for(int s=0;s<SUB;++s){
      const double L1=0.5*(vx*vx+vy*vy)-(gauge?Vp(x,y):0.0);
      vx+=0.5*h*(ax+fax); vy+=0.5*h*(ay+fay); x+=h*vx; y+=h*vy; ac(x,y,ax,ay); vx+=0.5*h*(ax+fax); vy+=0.5*h*(ay+fay);
      const double L2=0.5*(vx*vx+vy*vy)-(gauge?Vp(x,y):0.0); phis+=0.5*(L1+L2)*h; }
    t+=dt;
    if(x*x+y*y>(Lb/2-2.5)*(Lb/2-2.5)){ std::printf("# left grid region at t=%.1f\n",t); break; }
    const double vr=x*vx+y*vy;
    if(vrPrev<0&&vr>=0){ ++k; if(onep) srcOn=false;
      const double r=std::sqrt(x*x+y*y), E=0.5*(vx*vx+vy*vy)+Vp(x,y);
      if(!(E<0)){ std::printf("# unbound\n"); break; }
      const double L=x*vy-y*vx, Nn=1/std::sqrt(-2*E), Jr=Nn-std::sqrt(L*L+2*beta); (void)r;
      if(k%every==0){ double m=0; for(auto& c:psi) m=std::max(m,std::abs(c));
        const C Z=Zs/(double)ns; std::printf("%d %.4f %.9f %.9f %.6e %.6e %.6e %.3e\n",k,t,Jr,L,Z.real(),Z.imag(),Ns/ns,m); }
      Zs=0; Ns=0; ns=0; }
    vrPrev=vr; }
  return 0; }
int main(int argc,char** argv){ return std::getenv("JR0")? run743(argc,argv): main446(argc,argv); }
