// Audit 446: particle as an oscillating source of a matter-wave field (double-solution / HQFT, non-relativistic envelope).
// 2D, hbar = mu = k = 1.  Field: i dpsi/dt = -lap psi/2 + V psi - i Gamma psi - i absorb psi + source,
//   source = -i A G_sigma(x - x_p) exp(i phi_s),  phi_s = int K dt (de Broglie clock minus omega_C t; <K> = -E by virial).
// Particle: a = -grad V - gamma v + eps Re[exp(-i (phi_s + chi)) grad psi(x_p)].
// argv: J0 e0 orbits eps chi gamma A Gamma_orbits [N=128 Lbox=24]
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>
typedef std::complex<double> C;
static void fft1(std::vector<C>& a,bool inv){ const int n=a.size();
  for(int i=1,j=0;i<n;++i){ int bit=n>>1; for(;j&bit;bit>>=1) j^=bit; j^=bit; if(i<j) std::swap(a[i],a[j]); }
  for(int len=2;len<=n;len<<=1){ const double ang=2*M_PI/len*(inv?1:-1); const C wl(std::cos(ang),std::sin(ang));
    for(int i=0;i<n;i+=len){ C w(1); for(int j=0;j<len/2;++j){ const C u=a[i+j], v=a[i+j+len/2]*w; a[i+j]=u+v; a[i+j+len/2]=u-v; w*=wl; } } }
  if(inv) for(auto& x:a) x/=n; }
static void fft2(std::vector<C>& f,int N,bool inv){ std::vector<C> row(N);
  for(int y=0;y<N;++y){ for(int x=0;x<N;++x) row[x]=f[y*N+x]; fft1(row,inv); for(int x=0;x<N;++x) f[y*N+x]=row[x]; }
  for(int x=0;x<N;++x){ for(int y=0;y<N;++y) row[y]=f[y*N+x]; fft1(row,inv); for(int y=0;y<N;++y) f[y*N+x]=row[y]; } }
int main(int argc,char** argv){ setvbuf(stdout,nullptr,_IOLBF,0);
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
      if(count%(std::getenv("EVERY")?std::atoi(std::getenv("EVERY")):10)==0){ double m=0; for(auto& c:psi) m=std::max(m,std::abs(c)); std::printf("%d %.2f %.5f %.5f %.4f %.3e\n",count,t,sJ/nb,sL/nb,std::sqrt(std::max(0.0,1-(sL/nb)*(sL/nb)/((sJ/nb)*(sJ/nb)))),m); sJ=sL=0; nb=0; } } } }
