// Audit 451: two particles (e-, e+) each sourcing its own matter-wave field (envelope, Schrodinger-type, mass 1) that evolves in
// the partner's moving Coulomb potential; each particle is pushed by the gradient of its own and/or the partner's field, in the
// phase of its de Broglie clock phi_i = int K_i dt.  2D, CM frame, hbar = m_e = k = 1.
// argv: J0 e0 orbits eps chi w_own w_partner [A=1] [Gamma_orbits=3] [N=128] [L=24]
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
  const double eps=std::atof(argv[4]), chi=std::atof(argv[5]), wo=std::atof(argv[6]), wp=std::atof(argv[7]);
  const double A=argc>8?std::atof(argv[8]):1.0, Gorb=argc>9?std::atof(argv[9]):3.0; const int N=argc>10?std::atoi(argv[10]):128;
  const double Lb=argc>11?std::atof(argv[11]):24.0, dx=Lb/N, soft=0.1, sigma=0.3;
  // relative orbit with mu = 1/2: a = 2 J^2 ... for Ps (k = 1, mu = 1/2): E = -1/(4 J^2), a = 2 J^2, T = 4 pi J^3
  const double T0=4*M_PI*J0*J0*J0, Gam=1.0/(Gorb*T0), dt=T0/1000;
  std::vector<C> psi[2]={std::vector<C>(N*N,0.0),std::vector<C>(N*N,0.0)}; std::vector<C> kin(N*N); std::vector<double> absb(N*N), X(N*N), Y(N*N);
  for(int iy=0;iy<N;++iy) for(int ix=0;ix<N;++ix){ const int i=iy*N+ix; X[i]=-Lb/2+ix*dx; Y[i]=-Lb/2+iy*dx;
    const double d=std::min(std::min(X[i]+Lb/2,Lb/2-X[i]),std::min(Y[i]+Lb/2,Lb/2-Y[i])), w=2.0; absb[i]= d<w?0.5*std::pow((w-d)/w,2):0.0;
    const int kx=ix<N/2?ix:ix-N, ky=iy<N/2?iy:iy-N; const double qx=2*M_PI*kx/Lb, qy=2*M_PI*ky/Lb; kin[i]=std::exp(C(0,-0.5*(qx*qx+qy*qy)*dt)); }
  // relative orbit at periapsis: r = a(1-e), relative speed v from vis-viva with mu = 1/2: v^2 = 2 (2/r - 1/a)... k/mu = 2
  const double a=2*J0*J0, rp=a*(1-e0), vrel=std::sqrt(2*(2/rp-1/a));
  double x[2]={rp/2,-rp/2}, y[2]={0,0}, vx[2]={0,0}, vy[2]={vrel/2,-vrel/2}, ph[2]={0,0};
  auto sample=[&](const std::vector<C>& f,double px,double py,C& gx,C& gy){ const double fx=(px+Lb/2)/dx, fy=(py+Lb/2)/dx; const int ix=(int)std::floor(fx), iy=(int)std::floor(fy);
    const double ux=fx-ix, uy=fy-iy; auto P=[&](int i,int j){ i=(i%N+N)%N; j=(j%N+N)%N; return f[j*N+i]; };
    auto bil=[&](auto g){ return (1-ux)*(1-uy)*g(ix,iy)+ux*(1-uy)*g(ix+1,iy)+(1-ux)*uy*g(ix,iy+1)+ux*uy*g(ix+1,iy+1); };
    gx=bil([&](int i,int j){return (P(i+1,j)-P(i-1,j))/(2*dx);}); gy=bil([&](int i,int j){return (P(i,j+1)-P(i,j-1))/(2*dx);}); };
  auto accel=[&](int i,const double* xs,const double* ys,const double* phs,double& ax,double& ay){ const int j=1-i; const double dxp=xs[i]-xs[j], dyp=ys[i]-ys[j];
    const double r=std::sqrt(dxp*dxp+dyp*dyp+soft*soft); ax=-dxp/(r*r*r); ay=-dyp/(r*r*r);   // attraction (opposite charges, k = 1, mass 1)
    if(eps!=0){ C gx(0),gy(0),hx,hy; if(wo!=0){ sample(psi[i],xs[i],ys[i],hx,hy); gx+=wo*hx; gy+=wo*hy; } if(wp!=0){ sample(psi[j],xs[i],ys[i],hx,hy); gx+=wp*hx; gy+=wp*hy; }
      const C e=std::exp(C(0,-(phs[i]+chi))); ax+=eps*std::real(e*gx); ay+=eps*std::real(e*gy); } };
  static const bool sym=std::getenv("SYMMETRIZE")!=nullptr;
  auto accelAll=[&](const double* xs,const double* ys,const double* phs,double* AX,double* AY){
    for(int i=0;i<2;++i) accel(i,xs,ys,phs,AX[i],AY[i]);
    if(sym){ // remove the net (non-reciprocal) part of the wave forces: Coulomb parts are equal and opposite already
      const double mx=0.5*(AX[0]+AX[1]), my=0.5*(AY[0]+AY[1]); AX[0]-=mx; AX[1]-=mx; AY[0]-=my; AY[1]-=my; } };
  double ax[2],ay[2]; accelAll(x,y,ph,ax,ay);
  double angAcc=0, t=0; int count=0; double sJ=0,sL=0; int nb=0;
  std::printf("# orbit t <J_rel> <L_rel> e |P_cm| |psi|max\n");
  std::vector<C> pot(N*N);
  while(count<orbits){
    for(int f=0;f<2;++f){ const int j=1-f;   // field f evolves in the potential of particle j (the partner of its source)
      for(int i=0;i<N*N;++i){ const double dxp=X[i]-x[j], dyp=Y[i]-y[j]; const double V=-1/std::sqrt(dxp*dxp+dyp*dyp+soft*soft); pot[i]=std::exp(C(-(Gam+absb[i])*dt*0.5,-V*dt*0.5)); psi[f][i]*=pot[i]; }
      fft2(psi[f],N,false); for(int i=0;i<N*N;++i) psi[f][i]*=kin[i]; fft2(psi[f],N,true);
      for(int i=0;i<N*N;++i) psi[f][i]*=pot[i];
      if(A!=0){ const C s=C(0,-1)*A*std::exp(C(0,ph[f]))*dt; const int cx=(int)((x[f]+Lb/2)/dx), cy=(int)((y[f]+Lb/2)/dx), R=(int)(3*sigma/dx)+1;
        for(int jy=cy-R;jy<=cy+R;++jy) for(int jx=cx-R;jx<=cx+R;++jx){ if(jx<0||jy<0||jx>=N||jy>=N) continue; const double Xd=X[jy*N+jx]-x[f], Yd=Y[jy*N+jx]-y[f];
          psi[f][jy*N+jx]+=s*std::exp(-(Xd*Xd+Yd*Yd)/(2*sigma*sigma))/(2*M_PI*sigma*sigma); } } }
    double nx[2],ny[2],hvx[2],hvy[2];
    for(int i=0;i<2;++i){ hvx[i]=vx[i]+0.5*dt*ax[i]; hvy[i]=vy[i]+0.5*dt*ay[i]; nx[i]=x[i]+dt*hvx[i]; ny[i]=y[i]+dt*hvy[i]; ph[i]+=0.5*(vx[i]*vx[i]+vy[i]*vy[i])*dt; }
    const double relAng0=std::atan2(y[0]-y[1],x[0]-x[1]);
    accelAll(nx,ny,ph,ax,ay);
    for(int i=0;i<2;++i){ vx[i]=hvx[i]+0.5*dt*ax[i]; vy[i]=hvy[i]+0.5*dt*ay[i]; x[i]=nx[i]; y[i]=ny[i]; }
    double d=std::atan2(y[0]-y[1],x[0]-x[1])-relAng0; if(d>M_PI) d-=2*M_PI; if(d<-M_PI) d+=2*M_PI; angAcc+=d; t+=dt;
    if(std::hypot(x[0],y[0])>9.5||std::hypot(x[1],y[1])>9.5){ std::printf("# left grid region at t=%.1f\n",t); break; }
    if(angAcc>=2*M_PI){ angAcc-=2*M_PI; ++count;
      const double rx=x[0]-x[1], ry=y[0]-y[1], r=std::sqrt(rx*rx+ry*ry+soft*soft), prx=(vx[0]-vx[1])/2, pry=(vy[0]-vy[1])/2;   // p_rel = mu v_rel, mu = 1/2
      const double E=(prx*prx+pry*pry)-1/r;   // p^2/(2 mu) = p^2
      if(!(E<0)){ std::printf("# unbound\n"); break; } const double J=1/std::sqrt(-4*E), L=rx*pry-ry*prx; sJ+=J; sL+=L; ++nb;
      if(count%10==0){ double m=0; for(int f=0;f<2;++f) for(auto& c:psi[f]) m=std::max(m,std::abs(c));
        const double Pcm=std::hypot(vx[0]+vx[1],vy[0]+vy[1]); const double Jm=sJ/nb, Lm=sL/nb;
        std::printf("%d %.2f %.5f %.5f %.4f %.3e %.3e\n",count,t,Jm,Lm,std::sqrt(std::max(0.0,1-Lm*Lm/(Jm*Jm))),Pcm,m); sJ=sL=0; nb=0; } } } }
