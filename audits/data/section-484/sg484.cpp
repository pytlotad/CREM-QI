// Audit 484 (from 458): Stern-Gerlach with a sqrt(3)/2 hbar spin object.  argv: variant (0 rigid 0.866, 1 cone-averaged hbar/2 along
// the cone axis, 2 rigid 0.866 + cap-action attractor of 458, 3 cone held at 54.7 deg to B with the sign from the entry hemisphere)
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
int main(int argc,char** argv){ const int V=std::atoi(argv[1]); const double ell=(V==1)?0.5:std::sqrt(0.75); const bool attractor=(V==2);
  const double mu=firstMagneticMoment*(ell/0.5), g=std::abs(firstGFactor), M=107.8682*1.66053906660e-27, G=1000, B0=1.0, Lm=0.035, v=550, dt=1e-9;
  const double gamma=g*(1.602176634e-19)/(2*9.1093837015e-31);   // |gyromagnetic ratio| of the electron
  const double zmax=firstMagneticMoment*G*Lm*Lm/(2*M*v*v);   // 457 scale (full alignment of the model moment)
  const double q=std::exp(-1.0/3.0), lam=1e-4, gainF=lam*(1-q)*(1-q)/(2*M_PI*q);
  std::printf("# mu %.6e J/T, g %.8f, zmax(analytic) %.6e m\n",mu,g,zmax);
  for(int variant=0;variant<1;++variant){ int settled=0;
    std::mt19937_64 rng(457); std::uniform_real_distribution<double> U(0,1);
    const int N=2000, bins=30; std::vector<int> hist(bins,0); int centre=0, outer=0, up=0; double zmaxSim=0;
    for(int a=0;a<N;++a){
      const double ct=2*U(rng)-1, st=std::sqrt(1-ct*ct), ph=2*M_PI*U(rng);
      Vec3 m{mu*st*std::cos(ph),mu*st*std::sin(ph),mu*ct};          // moment direction (spin opposite for the electron; irrelevant here)
      double x=(U(rng)-0.5)*1e-4, z=0, vx=0, vz=0, t=0;
      if(V==3){ const Vec3 B{-G*x,0,B0+G*z}; const Vec3 bh=B*(1/B.norm()); const double c=dot(m,bh)/mu; const double cc=(c>=0?1:-1)/std::sqrt(3.0);
        Vec3 perp=m-bh*dot(m,bh); const double pn=perp.norm(); if(pn>1e-30){ perp=perp*(1/pn); m=(bh*cc+perp*std::sqrt(1-cc*cc))*mu; } }
      while(t<Lm/v){
        const Vec3 B{-G*x,0,B0+G*z};
        // force grad(m.B) = (-G m_x, 0, G m_z)
        const double ax=-G*m.x/M, az=G*m.z/M;
        vx+=ax*dt; vz+=az*dt; x+=vx*dt; z+=vz*dt; t+=dt;
        // exact precession about B over dt (Larmor; sense irrelevant for the projection)
        m=rotateDipoleByAngularVelocity(m,B*gamma,dt);
        // attractor on the cap action, applied per Larmor-averaged step: dA = -lambda_t gain F(2 pi (A - 1/2)) dt
        if(attractor){ const Vec3 bh=B*(1/B.norm()); double c=dot(m,bh)/mu; c=std::max(-1.0,std::min(1.0,c));
          const double A=ell*(1-c), th=2*M_PI*(A-0.5), F=q*std::sin(th)/(1-2*q*std::cos(th)+q*q);
          const double wL=gamma*B.norm(), lamT=gainF*wL/(2*M_PI);
          const double A2=std::max(0.0,std::min(2*ell,A-lamT*F*dt)); const double c2=1-A2/ell;
          // rotate m in the (m, B) plane to the new cosine
          Vec3 perp=m-bh*dot(m,bh); const double pn=perp.norm();
          if(pn>1e-30){ perp=perp*(1/pn); m=(bh*c2+perp*std::sqrt(std::max(0.0,1-c2*c2)))*mu; } } 
      }
      { const Vec3 B{-G*x,0,B0+G*z}; const double c=dot(m,B*(1/B.norm()))/mu; const double A=ell*(1-c); const double rung=std::floor(A)+0.5;
        if(std::abs((1-rung/ell)-c)<1e-3) ++settled; }
      int b=(int)((z/zmax+1.5)/3.0*bins); if(b<0) b=0; if(b>=bins) b=bins-1; ++hist[b];
      if(std::abs(z)<0.2*zmax) ++centre; if(std::abs(z)>0.9*zmax) ++outer; if(z>0) ++up; zmaxSim=std::max(zmaxSim,std::abs(z));
    }
    std::printf("== variant %d, l = %.4f: |z|<0.2 zmax %.4f, |z|>0.9 zmax %.4f, up fraction %.4f, max|z|/zmax %.4f, settled %.4f\n",V,ell,
      double(centre)/N,double(outer)/N,double(up)/N,zmaxSim/zmax,double(settled)/N);
    for(int b=0;b<bins;++b) std::printf("  z/zmax %+.2f  %4d\n",-1.5+(b+0.5)*3.0/bins,hist[b]);
  }
}
