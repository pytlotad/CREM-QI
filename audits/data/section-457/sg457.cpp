// Audit 457: Stern-Gerlach with the model's electron moment.  (A) classical spin dynamics only, (B) Born projection at entry.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
int main(){
  const double mu=firstMagneticMoment, g=std::abs(firstGFactor), M=107.8682*1.66053906660e-27, G=1000, B0=1.0, Lm=0.035, v=550, dt=1e-9;
  const double gamma=g*(1.602176634e-19)/(2*9.1093837015e-31);   // |gyromagnetic ratio| of the electron
  const double zmax=mu*G*Lm*Lm/(2*M*v*v);
  std::printf("# mu %.6e J/T, g %.8f, zmax(analytic) %.6e m\n",mu,g,zmax);
  for(int variant=0;variant<2;++variant){
    std::mt19937_64 rng(457); std::uniform_real_distribution<double> U(0,1);
    const int N=2000, bins=30; std::vector<int> hist(bins,0); int centre=0, outer=0, up=0; double zmaxSim=0;
    for(int a=0;a<N;++a){
      const double ct=2*U(rng)-1, st=std::sqrt(1-ct*ct), ph=2*M_PI*U(rng);
      Vec3 m{mu*st*std::cos(ph),mu*st*std::sin(ph),mu*ct};          // moment direction (spin opposite for the electron; irrelevant here)
      double x=(U(rng)-0.5)*1e-4, z=0, vx=0, vz=0, t=0;
      if(variant==1){ const Vec3 B{-G*x,0,B0+G*z}; const Vec3 bh=B*(1/B.norm()); const double c=dot(m,bh)/mu;
        m = (U(rng)<0.5*(1+c))? bh*mu : bh*(-mu); }
      while(t<Lm/v){
        const Vec3 B{-G*x,0,B0+G*z};
        // force grad(m.B) = (-G m_x, 0, G m_z)
        const double ax=-G*m.x/M, az=G*m.z/M;
        vx+=ax*dt; vz+=az*dt; x+=vx*dt; z+=vz*dt; t+=dt;
        // exact precession about B over dt (Larmor; sense irrelevant for the projection)
        m=rotateDipoleByAngularVelocity(m,B*gamma,dt);
      }
      int b=(int)((z/zmax+1.5)/3.0*bins); if(b<0) b=0; if(b>=bins) b=bins-1; ++hist[b];
      if(std::abs(z)<0.2*zmax) ++centre; if(std::abs(z)>0.9*zmax) ++outer; if(z>0) ++up; zmaxSim=std::max(zmaxSim,std::abs(z));
    }
    std::printf("== variant %s: |z|<0.2 zmax %.4f, |z|>0.9 zmax %.4f, up fraction %.4f, max|z|/zmax %.4f\n",variant?"(B) Born at entry":"(A) classical",
      double(centre)/N,double(outer)/N,double(up)/N,zmaxSim/zmax);
    for(int b=0;b<bins;++b) std::printf("  z/zmax %+.2f  %4d\n",-1.5+(b+0.5)*3.0/bins,hist[b]);
  }
}
