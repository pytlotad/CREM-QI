// Audit 645: Thomas-BMT of the engine in a g-2 storage-ring setting.  Free mu+ (antimuon species) in uniform B along z, optional radial electric field
// E = -E0 r^; motion by RK4 with the engine's lorentzForce / relativisticAcceleration, spin by the exact rotation omega = -(q/m) B_BMT
// (thomasBmtEffectiveField, midpoint velocity).  omega_a from a linear fit of the unwrapped spin-velocity angle over one 2 pi/omega_a; Eq. 16
// a = R/(lambda - R), R = omega_a/omega_p, omega_p = proton Larmor at rest (same BMT function), lambda = mu_mu/mu_p from the engine's moments.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <vector>
struct Result { double omegaA, steps; };
static Result run(double gammaValue,double B,double E0,int stepsPerTurn){
  const ParticleSpecies& s=antimuon; const double q=s.charge, m=s.mass, g=s.gFactor, a=0.5*(g-2);
  const double beta=std::sqrt(1-1/(gammaValue*gammaValue)), v=beta*c;
  const double R=gammaValue*m*v*v/(q*(v*B+E0));                 // circle with both forces inward
  Vec3 x{R,0,0}, u{0,-v,0};                                        // clockwise for mu+ in +z
  const Vec3 Bv{0,0,B};
  const auto field=[&](const Vec3& p){ const double rn=std::sqrt(p.x*p.x+p.y*p.y); return ElectromagneticField{Vec3{-E0*p.x/rn,-E0*p.y/rn,0},Bv}; };
  const auto acc=[&](const Vec3& p,const Vec3& w){ return relativisticAcceleration(w,lorentzForce(q,w,field(p)),m); };
  Vec3 spin=u*(1/u.norm());                                        // longitudinal start
  const double wc=q*(v*B+E0)/(gammaValue*m*v);                    // orbital angular speed
  const double wa_guess=q/m*std::abs(a*B+(a-1/(gammaValue*gammaValue-1))*beta*E0/c);
  const double dt=2*pi/wc/stepsPerTurn; const long N=(long)std::ceil(2*pi/wa_guess/dt);
  double sx=0,sy=0,sxx=0,sxy=0,prev=0,unwrap=0; long n=0;
  for(long k=0;k<=N;++k){
    if(k%16==0||k==N){ const Vec3 vh=u*(1/u.norm()); const double ang=std::atan2(dot(cross(vh,spin),Vec3{0,0,1}),dot(vh,spin));
      double d=ang-prev; while(d>pi) d-=2*pi; while(d<-pi) d+=2*pi; unwrap+=d; prev=ang; const double t=k*dt;
      sx+=t; sy+=unwrap; sxx+=t*t; sxy+=t*unwrap; ++n; }
    if(k==N) break;
    // RK4 for (x, u); spin by the BMT rotation at the midpoint velocity
    const Vec3 k1x=u, k1u=acc(x,u);
    const Vec3 k2x=u+k1u*(dt/2), k2u=acc(x+k1x*(dt/2),u+k1u*(dt/2));
    const Vec3 k3x=u+k2u*(dt/2), k3u=acc(x+k2x*(dt/2),u+k2u*(dt/2));
    const Vec3 k4x=u+k3u*dt, k4u=acc(x+k3x*dt,u+k3u*dt);
    const Vec3 xn=x+(k1x+k2x*2+k3x*2+k4x)*(dt/6), un=u+(k1u+k2u*2+k3u*2+k4u)*(dt/6);
    // midpoint state: mean direction at the mean speed (a stage velocity of RK4 is off |u| by O(dt^2), which moves gamma)
    Vec3 umid=u+un; umid=umid*(0.5*(u.norm()+un.norm())/umid.norm());
    Vec3 xmid=x+xn; xmid=xmid*(0.5*(x.norm()+xn.norm())/xmid.norm());
    const Vec3 omega=thomasBmtEffectiveField(umid,field(xmid),g)*(-q/m);
    spin=rotateDipoleByAngularVelocity(spin,omega,dt);
    x=xn; u=un;
  }
  const double slope=(n*sxy-sx*sy)/(n*sxx-sx*sx);
  return {std::abs(slope),(double)N};
}
int main(){
  const double B=1.45; const ParticleSpecies& mu=antimuon; const double a=0.5*(mu.gFactor-2);
  const double wp=(proton.charge/proton.mass*thomasBmtEffectiveField(Vec3{},ElectromagneticField{Vec3{},Vec3{0,0,B}},proton.gFactor)).norm();
  const double lambda=magneticMoment(mu)/magneticMoment(proton);
  std::printf("engine: a_mu = (g-2)/2 = %.14f, omega_p = %.9e rad/s (%.6f MHz), lambda = mu_mu/mu_p = %.12f; Fermilab 2023 a_mu 0.00116592059(22)\n",
    a,wp,wp/2/pi/1e6,lambda);
  for(double E0: {0.0,1.0e6}) for(double gm: {2.0,10.0,29.3,100.0}){
    const double beta=std::sqrt(1-1/(gm*gm));
    const double expect=mu.charge/mu.mass*std::abs(a*B+(a-1/(gm*gm-1))*beta*E0/c);
    for(int spt: {2000,4000,8000,16000}){
      const Result r=run(gm,B,E0,spt); const double Rr=r.omegaA/wp, aEq=Rr/(lambda-Rr);
      std::printf("E0 %.0e V/m gamma %6.1f steps/turn %d: omega_a %.12e rad/s, BMT formula %.12e (rel %+.2e); Eq.16 a = %.14f (rel to engine a %+.2e)\n",
        E0,gm,spt,r.omegaA,expect,r.omegaA/expect-1,aEq,aEq/a-1);
    }
  }
}
