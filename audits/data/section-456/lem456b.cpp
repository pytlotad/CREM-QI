// Audit 456: field angular momentum of a point charge q (at origin) and a point magnetic dipole m (at R), about the charge:
//   L = eps0 int r x (E x B) d^3r = (q/4pi) int [rhat (rhat.B) - B]/r d^3r   (E = q rhat/(4 pi eps0 r^2))
// B = (mu0/4pi)[3 n (n.m) - m]/s^3 + (2 mu0/3) m delta^3(s), s = r - R.  Units mu0/4pi = 1, q = 1, |R| = 1.
// Quadrature in shells around the dipole (Gauss-Legendre in cos(theta), uniform in phi), s mapped to (0, inf).
#include <cmath>
#include <cstdio>
#include <vector>
#include <cstdlib>
struct V{double x,y,z;}; V operator+(V a,V b){return{a.x+b.x,a.y+b.y,a.z+b.z};} V operator-(V a,V b){return{a.x-b.x,a.y-b.y,a.z-b.z};}
V operator*(V a,double s){return{a.x*s,a.y*s,a.z*s};} double dot(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;} double nrm(V a){return std::sqrt(dot(a,a));}
void gauleg(int n,std::vector<double>& x,std::vector<double>& w){ x.resize(n); w.resize(n);
  for(int i=0;i<n;++i){ double z=std::cos(M_PI*(i+0.75)/(n+0.5)),z1,pp; do{ double p1=1,p2=0; for(int j=0;j<n;++j){ double p3=p2; p2=p1; p1=((2*j+1)*z*p2-j*p3)/(j+1);} pp=n*(z*p1-p2)/(z*z-1); z1=z; z=z1-p1/pp; }while(std::fabs(z-z1)>1e-15); x[i]=z; w[i]=2/((1-z*z)*pp*pp);} }
V Lfield(V m){ const V R{0,0,1};
  std::vector<double> ct,wt; gauleg(192,ct,wt); const int NPHI=384; std::vector<double> u,wu; gauleg(800,u,wu);
  V L{0,0,0};
  for(size_t iu=0;iu<u.size();++iu){ // s = t/(1-t), t in (0,1): t = (u+1)/2
    const double t=0.5*(u[iu]+1), s=t/(1-t), ds=0.5*wu[iu]/((1-t)*(1-t));
    for(size_t ic=0;ic<ct.size();++ic){ const double c=ct[ic], sn=std::sqrt(1-c*c);
      for(int ip=0;ip<NPHI;++ip){ const double ph=2*M_PI*(ip+0.5)/NPHI; const V n{sn*std::cos(ph),sn*std::sin(ph),c};
        const V B=(n*(3*dot(n,m))-m)*(1/(s*s*s)); const V r=R+n*s; const double rr=nrm(r); const V rh=r*(1/rr);
        const V f=(rh*dot(rh,B)-B)*(1/rr); const double w=ds*s*s*wt[ic]*(2*M_PI/NPHI);
        L=L+f*w; } } }
  // contact term (2 mu0/3) m delta^3 at R: mu0 = 4 pi in these units -> (8 pi/3) m
  const V Bd=m*(8*M_PI/3); const V fd=(R*dot(R,Bd)-Bd)*(1.0/1.0);
  static const bool noContact=std::getenv("NO_CONTACT")!=nullptr; if(!noContact) L=L+fd;
  return L*(1/(4*M_PI)); }
int main(){
  V a=Lfield({0,0,1}), b=Lfield({1,0,0}), c=Lfield({0,1,0});
  std::printf("m || R: L = (%.6f, %.6f, %.6f)\nm _|_ R (x): L = (%.6f, %.6f, %.6f)\nm _|_ R (y): L = (%.6f, %.6f, %.6f)\n",a.x,a.y,a.z,b.x,b.y,b.z,c.x,c.y,c.z);
  std::printf("in units (mu0 q/(4 pi R)): L = a m + b (m.Rhat) Rhat with a = %.6f (from m _|_ R), a + b = %.6f (from m || R)\n",b.x,a.z);
}
