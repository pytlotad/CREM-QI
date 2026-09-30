// Audit 252 control: the same stencil on a field that satisfies
// Faraday's law EXACTLY, by closed form.
//
// A point charge in uniform motion has
//   E = q/(4 pi eps0) (1-b^2) Rhat / (R^2 (1-b^2 sin^2 th)^{3/2})
//   B = (v x E)/c^2
// with R the vector from the PRESENT position.  This is an exact
// Maxwell solution, so any residual the stencil reports is the
// stencil, not the field.  If it reproduces the model's beta^-2 law
// and its magnitude, the model's Faraday floor is exonerated.
#include "modules/crem_collapse.hpp"
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double q=1.602176634e-19;
    const double pre=q/(4.0*pi*epsilon0);
    const double r0=2.6462e-11;              // as in the model probe
    std::printf("%9s %11s %11s %11s %11s\n",
        "beta","|curlE|/c","|rotE+dB|/c","Far rel","Amp rel");
    for(double beta:{7.297e-3,1.032e-2,1.459e-2,2.063e-2,2.914e-2}){
        const Vec3 vd{1.0,0.2,-0.3};
        const Vec3 v=vd*(beta*c/vd.norm());
        const auto E=[&](const Vec3& x,double t){
            const Vec3 R=x-v*t;
            const double R2=R.squaredNorm();
            const double b2=beta*beta;
            const double ct=dot(R,v)/(std::sqrt(R2)*v.norm());
            const double s2=1.0-ct*ct;
            return R*(pre*(1.0-b2)
                /(R2*std::sqrt(R2)*std::pow(1.0-b2*s2,1.5)));
        };
        const auto field=[&](const Vec3& x,double t){
            const Vec3 e=E(x,t);
            return ElectromagneticField{e,cross(v,e)*(1.0/(c*c))};
        };
        const Vec3 xd{0.3,0.5,0.8};
        const Vec3 x=xd*(r0/xd.norm());
        const double hs=1.0e-4*r0, ht=hs/c;
        Vec3 curlB{},curlE{};
        for(int ax=0;ax<3;++ax){
            Vec3 d{}; (ax==0?d.x:ax==1?d.y:d.z)=hs;
            const ElectromagneticField p=field(x+d,0.0);
            const ElectromagneticField m=field(x-d,0.0);
            const Vec3 dB=(p.magnetic-m.magnetic)*(1.0/(2.0*hs));
            const Vec3 dE=(p.electric-m.electric)*(1.0/(2.0*hs));
            const int k2=(ax+1)%3,l=(ax+2)%3;
            double* cb[3]={&curlB.x,&curlB.y,&curlB.z};
            double* ce[3]={&curlE.x,&curlE.y,&curlE.z};
            const double dBk=(k2==0?dB.x:k2==1?dB.y:dB.z);
            const double dBl=(l==0?dB.x:l==1?dB.y:dB.z);
            const double dEk=(k2==0?dE.x:k2==1?dE.y:dE.z);
            const double dEl=(l==0?dE.x:l==1?dE.y:dE.z);
            *cb[l]+=dBk; *cb[k2]-=dBl;
            *ce[l]+=dEk; *ce[k2]-=dEl;
        }
        const ElectromagneticField tp=field(x,ht),tm=field(x,-ht);
        const Vec3 dEdt=(tp.electric-tm.electric)*(1.0/(2.0*ht));
        const Vec3 dBdt=(tp.magnetic-tm.magnetic)*(1.0/(2.0*ht));
        const double farAbs=(curlE+dBdt).norm()/c;
        const double ampAbs=(curlB-dEdt*(1.0/(c*c))).norm();
        const double sF=std::max(curlE.norm(),dBdt.norm());
        const double sA=std::max(curlB.norm(),dEdt.norm()/(c*c));
        std::printf("%9.3e %11.4e %11.4e %11.4e %11.4e\n",
            beta,curlE.norm()/c,farAbs,farAbs*c/sF,ampAbs/sA);
    }
}
