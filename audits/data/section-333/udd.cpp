// Audit 333: orbit-averaged dipole-dipole energy U = -m1.B2 (tensor = Plummer pole field, contact =
// magnetization term) from CREM_SPINVEC lines, with the engine's own field functions.
#include "modules/electrodynamics.hpp"
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
static Vec3 v3(const char* s){ Vec3 v; std::sscanf(s,"%lf,%lf,%lf",&v.x,&v.y,&v.z); return v; }
int main(int argc,char** argv){
    const int N=argc>1?std::atoi(argv[1]):20001; const double eps=magneticDipoleRadius();
    const double eV=1.602176634e-19;
    std::string line;
    std::printf("# t[ps]  a/a_pair  e   cos12   U_tensor[meV]  U_contact[meV]  U_total[meV]  binding[eV]\n");
    while(std::getline(std::cin,line)){
        if(line.rfind("CREM_SPINVEC",0)!=0) continue;
        double t,a,e; char m1s[200],m2s[200],Ls[200],Ps[200];
        if(std::sscanf(line.c_str(),"CREM_SPINVEC t=%lf a=%lf e=%lf m1=%199s m2=%199s L=%199s P=%199s",&t,&a,&e,m1s,m2s,Ls,Ps)!=7) continue;
        const Vec3 m1=v3(m1s), m2=v3(m2s); Vec3 L=v3(Ls), P=v3(Ps);
        L=L*(1.0/L.norm()); P=P-L*dot(P,L); P=P*(1.0/P.norm()); const Vec3 Q=cross(L,P);
        const double ec=std::sqrt(std::max(0.0,1-e*e));
        double Ut=0,Uc=0,W=0;
        for(int s=-1;s<=1;s+=2) for(int k=0;k<N;++k){
            const double u=(k+0.5)/N, E=s*pi*u*u*u*u, dE=pi*4*u*u*u/N, w=(1-e*std::cos(E))*dE;
            const Vec3 r=P*(a*(std::cos(E)-e))+Q*(a*ec*std::sin(E));
            Ut+=w*(-dot(m1,plummerDipoleField(r,m2,eps))); Uc+=w*(-dot(m1,plummerMagnetizationField(r,m2,eps))); W+=w;
        }
        Ut/=W; Uc/=W;
        const double binding=pairCoulombStrength/(2*a);
        std::printf("%9.4f  %.5f  %.5f  %+.5f  %+.5e  %+.5e  %+.5e  %.4f\n",t*1e12,a/pairBohrRadius(activePair),e,
            dot(m1,m2)/(m1.norm()*m2.norm()),Ut/eV*1e3,Uc/eV*1e3,(Ut+Uc)/eV*1e3,binding/eV);
    }
}
