#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 515: import 21 (CREM_SOMMERFELD_J=1) -- hydrogen level energy F at k = j + 1/2 for l >= 1, read through the engine's
// levelEnergyBeyondBohr over isotropic electron-spin directions; basins, identity check, fine-structure intervals vs measurements.
int main(){
  setenv("CREM_SOMMERFELD_J","1",1);
  applyPairFromOption("proton,electron");
  const double aB=pairBohrRadius(activePair), hP=2*pi*hbar, MHz=1e-6/hP;
  const double ae=0.5*(std::abs(secondGFactor)-2.0);
  const Vec3 z{0,0,1}, x{1,0,0};
  const Vec3 mp=z*(firstMagneticMoment);                     // proton moment, fixed (hyperfine is subtracted below)
  const int N=20000;
  double Ej[16][16][2]={};                                   // [n][l][k - l]
  std::printf("hydrogen, import 21 on; electron a_e = %.9f, 1 + 2 a_e = %.6f\n",ae,1+2*ae);
  for(int n:{2,3,4,8,12}){
    const double a=n*n*aB;
    for(int l=1;l<n;++l){
      const double e=std::sqrt(1.0-std::pow((l+0.5)/n,2)), L=(l+0.5)*hbar;
      long cnt[2]={0,0}; double sum[2]={0,0}, maxdev=0, cmin[2]={2,2}, cmax[2]={-2,-2};
      for(int i=0;i<N;++i){
        const double c=1.0-2.0*(i+0.5)/N, s=std::sqrt(1-c*c), ph=pi*(3.0-std::sqrt(5.0))*i;
        const Vec3 sh{s*std::cos(ph),s*std::sin(ph),c};
        const Vec3 me=sh*(-1.0*secondMagneticMoment);         // electron: s = sign(q) m^ -> m = -|mu| s
        const double k=sommerfeldLevelAngularMomentum(z,mp,me,L)/hbar;
        const int b=static_cast<int>(std::lround(k))-l;
        const double Flev=levelEnergyBeyondBohr(a,e,z,x,mp,me,L)-spinCouplingEnergy(a,e,z,x,mp,me,L);
        const double Fk=levelFineStructureEnergy(a,k*hbar);
        maxdev=std::max(maxdev,std::abs(Flev-Fk)/std::abs(Fk));
        cnt[b]++; sum[b]+=Flev; cmin[b]=std::min(cmin[b],c); cmax[b]=std::max(cmax[b],c);
      }
      for(int b=0;b<2;++b) Ej[n][l][b]=sum[b]/cnt[b];
      std::printf("n %2d l %2d: basin j = l - 1/2 (k = %d) %.4f [QM %.4f], cos(L,s) in [%+.4f, %+.4f]; j = l + 1/2 (k = %d) %.4f [QM %.4f], cos in [%+.4f, %+.4f]"
                  " (threshold -1/(2l+1) = %+.4f); |F_level - F(k)|/|F| max %.1e\n",
        n,l,l,cnt[0]/(double)N,l/(2.0*l+1),cmin[0],cmax[0],l+1,cnt[1]/(double)N,(l+1.0)/(2.0*l+1),cmin[1],cmax[1],-1.0/(2*l+1),maxdev);
      std::printf("          E(j = l - 1/2) %+.4f MHz, E(j = l + 1/2) %+.4f MHz, interval %.4f MHz; without import (Langer F) interval 0\n",
        Ej[n][l][0]*MHz,Ej[n][l][1]*MHz,(Ej[n][l][1]-Ej[n][l][0])*MHz);
    }
    const double FS=levelFineStructureEnergy(a,0.5*hbar), F1=levelFineStructureEnergy(a,1.0*hbar);
    std::printf("n %2d S state: F(Langer 1/2 + contact) %+.4f MHz, F(k = 1) %+.4f MHz, difference %.2e MHz\n",n,FS*MHz,F1*MHz,(FS-F1)*MHz);
    Ej[n][0][0]=FS;
  }
  struct M{const char* name; double model, meas, err;};
  const M m[]={
    {"2P3/2 - 2P1/2 (Hagley, Pipkin 1994 + Lundeen, Pipkin 1981)",(Ej[2][1][1]-Ej[2][1][0])*MHz,10969.045,0.015},
    {"4P3/2 - 4P1/2 (Berkeland, Hinds, Boshier 1995)",(Ej[4][1][1]-Ej[4][1][0])*MHz,1371.104,0.018},
    {"4D5/2 - 4P3/2 (Weitz et al. 1995 - Berkeland et al. 1995)",(Ej[4][2][1]-Ej[4][1][1])*MHz,454.771,0.026},
    {"8D5/2 - 8D3/2 (de Beauvoir et al. 1997)",(Ej[8][2][1]-Ej[8][2][0])*MHz,57.1342,0.0105},
    {"12D5/2 - 12D3/2 (Schwob et al. 1999)",(Ej[12][2][1]-Ej[12][2][0])*MHz,16.9310,0.0117},
    {"3D5/2 - 3S1/2 (Clark, Van Baak, Lundeen, Pipkin 1979)",(Ej[3][2][1]-Ej[3][0][0])*MHz,4013.204,0.078},
  };
  std::printf("intervals (model F with import 21 vs measurement, MHz):\n");
  for(const M& q:m) std::printf("  %-60s model %11.4f  meas %11.4f(%g)  residual meas - model %+9.4f  ratio %.6f\n",q.name,q.model,q.meas,q.err,q.meas-q.model,q.meas/q.model);
  std::printf("predictions without a measurement in hand: 3P3/2 - 3P1/2 %.4f MHz, 3D5/2 - 3D3/2 %.4f MHz, 3D3/2 - 3P3/2 %.4f MHz\n",
    (Ej[3][1][1]-Ej[3][1][0])*MHz,(Ej[3][2][1]-Ej[3][2][0])*MHz,(Ej[3][2][0]-Ej[3][1][1])*MHz);
}
