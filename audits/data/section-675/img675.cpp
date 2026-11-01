// Audit 675: the engine's classical 1S energy of the 661 image (I0) and of images with n_z = <sigma_z> (I1a antiparallel / I1b parallel transverse)
// against Breit-Rabi, for the photons of ../section-668/runs.  build: g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules img675.cpp $(root-config --libs)
// run (repository root): PAIR=antiproton,positron ./img675 < audits/data/section-675/photons675.txt
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
int main(){
  applyPairFromOption(std::getenv("PAIR"));
  const double reducedMass=pairReducedMass;
  gSingletChannelEligible=false;
  const double u=1.602176634e-25, Bm=1.0329, m1=firstMagneticMoment, m2=secondMagneticMoment;
  const double sg1=firstCharge<0?-1:1, sg2=secondCharge<0?-1:1, binding=pairBindingEnergy(activePair);
  const double attr=pairCoulombStrength/reducedMass, L1=0.5*hbar, e1=-binding/reducedMass;
  const double semi1=attr/(2*std::abs(e1)), ecc1=std::min(0.999999,std::sqrt(std::max(0.0,1+2*e1*(L1/reducedMass)*(L1/reducedMass)/(attr*attr))));
  std::string name; int sflag; double bx,by,bz,a1,a2,K,E1S,after,S1S,Zsa,PS,EA,LA; int m; Vec3 ax,pe,mA1,mA2;
  double mx[5]={0,0,0,0,0}; int n=0;
  std::printf("photon m P_S | E1S check | after check | S1S* engine, I0, I1a, I1b [ueV] | Zeeman part I0, I1a [ueV]\n");
  while(std::cin>>name>>sflag>>bx>>by>>bz>>a1>>a2>>K>>E1S>>after>>S1S>>Zsa>>m>>PS>>EA>>LA>>ax.x>>ax.y>>ax.z>>pe.x>>pe.y>>pe.z
                >>mA1.x>>mA1.y>>mA1.z>>mA2.x>>mA2.y>>mA2.z){
    const Vec3 b{bx,by,bz}, B=b*Bm;
    std::array<double,16> H{K+a1+a2,0,0,0, 0,-K+a1-a2,2*K,0, 0,2*K,-K-a1+a2,0, 0,0,0,K-a1-a2};
    std::array<double,4> w; std::array<double,16> V; psJacobiEigen<4>(H,w,V);
    int pick=-1; double best=1e9;
    for(int i=0;i<4;++i){ const double v0=V[i],v1=V[4+i],v2=V[8+i],v3=V[12+i]; const int mi=v0*v0>0.5?1:(v3*v3>0.5?-1:0);
      const double ps=0.5*(v1-v2)*(v1-v2); if(mi==m&&std::abs(ps-PS)<best){best=std::abs(ps-PS);pick=i;} }
    const double v0=V[pick],v1=V[4+pick],v2=V[8+pick],v3=V[12+pick];
    const double s1=v0*v0+v1*v1-v2*v2-v3*v3, s2=v0*v0-v1*v1+v2*v2-v3*v3;
    const double U1S=levelEnergyBeyondBohr(semi1,ecc1,ax,pe,b*(sg1*m1),b*(sg2*m2),L1);
    const double EBR=U1S-K+w[pick];
    const double eA=EA/reducedMass, lA=LA/reducedMass, semiA=attr/(2*std::abs(eA));
    const double eccA=std::min(0.999999,std::sqrt(std::max(0.0,1+2*eA*lA*lA/(attr*attr))));
    // the engine's singlet-channel flag at the photon (ledger "singlet"), as when it evaluates "after"
    const auto U=[&](const Vec3& n1,const Vec3& n2){ gSingletChannelEligible=sflag!=0;
      const double r=levelEnergyBeyondBohr(semiA,eccA,ax,pe,n1*(sg1*m1),n2*(sg2*m2),LA); gSingletChannelEligible=false; return r; };
    const auto Zs=[&](const Vec3& n1,const Vec3& n2){ return -dot(n1*(sg1*m1)+n2*(sg2*m2),B); };
    const Vec3 n1=mA1*(sg1/mA1.norm()), n2=mA2*(sg2/mA2.norm());
    Vec3 perp=n1-b*dot(n1,b); const bool hasPerp=perp.norm()>1e-9; if(hasPerp) perp=perp*(1.0/perp.norm());
    const auto img=[&](double sgnT){ if(!hasPerp||m!=0) return std::pair<Vec3,Vec3>{n1,n2};
      return std::pair<Vec3,Vec3>{perp*std::sqrt(std::max(0.0,1-s1*s1))+b*s1, perp*(sgnT*std::sqrt(std::max(0.0,1-s2*s2)))+b*s2}; };
    const double s0=(EBR-(U(n1,n2)+Zs(n1,n2)))/u;
    const auto ia=img(-1.0), ib=img(1.0);
    const double sa=(EBR-(U(ia.first,ia.second)+Zs(ia.first,ia.second)))/u, sb=(EBR-(U(ib.first,ib.second)+Zs(ib.first,ib.second)))/u;
    const double zBR=(a1*s1+a2*s2)/u;
    const double c1=std::abs(EBR-E1S)/u, c2=std::abs(U(n1,n2)-after)/u, c3=std::abs(s0-(S1S-Zsa)/u);
    mx[0]=std::max(mx[0],c1); mx[1]=std::max(mx[1],c2); mx[2]=std::max(mx[2],c3); mx[3]=std::max(mx[3],std::abs(sa)); mx[4]=std::max(mx[4],std::abs(sb)); ++n;
    std::printf("%s %+d %.4f | %.1e | %.1e | %+8.3f %+8.3f %+8.3f %+8.3f | %+8.3f %+8.3f\n",name.c_str(),m,PS,c1,c2,(S1S-Zsa)/u,s0,sa,sb,
                zBR-Zs(n1,n2)/u,zBR-Zs(ia.first,ia.second)/u);
  }
  std::printf("\n%d photons: max |E_BR - E1S printed| %.2e ueV, max |U(I0) - after printed| %.2e ueV, max |S1S*(I0) - engine| %.2e ueV -> H675v %s\n",
              n,mx[0],mx[1],mx[2],(mx[2]<=0.01&&mx[0]<=0.01)?"TAK":"NIE");
  std::printf("max |S1S*| I1a %.3f ueV -> H675a %s;  I1b %.3f ueV -> H675b %s\n",mx[3],mx[3]<=1?"TAK":"NIE",mx[4],mx[4]<=1?"TAK":"NIE");
}
