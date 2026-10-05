// Audit 331: orbital vs partner-dipole BMT precession, para vs ortho, with the model's
// own field and BMT functions (the per-node body of orbitAveragedBmtAngularVelocities).
#include "modules/secular_spin_orbit.hpp"
#include <cstdio>
#include <random>
struct Rates { Vec3 w1, w2; };
static double mu_red(){ return firstMass*secondMass/(firstMass+secondMass); }
static Rates ratesAt(double a,double e,double E,const Vec3& m1,const Vec3& m2){
    const double M=firstMass+secondMass, mu=mu_red();
    const double n=std::sqrt(pairCoulombStrength/(mu*a*a*a)), ec=std::sqrt(1-e*e);
    const double cs=std::cos(E), sn=std::sin(E), tw=1-e*cs;
    const Vec3 r={a*(cs-e),a*ec*sn,0}, v=Vec3{-a*n*sn,a*n*ec*cs,0}*(1.0/tw);
    State s{}; s.firstPosition=r*(secondMass/M); s.secondPosition=r*(-firstMass/M);
    s.firstVelocity=v*(secondMass/M); s.secondVelocity=v*(-firstMass/M);
    s.firstProperDipole=m1; s.secondProperDipole=m2; synchronizeCovariantDipoles(s);
    const StateHistory h{State{s}}; const LocalElectromagneticFields f=localRelativisticFields(s,h);
    return {thomasBmtEffectiveField(s.firstVelocity,f.atFirst,firstGFactor)*(-firstCharge/firstMass),
            thomasBmtEffectiveField(s.secondVelocity,f.atSecond,secondGFactor)*(-secondCharge/secondMass)};
}
int main(int argc,char**argv){
    const double mu=firstMagneticMoment, a=pairBohrRadius(activePair), rs=comptonBarrierRadius;
    const Vec3 zero{1e-300*0,0,0};
    // (1) profile along a near-radial orbit with periapsis = r*: e = 1 - r*/a
    const double e=1.0-rs/a;
    std::printf("# PROFILE a=a_pair, e=%.6f (periapsis = r*), moments along z (normal) and x (in plane)\n",e);
    std::printf("# r/r*   |w_orb|(rad/s)  para: w_partner.mhat2   ortho: w_partner.mhat2   [moment axis]\n");
    for(const char* ax: {"z","x"}){
      const Vec3 m=(ax[0]=='z')?Vec3{0,0,mu}:Vec3{mu,0,0};
      for(double E: {3.14159265,1.5,0.6,0.3,0.15,0.08,0.04,0.02,0.01,0.0}){
        const double r=a*(1-e*std::cos(E));
        const Rates orb=ratesAt(a,e,E,m,zero*0.0);           // partner moment zero
        const Rates par=ratesAt(a,e,E,m,m), ort=ratesAt(a,e,E,m,m*(-1.0));
        const Vec3 wpP=par.w1-orb.w1, wpO=ort.w1-orb.w1;
        std::printf("%9.3f  %.4e   %+.4e   %+.4e   [%s]\n",r/rs,orb.w1.norm(),dot(wpP,m)/mu,dot(wpO,m*(-1.0))/mu,ax);
      }
    }
    // (2) microcanonical ensemble, deterministic: q = J0^2 = 1 - e^2 on a geometric grid down to 1e-7
    //     (as in audit 324), time average over E = pi u^4 nodes (periapsis-clustered), isotropic moment
    //     average EXACT via the three axes (m.B(m) is quadratic in mhat).
    const int Q=argc>1?std::atoi(argv[1]):120, U=argc>2?std::atoi(argv[2]):4001;
    const double qmin=argc>3?std::atof(argv[3]):3e-4; std::vector<double> q; for(int i=0;i<Q;++i) q.push_back(qmin*std::pow(1.0/qmin,(double)i/(Q-1)));
    std::vector<double> fP,fOrb,fN; const double eps=0.96682*rs;
    auto nker=[&](double r){ return 3*eps*eps/(4*pi*std::pow(r*r+eps*eps,2.5)); };
    for(double qq: q){
        const double ee=std::sqrt(std::max(0.0,1.0-qq)); double num=0,orbn=0,den=0,nn=0,nden=0;
        for(int ax=0;ax<3;++ax){ Vec3 m{0,0,0}; if(ax==0) m.x=mu; else if(ax==1) m.y=mu; else m.z=mu;
          for(int k=0;k<U;++k){ const double u=(k+0.5)/U, E=pi*u*u*u*u, dE=pi*4*u*u*u/U, tw=1-ee*std::cos(E);
            const Rates orb=ratesAt(a,ee,E,m,zero*0.0), par=ratesAt(a,ee,E,m,m);
            num+=tw*dE*dot(par.w1-orb.w1,m)/mu; orbn+=tw*dE*orb.w1.norm(); den+=tw*dE; nn+=tw*dE*nker(a*tw); nden+=tw*dE; } }
        fP.push_back(num/den); fOrb.push_back(orbn/den); fN.push_back(nn/nden);
    }
    double avg=0,avgOrb=0,avgN=0; for(size_t i=1;i<q.size();++i){ avg+=0.5*(fP[i]+fP[i-1])*(q[i]-q[i-1]); avgOrb+=0.5*(fOrb[i]+fOrb[i-1])*(q[i]-q[i-1]); avgN+=0.5*(fN[i]+fN[i-1])*(q[i]-q[i-1]); }
    const double dE324=2*(2*mu0/3)*mu*mu*avgN/hbar;
    const double hfs=2*pi*203.3942e9;
    std::printf("# MICROCANONICAL n=1 deterministic: %d q-points, %d nodes, 3 axes\n",(int)q.size(),U);
    std::printf("Delta E/hbar = <mhat . w_partner(para)> = %.4e rad/s = %.2f GHz x 2pi = %.3f x measured 203.3942 GHz\n",avg,avg/2/pi/1e9,avg/hfs);
    std::printf("  classical channel energies: E_para = -Delta E/2, E_ortho = +Delta E/2 (QM: -3/4, +1/4 of Delta E)\n");
    std::printf("<|w_orbital|> = %.4e rad/s\n",avgOrb);
    std::printf("same q-range, audit-324 contact formula 2(2mu0/3)mu^2<n>/hbar = %.4e rad/s = %.2f GHz x 2pi; engine/formula = %.4f\n",dE324,dE324/2/pi/1e9,avg/dE324);
    std::printf("q range [%.1e, 1]: circle %.4e rad/s, q_min %.4e rad/s\n",qmin,fP.back(),fP[0]);
}
