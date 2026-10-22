// Is the third-order residual of audit 142 a defect, or the momentum the
// FIELD and the RADIATION carry (audit section 564)?
//
// The synthetic circular pair of audits 122-142 (1 r*, tilt 0, exact circular
// kinematics in the history, moments along the axis) is rigidly rotating:
// every field, the hidden momentum and the interaction field momentum turn
// with angular velocity omega z.  Momentum conservation for the mutual
// (cross) part of the fields then says, with the engine force on a moment
// already carrying -d p_hid/dt (audit 123),
//
//   F1 + F2 = oint T_cross . n dA - omega z x (P_V + p_hid),
//
// P_V = eps0 int_{|x|<R} (E1 x B2 + E2 x B1) dV the interaction field
// momentum inside the sphere R about the centre, T_cross the cross terms of
// the Maxwell stress on that sphere, p_hid = sum m_i x E_j(r_i)/c^2.  The
// right side is computed from the model's own field functions -- no new
// physics -- and must not depend on R.  If it equals the engine's F1 + F2,
// the residual is physical (the field's and the radiation's momentum) and
// nothing is missing; if not, the difference is a defect of the force law.
//
// Quadrature: a smooth partition of unity -- a C-infinity bump of radius rho
// about each particle (integrated on a grid centred on that particle, which
// absorbs the 1/r^2 of the point charge field), and the rest on a grid
// centred on the centre of mass, cut at R.
//
// Usage: rotating_pair_momentum_balance [beta] [history nodes] [span r/c]
//        [R in separations] [resolution multiplier] [separation in r*]
// Build from the repository root (OpenMP for the volume integral):
// g++ -std=c++20 -O2 -fopenmp -I . $(root-config --cflags)
//     tools/rotating_pair_momentum_balance.cpp -o /tmp/rotating
//     $(root-config --libs)
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {
struct Rule { std::vector<double> x,w; };
// Gauss-Legendre on [a,b].
Rule gauss(int n,double a,double b) {
    Rule r; r.x.resize(n); r.w.resize(n);
    for(int i=0;i<n;++i) {
        double x=std::cos(pi*(i+0.75)/(n+0.5)), d=0.0;
        for(int it=0;it<100;++it) {
            double p0=1.0,p1=0.0;
            for(int k=0;k<n;++k){ const double p2=p1; p1=p0;
                p0=((2.0*k+1.0)*x*p1-k*p2)/(k+1.0); }
            d=n*(x*p0-p1)/(x*x-1.0);
            const double s=p0/d; x-=s; if(std::abs(s)<1e-15) break;
        }
        double p0=1.0,p1=0.0;
        for(int k=0;k<n;++k){ const double p2=p1; p1=p0;
            p0=((2.0*k+1.0)*x*p1-k*p2)/(k+1.0); }
        d=n*(x*p0-p1)/(x*x-1.0);
        r.x[i]=0.5*(b-a)*x+0.5*(b+a);
        r.w[i]=0.5*(b-a)*2.0/((1.0-x*x)*d*d);
    }
    return r;
}
Rule composite(const std::vector<double>& edges,int perSegment) {
    Rule all;
    for(size_t s=0;s+1<edges.size();++s) {
        const Rule g=gauss(perSegment,edges[s],edges[s+1]);
        all.x.insert(all.x.end(),g.x.begin(),g.x.end());
        all.w.insert(all.w.end(),g.w.begin(),g.w.end());
    }
    return all;
}
// C-infinity step: 1 for u <= 1/2, 0 for u >= 1.
double bump(double u) {
    if(u<=0.5) return 1.0;
    if(u>=1.0) return 0.0;
    const double t=(u-0.5)/0.5;
    const auto f=[](double y){ return y>0.0?std::exp(-1.0/y):0.0; };
    return f(1.0-t)/(f(1.0-t)+f(t));
}
}

int main(int argc,char** argv) {
    const double beta=argc>1?atof(argv[1]):0.1;
    const int nodes=argc>2?atoi(argv[2]):65;
    const double span=argc>3?atof(argv[3]):16.0;
    const double outerFactor=argc>4?atof(argv[4]):5.0;
    const double resolution=argc>5?atof(argv[5]):1.0;
    // argument 6: separation in r* (audit 564, uzupelnienie 1); R and the
    // bumps scale with it, the history span is in separation/c already
    const double separation=(argc>6?atof(argv[6]):1.0)*comptonBarrierRadius;
    const Vec3 axis{0.0,0.0,1.0};
    const double armRadius=0.5*separation;
    const double angular=beta*c/armRadius;
    const auto build=[&](double time) {
        State s{};
        s.time=time;
        const double phase=angular*time;
        const Vec3 offset{armRadius*std::cos(phase),armRadius*std::sin(phase),0.0};
        const Vec3 velocity{-armRadius*angular*std::sin(phase),
                             armRadius*angular*std::cos(phase),0.0};
        const Vec3 acceleration=offset*(-angular*angular);
        s.firstPosition=offset; s.secondPosition=offset*-1.0;
        s.firstVelocity=velocity; s.secondVelocity=velocity*-1.0;
        s.firstAcceleration=acceleration; s.secondAcceleration=acceleration*-1.0;
        s.firstProperDipole=axis*firstMagneticMoment;
        s.secondProperDipole=axis*secondMagneticMoment;
        synchronizeCovariantDipoles(s);
        return s;
    };
    StateHistory history;
    const double lightTime=separation/c;
    for(int i=0;i<nodes;++i)
        history.push_back(build(-span*lightTime*(1.0-static_cast<double>(i)/nodes)));
    const State present=build(0.0);
    const MutualForces forces=retardedExternalForces(present,history);
    const Vec3 net=forces.first+forces.second;
    const Vec3 omega=axis*angular;
    // the engine's own hidden-momentum rate force, and p_hid with the field
    // the engine uses for it (matched softening, audit 126)
    const Vec3 rateForce=hiddenMomentumRateForce(present,history,true)
        +hiddenMomentumRateForce(present,history,false);
    const auto electricAt=[&](bool targetIsFirst) {
        return fieldFromOtherParticleAt(targetIsFirst?present.firstPosition
            :present.secondPosition,present.time,present,history,
            targetIsFirst,1.0e-5,matchedMomentSoftening()).electric;
    };
    const Vec3 hidden=(cross(present.firstDipole,electricAt(true))
        +cross(present.secondDipole,electricAt(false)))*(1.0/(c*c));
    // fields of each source at a point (charge LW + moment), point charge
    const auto sources=[&](const Vec3& p,ElectromagneticField& one,
                           ElectromagneticField& two) {
        one=fieldFromOtherParticleAt(p,present.time,present,history,false);
        two=fieldFromOtherParticleAt(p,present.time,present,history,true);
    };
    const auto density=[&](const Vec3& p) {
        ElectromagneticField one,two;
        sources(p,one,two);
        return (cross(one.electric,two.magnetic)
               +cross(two.electric,one.magnetic))*epsilon0;
    };
    const double outer=outerFactor*separation;
    const double rho=0.4*separation;
    const auto spherical=[&](const Vec3& centre,const Rule& radial,int nPolar,
                             int nAzimuth,auto weightOf) {
        const Rule polar=gauss(nPolar,-1.0,1.0);
        std::vector<Vec3> partial(radial.x.size());
        #pragma omp parallel for schedule(dynamic)
        for(int ir=0;ir<static_cast<int>(radial.x.size());++ir) {
            const double r=radial.x[ir];
            Vec3 sum;
            for(int ip=0;ip<nPolar;++ip) {
                const double ct=polar.x[ip], st=std::sqrt(std::max(0.0,1.0-ct*ct));
                for(int ia=0;ia<nAzimuth;++ia) {
                    const double ph=2.0*pi*(ia+0.5)/nAzimuth;
                    const Vec3 p=centre+Vec3{r*st*std::cos(ph),r*st*std::sin(ph),r*ct};
                    const double w=weightOf(p);
                    if(w==0.0) continue;
                    const Vec3 d=density(p);
                    if(isFinite(d)) sum+=d*(w*polar.w[ip]);
                }
            }
            partial[ir]=sum*(radial.w[ir]*r*r*2.0*pi/nAzimuth);
        }
        Vec3 total; for(const Vec3& v:partial) total+=v;
        return total;
    };
    const auto bumpOf=[&](const Vec3& p,const Vec3& centre){
        return bump((p-centre).norm()/rho); };
    const int nr=static_cast<int>(24*resolution), np=static_cast<int>(32*resolution),
              na=static_cast<int>(64*resolution);
    const Rule nearRadial=composite({0.0,0.25*rho,0.5*rho,rho},nr/2);
    const Vec3 nearOne=spherical(present.firstPosition,nearRadial,np,na,
        [&](const Vec3& p){ return bumpOf(p,present.firstPosition); });
    const Vec3 nearTwo=spherical(present.secondPosition,nearRadial,np,na,
        [&](const Vec3& p){ return bumpOf(p,present.secondPosition); });
    std::vector<double> edges{0.0,0.1*separation,0.3*separation,0.45*separation,
        0.55*separation,0.7*separation,0.9*separation,1.2*separation,
        2.0*separation};
    for(double e=3.0*separation;e<outer;e*=1.5) edges.push_back(e);
    edges.push_back(outer);
    const Rule farRadial=composite(edges,nr/2);
    const Vec3 farPart=spherical(Vec3{},farRadial,2*np,2*na,
        [&](const Vec3& p){ return 1.0-bumpOf(p,present.firstPosition)
                                     -bumpOf(p,present.secondPosition); });
    const Vec3 fieldMomentum=nearOne+nearTwo+farPart;
    // momentum flux of the cross terms through the sphere: oint T . n dA
    Vec3 flux;
    {
        const Rule polar=gauss(2*np,-1.0,1.0);
        const int nAz=2*na;
        for(int ip=0;ip<2*np;++ip) {
            const double ct=polar.x[ip], st=std::sqrt(std::max(0.0,1.0-ct*ct));
            for(int ia=0;ia<nAz;++ia) {
                const double ph=2.0*pi*(ia+0.5)/nAz;
                const Vec3 n{st*std::cos(ph),st*std::sin(ph),ct};
                ElectromagneticField one,two;
                sources(n*outer,one,two);
                const auto stress=[&](const Vec3& a,const Vec3& b,double k){
                    return (a*dot(b,n)+b*dot(a,n)-n*dot(a,b))*k; };
                const Vec3 t=stress(one.electric,two.electric,epsilon0)
                    +stress(one.magnetic,two.magnetic,1.0/mu0);
                flux+=t*(polar.w[ip]*outer*outer*2.0*pi/nAz);
            }
        }
    }
    const Vec3 expected=flux-cross(omega,fieldMomentum+hidden);
    std::printf("beta %.4f  nodes %d  span %.0f r/c  R %.2f r*  resolution %.2f"
                "  separation %.2f r*\n",beta,nodes,span,
                outer/comptonBarrierRadius,resolution,
                separation/comptonBarrierRadius);
    std::printf("  engine F1+F2              %+.6e %+.6e %+.6e  |.| %.6e N"
                "  (|F1| %.6e N, ratio %.4e)\n",net.x,net.y,net.z,net.norm(),
                forces.first.norm(),net.norm()/forces.first.norm());
    std::printf("  engine hidden rate force  %+.6e %+.6e %+.6e  |.| %.6e N\n",
                rateForce.x,rateForce.y,rateForce.z,rateForce.norm());
    const Vec3 rotHidden=cross(omega,hidden);
    std::printf("  -omega x p_hid            %+.6e %+.6e %+.6e  (rate-force check |diff|/|.| %.3e)\n",
                -rotHidden.x,-rotHidden.y,-rotHidden.z,
                (rateForce+rotHidden).norm()/std::max(rotHidden.norm(),1e-300));
    std::printf("  p_hid                     %+.6e %+.6e %+.6e\n",hidden.x,hidden.y,hidden.z);
    std::printf("  P_V (near1 near2 far)     %+.6e %+.6e %+.6e  (|near| %.3e %.3e |far| %.3e)\n",
                fieldMomentum.x,fieldMomentum.y,fieldMomentum.z,
                nearOne.norm(),nearTwo.norm(),farPart.norm());
    std::printf("  |P_V + p_hid|/|p_hid|     %.6e\n",
                (fieldMomentum+hidden).norm()/std::max(hidden.norm(),1e-300));
    std::printf("  -omega x (P_V + p_hid)    %+.6e %+.6e %+.6e\n",
                -cross(omega,fieldMomentum+hidden).x,-cross(omega,fieldMomentum+hidden).y,
                -cross(omega,fieldMomentum+hidden).z);
    std::printf("  flux oint T.n dA          %+.6e %+.6e %+.6e\n",flux.x,flux.y,flux.z);
    std::printf("  expected                  %+.6e %+.6e %+.6e  |.| %.6e N\n",
                expected.x,expected.y,expected.z,expected.norm());
    std::printf("  |F1+F2 - expected|/|F1+F2| %.4e\n",
                (net-expected).norm()/std::max(net.norm(),1e-300));
    // CANDIDATES for the missing third-order term (audit 566).  Delta is
    // what the balance asks the forces to add.
    const Vec3 delta=expected-net;
    const double h=1.0e-4*separation;
    const auto partner=[&](bool targetIsFirst,const Vec3& x){
        return fieldFromOtherParticleAt(x,present.time,present,history,
            targetIsFirst,1.0e-5,matchedMomentSoftening()); };
    Vec3 k1, k2, k2a, k2b;
    for(int which=0;which<2;++which) {
        const bool first=which==0;
        const Vec3 x0=first?present.firstPosition:present.secondPosition;
        const Vec3 v=first?present.firstVelocity:present.secondVelocity;
        const Vec3 p=first?present.firstElectricDipole:present.secondElectricDipole;
        const Vec3 mLab=first?present.firstDipole:present.secondDipole;
        const Vec3 m0=first?present.firstProperDipole:present.secondProperDipole;
        Vec3 dE[3], dB[3];   // d/dx_k of E and B
        const Vec3 e[3]={{h,0,0},{0,h,0},{0,0,h}};
        for(int k=0;k<3;++k){
            const ElectromagneticField plus=partner(first,x0+e[k]), minus=partner(first,x0-e[k]);
            dE[k]=(plus.electric-minus.electric)*(1.0/(2*h));
            dB[k]=(plus.magnetic-minus.magnetic)*(1.0/(2*h)); }
        const Vec3 curlE{dE[1].z-dE[2].y,dE[2].x-dE[0].z,dE[0].y-dE[1].x};
        const Vec3 pGradB=dB[0]*p.x+dB[1]*p.y+dB[2]*p.z;
        const ElectromagneticField here=partner(first,x0);
        const Vec3 pRate=cross(omega,p);
        k1+=cross(p,curlE)*-1.0+cross(pRate,here.magnetic)+cross(v,pGradB);
        const Vec3 hidLab=cross(mLab,here.electric)*(1.0/(c*c));
        const Vec3 hidCov=cross(m0,here.electric+cross(v,here.magnetic))*(1.0/(c*c));
        k2+=cross(omega,hidCov-hidLab)*-1.0;
        k2a+=cross(omega,cross(m0,cross(v,here.magnetic)))*(-1.0/(c*c));
        k2b+=cross(omega,cross(m0-mLab,here.electric))*(-1.0/(c*c));
    }
    const auto report=[&](const char* name,const Vec3& k){
        std::printf("  candidate %-6s %+.6e %+.6e %+.6e  |Delta-K|/|Delta| %.4f  K.x/Delta.x %.4f\n",
            name,k.x,k.y,k.z,(delta-k).norm()/std::max(delta.norm(),1e-300),
            delta.x!=0.0?k.x/delta.x:0.0); };
    std::printf("  Delta = expected - F1+F2   %+.6e %+.6e %+.6e\n",delta.x,delta.y,delta.z);
    report("K1",k1); report("K2",k2); report("K1+K2",k1+k2);
    report("K2a",k2a); report("K2b",k2b);
    // the Thomas back-reaction's own net force (audit 567): L_T depends on
    // r_p - r_o, so F_p + F_o = -dA_T/dt, not zero
    {
        const MutualForces th=thomasBackReactionForces(present,
            thomasBmtDipoleDerivatives(present,history));
        const Vec3 tn=th.first+th.second;
        std::printf("  Thomas back-reaction net  %+.6e %+.6e %+.6e  |.|/|F1+F2| %.4f\n",
            tn.x,tn.y,tn.z,tn.norm()/std::max(net.norm(),1e-300));
    }
    std::printf("  |m_lab|/|m_0| - 1  %.6e   v.B check |v||B_partner|/(|E|) %.6e\n",
        present.firstDipole.norm()/present.firstProperDipole.norm()-1.0,
        present.firstVelocity.norm()*partner(true,present.firstPosition).magnetic.norm()
            /partner(true,present.firstPosition).electric.norm());
}
