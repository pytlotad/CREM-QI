// Audit 263: the round trip.  Does the rotation at the emission site
// put the pattern where it belongs in the lab?
//
// 262d closed the ACCUMULATION side: the m=2 phase reproduces the
// Kepler true anomaly, so the moments really are in the orbital frame
// the bookkeeping claims.  What was never checked is the other end --
// the rotation out of that frame when a photon is drawn.  A handedness
// error there would mirror the azimuth, and nothing downstream would
// notice: the amplitude, the power and the polar part would all be
// untouched.
//
// The test calls drawLabDirectionFromOrbitalPattern -- the same
// function the estimator calls -- and checks the SAMPLE against an
// analytic tensor rotation.  For p(n) ~ c0 + n.Q.n,
//   <n_i n_j> = delta_ij/3 + 2 Q_ij/(15 c0),
// so Q is recoverable from the sample's second moment as
//   Q_ij = (15 c0/2)(<n_i n_j> - delta_ij/3),
// and the prediction is Q_lab = R Q_orb R^T with R the frame's columns.
// Two independent objects, not a retyping of the code under test.
//
// The mirrored column is the meta-test: it predicts with the OPPOSITE
// handedness.  If that residual were also small, this test would have
// no teeth and would prove nothing.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <random>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    // A pattern with ALL FIVE moments nonzero, so a swapped or
    // sign-flipped basis vector cannot hide.
    const double c0=1.0;
    const double cm0=0.40,cm2=0.30,cxy=0.25,cyz=0.20,cxz=0.15;
    AngularPatternMoments m;
    m.a00=4.0*pi*c0;
    m.a2[0]=cxy*(4.0*pi/15.0);
    m.a2[1]=cyz*(4.0*pi/15.0);
    m.a2[2]=cm0*(4.0*pi/5.0);
    m.a2[3]=cxz*(4.0*pi/15.0);
    m.a2[4]=cm2*(4.0*pi/15.0);
    double Qorb[3][3]={
        {-0.5*cm0+0.5*cm2, 0.5*cxy,          0.5*cxz},
        { 0.5*cxy,         -0.5*cm0-0.5*cm2, 0.5*cyz},
        { 0.5*cxz,          0.5*cyz,         cm0}};
    // A frame aligned with nothing, and a periapsis that needs
    // orthogonalizing -- so that line is exercised too.
    const Vec3 rawThird{1.0,2.0,3.0};
    const Vec3 third=rawThird*(1.0/rawThird.norm());
    const Vec3 rawPeriapsis{3.0,-1.0,1.0};
    Vec3 first=rawPeriapsis-third*dot(rawPeriapsis,third);
    first=first*(1.0/first.norm());
    const Vec3 second=cross(third,first);
    const Vec3 mirrored=cross(first,third);   // opposite handedness
    const auto rotate=[&](const Vec3& e1,const Vec3& e2,const Vec3& e3,
                          double out[3][3]){
        const Vec3 cols[3]={e1,e2,e3};
        for(int i=0;i<3;++i) for(int j=0;j<3;++j){
            double sum=0.0;
            for(int p=0;p<3;++p) for(int q=0;q<3;++q){
                const double Rip=(p==0?cols[0]:p==1?cols[1]:cols[2])
                    .x*0.0; (void)Rip;
                const double a=(i==0?cols[p].x:i==1?cols[p].y:cols[p].z);
                const double b=(j==0?cols[q].x:j==1?cols[q].y:cols[q].z);
                sum+=a*Qorb[p][q]*b;
            }
            out[i][j]=sum;
        }
    };
    double Qgood[3][3],Qbad[3][3];
    rotate(first,second,third,Qgood);
    rotate(first,mirrored,third,Qbad);
    std::mt19937_64 rng(98765);
    std::uniform_real_distribution<double> u(0.0,1.0);
    const auto draw=[&]{ return u(rng); };
    const int n=4000000;
    double acc[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    long long got=0,failed=0;
    for(int i=0;i<n;++i){
        const Vec3 d=drawLabDirectionFromOrbitalPattern(
            m,third,rawPeriapsis,draw);
        if(d.squaredNorm()==0.0){ ++failed; continue; }
        ++got;
        const double v[3]={d.x,d.y,d.z};
        for(int a=0;a<3;++a) for(int b=0;b<3;++b) acc[a][b]+=v[a]*v[b];
    }
    double Qhat[3][3];
    for(int a=0;a<3;++a) for(int b=0;b<3;++b)
        Qhat[a][b]=(15.0*c0/2.0)*(acc[a][b]/got-(a==b?1.0/3.0:0.0));
    double worstGood=0.0,worstBad=0.0;
    for(int a=0;a<3;++a) for(int b=0;b<3;++b){
        worstGood=std::max(worstGood,std::abs(Qhat[a][b]-Qgood[a][b]));
        worstBad=std::max(worstBad,std::abs(Qhat[a][b]-Qbad[a][b]));
    }
    std::printf("# probek=%lld nieudanych=%lld\n",got,failed);
    std::printf("%10s %14s %14s %14s\n","element","z probki","przewidziane",
        "lustrzane");
    const char* nm[3][3]={{"xx","xy","xz"},{"yx","yy","yz"},{"zx","zy","zz"}};
    for(int a=0;a<3;++a) for(int b=a;b<3;++b)
        std::printf("%10s %14.6f %14.6f %14.6f\n",nm[a][b],
            Qhat[a][b],Qgood[a][b],Qbad[a][b]);
    std::printf("\n# najgorsza roznica wzgledem PRZEWIDZIANEGO:  %.3e\n",
        worstGood);
    std::printf("# najgorsza roznica wzgledem LUSTRZANEGO:      %.3e\n",
        worstBad);
    std::printf("# test ma zeby, jesli druga liczba jest duzo wieksza\n");
}
