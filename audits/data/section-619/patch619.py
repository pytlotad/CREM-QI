# Audit 619: patch for crem_collapse.hpp -- Born at formation (J, j, F without a field; 1S eigenstates in a field).
p='modules/crem_collapse.hpp'; s=open(p).read()
def rep(old,new,once=False):
    global s
    n=s.count(old)
    assert n==1 or (once and n>=1), (n,old[:60]); s=s.replace(old,new,1)
helpers=r'''// BORN AT FORMATION (audit 619, the author's proposal; default,
// CREM_FORMATION_BASINS=1 restores rule Q's basins at formation).  A basin
// depends on which configurations make up the ensemble, not on its density
// matrix (496); at formation the label (J of the Ps triplet P, j and F of
// hydrogen / muonium) is drawn with the Born probability of the coherent state
// of the classical configuration, and the configuration is turned into the
// centre of the drawn label's basin, so that the on-the-fly rounding (462, 515,
// 520) keeps returning it.  Isotropic weights are those of the basins exactly
// (born619.txt); a single configuration differs (mean 0.27-0.42).
inline bool formationBasins() {
    static const bool on=std::getenv("CREM_FORMATION_BASINS")!=nullptr;
    return on;
}
inline double formationCg(double j1,double m1,double j2,double m2,double J,double M) {
    if(std::abs(m1+m2-M)>1e-9||J<std::abs(j1-j2)-1e-9||J>j1+j2+1e-9||std::abs(M)>J+1e-9) return 0.0;
    const auto f=[](double x){ return std::tgamma(std::round(x)+1.0); };
    const double pre=std::sqrt((2*J+1)*f(J+j1-j2)*f(J-j1+j2)*f(j1+j2-J)/f(j1+j2+J+1))
        *std::sqrt(f(J+M)*f(J-M)*f(j1-m1)*f(j1+m1)*f(j2-m2)*f(j2+m2));
    double sum=0.0;
    for(int k=0;k<200;++k) {
        const double a[6]={double(k),j1+j2-J-k,j1-m1-k,j2+m2-k,J-j2+m1+k,J-j1-m2+k};
        bool ok=true; for(double x: a) if(x<-1e-9) ok=false;
        if(!ok) { if(a[1]<-1e-9||a[2]<-1e-9||a[3]<-1e-9) break; continue; }
        double den=1.0; for(double x: a) den*=f(x);
        sum+=(k%2?-1.0:1.0)/den;
    }
    return pre*sum;
}
// Born probabilities of |j1, axis> (x) |j2, theta from the axis> on |J>, J = |j1 - j2| .. j1 + j2.
inline std::vector<double> formationBornWeights(double j1,double j2,double theta) {
    std::vector<double> p;
    const double c=std::cos(0.5*theta), sn=std::sin(0.5*theta);
    for(double J=std::abs(j1-j2);J<=j1+j2+1e-9;J+=1.0) {
        double w=0.0;
        for(double m2=-j2;m2<=j2+1e-9;m2+=1.0) {
            const double binom=std::tgamma(2*j2+1)/(std::tgamma(j2-m2+1)*std::tgamma(j2+m2+1));
            const double d2=binom*std::pow(c,2*(j2+m2))*std::pow(sn,2*(j2-m2));
            const double g=formationCg(j1,j1,j2,m2,J,j1+m2); w+=d2*g*g;
        }
        p.push_back(w);
    }
    return p;
}
inline int formationDraw(const std::vector<double>& p,double u) {
    double t=0.0; for(double x: p) t+=x;
    double acc=0.0;
    for(size_t i=0;i<p.size();++i) { acc+=p[i]/t; if(u<acc) return static_cast<int>(i); }
    return static_cast<int>(p.size())-1;
}
// Turn vv about k^ = a^ x v^ by the angle that brings angle(a^, v^) to target.
inline Vec3 formationTurn(const Vec3& v,const Vec3& a,double target,const Vec3& vv) {
    const Vec3 vh=v*(1.0/v.norm());
    Vec3 k=cross(a,vh); double sk=k.norm();
    if(!(sk>1e-12)) { k=cross(a,orbitPlaneDirection(a,{})); sk=k.norm(); }
    k=k*(1.0/sk);
    const double now=std::acos(std::clamp(dot(a,vh),-1.0,1.0)), ang=target-now;
    const double cs=std::cos(ang), sn=std::sin(ang);
    return vv*cs+cross(k,vv)*sn+k*(dot(k,vv)*(1.0-cs));
}
inline CremCollapseEstimate estimateCremCollapse(std::uint64_t seed,'''
rep("inline CremCollapseEstimate estimateCremCollapse(std::uint64_t seed,", helpers, once=True)
rep("""        if(bornTripletInField&&reducedMass*elements.specificAngularMomentum
                                <hbar*(1.0-1.0e-9)) {
            const Vec3 bHat=gExternalMagneticField*(1.0/gExternalMagneticField.norm());
            const Vec3 axisNow=firstDipole*(firstSign/firstMoment);""","""        if(bornTripletInField&&reducedMass*elements.specificAngularMomentum
                                <hbar*(1.0-1.0e-9)&&!formationBasins()) {
            bornTriplet1SInField=true;   // audit 619: Born on the 1S field eigenstates below
        } else if(bornTripletInField&&reducedMass*elements.specificAngularMomentum
                                <hbar*(1.0-1.0e-9)) {
            const Vec3 bHat=gExternalMagneticField*(1.0/gExternalMagneticField.norm());
            const Vec3 axisNow=firstDipole*(firstSign/firstMoment);""")
rep("    bool bornTripletInField=false;\n","    bool bornTripletInField=false;\n    bool bornTriplet1SInField=false;   // audit 619\n")
rep("    if(!(periapsisDirection.norm()>0.0)) return result;\n    // RULE Q BASIN ON THE CONSERVED AXIS AT l >= 1 (audit 544;",
r"""    if(!(periapsisDirection.norm()>0.0)) return result;
    // Audit 619 (a): a triplet formed in the s band in a field -- Born weights of
    // |s^ s^> on the 1S field eigenstates (H = (Delta E/4) sigma1.sigma2 + Zeeman)
    // and the classical image of the drawn state (539, 559), as at the photon of 547.
    if(bornTriplet1SInField&&elements.specificEnergy<0.0) {
        const double orbital=reducedMass*elements.specificAngularMomentum;
        const double semi=-attractionParameter/(2.0*elements.specificEnergy);
        const double ecc=std::min(0.999999,std::sqrt(std::max(0.0,1.0
            +2.0*elements.specificEnergy*elements.specificAngularMomentum
                *elements.specificAngularMomentum
                /(attractionParameter*attractionParameter))));
        const double m1=firstDipole.norm(), m2=secondDipole.norm();
        const double sgn1=firstCharge<0.0?-1.0:1.0, sgn2=secondCharge<0.0?-1.0:1.0;
        const Vec3 bHat=gExternalMagneticField*(1.0/gExternalMagneticField.norm());
        const double triplet=spinCouplingEnergy(semi,ecc,angularMomentumDirection,
            periapsisDirection,bHat*(sgn1*m1),bHat*(sgn2*m2),orbital);
        const double singletE=spinCouplingEnergy(semi,ecc,angularMomentumDirection,
            periapsisDirection,bHat*(sgn1*m1),bHat*(-sgn2*m2),orbital);
        const double Bm=gExternalMagneticField.norm();
        PsFieldEigenState saved=gPsFieldEigen, st;
        psFieldFrame(gExternalMagneticField,st.e1,st.e2,st.b);
        st.a1=-sgn1*m1*Bm; st.a2=-sgn2*m2*Bm;
        const Vec3 sh=firstDipole*(sgn1/m1);
        const Vec3 inF{dot(sh,st.e1),dot(sh,st.e2),dot(sh,st.b)};
        std::complex<double> ch[2]; psCoherentHalf(inF,ch);
        std::complex<double> psi4[4];
        for(int a=0;a<2;++a) for(int b=0;b<2;++b) psi4[2*a+b]=ch[a]*ch[b];
        for(int a=0;a<4;++a) for(int c=0;c<4;++c) st.rho[a][c]=std::real(psi4[a]*std::conj(psi4[c]));
        gPsFieldEigen=st;
        std::uint64_t bornStream=splitMix64(seed^0x619b0e1fULL);
        const double u=drawUniformUnit(bornStream), psi=2.0*pi*drawUniformUnit(bornStream);
        const PsFieldEigenImage im=psFieldEigenAtPhoton(triplet-singletE,u,psi);
        gPsFieldEigen=saved;
        firstDipole=im.n1*(sgn1*m1); secondDipole=im.n2*(sgn2*m2);
        if(im.singlet&&!psSingletFieldImage()) {
            secondDipole=im.n1*(-sgn2*m2);
            gSingletSpinTransport=std::getenv("CREM_NO_SINGLET_TRANSPORT")==nullptr;
        }
        if(std::getenv("CREM_DEBUG_BORN"))
            std::fprintf(stderr,"BORN 1S field eigen m=%d P_S=%.6f singlet=%d\n",im.m,im.singletWeight,int(im.singlet));
    }
    // RULE Q BASIN ON THE CONSERVED AXIS AT l >= 1 (audit 544;""")
rep("""    // Osculating elements of the prepared orbit and the closed-form classical
    // prediction they imply.""",r"""    // Audit 619 (b), (c): Born draw of J (Ps triplet P) and of j, F (hydrogen,
    // muonium) at formation without a field; the configuration is turned into the
    // centre of the drawn label's basin (rule Q lengths).
    if(!formationBasins()&&!(gExternalMagneticField.norm()>0.0)&&elements.specificEnergy<0.0
       &&reducedMass*elements.specificAngularMomentum>=hbar*(1.0-1.0e-9)
       &&firstDipole.norm()>0.0&&secondDipole.norm()>0.0) {
        const double orbital=reducedMass*elements.specificAngularMomentum;
        const int l=std::max(1,static_cast<int>(std::lround(orbital/hbar-0.5)));
        const double LQ=l+0.5;
        const Vec3 Lh=angularMomentumDirection;
        const double sgn1=firstCharge<0.0?-1.0:1.0, sgn2=secondCharge<0.0?-1.0:1.0;
        std::uint64_t bornStream=splitMix64(seed^0x619f0a71ULL);
        if(particleAntiparticlePair()&&gSingletChannelEligible
           &&!spinsAntiparallel(firstDipole,secondDipole)) {
            const Vec3 Sh=firstDipole*(sgn1/firstDipole.norm());
            const double th=std::acos(std::clamp(dot(Lh,Sh),-1.0,1.0));
            const std::vector<double> p=formationBornWeights(l,1.0,th);
            const int J=std::abs(l-1)+formationDraw(p,drawUniformUnit(bornStream));
            const double target=std::acos(std::clamp(((J+0.5)*(J+0.5)-LQ*LQ-2.25)/(2.0*LQ*1.5),-1.0,1.0));
            const Vec3 d1=formationTurn(Sh,Lh,target,firstDipole), d2=formationTurn(Sh,Lh,target,secondDipole);
            firstDipole=d1; secondDipole=d2;
            if(std::getenv("CREM_DEBUG_BORN")) std::fprintf(stderr,"BORN formation J=%d (l=%d)\n",J,l);
        } else if(!particleAntiparticlePair()) {
            const bool firstLight=firstMass<secondMass;
            Vec3& light=firstLight?firstDipole:secondDipole; Vec3& heavy=firstLight?secondDipole:firstDipole;
            const double sl=(firstLight?sgn1:sgn2), shv=(firstLight?sgn2:sgn1);
            const Vec3 sHat=light*(sl/light.norm());
            const double th=std::acos(std::clamp(dot(Lh,sHat),-1.0,1.0));
            const int jIndex=formationDraw(formationBornWeights(l,0.5,th),drawUniformUnit(bornStream));
            const int k=l+jIndex; const double j=k-0.5;
            const double target=std::acos(std::clamp((double(k)*k-LQ*LQ-1.0)/(2.0*LQ),-1.0,1.0));
            light=formationTurn(sHat,Lh,target,light);
            const Vec3 sNew=light*(sl/light.norm()); const Vec3 Jv=Lh*LQ+sNew; const Vec3 Jh=Jv*(1.0/Jv.norm());
            const Vec3 iHat=heavy*(shv/heavy.norm());
            const double th2=std::acos(std::clamp(dot(Jh,iHat),-1.0,1.0));
            const int fIndex=formationDraw(formationBornWeights(j,0.5,th2),drawUniformUnit(bornStream));
            const double F=std::abs(j-0.5)+fIndex;
            const double target2=std::acos(std::clamp(((F+0.5)*(F+0.5)-double(k)*k-1.0)/(2.0*k),-1.0,1.0));
            heavy=formationTurn(iHat,Jh,target2,heavy);
            if(std::getenv("CREM_DEBUG_BORN")) std::fprintf(stderr,"BORN formation j=%.1f F=%.1f (l=%d)\n",j,F,l);
        }
    }

    // Osculating elements of the prepared orbit and the closed-form classical
    // prediction they imply.""")
open(p,'w').write(s)
print("patched")
