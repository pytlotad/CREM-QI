#pragma once
// IMPORT 23 (audit 547): FIELD EIGENSTATES OF Ps AT FORMATION IN 2P AND AT THE
// PHOTON TO 1S.  Default with the field rules; CREM_NO_PS_FIELD_EIGEN=1
// restores audit 544 (rule-Q basin on J, classical S/T mixing in 2P),
// CREM_NO_PS_FIELD_RULES=1 removes all field rules.
//
// Why an import.  In a field the 2P eigenstates of Ps change character singlet
// <-> triplet through avoided crossings within each M block, and |up down>
// (P_S = 1/2) has no classical stationary image; neither adiabatic switching
// nor any classical basin reproduces the singlet admixture above ~0.1 T
// (audits 543, 545).  Born's rule at the formation -- an event that resolves
// energy (audit 422) -- does: the coherent state of the classical
// configuration, |L^>_{l=1} (x) |s^ s^> (triplet) or |L^> (x) |S> (singlet),
// projected on the eigenstates of the model's own 2P Hamiltonian
// (levels E_J of import 18, the singlet at 0 by rule B, Zeeman), averages to
// the QM sudden-start result within 0.3 % at 0.01-3 T (audit 546).  The E1
// photon to 1S does not act on the spins, so the spin state after it is
// Tr_L |k><k|; a second Born projection there, on the 1S field eigenstates
// (the model's 1S splitting, Zeeman), and the classical stationary image of
// the chosen state (audit 539: spins along +-B, or the m = 0 state tilted to
// its Breit-Rabi P_S; the mostly-singlet state as the channel singlet) close
// the chain without the early loss of a nutating start (540).
// Numbers: only the engine's levels and moments; the import is the
// eigenstate structure (coherences) of QM.  l = 1 only (12 states); l >= 2
// keeps 544.
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>

struct PsFieldEigenState {
    bool active=false;
    double rho[4][4]={};      // Tr_L |k><k| in |m1 m2> along B^ (up, down)
    Vec3 e1, e2, b;           // the field frame
    double a1=0.0, a2=0.0;    // Zeeman energy per sigma_iz
    std::uint64_t stream=0;   // own random stream
    int k=-1; double singletWeight=0.0, probability=0.0;
    double energy=0.0;        // w_k, eigenvalue of the 2P Hamiltonian (audit 651)
};
inline thread_local PsFieldEigenState gPsFieldEigen;

// S' AS THE SECOND FIXED POINT (audit 559; with the field rules, default;
// CREM_NO_PS_SINGLET_FIELD_IMAGE=1 restores the pure channel singlet with the
// averaged singlet transport, the field part of import 8).  In the s band in a
// field the singlet-like eigenstate S' of {S, T0} is the second smooth fixed
// point of the energy of audit 539 (tan 2 beta = x, the branch beta - pi/2),
// w = P_S(S') exactly (audit 558); the pair then follows the rates of 541 like
// the triplet.  A pure singlet decays at Gamma_p, 1.8 % (1 T) and 13 % (3 T)
// too fast, and its survival after 0.5 ns was 7 % / 37 % low (audit 558).
inline bool psSingletFieldImage() {
    static const bool on=std::getenv("CREM_NO_PS_SINGLET_FIELD_IMAGE")==nullptr
        &&std::getenv("CREM_NO_PS_FIELD_RULES")==nullptr;
    return on;
}
inline bool psFieldEigenEnabled() {
    static const bool on=std::getenv("CREM_NO_PS_FIELD_EIGEN")==nullptr
        &&std::getenv("CREM_NO_PS_FIELD_RULES")==nullptr;
    return on;
}
// PHOTON ENERGY FROM THE EIGENVALUES (audit 651; default, CREM_PS_CLASSICAL_FIELD_PHOTON=1
// restores the classical spin balance of audits 378-650).
// The photon that closes the chain 2P k -> 1S q carries the spin part
// [U(2^1P_1) + w_k] - [U(1^3S) - Delta E/4 + w_q] instead of the classical
// U_before - U_after of audit 378: w_k, w_q are the eigenvalues of the two
// Hamiltonians the labels were drawn from (Zeeman included), U the engine's
// level energy with spins along B, antiparallel (the zero of the 2P
// Hamiltonian, rule B) / parallel (1^3S = zero + Delta E/4).  Gap, anchors
// and Delta E are taken on the rung (n, L = (l + 1/2) hbar), not on the orbit
// the classical spin-orbit exchange has pushed off it.  The orbit lands on the
// ladder; the classical spin energy still drops by the classical jump, the
// difference is the energy of the import.
inline bool psEigenPhotonEnergy() {
    static const bool on=std::getenv("CREM_PS_CLASSICAL_FIELD_PHOTON")==nullptr;
    return on&&psFieldEigenEnabled();
}

// Cyclic Jacobi for a real symmetric n x n matrix (row-major): eigenvalues in
// w, eigenvectors in the COLUMNS of V.
template<int n>
inline void psJacobiEigen(std::array<double,n*n> A,std::array<double,n>& w,
                          std::array<double,n*n>& V) {
    V.fill(0.0);
    for(int i=0;i<n;++i) V[i*n+i]=1.0;
    for(int sweep=0;sweep<100;++sweep) {
        double off=0.0, scale=0.0;
        for(int i=0;i<n;++i) for(int j=0;j<n;++j) {
            if(i!=j) off+=A[i*n+j]*A[i*n+j];
            scale+=A[i*n+j]*A[i*n+j];
        }
        if(!(off>1.0e-30*scale)) break;
        for(int p=0;p<n-1;++p) for(int q=p+1;q<n;++q) {
            const double apq=A[p*n+q];
            if(std::abs(apq)<=1.0e-300) continue;
            const double theta=(A[q*n+q]-A[p*n+p])/(2.0*apq);
            const double t=(theta>=0.0?1.0:-1.0)
                /(std::abs(theta)+std::sqrt(theta*theta+1.0));
            const double cs=1.0/std::sqrt(t*t+1.0), sn=t*cs;
            for(int k=0;k<n;++k) {
                const double akp=A[k*n+p], akq=A[k*n+q];
                A[k*n+p]=cs*akp-sn*akq; A[k*n+q]=sn*akp+cs*akq;
            }
            for(int k=0;k<n;++k) {
                const double apk=A[p*n+k], aqk=A[q*n+k];
                A[p*n+k]=cs*apk-sn*aqk; A[q*n+k]=sn*apk+cs*aqk;
            }
            for(int k=0;k<n;++k) {
                const double vkp=V[k*n+p], vkq=V[k*n+q];
                V[k*n+p]=cs*vkp-sn*vkq; V[k*n+q]=sn*vkp+cs*vkq;
            }
        }
    }
    for(int i=0;i<n;++i) w[i]=A[i*n+i];
}

// 2P Hamiltonian (l = 1) in |m_L>|m1>|m2>, index 4 i_L + 2 i1 + i2, i_L = 0, 1, 2
// for m_L = +1, 0, -1, i = 0 (up), 1 (down); all operators real (L_y = -i D).
inline std::array<double,144> psTwoPHamiltonian(const double E[3],double a1,
                                                double a2) {
    using M3=std::array<double,9>; using M2=std::array<double,4>;
    using M12=std::array<double,144>;
    const double r=std::sqrt(0.5);
    const M3 lx{0,r,0, r,0,r, 0,r,0}, ld{0,r,0, -r,0,r, 0,-r,0}, lz{1,0,0, 0,0,0, 0,0,-1};
    const M3 l1{1,0,0, 0,1,0, 0,0,1};
    const M2 sx{0,0.5,0.5,0}, sd{0,0.5,-0.5,0}, sz{0.5,0,0,-0.5}, s1{1,0,0,1};
    const auto kron=[](const M3& a,const M2& b,const M2& c) {
        M12 m{};
        for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        for(int p=0;p<2;++p) for(int q=0;q<2;++q)
        for(int u=0;u<2;++u) for(int v=0;v<2;++v)
            m[(4*i+2*p+u)*12+(4*j+2*q+v)]=a[3*i+j]*b[2*p+q]*c[2*u+v];
        return m;
    };
    const auto add=[](M12 a,const M12& b,double f) {
        for(int i=0;i<144;++i) a[i]+=f*b[i];
        return a;
    };
    const auto mul=[](const M12& a,const M12& b) {
        M12 m{};
        for(int i=0;i<12;++i) for(int k=0;k<12;++k) {
            const double aik=a[i*12+k];
            if(aik==0.0) continue;
            for(int j=0;j<12;++j) m[i*12+j]+=aik*b[k*12+j];
        }
        return m;
    };
    M12 I{}; for(int i=0;i<12;++i) I[i*12+i]=1.0;
    const M12 S1x=kron(l1,sx,s1), S2x=kron(l1,s1,sx), S1d=kron(l1,sd,s1),
        S2d=kron(l1,s1,sd), S1z=kron(l1,sz,s1), S2z=kron(l1,s1,sz);
    const M12 Sx=add(S1x,S2x,1.0), Sd=add(S1d,S2d,1.0), Sz=add(S1z,S2z,1.0);
    const M12 Jx=add(kron(lx,s1,s1),Sx,1.0), Jd=add(kron(ld,s1,s1),Sd,1.0),
        Jz=add(kron(lz,s1,s1),Sz,1.0);
    const M12 S2=add(add(mul(Sx,Sx),mul(Sd,Sd),-1.0),mul(Sz,Sz),1.0);
    const M12 J2=add(add(mul(Jx,Jx),mul(Jd,Jd),-1.0),mul(Jz,Jz),1.0);
    M12 PT=S2; for(double& x:PT) x*=0.5;
    M12 H{};
    for(int J=0;J<=2;++J) {
        M12 P=PT;
        for(int K=0;K<=2;++K) if(K!=J) {
            M12 f=add(J2,I,-K*(K+1.0));
            for(double& x:f) x/=(J*(J+1.0)-K*(K+1.0));
            P=mul(P,f);
        }
        H=add(H,P,E[J]);
    }
    H=add(H,S1z,2.0*a1);
    H=add(H,S2z,2.0*a2);
    return H;
}

// Coherent states in the field frame: l = 1 along u (m = +1, 0, -1) and spin
// 1/2 up along u (up, down).
inline void psCoherentL1(const Vec3& u,std::complex<double> out[3]) {
    const double ct=std::clamp(u.z,-1.0,1.0), ph=std::atan2(u.y,u.x);
    const std::complex<double> e(std::cos(ph),-std::sin(ph));
    out[0]=e*(1.0+ct)*0.5; out[1]=std::sqrt(std::max(0.0,1.0-ct*ct))/std::sqrt(2.0);
    out[2]=std::conj(e)*(1.0-ct)*0.5;
}
inline void psCoherentHalf(const Vec3& u,std::complex<double> out[2]) {
    const double th=std::acos(std::clamp(u.z,-1.0,1.0)), ph=std::atan2(u.y,u.x);
    out[0]=std::polar(std::cos(0.5*th),-0.5*ph);
    out[1]=std::polar(std::sin(0.5*th),0.5*ph);
}
inline void psFieldFrame(const Vec3& field,Vec3& e1,Vec3& e2,Vec3& b) {
    b=field*(1.0/field.norm());
    Vec3 t=std::abs(b.x)<0.9?Vec3{1.0,0.0,0.0}:Vec3{0.0,1.0,0.0};
    e1=t-b*dot(t,b); e1=e1*(1.0/e1.norm()); e2=cross(b,e1);
}

// Formation: project the classical configuration on the 2P field eigenstates,
// draw k with u, keep Tr_L |k><k|.  spinHat is the triplet's spin axis;
// singlet chooses |S>.  a_i: Zeeman energy per sigma_iz.
inline void psFieldEigenAtFormation(const double E[3],double a1,double a2,
        const Vec3& field,const Vec3& orbitalHat,const Vec3& spinHat,
        bool singlet,double u,std::uint64_t stream) {
    PsFieldEigenState st;
    psFieldFrame(field,st.e1,st.e2,st.b);
    st.a1=a1; st.a2=a2; st.stream=stream;
    const auto inFrame=[&](const Vec3& v) {
        return Vec3{dot(v,st.e1),dot(v,st.e2),dot(v,st.b)};
    };
    std::complex<double> lc[3], sc[2], psi[12];
    psCoherentL1(inFrame(orbitalHat),lc);
    psCoherentHalf(inFrame(spinHat),sc);
    const double r=std::sqrt(0.5);
    for(int i=0;i<3;++i) for(int p=0;p<2;++p) for(int q=0;q<2;++q)
        psi[4*i+2*p+q]=singlet
            ?lc[i]*((p==0&&q==1)?r:(p==1&&q==0)?-r:0.0)
            :lc[i]*sc[p]*sc[q];
    std::array<double,12> w; std::array<double,144> V;
    psJacobiEigen<12>(psTwoPHamiltonian(E,a1,a2),w,V);
    double prob[12], total=0.0;
    for(int k=0;k<12;++k) {
        std::complex<double> o=0.0;
        for(int j=0;j<12;++j) o+=V[j*12+k]*psi[j];
        prob[k]=std::norm(o); total+=prob[k];
    }
    int k=11; double acc=0.0;
    for(int i=0;i<12;++i) { acc+=prob[i]/total; if(u<acc) { k=i; break; } }
    for(int a=0;a<4;++a) for(int c=0;c<4;++c) {
        double s=0.0;
        for(int i=0;i<3;++i) s+=V[(4*i+a)*12+k]*V[(4*i+c)*12+k];
        st.rho[a][c]=s;
    }
    st.k=k; st.probability=prob[k]/total; st.energy=w[k];
    st.singletWeight=0.5*(st.rho[1][1]+st.rho[2][2]-st.rho[1][2]-st.rho[2][1]);
    st.active=true;
    gPsFieldEigen=st;
}

// The photon to l = 0: the 1S field eigenstates of H = (Delta E/4) sigma1.sigma2
// + a1 sigma1z + a2 sigma2z, Born weights <i|rho|i>, the classical stationary
// image of the drawn state as unit spin directions; singlet = the channel
// singlet.  u: uniform for the draw, psi: azimuth in the plane _|_ B.
struct PsFieldEigenImage { Vec3 n1, n2; bool singlet=false; int m=0; double singletWeight=0.0, q=0.0, energy=0.0; };
inline PsFieldEigenImage psFieldEigenAtPhoton(double deltaE,double u,double psi) {
    const PsFieldEigenState& st=gPsFieldEigen;
    const double K=0.25*deltaE;
    std::array<double,16> H{K+st.a1+st.a2,0,0,0,
                            0,-K+st.a1-st.a2,2*K,0,
                            0,2*K,-K-st.a1+st.a2,0,
                            0,0,0,K-st.a1-st.a2};
    std::array<double,4> w; std::array<double,16> V;
    psJacobiEigen<4>(H,w,V);
    double q[4], total=0.0;
    for(int i=0;i<4;++i) {
        double s=0.0;
        for(int a=0;a<4;++a) for(int c=0;c<4;++c) s+=V[a*4+i]*st.rho[a][c]*V[c*4+i];
        q[i]=std::max(0.0,s); total+=q[i];
    }
    int pick=3; double acc=0.0;
    for(int i=0;i<4;++i) { acc+=q[i]/total; if(u<acc) { pick=i; break; } }
    const double v0=V[0*4+pick], v1=V[1*4+pick], v2=V[2*4+pick], v3=V[3*4+pick];
    PsFieldEigenImage im;
    im.q=q[pick]/total; im.energy=w[pick];
    im.singletWeight=0.5*(v1-v2)*(v1-v2);
    const Vec3 perp=st.e1*std::cos(psi)+st.e2*std::sin(psi);
    if(v0*v0>0.5) { im.m=1; im.n1=st.b; im.n2=st.b; }
    else if(v3*v3>0.5) { im.m=-1; im.n1=st.b*-1.0; im.n2=st.b*-1.0; }
    else if(im.singletWeight>=0.5&&!psSingletFieldImage()) { im.singlet=true; im.n1=perp; im.n2=perp*-1.0; }
    else {
        // m = 0: the stationary image of audit 539 for T0' and, since audit 559,
        // also for S' -- the second fixed point of the same energy (audit 558:
        // w = P_S of S' to six digits; the transport with the rules of 541
        // reproduces its QM survival to 1e-4 over the first nanoseconds).
        im.singlet=im.singletWeight>=0.5;
        const double sigma1=v0*v0+v1*v1-v2*v2-v3*v3;
        const double sinPhi=std::sqrt(im.singletWeight), cosPhi=std::sqrt(1.0-im.singletWeight);
        const double sgn=sigma1>=0.0?1.0:-1.0;
        im.n1=perp*cosPhi+st.b*(sgn*sinPhi);
        im.n2=perp*cosPhi-st.b*(sgn*sinPhi);
    }
    return im;
}

// The eigenvalue of H_1S = K sigma1.sigma2 + a1 sigma1z + a2 sigma2z for the
// image the photon landed on (audit 651): m = +-1 -> K +- (a1 + a2); m = 0 ->
// -K -+ sqrt(4K^2 + (a1 - a2)^2), the lower root for S' (singlet weight >= 1/2).
// T0' - S' = Delta E sqrt(1 + x^2), Breit-Rabi with Delta E = 4K.
inline double psOneSEigenvalue(double K,double a1,double a2,int m,bool singletLike) {
    if(m==1) return K+a1+a2;
    if(m==-1) return K-a1-a2;
    return -K+(singletLike?-1.0:1.0)*std::sqrt(4.0*K*K+(a1-a2)*(a1-a2));
}

// The S' image for a singlet formed in the s band in a field (audit 559):
// the m = 0 eigenstate of H_1S = (Delta E/4) sigma1.sigma2 + a1 sigma1z +
// a2 sigma2z with the larger singlet weight, as in psFieldEigenAtPhoton.
inline PsFieldEigenImage psSingletImageInField(double deltaE,double a1,double a2,
        const Vec3& field,double psi) {
    PsFieldEigenState saved=gPsFieldEigen;
    PsFieldEigenState st;
    psFieldFrame(field,st.e1,st.e2,st.b);
    st.a1=a1; st.a2=a2;
    const double r=std::sqrt(0.5);   // rho of the pure singlet: its Born weights pick S' (and T0')
    const double v[4]={0.0,r,-r,0.0};
    for(int a=0;a<4;++a) for(int c=0;c<4;++c) st.rho[a][c]=v[a]*v[c];
    gPsFieldEigen=st;
    // pick S' deterministically: u just below its cumulative weight is not
    // needed -- take the state of the larger singlet weight among m = 0.
    PsFieldEigenImage best; best.singletWeight=-1.0;
    for(int trial=0;trial<2;++trial) {
        const PsFieldEigenImage im=psFieldEigenAtPhoton(deltaE,trial==0?0.0:0.999999,psi);
        if(im.m==0&&im.singletWeight>best.singletWeight) best=im;
    }
    gPsFieldEigen=saved;
    return best;
}
