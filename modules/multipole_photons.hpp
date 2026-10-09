#pragma once

// MULTIPOLE HIERARCHY OF THE PHOTONS (audit 396; the author: "Zaimplementuj
// prawidłową hierarchię typów multipolowych fotonów.  Wszystkie reguły fotonu
// są napisane dla E1 i tylko dla E1").
//
// Each photon has a type, and the type -- not E1 by default -- decides what
// it carries:
//
//   type  j  parity    orbital Delta l          source in the model
//   E1    1  (-1)^1    +-1                      orbiting charges (pair dipole)
//   M1    1  (-1)^2    0 (spin flip)            precessing intrinsic moments
//   E2    2  (-1)^2    0, +-2 (0 -> 0 forbidden) orbiting charges (quadrupole)
//
// E2 exists only for pairs with a charge quadrupole, kappa_2 =
// (q1 m2^2 + q2 m1^2)/M^2 (zero for e+e-: the orbital multipoles of a
// mass-symmetric neutral pair are the odd ones only).  The next rungs, M2 and
// E3, sit at beta^4 relative to E1 (2.8e-9 at the pair Bohr radius) and are
// not carried; see the ladder comment in crem_trajectory.hpp.
//
// E2 of a Kepler orbit.  The traceless quadrupole of the relative coordinate,
// Q_ij = kappa_2 (3 x_i x_j - r^2 delta_ij), has in the orbit plane (z along L)
// only the projections m = +-2 (from zeta^2, zeta = x + i y) and m = 0 (from
// Q_zz = -kappa_2 r^2); m = +-1 vanish for a planar orbit.  The radiated power
//     P_E2 = <sum_ij (d^3 Q_ij/dt^3)^2> / (1440 pi eps0 c^5)
// splits as sum_ij Q_ij^2 = (3/2) Q_zz^2 + (1/2) |3 zeta^2|^2 into those
// projections, and each projection into the harmonics k of the mean anomaly.
// The orbit average has the same closed form as Peters' gravitational
// quadrupole (Peters 1964, with G M -> k/mu):
//     <sum (d^3 Q/dt^3)^2> = 288 kappa_2^2 (k/mu)^3 f(e) / a^5,
//     f(e) = (1 + 73/24 e^2 + 37/96 e^4) / (1 - e^2)^{7/2},
// so P_E2 = kappa_2^2 k^3 f(e) / (5 pi eps0 c^5 mu^3 a^5).  The harmonic table
// below comes from an FFT of the orbit itself and is checked against f(e).
//
// A photon of harmonic k and projection m takes k quanta of action (the
// ladder step of the action rule) and m hbar of orbital angular momentum
// about L.  The Delta l rule then fixes what a projection cannot do: an E2
// photon with m = +2 from l < 2, or m = 0 from l = 0 (s -> s is forbidden),
// raises l by 2 instead (the same correspondence the E1 rule uses below
// |L| = hbar, audit 361).

#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <mutex>
#include <vector>

enum class PhotonMultipole { E1, M1, E2 };

inline int photonMultipoleOrder(PhotonMultipole type) {
    return type==PhotonMultipole::E2?2:1;
}
inline const char* photonMultipoleName(PhotonMultipole type) {
    switch(type) {
    case PhotonMultipole::E1: return "E1";
    case PhotonMultipole::M1: return "M1";
    case PhotonMultipole::E2: return "E2";
    }
    return "?";
}

// Charge quadrupole factor of the relative coordinate (zero for e+e-).
inline double pairQuadrupoleCharge(double q1,double m1,double q2,double m2) {
    const double total=m1+m2;
    return (q1*m2*m2+q2*m1*m1)/(total*total);
}

inline double petersEccentricityFactor(double e) {
    const double e2=e*e;
    return (1.0+73.0/24.0*e2+37.0/96.0*e2*e2)/std::pow(1.0-e2,3.5);
}

// Orbit-averaged E2 power of a Kepler orbit (closed form, see above).
inline double keplerQuadrupolePower(double kappa2,double coulombStrength,
        double reducedMass,double semiMajorAxis,double eccentricity,
        double eps0,double lightSpeed) {
    if(!(semiMajorAxis>0.0)||kappa2==0.0) return 0.0;
    return kappa2*kappa2*std::pow(coulombStrength,3)
        *petersEccentricityFactor(std::min(eccentricity,0.99))
        /(5.0*M_PI*eps0*std::pow(lightSpeed,5)*std::pow(reducedMass,3)
          *std::pow(semiMajorAxis,5));
}

// One harmonic/projection entry of the E2 spectrum and the cached spectrum
// per eccentricity bin.
struct QuadrupoleHarmonic { int k; int m; double power; };
struct QuadrupoleSpectrum {
    std::vector<QuadrupoleHarmonic> entries;
    std::vector<double> countCdf;   // over entries, weights power/k
    double totalPower=0.0;          // units a = 1, omega = 1, kappa = 1
    double countFactor=1.0;         // (sum power/k) / (sum power)
};

namespace multipole_detail {
inline void fft(std::vector<std::complex<double>>& a) {
    const int n=static_cast<int>(a.size());
    for(int i=1,j=0;i<n;++i) {
        int bit=n>>1;
        for(;j&bit;bit>>=1) j^=bit;
        j^=bit;
        if(i<j) std::swap(a[i],a[j]);
    }
    for(int len=2;len<=n;len<<=1) {
        const std::complex<double> w=std::polar(1.0,-2.0*M_PI/len);
        for(int i=0;i<n;i+=len) {
            std::complex<double> v=1.0;
            for(int j=0;j<len/2;++j) {
                const auto x=a[i+j], y=a[i+j+len/2]*v;
                a[i+j]=x+y; a[i+j+len/2]=x-y; v*=w;
            }
        }
    }
}
}

// E2 spectrum of a Kepler orbit, eccentricity from e^2 bins (4096 bins,
// e^2 <= 0.98: at higher e the k^6 tail needs more than the 2^14-point FFT;
// such orbits use the e^2 = 0.98 spectrum, stated in audit 396).
inline const QuadrupoleSpectrum& keplerQuadrupoleSpectrum(double eccentricity) {
    constexpr int bins=4096, samples=1<<14, kmax=samples/4;
    const double e2=std::clamp(eccentricity*eccentricity,0.0,0.98);
    const int bin=static_cast<int>(std::lround(e2*bins));
    static std::mutex cacheMutex;
    static std::map<int,QuadrupoleSpectrum> cache;
    std::lock_guard<std::mutex> lock(cacheMutex);
    auto it=cache.find(bin);
    if(it!=cache.end()) return it->second;
    const double e=std::sqrt(static_cast<double>(bin)/bins);
    const double u=std::sqrt(std::max(0.0,1.0-e*e));
    std::vector<std::complex<double>> radial(samples), rotating(samples);
    double E=0.0;
    for(int j=0;j<samples;++j) {
        const double M=2.0*M_PI*j/samples;
        double Ej=j==0?0.0:E;
        for(int it2=0;it2<80;++it2) {
            const double f=Ej-e*std::sin(Ej)-M, d=1.0-e*std::cos(Ej);
            const double step=f/d; Ej-=step;
            if(std::abs(step)<1e-15) break;
        }
        E=Ej;
        const double x=std::cos(E)-e, y=u*std::sin(E);
        radial[j]=x*x+y*y;                               // r^2 (m = 0)
        rotating[j]=3.0*std::complex<double>(x,y)*std::complex<double>(x,y); // 3 zeta^2
    }
    multipole_detail::fft(radial);
    multipole_detail::fft(rotating);
    QuadrupoleSpectrum s;
    const double norm=1.0/(static_cast<double>(samples)*samples);
    for(int k=1;k<=kmax;++k) {
        const double k6=std::pow(static_cast<double>(k),6);
        // m = 0: real signal, harmonics +k and -k together (factor 2), weight 3/2.
        const double p0=1.5*2.0*std::norm(radial[k])*norm*k6;
        const double pp=0.5*std::norm(rotating[k])*norm*k6;            // m = +2
        const double pm=0.5*std::norm(rotating[samples-k])*norm*k6;    // m = -2
        if(p0>0.0) s.entries.push_back({k,0,p0});
        if(pp>0.0) s.entries.push_back({k,2,pp});
        if(pm>0.0) s.entries.push_back({k,-2,pm});
    }
    // Drop round-off entries (a circle has exactly one physical one, k = 2,
    // m = +2; the rest of its FFT is ~1e-30).
    double peak=0.0;
    for(const auto& h: s.entries) peak=std::max(peak,h.power);
    std::erase_if(s.entries,[peak](const QuadrupoleHarmonic& h) {
        return h.power<1.0e-14*peak; });
    double count=0.0;
    for(const auto& h: s.entries) { s.totalPower+=h.power; count+=h.power/h.k; }
    s.countFactor=s.totalPower>0.0?count/s.totalPower:1.0;
    double running=0.0;
    for(const auto& h: s.entries) { running+=h.power/h.k; s.countCdf.push_back(running/count); }
    return cache.emplace(bin,std::move(s)).first->second;
}

inline QuadrupoleHarmonic drawQuadrupoleHarmonic(double eccentricity,double uniformDraw) {
    const QuadrupoleSpectrum& s=keplerQuadrupoleSpectrum(eccentricity);
    if(s.entries.empty()) return {2,2,0.0};
    const auto pos=std::lower_bound(s.countCdf.begin(),s.countCdf.end(),
                                    std::clamp(uniformDraw,0.0,1.0));
    const auto index=std::min<std::ptrdiff_t>(pos-s.countCdf.begin(),
        static_cast<std::ptrdiff_t>(s.entries.size())-1);
    return s.entries[static_cast<size_t>(index)];
}

// Orbital angular momentum the photon removes, in hbar along L (signed:
// negative raises |L|), after the Delta l rule of its type.  E1 keeps the
// existing rule (axialPhotonSpinSign, passed in).  angularHbar = |L|/hbar.
inline double photonOrbitalTransferHbar(PhotonMultipole type,int projection,
                                        double angularHbar,double e1Sign) {
    switch(type) {
    case PhotonMultipole::E1: return e1Sign;
    case PhotonMultipole::M1: return 0.0;
    case PhotonMultipole::E2: {
        // Langer: |L| = (l + 1/2) hbar.
        const int l=std::max(0,static_cast<int>(std::lround(angularHbar-0.5)));
        if(projection==2) return l>=2?2.0:-2.0;
        if(projection==0) return l==0?-2.0:0.0;
        return -2.0;
    }
    }
    return e1Sign;
}

// cos(theta) about the emission axis for a (j, |m|) photon by inverting the
// angular CDF with ONE uniform (the draw count of the E1 sampler is kept).
//   j = 2, |m| = 2:  (1 - c^4);   j = 2, m = 0:  c^2 (1 - c^2).
inline double sampleQuadrupoleCosTheta(int projection,double uniformDraw) {
    const double u=std::clamp(uniformDraw,0.0,1.0);
    const auto cdf=[projection](double c) {
        if(std::abs(projection)==2)          // int (1 - c^4) from -1, / (8/5)
            return (c-std::pow(c,5)/5.0+0.8)/1.6;
        return (std::pow(c,3)/3.0-std::pow(c,5)/5.0+2.0/15.0)/(4.0/15.0);
    };
    double lo=-1.0, hi=1.0;
    for(int it=0;it<60;++it) {
        const double mid=0.5*(lo+hi);
        if(cdf(mid)<u) lo=mid; else hi=mid;
    }
    return 0.5*(lo+hi);
}

// Unit vector at cos(theta) from a unit axis, azimuth phi.
template<class V>
inline V multipoleDirectionFromAxis(const V& axis,double cosTheta,double phi) {
    const V trial=std::abs(axis.z)<0.9?V{0.0,0.0,1.0}:V{1.0,0.0,0.0};
    V e1=cross(axis,trial);
    e1=e1*(1.0/e1.norm());
    const V e2=cross(axis,e1);
    const double s=std::sqrt(std::max(0.0,1.0-cosTheta*cosTheta));
    return axis*cosTheta+(e1*std::cos(phi)+e2*std::sin(phi))*s;
}
