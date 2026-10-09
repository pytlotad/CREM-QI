#pragma once

// LABORATORY SOURCE OF POSITRONIUM (audit 390, --ps-source).
//
// The model prepares the pair at rest; its centre of mass moves only through
// photon recoil.  What a laboratory measures -- the Doppler width of the
// 511 keV line, the angular correlation of the two photons -- depends on how
// the positronium was moving when it decayed, i.e. on its source.  Three
// sources:
//   rest            the model's own frame (default; photons as emitted);
//   thermal[:T]     Maxwell-Boltzmann at T kelvin (default 300 K) for the
//                   total mass M = m1 + m2, isotropic;
//   beam:<E_eV>     kinetic energy E along +z (relativistic).
// The source velocity is composed relativistically with the recoil velocity
// the cascade left (CremCollapseEstimate::centreOfMassVelocityAtStop), and
// every photon is boosted from the pair's rest frame to the laboratory.  No
// detector response (resolution, efficiency, geometry) is applied: the
// distributions are ideal (README, plan of the distribution redesign).

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>

enum class PsSourceKind { Rest, Thermal, Beam };

struct PsSource {
    PsSourceKind kind=PsSourceKind::Rest;
    double temperatureKelvin=300.0;
    double beamKineticEnergyEv=0.0;
};

inline PsSource gPsSource;

inline PsSource parsePsSource(const std::string& text) {
    PsSource source;
    const auto colon=text.find(':');
    const std::string head=text.substr(0,colon);
    const std::string tail=colon==std::string::npos?"":text.substr(colon+1);
    if(head=="rest"&&tail.empty()) { source.kind=PsSourceKind::Rest; return source; }
    if(head=="thermal") {
        source.kind=PsSourceKind::Thermal;
        if(!tail.empty()) source.temperatureKelvin=std::stod(tail);
        if(!(source.temperatureKelvin>0.0)||!std::isfinite(source.temperatureKelvin))
            throw std::invalid_argument("--ps-source thermal:T needs T > 0 (kelvin)");
        return source;
    }
    if(head=="beam"&&!tail.empty()) {
        source.kind=PsSourceKind::Beam;
        source.beamKineticEnergyEv=std::stod(tail);
        if(!(source.beamKineticEnergyEv>=0.0)||!std::isfinite(source.beamKineticEnergyEv))
            throw std::invalid_argument("--ps-source beam:E needs E >= 0 (eV)");
        return source;
    }
    throw std::invalid_argument(
        "--ps-source must be rest, thermal[:T_K] or beam:<E_kin_eV>");
}

inline std::string describePsSource(const PsSource& source) {
    char buffer[96];
    switch(source.kind) {
    case PsSourceKind::Rest: return "at rest (the model's frame)";
    case PsSourceKind::Thermal:
        std::snprintf(buffer,sizeof buffer,"thermal, T = %g K",source.temperatureKelvin);
        return buffer;
    case PsSourceKind::Beam:
        std::snprintf(buffer,sizeof buffer,"beam along +z, E_kin = %g eV",
                      source.beamKineticEnergyEv);
        return buffer;
    }
    return "?";
}

// Standard normal from the trajectory stream (Box-Muller, one value).
inline double drawPsSourceNormal(std::uint64_t& stream) {
    double u=drawUniformUnit(stream);
    while(!(u>0.0)) u=drawUniformUnit(stream);
    const double v=drawUniformUnit(stream);
    return std::sqrt(-2.0*std::log(u))*std::cos(2.0*pi*v);
}

// Laboratory velocity of the positronium centre of mass from the source.
inline Vec3 drawPsSourceVelocity(const PsSource& source,std::uint64_t& stream) {
    const double totalMass=firstMass+secondMass;
    switch(source.kind) {
    case PsSourceKind::Rest: return Vec3{0.0,0.0,0.0};
    case PsSourceKind::Thermal: {
        // Non-relativistic Maxwell-Boltzmann: each component N(0, kT/M);
        // at 300 K beta ~ 2.7e-6, so the relativistic correction is ~1e-11.
        const double sigma=std::sqrt(boltzmannConstant*source.temperatureKelvin
                                     /totalMass);
        return Vec3{drawPsSourceNormal(stream)*sigma,
                    drawPsSourceNormal(stream)*sigma,
                    drawPsSourceNormal(stream)*sigma};
    }
    case PsSourceKind::Beam: {
        const double restEnergy=totalMass*c*c;
        const double gammaFactor=1.0+source.beamKineticEnergyEv*eCharge/restEnergy;
        const double beta=std::sqrt(std::max(0.0,1.0-1.0/(gammaFactor*gammaFactor)));
        return Vec3{0.0,0.0,beta*c};
    }
    }
    return Vec3{0.0,0.0,0.0};
}

// Relativistic composition: the laboratory velocity of a body moving with v
// in a frame that itself moves with u.
inline Vec3 composePsVelocities(const Vec3& u,const Vec3& v) {
    const double u2=u.squaredNorm();
    if(!(u2>0.0)) return v;
    const double gammaU=1.0/std::sqrt(1.0-u2/(c*c));
    const double uv=dot(u,v);
    const Vec3 vParallel=u*(uv/u2);
    const Vec3 vPerp=v-vParallel;
    return (u+vParallel+vPerp*(1.0/gammaU))*(1.0/(1.0+uv/(c*c)));
}

struct LabPhoton { double energyJoules=0.0; Vec3 direction; };

// Boost of a photon (energy, unit direction in the pair's rest frame) to the
// laboratory where the pair moves with velocity v.
inline LabPhoton boostPhotonToLab(double energyJoules,const Vec3& direction,
                                  const Vec3& velocity) {
    LabPhoton out{energyJoules,direction};
    const double speed=velocity.norm();
    if(!(speed>0.0)) return out;
    const double beta=speed/c, gammaFactor=1.0/std::sqrt(1.0-beta*beta);
    const Vec3 hat=velocity*(1.0/speed);
    const double cosine=dot(direction,hat);
    out.energyJoules=gammaFactor*energyJoules*(1.0+beta*cosine);
    const Vec3 momentum=hat*(gammaFactor*(cosine+beta))+(direction-hat*cosine);
    out.direction=momentum*(1.0/momentum.norm());
    return out;
}
