#pragma once

#include <deque>

#include "vector3.hpp"
#include "dipole_tensor.hpp"

namespace positronium::objects {

// Even-l real spherical-harmonic moments of the radiated angular power
// pattern dP/dOmega, in the FIXED laboratory basis, to l=2.
//
// These are what a photon's direction has to be drawn from.  Production
// draws it instead from a prescribed (1+cos^2 theta) about the orbital
// normal with a UNIFORM azimuth (crem_collapse.hpp:4750-4803, identified
// in audit 256a).  257c measured the instantaneous pattern carrying an
// m=2 azimuthal harmonic at 0.2475 of the power, which the uniform draw
// discards entirely; 258a then measured that the same harmonic averages
// down to 5e-03 over a FULL orbit, so what actually survives into an
// emitted photon is the fractional-window remainder, a few percent.
//
// l=0 and l=2 together span the whole E1 pattern, azimuthal structure
// included: (1+cos^2) is a00 plus a2[2], and the m=2 harmonic is
// a2[0] and a2[4].  ODD l is deliberately absent.  It carries the
// front-back asymmetry and hence the net momentum, and 256d measured
// that channel as capped by the history grid and swinging through
// order-unity angles as the integrator tolerance tightens -- including
// it would import that noise straight into the photon direction.
// Lebedev 302 integrates to degree 29, so the projection onto l<=2 is
// exact rather than a fit.
//
// Basis functions, in a2[] order: xy, yz, (3z^2-1)/2, xz, (x^2-y^2)/2.
// Their squared norms over the sphere are 4pi/15 except the third,
// which is 4pi/5; a00's is 4pi.  a00 must equal the total power
// identically, since its basis function is 1 -- a free invariant, and
// the first thing to check if anything here looks wrong.
struct AngularPatternMoments {
    double a00=0.0;
    double a2[5]={0.0,0.0,0.0,0.0,0.0};
};

struct State {
    Vec3 firstPosition,secondPosition;
    Vec3 firstVelocity,secondVelocity;
    Vec3 firstAcceleration,secondAcceleration;
    Vec3 firstDipole,secondDipole;
    Vec3 firstElectricDipole,secondElectricDipole;
    // Magnetic moments in the instantaneous particle rest frames.  These are
    // the independent degrees of freedom; the laboratory p and mu above are
    // reconstructed from the covariant tensor after every update.
    Vec3 firstProperDipole,secondProperDipole;
    double time=0.0;
    double radiatedEnergy=0.0;
    // Same accumulation as radiatedEnergy, but the E1 (charge-charge,
    // orbital) far-zone flux only, with the M1 (magnetic-dipole/spin-
    // precession) channel excluded rather than merged in.  Exists because
    // radiatedEnergy's merge is exactly what makes it useless as an
    // orbital-decay diagnostic near/below the Compton barrier: M1 dominates
    // it by orders of magnitude there (does not recoil the orbit at all),
    // while orbitalRadiatedEnergy is computed from the exact retarded far-
    // zone Poynting integral and does not share conservativeParticleEnergy's
    // Darwin/near-field approximation, so it stays meaningful exactly where
    // that approximation (period comparable to the light-crossing time)
    // breaks down.
    //
    // NOT diagnostic-only, whatever this comment used to say: it is THE
    // production energy channel.  crem_collapse.hpp's secular estimator sets
    // deltaEnergyPerOrbit = -finalState.orbitalRadiatedEnergy/reducedMass,
    // so the M1 exclusion above is not a bookkeeping nicety -- it is the
    // reason the magnetic-dipole channel never tightens the orbit through
    // the secular path in ANY reaction model.
    //
    // One exception, and it is the only route M1 has into the dynamics: in
    // the quantized mode crem_trajectory.hpp's quantizedPower sums E1, M1
    // and E2 before dividing by hbar*omega, so M1 does raise the in-orbit
    // photon hazard there.  "Does not recoil the orbit at all" above is
    // therefore true of the secular ledger and false of that hazard.
    double orbitalRadiatedEnergy=0.0;
    double dipoleConstraintEnergy=0.0;
    // Accumulated phase of the orbit-following zero-point field, integral of
    // the osculating orbital angular frequency.  Carried in the state because
    // a band that tracks the orbit has a CHANGING mode frequency, and a phase
    // written as -omega(t)*t would jump whenever omega moved, injecting energy
    // out of nowhere.  Integrating d(phase)/dt = omega keeps it continuous.
    // Zero and unused unless --zpf is on.
    double zeroPointPhase=0.0;
    Vec3 radiatedMomentum,radiatedAngularMomentum;
    // Accumulated angular pattern of the radiated energy, same trapezoidal
    // accumulation as orbitalRadiatedEnergy beside it.  Accumulated in the
    // LAB basis, not one tied to the instantaneous separation: 258c measured
    // that a separation-locked accumulation keeps 0.2449 and is simply wrong,
    // because it would hand the photon a structure the emitted energy does
    // not have.
    AngularPatternMoments radiatedPattern;
    // Reconstructed bound/interference field reservoir required to close the
    // particle-plus-field conservation laws on the control surface.
    double boundFieldEnergy=0.0;
    Vec3 boundFieldMomentum,boundFieldAngularMomentum;
    // Independent integral of the mismatch between individual LL reaction
    // and the coherent charge-sector far-field flux.
    double reactionEnergyMismatch=0.0;
    // Previous COMMITTED step length and the rates sampled at its start.
    // They exist only to give the radiation bookkeeping trapezoidal weights,
    // and they live in the State so a rejected trial step discards them with
    // everything else.
    double previousStepDt=0.0;
    bool hasPreviousRates=false;
    double previousFluxEnergy=0.0,previousDipoleFluxEnergy=0.0;
    AngularPatternMoments previousPatternRate;
    Vec3 previousFluxMomentum,previousFluxAngularMomentum;
    double previousMismatchEnergy=0.0;
    Vec3 previousMismatchMomentum,previousMismatchAngularMomentum;
    Vec3 reactionMomentumMismatch,reactionAngularMomentumMismatch;
};

using StateHistory=std::deque<State>;

} // namespace positronium::objects
