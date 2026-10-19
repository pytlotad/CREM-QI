#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 513: Sommerfeld's angular quantum number assigned to the TOTAL angular momentum, |J| = k hbar with k = j + 1/2 -- hydrogen n = 2 fine
// structure from the engine's classical relativistic level energy levelFineStructureEnergy(a, L).
int main(){
  applyPairFromOption("proton,electron");
  const double aB=pairBohrRadius(activePair), a=4*aB, h=2*pi*hbar, MHz=1e-6/h;
  const double F1=levelFineStructureEnergy(a,1.0*hbar), F2=levelFineStructureEnergy(a,2.0*hbar), FL=levelFineStructureEnergy(a,1.5*hbar), FS=levelFineStructureEnergy(a,0.5*hbar);
  std::printf("hydrogen n = 2, a = 4 a_pair = %.6e m\n",a);
  std::printf("  F(L = 1 hbar, k = 1: 2p1/2, 2s1/2) = %+.4f MHz\n  F(L = 2 hbar, k = 2: 2p3/2)        = %+.4f MHz\n",F1*MHz,F2*MHz);
  std::printf("  F(L = 1.5 hbar, Langer 2p, no j)   = %+.4f MHz\n  F(L = 0.5 hbar, Langer 2s, contact)= %+.4f MHz\n",FL*MHz,FS*MHz);
  std::printf("  QR contact density at L = 1 hbar: %.3e m^-3 (contact term only below hbar)\n",quiggaRosnerContactDensity(a,1.0*hbar));
  const double meas=9911.200+1057.845;
  std::printf("Sommerfeld k = j + 1/2: 2P3/2 - 2P1/2 = F(2) - F(1) = %.3f MHz; measured 9911.200 + 1057.845 = %.3f MHz; ratio %.6f\n",(F2-F1)*MHz,meas,(F2-F1)*MHz/meas);
  std::printf("  2S1/2 - 2P3/2 with 2s at k = 1: %.3f MHz (measured 9911.200; Lamb shift 1057.845 absent); with engine's Langer 2s (L = 1/2 + contact): %.3f MHz\n",(F2-F1)*MHz,(F2-FS)*MHz);
  // Dirac with reduced mass, leading order: E = -mu c^2 alpha^2/(2n^2) [1 + alpha^2/n^2 (n/(j+1/2) - 3/4)]
  const double al=fineStructureConstant, mu=firstMass*secondMass/(firstMass+secondMass);
  const auto dirac=[&](double k){ return -mu*c*c*al*al/8.0*(al*al/4.0*(2.0/k-0.75)); };
  std::printf("  reference (theory, not measurement): leading Dirac with reduced mass 2P3/2 - 2P1/2 = %.3f MHz\n",(dirac(2)-dirac(1))*MHz);
}
