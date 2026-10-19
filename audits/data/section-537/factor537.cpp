#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){ applyPairFromOption("electron,positron"); gSingletChannelEligible=true; std::printf("Ps contact torque factor f = %.6f (QM 7/4)\n",sStateContactTorqueFactor());
  applyPairFromOption("proton,electron"); std::printf("H: f = %.6f (not a particle-antiparticle pair)\n",sStateContactTorqueFactor()); }
