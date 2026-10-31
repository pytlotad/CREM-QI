// Audit 663 (post hoc, uzupelnienie 2): offline replay of the delta stream splitMix64(splitMix64(623000+i) ^ 0x663 ^ photon count 0).
#include "crem_collapse.hpp"
#include <cstdio>
int main(){ for(unsigned c=0;c<4;++c){ std::uint64_t st=splitMix64(splitMix64(623000ULL)^0x663ULL^c); std::printf("count %u x %.6f\n",c,0.5*std::tan(M_PI*(drawUniformUnit(st)-0.5)));}
 int k=0; for(int sd=0;sd<64;++sd){ std::uint64_t st=splitMix64(splitMix64(623000ULL+sd)^0x663ULL^0ULL); double x=0.5*std::tan(M_PI*(drawUniformUnit(st)-0.5)); k+=std::abs(x)<0.5;} std::printf("count0 frac64 %.4f\n",k/64.0);}
