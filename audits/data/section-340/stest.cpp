#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  const double keys[]={0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.93,0.95,0.97,0.99};
  const double vals[]={1.0,0.995009,0.980140,0.955704,0.922208,0.880335,0.830924,0.774946,0.713472,0.647646,0.578663,0.507742,0.436110,0.364990,0.295604,0.229184,0.167007,0.110485,0.061352,0.036433,0.022231,0.010466,0.000853};
  double worst=0; for(int i=0;i<23;++i) worst=std::max(worst,std::abs(eccentricOrbitHazardSuppression(keys[i])-vals[i]));
  std::printf("worst node deviation %.2e\n",worst);
  for(double e: {0.002,0.005,0.009,0.02,0.035,0.05,0.5,0.95}){ double f=std::pow(1-e*e,1.5)/(1+e*e/2);
    std::printf("e=%.3f  1-S=%.4e  exact 1-f=%.4e  ratio %.3f\n",e,1-eccentricOrbitHazardSuppression(e),1-f,(1-eccentricOrbitHazardSuppression(e))/(1-f)); }
}
