#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
// Audit 357: advance the post-landing o-Ps state (356) by the secular transport for 1 ns in 60 ps steps; print |L|.
static Vec3 parse(const std::string& s){ Vec3 v; std::sscanf(s.c_str(),"%lf,%lf,%lf",&v.x,&v.y,&v.z); return v; }
int main(int argc,char** argv){
  std::ifstream f(argv[1]); std::map<std::string,std::string> kv; std::string tok;
  while(f>>tok){ auto p=tok.find('='); if(p!=std::string::npos) kv[tok.substr(0,p)]=tok.substr(p+1); }
  const double a=std::stod(kv["a"]), mred=firstMass*secondMass/(firstMass+secondMass);
  Vec3 Lhat=parse(kv["L"]); Lhat=Lhat/Lhat.norm(); SecularSpinOrbitState st{Lhat*(0.50011*hbar),parse(kv["m1"]),parse(kv["m2"]),0.0,parse(kv["P"])};
  std::printf("start |L| = %.6f hbar\n",st.orbitalAngularMomentum.norm()/hbar);
  for(int k=1;k<=16;++k){
    auto adv=advanceCoupledSecularSpinOrbit(st,a,mred,6.0e-11);
    if(!adv.completed){ std::printf("failed\n"); return 1; }
    st=adv.state;
    if(k%4==0) std::printf("t = %.2f ns  |L| = %.6f hbar  substeps %d\n",k*0.06,st.orbitalAngularMomentum.norm()/hbar,adv.substeps);
  }
}
