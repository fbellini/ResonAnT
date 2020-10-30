//#include "/Users/fbellini/alice/macros/ResonAnT/phiXeXe/projectPhiXeXe.C"
#include "/Users/fbellini/alice/macros/ResonAnT/phiXeXe/FitSpectrum.C"
void runFinalAnalysis()
{
  //projectPhiXeXe("20200615_RsnOut_B3.root", "phi", "RsnOut_tpc2sPtDep_tof3sveto5smism", "B3", 0, 0);
  FitSpectrum(-1, "bgbw", 0.5, 10.0, "30oct20");
  FitSpectrum(-1, "boltz", 0.5, 10.0, "30oct20");
  FitSpectrum(-1, "bose", 0.5, 10.0, "30oct20");
  FitSpectrum(-1, "levy", 0.5, 10.0, "30oct20");
  FitSpectrum(-1, "mtexp", 0.5, 10.0, "30oct20");
  return;
}
