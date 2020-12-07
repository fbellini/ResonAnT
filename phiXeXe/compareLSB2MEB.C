#include "/Users/fbellini/alice/macros/utils/GetPlotRatio.C"

void compareLSB2MEB(Int_t centrality=-1)
{

  TString lsbName = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/final/analysis/phifinal_default_LowBdca/norm1.07-1.10/fit_Like_VOIGTpoly2_fixW/fit_r0.992-1.100/RAW_fitResult.root";
  //"/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiA3_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Like_VOIGTpoly1_fixW/fit_r0.998-1.046/RAW_fitResult.root";

  TString mebName = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/final/analysis/phifinal_default_LowBdca/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070/RAW_fitResult.root";
  //"/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiA3_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070/RAW_fitResult.root";
  //"/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0221/phiC3_tpc2s_tof3sveto/norm1.05-1.10/fit_Mixing_VOIGTpoly1_FixRes/fit_r0.990-1.070/RAW_fitResult.root";

  
  TFile * finlsb = TFile::Open(lsbName.Data());
  if (!finlsb) return;

  TFile * finmeb = TFile::Open(mebName.Data());
  if (!finmeb) return;

  TH1D * hstat[5][2];
  Int_t centEdges[] = {0, 10, 30, 50, 70, 90};
  
  for (int ic =0; ic<5; ic++){
    if (centrality>=0 && ic!=centrality) continue;
    hstat[ic][0] = (TH1D *) finlsb->Get(Form("hRawYieldVsPt_%i", ic));
    if (!hstat[ic][0]) return;
    hstat[ic][1] = (TH1D *) finmeb->Get(Form("hRawYieldVsPt_%i", ic));
    if (!hstat[ic][1]) return;
    GetPlotRatio(hstat[ic][0], hstat[ic][1], 1, Form("LSB2MEB_c%i.eps", ic),"LSB", "MEB", Form("%i-%i%%", centEdges[ic], centEdges[ic+1]), "d^{2}#it{N}_{raw}/(d#it{y}d#it{p}_{T})");
  }

}
