/*
effType: 
match1 for single track matching efficiency
match2 for two-tracks (rsn) matching efficiency
pid for PID efficiency
ratio for ratio of any 2 effs
*/
#define DOPAP 1

Color_t color[5]={kRed+1, kOrange+1, kGreen+2, kBlue+1, kMagenta+1};

Int_t GetEfficiency(Int_t nsigma = 2)
{
  for (Int_t i =0;i<5;i++){ 
    GetEfficiency(i, nsigma);
  }
  return 0;
}
Int_t GetEfficiency(
		    Int_t centbin,
		    Int_t nsigma,
		    TString fpidname = "efficiency_RsnOut_tof2s_centBin00.root", 
		    TString fmatchname = "efficiency_RsnOut_match2011_centBin00.root", 
		    TString fqualname = "efficiency_RsnOut_quality2011_centBin00.root",
		    TString legendOpt = "(0-20%)")
{
  //load macros and set style
  gROOT->LoadMacro("$ASD/GetPlotRatio.C");
  gROOT->LoadMacro("$ASD/SetGraphicStyle.C");
  SetGraphicStyle(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TString effsumname = "hEffVsPt";
  TString particleName = "Kstar";
  
  if (centbin>0){
    fpidname.ReplaceAll("centBin00",Form("centBin0%i",centbin));
    fmatchname.ReplaceAll("centBin00",Form("centBin0%i",centbin));
    fqualname.ReplaceAll("centBin00",Form("centBin0%i",centbin));
    legendOpt.ReplaceAll("(0-20%)",Form("(%i-%i%)", 20*centbin, 20*(centbin+1)));
  }

  if (nsigma!=2){
    fpidname.ReplaceAll("2s",Form("%is",nsigma));
    fmatchname.ReplaceAll("2s",Form("%is",nsigma));
    fqualname.ReplaceAll("2s",Form("%is",nsigma));
  }

  TFile * fpid = TFile::Open(fpidname.Data());
  TFile * fmatch = TFile::Open(fmatchname.Data());
  TFile * fquality = TFile::Open(fqualname.Data());
  
  if (!fpid) { Printf("Error opening PID eff file - check that %s exists.", fpidname.Data()); return 1;}
  else Printf(" ===== PID eff. numerator from %s", fpid->GetName());
  if (!fmatch) { Printf("Error opening matching eff file - check %s that exists.", fmatchname.Data()); return 2;}
  else Printf(" ===== Match eff. numerator = PID eff denominator from %s", fmatch->GetName());
  if (!fquality) { Printf("Error opening TPC eff file - check that %s exists.", fqualname.Data()); return 3;}
  else Printf(" ===== Match eff. denominator from %s", fquality->GetName());
  
  TString effname(effsumname);
  TH1D * hpid = (TH1D*) fpid->Get(effname.Data());
  TH1D * hmatch = (TH1D*) fmatch->Get(effname.Data());
  TH1D * hquality = (TH1D*) fquality->Get(effname.Data());
  
  TH1D * hpurepid = (TH1D*) GetPlotRatio(hpid,hmatch,0,"TPC & match & PID","TPC & match","efficiency");
  hpurepid->SetNameTitle("pid_eff", Form("PID eff. %s",legendOpt.Data()));
  hpurepid->GetYaxis()->SetRangeUser(0.5,1.2);
  hpurepid->SetMarkerStyle(20);
  TH1D * hpurematch = (TH1D*) GetPlotRatio(hmatch,hquality, 0, "TPC & match","TPC", "efficiency");
  hpurematch->SetNameTitle("match_eff", Form("Match eff. %s",legendOpt.Data()));
  hpurematch->GetYaxis()->SetRangeUser(0.,1.2);
  hpurematch->SetMarkerStyle(25);

  TFile * feff = new TFile(Form("ratioEffs_%is_0%i_08aug2013.root",nsigma,centbin),"recreate");
  feff->cd();
  // hpurepid->Write();
  // hpurematch->Write();

#if DOPAP

  //get efficiency for particle
  TH1D * hpid_p = (TH1D*) fpid->Get(Form("%s%s", effname.Data(), particleName.Data()));
  TH1D * hmatch_p = (TH1D*) fmatch->Get(Form("%s%s", effname.Data(), particleName.Data()));
  TH1D * hquality_p = (TH1D*) fquality->Get(Form("%s%s", effname.Data(), particleName.Data()));
  
  TH1D * hpurepid_p = (TH1D*) GetPlotRatio(hpid_p,hmatch_p,0,"TPC & match & PID","TPC & match","efficiency");
  hpurepid_p->SetNameTitle("pid_eff_p", Form("%s: PID eff. %s", particleName.Data(),legendOpt.Data()));
  hpurepid_p->GetYaxis()->SetRangeUser(0.5,1.2);
  hpurepid_p->SetMarkerStyle(20);
  TH1D * hpurematch_p = (TH1D*) GetPlotRatio(hmatch_p,hquality_p, 0, "TPC & match","TPC", "efficiency");
  hpurematch_p->SetNameTitle("match_eff_p", Form("%s: Match eff. %s",particleName.Data(), legendOpt.Data()));
  hpurematch_p->GetYaxis()->SetRangeUser(0.,1.2);
  hpurematch_p->SetMarkerStyle(23);

    //get efficiency for particle
  TH1D * hpid_ap = (TH1D*) fpid->Get(Form("%sAnti%s", effname.Data(), particleName.Data()));
  if (!hpid_ap) Printf("ciccio");
  TH1D * hmatch_ap = (TH1D*) fmatch->Get(Form("%sAnti%s", effname.Data(), particleName.Data()));
  if (!hmatch_ap) Printf("ciccio2");
  TH1D * hquality_ap = (TH1D*) fquality->Get(Form("%sAnti%s", effname.Data(), particleName.Data()));
  if (!hquality_ap) Printf("ciccio3");
  
  TH1D * hpurepid_ap = (TH1D*) GetPlotRatio(hpid_ap,hmatch_ap,0,"TPC & match & PID","TPC & match","efficiency");
  hpurepid_ap->SetNameTitle("pid_eff_ap", Form("Anti%s: %i#sigma PID eff. %s", particleName.Data(), nsigma, legendOpt.Data()));
  hpurepid_ap->GetYaxis()->SetRangeUser(0.5,1.2);
  hpurepid_ap->SetMarkerStyle(24);
  TH1D * hpurematch_ap = (TH1D*) GetPlotRatio(hmatch_ap,hquality_ap, 0, "TPC & match","TPC", "efficiency");
  hpurematch_ap->SetNameTitle("match_eff_ap", Form("Anti%s: Match eff. %s",particleName.Data(), legendOpt.Data()));
  hpurematch_ap->GetYaxis()->SetRangeUser(0.,1.2);
  hpurematch_ap->SetMarkerStyle(27);

  TH1D * hpurepid_ap2p = (TH1D*) GetPlotRatio(hpurepid_ap,hpurepid_p,0,"#epsilon_{PID} anti-particle", "#epsilon_{PID} particle", "efficiency");
  hpurepid_ap2p->SetNameTitle("hpurepid_ap2p",Form("%i#sigma PID eff. anti%s/%s %s", nsigma, particleName.Data(), particleName.Data(), legendOpt.Data()));
  hpurepid_ap2p->GetYaxis()->SetRangeUser(0.7,1.3);
  
  TH1D * hpurematch_ap2p = (TH1D*) GetPlotRatio(hpurematch_ap,hpurematch_p,0,"#epsilon_{match} anti-particle", "#epsilon_{match} particle", "efficiency");
  hpurematch_ap2p->SetNameTitle("hpurematch_ap2p",Form("Match eff. anti%s/%s %s", particleName.Data(), particleName.Data(), legendOpt.Data()));
  hpurematch_ap2p->GetYaxis()->SetRangeUser(0.7,1.3);

  TH1D * hqual_ap2p = (TH1D*) GetPlotRatio(hquality_ap,hquality_p,0,"#epsilon_{TPC} anti-particle", "#epsilon_{TPC} particle", "efficiency");
  hqual_ap2p->SetNameTitle("hqual_ap2p",Form("TPC eff. anti%s/%s %s", particleName.Data(), particleName.Data(), legendOpt.Data()));
  hqual_ap2p->GetYaxis()->SetRangeUser(0.7,1.3);

  hpurepid_p->SetLineColor(color[centbin]);
  hpurepid_p->SetMarkerColor(color[centbin]);
  hpurepid_ap->SetLineColor(color[centbin]);
  hpurepid_ap->SetMarkerColor(color[centbin]);
  hpurepid_ap2p->SetLineColor(color[centbin]);
  hpurepid_ap2p->SetMarkerColor(color[centbin]);  
  hpurematch_p->SetLineColor(color[centbin]);
  hpurematch_p->SetMarkerColor(color[centbin]);
  hpurematch_ap->SetLineColor(color[centbin]);
  hpurematch_ap->SetMarkerColor(color[centbin]);
  hpurematch_ap2p->SetLineColor(color[centbin]);
  hpurematch_ap2p->SetMarkerColor(color[centbin]);
  hqual_ap2p->SetLineColor(color[centbin]);
  hqual_ap2p->SetMarkerColor(color[centbin]);

  //save in file
  feff->cd();
  hpurepid_p->Write();
  hpurepid_ap->Write();
  hpurepid_ap2p->Write();
  hpurematch_p->Write();
  hpurematch_ap->Write();
  hpurematch_ap2p->Write();
  hqual_ap2p->Write();
#endif

  feff->Close();
  return 0;
}
