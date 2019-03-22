void comparePrelimWithNewEff2Prelim(bool useSys=1,
				    TString final = "$HOME/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/final_sept14/finalCand_kstar_pPb_smoothSys.root",
				    TString prelim = "$HOME/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/preliminary/prelim_kstar_pPb_smoothSys.root")
{
  TString uncert="";

  if (useSys==1) uncert="_sys"; 
  gROOT->LoadMacro("$ASD/GetPlotRatio.C");
  TFile *_file0 = TFile::Open(final.Data());
  if (!_file0) {Printf("Problem with file %s", final.Data()); return;}

  TString prefixfinal = "hKstar"; 
  TString prefixprelim = "hKstar"; 
  
  TH1F * h0p = (TH1F*) _file0->Get(Form("%s_0%s", prefixfinal.Data(), uncert.Data()));
  TH1F * h1p = (TH1F*) _file0->Get(Form("%s_1%s", prefixfinal.Data(), uncert.Data()));
  TH1F * h2p = (TH1F*) _file0->Get(Form("%s_2%s", prefixfinal.Data(), uncert.Data()));
  TH1F * h3p = (TH1F*) _file0->Get(Form("%s_3%s", prefixfinal.Data(), uncert.Data()));
  TH1F * h4p = (TH1F*) _file0->Get(Form("%s_4%s", prefixfinal.Data(), uncert.Data()));

  TFile *_file1 = TFile::Open(prelim.Data());//
  if (!_file1) {Printf("Problem with file %s", prelim.Data()); return;}

  TH1F * h0old = (TH1F*) _file1->Get(Form("%s_0%s",prefixprelim.Data(),uncert.Data()));
  TH1F * h1old = (TH1F*) _file1->Get(Form("%s_1%s",prefixprelim.Data(),uncert.Data()));
  TH1F * h2old = (TH1F*) _file1->Get(Form("%s_2%s",prefixprelim.Data(),uncert.Data()));
  TH1F * h3old = (TH1F*) _file1->Get(Form("%s_3%s",prefixprelim.Data(),uncert.Data()));
  TH1F * h4old = (TH1F*) _file1->Get(Form("%s_4%s",prefixprelim.Data(),uncert.Data()));
  
  TH1D * r0 = GetPlotRatio(h0p, h0old, 1, "MB efficiency","multip. dep. efficiency (prelim.)","yields",0,0,Form("final2prel_0_num%sUncert_01sep14.png",(useSys?"":"Stat")),0.0,14.9,"num"); 
  TH1D * r1 = GetPlotRatio(h1p, h1old, 1, "MB efficiency","multip. dep. efficiency (prelim.)","yields",0,0,Form("final2prel_1_num%sUncert_01sep14.png",(useSys?"":"Stat")),0.0,14.9,"num");
  TH1D * r2 = GetPlotRatio(h2p, h2old, 1, "MB efficiency","multip. dep. efficiency (prelim.)","yields",0,0,Form("final2prel_2_num%sUncert_01sep14.png",(useSys?"":"Stat")),0.0,14.9,"num") ;
  TH1D * r3 = GetPlotRatio(h3p, h3old, 1, "MB efficiency with param.","multip. dep. efficiency (prelim.)","yields",0,0,Form("final2prel_3_num%sUncert_01sep14.png",(useSys?"":"Stat")),0.0,14.9,"num") ;
  TH1D * r4 = GetPlotRatio(h4p, h4old, 1, "MB efficiency with param.","multip. dep. efficiency (prelim.)","yields",0,0,Form("final2prel_4_num%sUncert_01sep14.png",(useSys?"":"Stat")),0.0,5.9,"num") ;
  return;
}
