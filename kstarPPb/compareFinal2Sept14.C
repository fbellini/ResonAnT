void compareFinal2Sept14(bool useSys=1,		      
			 Bool_t compare2prelim = 0,
			 TString final = "$HOME/alice/resonances/kstar_pA5.02TeV/output_LF5455/multi/finalCand_kstar_pPb_smoothSys.root")
{
  TString prelim;
  if (compare2prelim) prelim = "$HOME/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/preliminary/prelim_kstar_pPb_smoothSys.root";   /*PRELIMINARY QM 2014*/
  else prelim = "$HOME/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/final_sept14/finalCand_kstar_pPb_smoothSys.root"; /* SEPT. 2014*/
  
  TString uncert="";

  if (useSys==1) uncert="_sys"; 
  gROOT->LoadMacro("$ASD/GetPlotRatio.C");
  //TFile *_file0 = TFile::Open("finalCand_kstar_pPb_smoothSys.root");
  TFile *_file0 = TFile::Open(final.Data());
  TString prefixfinal = "hKstar";
  TH1D * h100p = (TH1D*) _file0->Get(Form("%s_100%s", prefixfinal.Data(), uncert.Data()));
  TH1D * h0p = (TH1D*) _file0->Get(Form("%s_0%s", prefixfinal.Data(), uncert.Data()));
  TH1D * h1p = (TH1D*) _file0->Get(Form("%s_1%s", prefixfinal.Data(), uncert.Data()));
  TH1D * h2p = (TH1D*) _file0->Get(Form("%s_2%s", prefixfinal.Data(), uncert.Data()));
  TH1D * h3p = (TH1D*) _file0->Get(Form("%s_3%s", prefixfinal.Data(), uncert.Data()));
  TH1D * h4p = (TH1D*) _file0->Get(Form("%s_4%s", prefixfinal.Data(), uncert.Data()));

  //TFile *_file1 = TFile::Open("../preliminary/prelim_kstar_pPb_smoothSys.root");
  TFile *_file1 = TFile::Open(prelim.Data());//
  TH1D * h100old = (TH1D*) _file1->Get(Form("hKstar_100%s",uncert.Data()));
  if (!compare2prelim) { h100old->SetLineColor(kCyan+2); h100old->SetMarkerColor(kCyan+2);}
  else {h100old->SetLineColor(kMagenta+1); h100old->SetMarkerColor(kMagenta+1);}
  TH1D * h0old = (TH1D*) _file1->Get(Form("%s_0%s",prefixfinal.Data(),uncert.Data()));
  TH1D * h1old = (TH1D*) _file1->Get(Form("%s_1%s",prefixfinal.Data(),uncert.Data()));
  TH1D * h2old = (TH1D*) _file1->Get(Form("%s_2%s",prefixfinal.Data(),uncert.Data()));
  TH1D * h3old = (TH1D*) _file1->Get(Form("%s_3%s",prefixfinal.Data(),uncert.Data()));
  TH1D * h4old = (TH1D*) _file1->Get(Form("%s_4%s",prefixfinal.Data(),uncert.Data()));
  
  TString oldname = (compare2prelim? "preliminary" : "sept14" );
  TString ytitle = "1/N_{evt}d^{2}#it{N}/(d#it{p}_{T}d#it{y} [(GeV/c)^{-2}]";
  TH1D * r100 = GetPlotRatio(h100p, h100old, 1, "final", oldname.Data(), ytitle.Data(),0,0, Form("final2%s_100_num%sUncert_08jul15.png",(compare2prelim?"prelim":"prev"),(useSys?"":"Stat")),0.0,14.9,"num"); 
  TH1D * r0 = GetPlotRatio(h0p, h0old, 1, "final", oldname.Data(), ytitle.Data(),0,0, Form("final2%s_0_num%sUncert_08jul15.png",(compare2prelim?"prelim":"prev"), (useSys?"":"Stat")),0.0,14.9,"num"); 
  TH1D * r1 = GetPlotRatio(h1p, h1old, 1, "final", oldname.Data(), ytitle.Data(),0,0, Form("final2%s_1_num%sUncert_08jul15.png",(compare2prelim?"prelim":"prev"),(useSys?"":"Stat")),0.0,14.9,"num");
  TH1D * r2 = GetPlotRatio(h2p, h2old, 1, "final", oldname.Data(), ytitle.Data(),0,0, Form("final2%s_2_num%sUncert_08jul15.png",(compare2prelim?"prelim":"prev"),(useSys?"":"Stat")),0.0,14.9,"num") ;
  TH1D * r3 = GetPlotRatio(h3p, h3old, 1, "final", oldname.Data(), ytitle.Data(),0,0, Form("final2%s_3_num%sUncert_08jul15.png",(compare2prelim?"prelim":"prev"),(useSys?"":"Stat")),0.0,14.9,"num") ;
  TH1D * r4 = GetPlotRatio(h4p, h4old, 1, "final", oldname.Data(), ytitle.Data(),0,0, Form("final2%s_4_num%sUncert_08jul15.png", (compare2prelim?"prelim":"prev"),(useSys?"":"Stat")),0.0,5.9,"num") ;
  return;
}
