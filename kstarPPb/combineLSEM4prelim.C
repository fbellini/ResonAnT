void combineLSEM4final(Bool_t useSmoothed=0)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  Color_t color[8]={ kRed+1, kPink+6, kGreen+1, kAzure+1, kBlue+3, kBlack, kMagenta, kOrange};  
  Color_t colorsys[9]={ kRed, kBlue+1, kGreen+1, kMagenta+1, kOrange, kCyan+2, kPink+2, kGray, kViolet+1};  
  
  //get bins
  TString filebins = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/proj_2424_tpc2s_tof3sveto.root";
  Int_t npt_axis = 0, ncent_axis=0; 
  TFile *f=TFile::Open(filebins.Data());
  if (!ptbins || !centbins) return;
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  npt_axis = ptbins->GetNbins();  
  TAxis *centbins = (TAxis*)f->Get("centbins");
  ncent_axis = centbins->GetNbins();
  f->Close();
  const Int_t npt = npt_axis;
  const Int_t ncent = ncent_axis;
  Double_t pt[npt+1];
  for (Int_t k=0; k<npt+1;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
  }
  Double_t cent[ncent+1]; 
  for (Int_t k=0; k<ncent+1;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
  }
 
  TCanvas * ctmp = new TCanvas("ctmp","tmp", 800, 800); 

  //define histogram names
  TString histCorrName = "hCorrected";
  
  //access files with min bias spectra EM
  TFile * fin100EM;
  if (useSmoothed) 
    fin100EM = TFile::Open("/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/0to100_binB/preliminary/finalWsyst_smooth2_EM_28apr14_0.root");
  else 
    fin100EM = TFile::Open("/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/0to100_binB/preliminary/finalWsyst_EM_28apr14_0.root");
  if (!fin100EM) Printf("cannot open file min bias EM.");  
  TH1D * hCorrMinBiasEM = (TH1D*) fin100EM->Get(Form("%s_0", histCorrName.Data()));
  if (!hCorrMinBiasEM) Printf("cannot find min bias spectrum 0-100 EM");
  ApplyCosmetics2Spectrum(hCorrMinBiasEM, color[5], 20); 
  hCorrMinBiasEM->SetTitle(Form("Min. bias"));
  
  TH1D * hCorrMinBiasEMsys = (TH1D*) fin100EM->Get(Form("%s_00_syst", histCorrName.Data()));
  if (!hCorrMinBiasEMsys) Printf("cannot find min bias spectrum with sys 0-100 EM");
  ApplyCosmetics2Spectrum(hCorrMinBiasEMsys, color[5], 0); 
  hCorrMinBiasEMsys->SetTitle(Form("Min. bias"));
  ctmp->cd();
  hCorrMinBiasEMsys->Draw("E2");
  hCorrMinBiasEM->Draw("same");

  TH1D * hCorrMinBiasEMuncorr = (TH1D*) fin100EM->Get(Form("%s_00_syst_uncorr", histCorrName.Data()));
  ApplyCosmetics2Spectrum(hCorrMinBiasEMuncorr, color[5], 0); 
  TH1D * hCorrMinBiasEMuncorrPi = (TH1D*) fin100EM->Get(Form("%s_00_syst_uncorrPi", histCorrName.Data()));
  ApplyCosmetics2Spectrum(hCorrMinBiasEMuncorrPi, color[5], 0); 
  TH1D * hCorrMinBiasEMuncorrKa = (TH1D*) fin100EM->Get(Form("%s_00_syst_uncorrKa", histCorrName.Data()));
  ApplyCosmetics2Spectrum(hCorrMinBiasEMuncorrKa, color[5], 0); 
  TH1D * hCorrMinBiasEMuncorrPro = (TH1D*) fin100EM->Get(Form("%s_00_syst_uncorrPro", histCorrName.Data()));
  ApplyCosmetics2Spectrum(hCorrMinBiasEMuncorrPro, color[5], 0); 
  

  //access files with min bias spectra LS
  TFile * fin100LS;
  if (useSmoothed) 
    fin100LS = TFile::Open("/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/0to100_binB/preliminary/finalWsyst_smooth2_LS_27apr14_0.root");
  else 
    fin100LS = TFile::Open("/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/0to100_binB/preliminary/finalWsyst_LS_27apr14_0.root");
  if (!fin100LS) Printf("cannot open file min bias LS.");
  TH1D * hCorrMinBiasLS = (TH1D*) fin100LS->Get(Form("%s_0", histCorrName.Data()));
  if (!hCorrMinBiasLS) Printf("cannot find min bias spectrum 0-100 LS");
  ApplyCosmetics2Spectrum(hCorrMinBiasLS, color[5], 20); 
  hCorrMinBiasLS->SetTitle(Form("Min. bias"));

  TH1D * hCorrMinBiasLSsys = (TH1D*) fin100LS->Get(Form("%s_00_syst", histCorrName.Data()));
  if (!hCorrMinBiasLSsys) Printf("cannot find min bias spectrum with sys 0-100 LS");
  ApplyCosmetics2Spectrum(hCorrMinBiasLSsys, color[5], 0); 
  hCorrMinBiasLSsys->SetTitle(Form("Min. bias"));

  TH1D * hCorrMinBiasLSuncorr = (TH1D*) fin100LS->Get(Form("%s_00_syst_uncorr", histCorrName.Data()));
  ApplyCosmetics2Spectrum(hCorrMinBiasLSuncorr, color[5]); 
  TH1D * hCorrMinBiasLSuncorrPi = (TH1D*) fin100LS->Get(Form("%s_00_syst_uncorrPi", histCorrName.Data()));
  ApplyCosmetics2Spectrum(hCorrMinBiasLSuncorrPi, color[5]); 
  TH1D * hCorrMinBiasLSuncorrKa = (TH1D*) fin100LS->Get(Form("%s_00_syst_uncorrKa", histCorrName.Data()));
  ApplyCosmetics2Spectrum(hCorrMinBiasLSuncorrKa, color[5]); 
  TH1D * hCorrMinBiasLSuncorrPro = (TH1D*) fin100LS->Get(Form("%s_00_syst_uncorrPro", histCorrName.Data()));
  ApplyCosmetics2Spectrum(hCorrMinBiasLSuncorrPro, color[5]); 

  ctmp->cd();
  hCorrMinBiasEMsys->Draw("E2same");
  hCorrMinBiasEM->Draw("same");
  
  //define histos
  TFile * finMultiCorrEM[5]; 
  TH1D * hCorrMultiEM[5];
  TH1D * hCorrMultiEMsys[5];
  TH1D * hCorrEMuncorr[5];
  TH1D * hCorrEMuncorrPi[5];
  TH1D * hCorrEMuncorrKa[5];
  TH1D * hCorrEMuncorrPro[5];
  
  TFile * finMultiCorrLS[5];
  TH1D * hCorrMultiLS[5];
  TH1D * hCorrMultiLSsys[5];
  TH1D * hCorrLSuncorr[5];
  TH1D * hCorrLSuncorrPi[5];
  TH1D * hCorrLSuncorrKa[5];
  TH1D * hCorrLSuncorrPro[5];

  for (Int_t ic=0;ic<5;ic++) {
    //access files with multi-dependent spectra EM
    if (useSmoothed) 
      finMultiCorrEM[ic] = TFile::Open(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/preliminary/finalWsyst_smooth2_EM_27apr14_%i.root", ic));
    else    
      finMultiCorrEM[ic] = TFile::Open(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/preliminary/finalWsyst_EM_27apr14_%i.root", ic));

    if (!finMultiCorrEM[ic]) Printf("cannot open file multi EM.");
    
    //access files with multi-dependent spectra LS
    if (useSmoothed) 
      finMultiCorrLS[ic] = TFile::Open(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/preliminary/finalWsyst_smooth2_LS_27apr14_%i.root", ic));
    else    
      finMultiCorrLS[ic] = TFile::Open(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/preliminary/finalWsyst_LS_27apr14_%i.root", ic));
    if (!finMultiCorrLS[ic]) Printf("cannot open file multi EM.");
    
    //get histo multi dependent
    //EM
    hCorrMultiEM[ic] = (TH1D*) finMultiCorrEM[ic]->Get(Form("%s_%i", histCorrName.Data(), ic));
    if (!hCorrMultiEM[ic]) Printf("cannot find EM spectrum for multi bin %i.",ic);  
    ApplyCosmetics2Spectrum(hCorrMultiEM[ic], color[ic]);
    hCorrMultiEM[ic]->SetTitle(Form("%i-%i%%", ic*20, (ic+1)*20));
    
    hCorrMultiEMsys[ic] = (TH1D*) finMultiCorrEM[ic]->Get(Form("%s_%i%i_syst", histCorrName.Data(), ic, ic));
    if (!hCorrMultiEMsys[ic]) Printf("cannot find EM with sys spectrum for multi bin %i.",ic); 
    ApplyCosmetics2Spectrum(hCorrMultiEMsys[ic], color[ic], 0); 
    hCorrMultiEMsys[ic]->SetTitle(Form("%i-%i%%", ic*20, (ic+1)*20));

    hCorrEMuncorr[ic] = (TH1D*) finMultiCorrEM[ic]->Get(Form("%s_%i%i_syst_uncorr", histCorrName.Data(), ic, ic));
    ApplyCosmetics2Spectrum(hCorrEMuncorr[ic], color[ic]); 
    hCorrEMuncorrPi[ic] = (TH1D*) finMultiCorrEM[ic]->Get(Form("%s_%i%i_syst_uncorrPi", histCorrName.Data(), ic, ic));
    ApplyCosmetics2Spectrum(hCorrEMuncorrPi[ic], color[ic]); 
    hCorrEMuncorrKa[ic] = (TH1D*) finMultiCorrEM[ic]->Get(Form("%s_%i%i_syst_uncorrKa", histCorrName.Data(), ic, ic));
    ApplyCosmetics2Spectrum(hCorrEMuncorrKa[ic], color[ic]); 
    hCorrEMuncorrPro[ic] = (TH1D*) finMultiCorrEM[ic]->Get(Form("%s_%i%i_syst_uncorrPro", histCorrName.Data(), ic, ic));
    ApplyCosmetics2Spectrum(hCorrEMuncorrPro[ic], color[ic]); 
    
    //LS
    hCorrMultiLS[ic] = (TH1D*) finMultiCorrLS[ic]->Get(Form("%s_%i", histCorrName.Data(), ic));
    if (!hCorrMultiLS[ic]) Printf("cannot find LS spectrum for multi bin %i.",ic);  
    ApplyCosmetics2Spectrum(hCorrMultiLS[ic], color[ic], 0); 
    hCorrMultiLS[ic]->SetTitle(Form("%i-%i%%", ic*20, (ic+1)*20));

    hCorrMultiLSsys[ic] = (TH1D*) finMultiCorrLS[ic]->Get(Form("%s_%i%i_syst", histCorrName.Data(), ic, ic));
    if (!hCorrMultiLSsys[ic]) Printf("cannot find LS with sys spectrum for multi bin %i.",ic); 
    ApplyCosmetics2Spectrum(hCorrMultiLSsys[ic], color[ic], 0); 
    hCorrMultiLSsys[ic]->SetTitle(Form("%i-%i%%", ic*20, (ic+1)*20));

    hCorrLSuncorr[ic] = (TH1D*) finMultiCorrLS[ic]->Get(Form("%s_%i%i_syst_uncorr", histCorrName.Data(), ic, ic));
    ApplyCosmetics2Spectrum(hCorrLSuncorr[ic], color[ic]); 
    hCorrLSuncorrPi[ic] = (TH1D*) finMultiCorrLS[ic]->Get(Form("%s_%i%i_syst_uncorrPi", histCorrName.Data(), ic, ic));
    ApplyCosmetics2Spectrum(hCorrLSuncorrPi[ic], color[ic]); 
    hCorrLSuncorrKa[ic] = (TH1D*) finMultiCorrLS[ic]->Get(Form("%s_%i%i_syst_uncorrKa", histCorrName.Data(), ic, ic));
    ApplyCosmetics2Spectrum(hCorrLSuncorrKa[ic], color[ic]); 
    hCorrLSuncorrPro[ic] = (TH1D*) finMultiCorrLS[ic]->Get(Form("%s_%i%i_syst_uncorrPro", histCorrName.Data(), ic, ic));
    ApplyCosmetics2Spectrum(hCorrLSuncorrPro[ic], color[ic]); 

    //display original plots
    ctmp->cd();
    hCorrMultiEM[ic]->Draw("same");
    hCorrMultiLS[ic]->Draw("same");
  }

  //use first two points from LS and remaining from EM
  TH1D * hFinalMinBias  = (TH1D*) hCorrMinBiasEM->Clone("hKstar_100");//new TH1D("hKstar_100","(K^{0}*+#bar{K^{0}*})/2, 0-100%", npt, pt);
  hFinalMinBias->SetBinContent(1, hCorrMinBiasLS->GetBinContent(1));
  hFinalMinBias->SetBinError(1, hCorrMinBiasLS->GetBinError(1));

  hFinalMinBias->SetBinContent(2, hCorrMinBiasLS->GetBinContent(2));
  hFinalMinBias->SetBinError(2, hCorrMinBiasLS->GetBinError(2));

  TH1D * hFinalMinBiasSys  = (TH1D*) hCorrMinBiasEMsys->Clone("hKstar_100_sys");//new TH1D("hKstar_100","(K^{0}*+#bar{K^{0}*})/2, 0-100%", npt, pt);
  hFinalMinBiasSys->SetBinContent(1, hCorrMinBiasLSsys->GetBinContent(1));
  hFinalMinBiasSys->SetBinError(1, hCorrMinBiasLSsys->GetBinError(1));
  
  hFinalMinBiasSys->SetBinContent(2, hCorrMinBiasLSsys->GetBinContent(2));
  hFinalMinBiasSys->SetBinError(2, hCorrMinBiasLSsys->GetBinError(2));
  
  TH1D * hFinalMinBiasSysUncorr  = (TH1D*) hCorrMinBiasEMuncorr->Clone("hKstar_100_sys_uncorr");//new TH1D("hKstar_100","(K^{0}*+#bar{K^{0}*})/2, 0-100%", npt, pt);
  hFinalMinBiasSysUncorr->SetBinContent(1, hCorrMinBiasLSuncorr->GetBinContent(1));
  hFinalMinBiasSysUncorr->SetBinError(1, hCorrMinBiasLSuncorr->GetBinError(1)); 
  hFinalMinBiasSysUncorr->SetBinContent(2, hCorrMinBiasLSuncorr->GetBinContent(2));
  hFinalMinBiasSysUncorr->SetBinError(2, hCorrMinBiasLSuncorr->GetBinError(2));

  TH1D * hFinalMinBiasSysUncorrPi  = (TH1D*) hCorrMinBiasEMuncorrPi->Clone("hKstar_100_sys_uncorrPi");//new TH1D("hKstar_100","(K^{0}*+#bar{K^{0}*})/2, 0-100%", npt, pt);
  hFinalMinBiasSysUncorrPi->SetBinContent(1, hCorrMinBiasLSuncorrPi->GetBinContent(1));
  hFinalMinBiasSysUncorrPi->SetBinError(1, hCorrMinBiasLSuncorrPi->GetBinError(1)); 
  hFinalMinBiasSysUncorrPi->SetBinContent(2, hCorrMinBiasLSuncorrPi->GetBinContent(2));
  hFinalMinBiasSysUncorrPi->SetBinError(2, hCorrMinBiasLSuncorrPi->GetBinError(2));

  TH1D * hFinalMinBiasSysUncorrKa  = (TH1D*) hCorrMinBiasEMuncorrKa->Clone("hKstar_100_sys_uncorrKa");//new TH1D("hKstar_100","(K^{0}*+#bar{K^{0}*})/2, 0-100%", npt, pt);
  hFinalMinBiasSysUncorrKa->SetBinContent(1, hCorrMinBiasLSuncorrKa->GetBinContent(1));
  hFinalMinBiasSysUncorrKa->SetBinError(1, hCorrMinBiasLSuncorrKa->GetBinError(1));
  hFinalMinBiasSysUncorrKa->SetBinContent(2, hCorrMinBiasLSuncorrKa->GetBinContent(2));
  hFinalMinBiasSysUncorrKa->SetBinError(2, hCorrMinBiasLSuncorrKa->GetBinError(2));

  TH1D * hFinalMinBiasSysUncorrPro  = (TH1D*) hCorrMinBiasEMuncorrPro->Clone("hKstar_100_sys_uncorrPro");//new TH1D("hKstar_100","(K^{0}*+#bar{K^{0}*})/2, 0-100%", npt, pt);
  hFinalMinBiasSysUncorrPro->SetBinContent(1, hCorrMinBiasLSuncorrPro->GetBinContent(1));
  hFinalMinBiasSysUncorrPro->SetBinError(1, hCorrMinBiasLSuncorrPro->GetBinError(1));
  hFinalMinBiasSysUncorrPro->SetBinContent(2, hCorrMinBiasLSuncorrPro->GetBinContent(2));
  hFinalMinBiasSysUncorrPro->SetBinError(2, hCorrMinBiasLSuncorrPro->GetBinError(2));
  
  TCanvas * cfinal = new TCanvas("cfinal","cfinal", 700, 900); 
  cfinal->cd();
  gPad->SetLogy();
  hFinalMinBiasSys->Draw("E2");
  hFinalMinBias->Draw("same");
  AddPaveTextErrors();
  AddPaveText_KStar_pPb("tr");

  TCanvas * cfinalMulti = new TCanvas("cfinalMulti","cfinalMulti", 700, 900); 
  cfinalMulti->cd();
  
  TH1D * hFinalMulti[5];
  TH1D * hFinalMultiSys[5];
  TH1D * hFinalMultiUncorr[5];
  TH1D * hFinalMultiUncorrPi[5];
  TH1D * hFinalMultiUncorrKa[5];
  TH1D * hFinalMultiUncorrPro[5];
 for (Int_t ic =0; ic<ncent; ic++) {
    //stat uncert only
    hFinalMulti[ic] = (TH1D*) hCorrMultiEM[ic]->Clone(Form("hKstar_%i",ic)); 
    hFinalMulti[ic]->SetBinContent(1, hCorrMultiLS[ic]->GetBinContent(1));
    hFinalMulti[ic]->SetBinError(1, hCorrMultiLS[ic]->GetBinError(1));
    hFinalMulti[ic]->SetBinContent(2, hCorrMultiLS[ic]->GetBinContent(2));
    hFinalMulti[ic]->SetBinError(2, hCorrMultiLS[ic]->GetBinError(2));
    
    //sys uncert
    hFinalMultiSys[ic] = (TH1D*) hCorrMultiEMsys[ic]->Clone(Form("hKstar_%i_sys",ic)); 
    hFinalMultiSys[ic]->SetBinContent(1, hCorrMultiLSsys[ic]->GetBinContent(1));
    hFinalMultiSys[ic]->SetBinError(1, hCorrMultiLSsys[ic]->GetBinError(1));
    hFinalMultiSys[ic]->SetBinContent(2, hCorrMultiLSsys[ic]->GetBinContent(2));
    hFinalMultiSys[ic]->SetBinError(2, hCorrMultiLSsys[ic]->GetBinError(2));
    
    //sys uncert uncorr pt
    hFinalMultiUncorr[ic] = (TH1D*) hCorrEMuncorr[ic]->Clone(Form("hKstar_%i_Uncorr",ic)); 
    hFinalMultiUncorr[ic]->SetBinContent(1, hCorrLSuncorr[ic]->GetBinContent(1));
    hFinalMultiUncorr[ic]->SetBinError(1, hCorrLSuncorr[ic]->GetBinError(1));
    hFinalMultiUncorr[ic]->SetBinContent(2, hCorrLSuncorr[ic]->GetBinContent(2));
    hFinalMultiUncorr[ic]->SetBinError(2, hCorrLSuncorr[ic]->GetBinError(2));
    //sys uncert uncorr hadrons
    hFinalMultiUncorrPi[ic] = (TH1D*) hCorrEMuncorrPi[ic]->Clone(Form("hKstar_%i_UncorrPi",ic)); 
    hFinalMultiUncorrPi[ic]->SetBinContent(1, hCorrLSuncorrPi[ic]->GetBinContent(1));
    hFinalMultiUncorrPi[ic]->SetBinError(1, hCorrLSuncorrPi[ic]->GetBinError(1));
    hFinalMultiUncorrPi[ic]->SetBinContent(2, hCorrLSuncorrPi[ic]->GetBinContent(2));
    hFinalMultiUncorrPi[ic]->SetBinError(2, hCorrLSuncorrPi[ic]->GetBinError(2));
   //sys uncert uncorr hadrons
    hFinalMultiUncorrKa[ic] = (TH1D*) hCorrEMuncorrKa[ic]->Clone(Form("hKstar_%i_UncorrKa",ic)); 
    hFinalMultiUncorrKa[ic]->SetBinContent(1, hCorrLSuncorrKa[ic]->GetBinContent(1));
    hFinalMultiUncorrKa[ic]->SetBinError(1, hCorrLSuncorrKa[ic]->GetBinError(1));
    hFinalMultiUncorrKa[ic]->SetBinContent(2, hCorrLSuncorrKa[ic]->GetBinContent(2));
    hFinalMultiUncorrKa[ic]->SetBinError(2, hCorrLSuncorrKa[ic]->GetBinError(2));
   //sys uncert uncorr hadrons
    hFinalMultiUncorrPro[ic] = (TH1D*) hCorrEMuncorrPro[ic]->Clone(Form("hKstar_%i_UncorrPro",ic)); 
    hFinalMultiUncorrPro[ic]->SetBinContent(1, hCorrLSuncorrPro[ic]->GetBinContent(1));
    hFinalMultiUncorrPro[ic]->SetBinError(1, hCorrLSuncorrPro[ic]->GetBinError(1));
    hFinalMultiUncorrPro[ic]->SetBinContent(2, hCorrLSuncorrPro[ic]->GetBinContent(2));
    hFinalMultiUncorrPro[ic]->SetBinError(2, hCorrLSuncorrPro[ic]->GetBinError(2));


    cfinalMulti->cd();
    gPad->SetLogy();
    hFinalMultiSys[ic]->Draw((ic==0?"E2":"E2same"));
    hFinalMulti[ic]->Draw("same");
  }
  cfinalMulti->cd();
  AddPaveTextErrors();
  AddPaveText_KStar_pPb("tr");
  TLegend * leg = new TLegend(0.55,0.5,0.89,0.8,"V0A Multiplicity Classes (Pb side)");
  for (Int_t ic=0;ic<5;ic++){ leg->AddEntry(hFinalMulti[ic], hFinalMulti[ic]->GetTitle(), "lp");}
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->Draw();
 
  //save into filehCorrected_44_syst_uncorrhCorrected_44_syst_uncorr
  TFile * fout = new TFile(Form("prelim_kstar_pPb%s.root", useSmoothed?"_smoothSys":""),"recreate");
  fout->cd(); 
  hFinalMinBiasSys->Write();
  hFinalMinBias->Write();
  hFinalMinBiasSysUncorr->Write();
  hFinalMinBiasSysUncorrPi->Write();
  hFinalMinBiasSysUncorrKa->Write();
  hFinalMinBiasSysUncorrPro->Write();
  for (Int_t ic=0;ic<5;ic++){
    hFinalMulti[ic]->Write();
    hFinalMultiSys[ic]->Write();
    hFinalMultiUncorr[ic]->Write();
    hFinalMultiUncorrPi[ic]->Write();
    hFinalMultiUncorrKa[ic]->Write();
    hFinalMultiUncorrPro[ic]->Write();
  }
  
  cfinal->Write();  
  cfinalMulti->Write();
  return;
}


void ApplyCosmetics2Spectrum(TH1D* h, Color_t color, Int_t marker=-1, Int_t lineWidth = 1)
{
  if (!h) return;
  h->GetXaxis()->SetRangeUser(0.0,14.9);
  h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  h->GetYaxis()->SetTitle("d^{2}#it{N}/(d#it{p}_{T}d#it{y}) (GeV/#it{c})^{-1}");
  h->GetYaxis()->SetTitleOffset(1.7);
  h->SetLineColor(color);
  h->SetLineWidth(1);
  h->SetFillColor(kGray);
  h->SetFillStyle(0);
  h->SetMarkerColor(color);
  if (marker>=0) h->SetMarkerStyle(marker);
  TString name = (h->GetName());
  if (name.Contains("sys")) h->SetDrawOption("E2");
  return;
}

TPaveText * AddPaveText_KStar_pPb(TString position = "bl"){
  TString text="ALICE, p-Pb #sqrt{s_{NN}} = 5.02 TeV";
  TString textks = "K*^{0}#rightarrow K#pi, 0 < y < 0.5";
  TString textcent="Centrality 0-100%";
  TPaveText * pave;
  if (position.Contains("tl")) pave = new TPaveText(0.12,0.75,0.45,0.89,"NDC");
  else
    if (position.Contains("bl")) pave = new TPaveText(0.22,0.22,0.55,0.36,"NDC");
    else 
      if (position.Contains("tr")) pave = new TPaveText(0.50,0.80,0.89,0.89,"NDC");
      else
	if (position.Contains("br")) pave = new TPaveText(0.55,0.12,0.89,0.26,"NDC");
  pave->SetTextColor(kBlack);
  pave->SetTextAlign(12);
  pave->SetTextFont(42);
  //pave->SetTextSize(18);
  pave->SetBorderSize(0);
  pave->SetFillColor(kWhite);
  pave->SetFillStyle(0);
  pave->InsertText(text.Data());
  pave->InsertText(textks.Data());
  //pave->InsertText(textcent.Data());
  pave->Draw();
  return pave;
}


TPaveText * AddPaveTextErrors(TString text="Uncertainties: stat.(bars), syst.(boxes)"){
  TPaveText * pave = new TPaveText(0.12,0.13,0.63,0.18,"NDC");
  pave->SetBorderSize(0);
  pave->SetFillColor(kWhite);
  pave->SetFillStyle(0);
  pave->SetTextColor(kBlack);
  pave->SetTextFont(42);
  pave->InsertText(text.Data());
  pave->Draw();
  return pave;
}
