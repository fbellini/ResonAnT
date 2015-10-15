void steerFit(TString func, Bool_t useofficial = 0, Bool_t refit = 1)
{
  
  gROOT->LoadMacro("$ASD/spectraTools/mySpectraUtils.C"); // pPb502.fullpT.SPECTRA.20141121.root only average p+ap

  //pion spectrum in pPb 
  TH1D * pi_020_stat = (TH1D *) getIdSpectrum20("stat", "pion", 0, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-PiKaPr.root"); Printf("pi");
  TH1D * pi_020_uncor = (TH1D *) getIdSpectrum20("uncor", "pion", 0, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-PiKaPr.root");
  TH1D * pi_020_sys = (TH1D *) getIdSpectrum20("sys", "pion", 0, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-PiKaPr.root");
  
  TH1D * ka_020_stat = (TH1D *) getIdSpectrum20("stat", "kaon", 0, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-PiKaPr.root");Printf("K");
  TH1D * ka_020_uncor = (TH1D *) getIdSpectrum20("uncor", "kaon", 0, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-PiKaPr.root");
  TH1D * ka_020_sys = (TH1D *) getIdSpectrum20("sys", "kaon", 0, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-PiKaPr.root");

  TH1D * prot_020_stat = (TH1D *) getIdSpectrum20("stat", "proton", 0, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-PiKaPr.root");Printf("proton");
  TH1D * prot_020_uncor = (TH1D *) getIdSpectrum20("uncor", "proton", 0, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-PiKaPr.root");
  TH1D * prot_020_sys = (TH1D *) getIdSpectrum20("sys", "proton", 0, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-PiKaPr.root");

  TH1D * lambda_020_stat = (TH1D *) getIdSpectrum20("Stat", "Lambda", 1, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-Lambda.root");Printf("Lambda");
  //TH1D * lambda_020_uncor = (TH1D *) getIdSpectrum20("Uncor", "Lambda", 1, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-Lambda.root");
  TH1D * lambda_020_sys = (TH1D *) getIdSpectrum20("Syst", "Lambda", 1, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-Lambda.root");

  TH1D * k0s_020_stat = (TH1D *) getIdSpectrum20("Stat", "K0Short", 1, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-K0s.root");Printf("K0s");
  //TH1D * k0s_020_uncor = (TH1D *) getIdSpectrum20("Uncor", "K0s", 1, "pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-K0s.root");
  TH1D * k0s_020_sys = (TH1D *) getIdSpectrum20("Syst", "K0Short", 1,"pPb", "$HOME/alice/physics/pPb5.02TeV/pPb-K0s.root");

  //Lambda in PbPb
  TH1D * lambda_PbPb_020_stat = (TH1D *) getIdSpectrum20("stat", "Lambda", 1, "PbPb", "$HOME/alice/physics/PbPb2.76TeV/k0s_lambda_final_spectra.root");
  TH1D * lambda_PbPb_020_sys = (TH1D *) getIdSpectrum20("syst", "Lambda", 1, "PbPb", "$HOME/alice/physics/PbPb2.76TeV/k0s_lambda_final_spectra.root");
 
 //Lambda in PbPb
  TH1D * lambda_PbPb_005_stat = (TH1D *) getIdSpectrum(0, "stat", "Lambda", 1, "PbPb", "$HOME/alice/physics/PbPb2.76TeV/k0s_lambda_final_spectra.root");
  TH1D * lambda_PbPb_005_sys = (TH1D *) getIdSpectrum(0, "syst", "Lambda", 1, "PbPb", "$HOME/alice/physics/PbPb2.76TeV/k0s_lambda_final_spectra.root");
 
 //K0s in PbPb
  TH1D * k0s_PbPb_020_stat = (TH1D *) getIdSpectrum20("stat", "K0s", 1, "PbPb","$HOME/alice/physics/PbPb2.76TeV/k0s_lambda_final_spectra.root");
  TH1D * k0s_PbPb_020_sys = (TH1D *) getIdSpectrum20("syst", "K0s", 1, "PbPb", "$HOME/alice/physics/PbPb2.76TeV/k0s_lambda_final_spectra.root");
    
  //Fit ranges
  Double_t pi_fitRange[2] = {0.1, 3.0};
  Double_t ka_fitRange[2] = {0.2, 2.5};
  Double_t prot_fitRange[2] = {0.3, 4.0};
  Double_t lambda_fitRange[2] = {0.4, 8.0};
  Double_t k0s_fitRange[2] = {0.3, 8.0};


  gROOT->LoadMacro("$HOME/alice/physics/meanpt/FitSpectrum.C");
  // for (Int_t j = 0; j<10; j++){
  //   for (Int_t k = 0; k<4; k++){
  //     FitSpectrum(pi_020_stat, pi_020_sys, "pion", "pPb", func.Data(), pi_fitRange[0]+0.1*k, pi_fitRange[1]-j*0.2, refit, useofficial);
  //     FitSpectrum(ka_020_stat, ka_020_sys, "kaon", "pPb",func.Data(), ka_fitRange[0]+0.1*k, ka_fitRange[1]-j*0.1, refit, useofficial);
  //     FitSpectrum(prot_020_stat, prot_020_sys, "proton", "pPb",func.Data(), prot_fitRange[0]+0.1*k, prot_fitRange[1]-j*0.2, refit, useofficial);
  //     FitSpectrum(k0s_020_stat, k0s_020_sys, "K0s", "pPb",func.Data(), k0s_fitRange[0]+0.1*k, k0s_fitRange[1]-j*0.5, refit, useofficial);
  //     FitSpectrum(lambda_020_stat, lambda_020_sys, "Lambda","pPb", func.Data(), lambda_fitRange[0]+0.1*k, lambda_fitRange[1]-j*0.5, refit, useofficial);
  //   }
  // }

  //PbPb
  // https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/belikov/2013-Jul-02-analysis_note-analysis_note_lk_2013.pdf (page 12)
  // It seems that the parameters (not yield obviously!) for all three bins are rather similar so I would hope that the fit could converge by using these to start from: 
  // Lambda: T=0.071, beta_s=0.583*1.5, n=0.68
  // K0: T=0.075, beta_s=0.61*1.5, n=0.88
  Double_t lambda_fitRangePbPb[2] = {0.6, 2.5};
  Double_t k0s_fitRangePbPb[2] = {0.4, 1.6};

  Float_t steprange = 0.2;
  for (Int_t j = 0; j<5; j++){
    for (Int_t k = 0; k<2; k++){
      // FitSpectrum(k0s_PbPb_020_stat, k0s_PbPb_020_sys, "K0s", "PbPb_020",func.Data(), k0s_fitRangePbPb[0]+k*steprange, k0s_fitRangePbPb[1]-j*steprange, refit, useofficial);
      // FitSpectrum(lambda_PbPb_020_stat, lambda_PbPb_020_sys, "Lambda", "PbPb_020", func.Data(), lambda_fitRangePbPb[0]+k*steprange, lambda_fitRangePbPb[1]-j*steprange, refit, useofficial);
      FitSpectrum(lambda_PbPb_005_stat, lambda_PbPb_005_sys, "Lambda", "PbPb_005", func.Data(), lambda_fitRangePbPb[0]+k*steprange, lambda_fitRangePbPb[1]-j*steprange, refit, useofficial);
    }
  }
  

  return;

}
