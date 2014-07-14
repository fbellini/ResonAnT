void compareCent(){

  TString massName = "hMassVsPt_";
  TString widthName = "hWidthVsPt_";
  TString dataDir = "/Users/bellini/alice/resonances/myKstar/pwglf_train_out/data/aod049";
  TFile *file0 = TFile::Open(Form("%s/roofit/fitEM_BW_EXP_sub_proj_LHC10h_AOD049_merged_AnalysisResults.root", dataDir.Data()));
  TFile *file1 = TFile::Open(Form("%s/roofit/fitEM_BW_P0LY1_sub_proj_LHC10h_AOD049_merged_AnalysisResults.root", dataDir.Data()));
  TFile *file2 = TFile::Open(Form("%s/roofit/fitEM_BW_P0LY1_sub_proj_LHC10h_AOD049_merged_AnalysisResults.root", dataDir.Data()));
  TFile *file3 = TFile::Open(Form("%s/roofit/fitEM_BW_P0LY1_sub_proj_LHC10h_AOD049_merged_AnalysisResults.root", dataDir.Data()));
  
  for (Int_t icent = 0; icent<4;icent++){
    massName.Append(Form("%i",icent));
    
  }

  }

