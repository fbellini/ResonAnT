void applyEff_pp276(Bool_t isRoofit=1, Bool_t isLS = 1){
  
  TString spectraFileName;
  if (isRoofit) {
    (isLS)? spectraFileName = "/Users/bellini/alice/resonances/myKstar/ESD_pp2011/data/ptBins/rawYields_fitLS_sub_proj_AnalysisResults.root" :  spectraFileName ="/Users/bellini/alice/resonances/myKstar/ESD_pp2011/data/ptBins/rawYields_fitEM_sub_proj_AnalysisResults.root"; 
  }
  else {
    (isLS)? spectraFileName = "/Users/bellini/alice/resonances/myKstar/ESD_pp2011/data/ptBins/uncorrectedYields_EM.root":spectraFileName = "/Users/bellini/alice/resonances/myKstar/ESD_pp2011/data/ptBins/uncorrectedYields_LS.root";
  }
  const TString effFileName = "/Users/bellini/alice/resonances/myKstar/ESD_pp2011/mc/efficiency_centBin00.root";
  const TString spectraHistName = "hRawYieldVsPt_0";
  const TString effHistName = "hEffVsPt";
  Float_t nEvents = 1./19101796.; 
  Float_t branchingRatio = 1./1.;
  
  TFile *fraw = TFile::Open(spectraFileName.Data());
  if (!fraw) return;
  TFile *feff = TFile::Open(effFileName.Data());
  if (!feff) return;
  TH1F * hraw = (TH1F*) fraw->Get(spectraHistName.Data());
  if (!hraw) return;
  TH1F * heff = (TH1F*) feff->Get(effHistName.Data());
  if (!heff) return;

  hraw->Divide(heff);
  hraw->Scale(nEvents);
  hraw->Scale(branchingRatio);
  hraw->SetTitle("Corrected p_{t} spectrum (LS bg)");
  hraw->GetYaxis()->SetTitle("1/N_{evt}* dN/dp_{t} * 1/#epsilon (|y|<0.5)");
  hraw->GetYaxis()->SetRangeUser(2.e-5,1.e-1);
  hraw->GetXaxis()->SetRangeUser(0.0,7.);
  hraw->SetLineColor(kBlue+1);
  hraw->SetMarkerColor(kBlue+1);

  TCanvas * c1 = new TCanvas("corrected","corrected spectra",600,500);
  c1->cd();
  hraw->Draw();
  TFile * fout = new TFile(Form("correctedYields_%s%s.root", isLS? "LS":"EM",isRoofit? "_roofit" : ""),"recreate");
  hraw->SetLineWidth(2);
  if (isRoofit) {
    hraw->SetLineColor(kGreen+2);
    hraw->SetMarkerColor(kGreen+2);
    hraw->SetMarkerStyle(24);
  }
  hraw->Write();
  fout->Close();
  
  return;
}
