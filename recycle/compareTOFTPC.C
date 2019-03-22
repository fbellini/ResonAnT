void compareTOFTPC(){

  TFile* ftpc = TFile::Open("norm_TPCraw.root","read");
  TFile* ftof = TFile::Open("norm_TOFraw.root","read");
  
  const Double_t tpcBins[]={0., 1., 1.5, 2., 2.5, 3., 4., 5.};
  TString tofhname = "hRawYieldVsPt_0"; 
  TString tpchname = "hraw0"; 
  
  TH1F* htpc =(TH1F*)ftpc->Get(tpchname.Data());
  TH1F* htof = (TH1F*)ftof->Get(tofhname.Data());
 
  TH1F* hratio = new TH1F("rawTOFoverTPC","rawTOFoverTPC", 7, tpcBins);//htof->Clone("TOFoverTPC");
  hratio->SetBinContent(1, 0.0);
  hratio->SetBinContent(2, htof->GetBinContent(3));
  hratio->SetBinContent(3, htof->GetBinContent(4));
  hratio->SetBinContent(4, htof->GetBinContent(5));
  hratio->SetBinContent(5, htof->GetBinContent(6));
  hratio->SetBinContent(6, htof->GetBinContent(7)+htof->GetBinContent(8));
  hratio->SetBinContent(7, htof->GetBinContent(9)+htof->GetBinContent(10));

  hratio->SetBinError(1, 0.0);
  hratio->SetBinError(2, htof->GetBinError(3));
  hratio->SetBinError(3, htof->GetBinError(4));
  hratio->SetBinError(4, htof->GetBinError(5));
  hratio->SetBinError(5, htof->GetBinError(6));
  hratio->SetBinError(6, htof->GetBinError(7)+htof->GetBinError(8));
  hratio->SetBinError(7, htof->GetBinError(9)+htof->GetBinError(10));
  
  TH1F* hcorrtpc =(TH1F*)htpc->Clone("corrTPC");
  TH1F* hcorrtof = (TH1F*)hratio->Clone("corrTOF");
  
  if (!hratio) return;
  TH1F* hratioraw = new TH1F("hratioraw","hratioraw",7,tpcBins);
  for (Int_t j=0;j<7;j++){
    hratioraw->SetBinContent(j, hratio->GetBinContent(j)/htpc->GetBinContent(j));
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  hratioraw->GetYaxis()->SetRangeUser(1.e-2,1.);
  hratioraw->Draw();
  
  TFile* fefftpc = TFile::Open("eff/efficiency_RsnOut_quality_centBin00.root","read");
  TFile* fefftof = TFile::Open("eff/efficiency_RsnOut_TofMatch_centBin00.root","read");
  TH1F* hefftpc =(TH1F*)fefftpc->Get("hEffVsPt");
  TH1F* hefftof = (TH1F*)fefftof->Get("hEffVsPt");

  hcorrtpc->Divide(hefftpc);
  hcorrtof->Divide(hefftof);
  hcorrtpc->SetTitle("1/N_{evt}^{TPC} * dN/dydp_{t} * 1/#epsilon TPC");
  hcorrtof->SetTitle("1/N_{evt}^{TOF} * dN/dydp_{t} * 1/#epsilon TOF");
  Float_t totScaleFactor = 1./(0.666*0.911); //BR * PID 2sigma factor
  hcorrtpc->Scale(totScaleFactor);
  hcorrtof->Scale(totScaleFactor);
  hcorrtof->Scale(1./(0.9999*0.9999));

  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  hcorrtof->SetLineColor(kRed);
  hcorrtof->SetMarkerColor(kRed);
  hcorrtpc->SetLineColor(kBlack);
  hcorrtpc->SetMarkerColor(kBlack);
  hcorrtof->GetYaxis()->SetRangeUser(0.001,10.);
  hcorrtof->Draw();
  hcorrtpc->Draw("same");

  TH1F* hratiocorr = new TH1F("hratiocorr","hratiocorr",7,tpcBins);
  for (Int_t j=0;j<7;j++){
    hratiocorr->SetBinContent(j, hcorrtof->GetBinContent(j)/hcorrtpc->GetBinContent(j));
  }
  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  hratiocorr->Draw();
  return;
}

