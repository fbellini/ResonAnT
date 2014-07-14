void compareTOFTPCcorrectedSpectra(TString fileTOF, TString fileTPC, Int_t centBin=0){

  TFile* ftpc = TFile::Open(fileTPC.Data(),"read");
  TFile* ftof = TFile::Open(fileTOF.Data(),"read");

  Color_t color[]={kRed, kOrange, kGreen, kBlue};
  Int_t marker[]={20, 21, 28, 22};

  
  const Double_t tpcBins[]={0.3, 0.8, 1.2, 1.6, 2., 2.5, 3., 4., 5.};
  TString tofhname = Form("sum_corr_%i", centBin); //Form("hRawYieldVsPt_%i",centBin);  
  TString tpchname = Form("htpc_%i",centBin);
  
  TH1F* htof = (TH1F*)ftof->Get(tofhname.Data());
  TH1F* htpc =(TH1F*)ftpc->Get(tpchname.Data());
 
  if (!htof){
    Printf("Tof data unaccessible");
    return;
  }
  if (!htpc){
    Printf("Tpc data unaccessible");
    return;
  }
  TH1F* hratio = new TH1F("rawTOFoverTPC","rawTOFoverTPC", 8, tpcBins);//htof->Clone("TOFoverTPC");
  hratio->SetBinContent(1, 0.0);
  hratio->SetBinContent(2, htof->GetBinContent(2));
  hratio->SetBinContent(3, htof->GetBinContent(3));
  hratio->SetBinContent(4, htof->GetBinContent(4));
  hratio->SetBinContent(5, htof->GetBinContent(5));
  hratio->SetBinContent(6, htof->GetBinContent(6));
  hratio->SetBinContent(7, (htof->GetBinContent(7)+htof->GetBinContent(8))/2.);
  hratio->SetBinContent(8, htof->GetBinContent(9));
  // hratio->SetBinContent(9, htof->GetBinContent(10));
  // hratio->SetBinContent(10, htof->GetBinContent(11));

  hratio->SetBinError(1, 0.0);
  hratio->SetBinError(2, htof->GetBinError(2));
  hratio->SetBinError(3, htof->GetBinError(3));
  hratio->SetBinError(4, htof->GetBinError(4));
  hratio->SetBinError(5, htof->GetBinError(5));
  hratio->SetBinError(6, htof->GetBinError(6));
  hratio->SetBinError(7, (htof->GetBinError(7)+htof->GetBinError(8))/2);
  hratio->SetBinError(8, htof->GetBinError(9));
  // hratio->SetBinError(9, htof->GetBinError(10));
  // hratio->SetBinError(10, htof->GetBinError(11));
  
  /*
   TH1F* trebtpc = new TH1F("rebinnedTPC","rebinned", 10, tpcBins);//htpc->Clone("TOFoverTPC");
  trebtpc->SetBinContent(1, 0.0);
  trebtpc->SetBinContent(2, htpc->GetBinContent(1));
  trebtpc->SetBinContent(3, htpc->GetBinContent(2));
  trebtpc->SetBinContent(4, htpc->GetBinContent(3));
  trebtpc->SetBinContent(5, htpc->GetBinContent(4));
  trebtpc->SetBinContent(6, htpc->GetBinContent(5));
  trebtpc->SetBinContent(7, htpc->GetBinContent(6));
  trebtpc->SetBinContent(8, htpc->GetBinContent(7));
  trebtpc->SetBinContent(9, 0.0);
  trebtpc->SetBinContent(10, 0.0);

  trebtpc->SetBinError(1, 0.0);
  trebtpc->SetBinError(2, htpc->GetBinError(1));
  trebtpc->SetBinError(3, htpc->GetBinError(2));
  trebtpc->SetBinError(4, htpc->GetBinError(3));
  trebtpc->SetBinError(5, htpc->GetBinError(4));
  trebtpc->SetBinError(6, htpc->GetBinError(5));
  trebtpc->SetBinError(7, htpc->GetBinError(6));
  trebtpc->SetBinError(8, htpc->GetBinError(7));
  trebtpc->SetBinError(9, 0.0);
  trebtpc->SetBinError(10, 0.0);
  */
  TH1F* hcorrtpc =(TH1F*)htpc->Clone("corrTPC");
  TH1F* hcorrtof = (TH1F*)hratio->Clone("corrTOF");
  hcorrtpc->SetLineColor(color[centBin]);
  hcorrtpc->SetMarkerColor(color[centBin]);
  hcorrtpc->SetMarkerStyle(marker[centBin]+4);
  
  hcorrtof->SetLineColor(color[centBin]);
  hcorrtof->SetMarkerColor(color[centBin]);
  hcorrtof->SetMarkerStyle(marker[centBin]);
  

  gROOT->LoadMacro("${ASD}/GetPlotRatio.C");
  GetPlotRatio(hcorrtof,hcorrtpc, kTRUE,"TOF","TPC","corrected spectra");
  // if (!hratio) return;
  // TH1F* hratioraw = new TH1F("hratioraw","hratioraw",7,tpcBins);
  // for (Int_t j=0;j<7;j++){
  //   hratioraw->SetBinContent(j, hratio->GetBinContent(j)/htpc->GetBinContent(j));
  // }
  
  // TCanvas *c1 = new TCanvas("c1","c1",600,600);
  // hratioraw->GetYaxis()->SetRangeUser(1.e-2,1.);
  // hratioraw->Draw();
  
  // TCanvas *c2 = new TCanvas("c2","c2",600,600);
  // hcorrtof->SetLineColor(kRed);
  // hcorrtof->SetMarkerColor(kRed);
  // hcorrtpc->SetLineColor(kBlack);
  // hcorrtpc->SetMarkerColor(kBlack);
  // hcorrtof->GetYaxis()->SetRangeUser(0.001,10.);
  // hcorrtof->Draw();
  // hcorrtpc->Draw("same");

  // TH1F* hratiocorr = new TH1F("hratiocorr","hratiocorr",7,tpcBins);
  // for (Int_t j=0;j<7;j++){
  //   hratiocorr->SetBinContent(j, hcorrtof->GetBinContent(j)/hcorrtpc->GetBinContent(j));
  // }
  // TCanvas *c3 = new TCanvas("c3","c3",600,600);
  // hratiocorr->Draw();
  return;
}

