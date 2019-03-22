void normalizeByEvtN(TString spectraFileName, Int_t aodN = 49, Bool_t isTPC=0)
{
  Float_t nEvents49[4] = {3.069377e+06, 3.055345e+06, 3.059987e+06, 3.073845e+06}; //aod049
  Float_t nEvents49tpc[4] = { 2.05506e+06, 2.05308e+06, 2.05239e+06, 2.06007e+06}; //aod49 TPC
  Float_t nEvents86[4] = { 3.248661e+06, 3.258249e+06, 3.255904e+06, 3.258934e+06}; //aod086
  Float_t nEvents95[4] = { 2.274366e+07, 9.840848e+06, 5.453893e+06, 1.098727e+06}; //aod095
  Float_t nEvents115[4] = { 2.240713e+07, 9.710764e+06, 5.388053e+06, 1.081009e+06}; //aod115
  Float_t nEvents49t47[4] = {2.992844e+06, 2.982032e+06, 2.985484e+06, 2.999089e+06}; //aod049 train47
  Float_t nEvents49tpct36[4] = { 2.707814e+06, 2.696713e+06, 2.700422e+06, 2.713401e+06}; //aod49 train 36 TPC
  Float_t nEvents4990[4] = { 2.732529e+06, 2.720139e+06, 2.723456e+06, 2.735298e+06}; //aod49 train 90
  
  Color_t color[]={kRed, kOrange, kGreen+2, kBlue, kBlack};
  Int_t marker[]={20, 21, 28, 22, 23};
  
  TCanvas * c1 = new TCanvas(Form("norm_%i",0),"normalized spectra",600,500);
  TFile *fraw = TFile::Open(spectraFileName.Data());
  if (!fraw){
    Printf("No input file");
    return;
  }
  TString newName;
  if (spectraFileName.Contains("rawYields")) newName = spectraFileName.ReplaceAll("rawYields","normYields");
  else newName = spectraFileName.ReplaceAll("raw","norm");
  TFile * fout = new TFile(newName.Data(),"recreate");

  for (Int_t ic = 0; ic<4;ic++){
    TString spectraHistName = Form("hRawYieldVsPt_%i",ic);
    if (isTPC) spectraHistName.ReplaceAll("hRawYieldVsPt_","h");
    TH1F * hraw = (TH1F*) fraw->Get(spectraHistName.Data());
    if (!hraw) continue;
    else Printf("--> file %s: normalising %s",spectraFileName.Data(), spectraHistName.Data() );
    
    if (aodN==49){
      if (isTPC) hraw->Scale(1./nEvents49tpc[ic]);
      else hraw->Scale(1./nEvents49[ic]);
    }

    if (aodN==4936) hraw->Scale(1./nEvents49tpct36[ic]);
    if (aodN==4947) hraw->Scale(1./nEvents49t47[ic]);
    if (aodN==86) hraw->Scale(1./nEvents86[ic]);
    if (aodN==95) hraw->Scale(1./nEvents95[ic]);
    if (aodN==115) hraw->Scale(1./nEvents115[ic]);
    if (aodN==4990) hraw->Scale(1./nEvents4990[ic]);

    
    hraw->SetTitle("Normalized raw yields (EM bg)");
    hraw->GetYaxis()->SetTitle("1/N_{evt}* dN/dydp_{t} (|y|<0.5)");
    hraw->GetYaxis()->SetRangeUser(5.e-7,1.2e1);
    hraw->GetXaxis()->SetRangeUser(0.0,10.0);
    hraw->SetLineColor(color[ic]);
    hraw->SetMarkerColor(color[ic]);
    hraw->SetMarkerStyle(marker[ic]);
    if (isTPC) hraw->SetMarkerStyle(marker[ic]+4);
    hraw->SetLineWidth(2);

    c1->cd();
    gPad->SetLogy();
    // if (ic==0)  hraw->Draw();
    // else     hraw->Draw("same"); 
    hraw->Draw();
    hraw->Write();
    TString pngname(newName.Data());
    pngname.ReplaceAll("normYields/","normYields/img/");
    pngname.ReplaceAll(".root",".png");
    c1->SaveAs(pngname.Data());
  }
  
  Printf("DONE");
  return;
}


