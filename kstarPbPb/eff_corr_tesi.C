void eff_corr_tesi(TString spectraFileName,Bool_t applyEff=1,Bool_t isRoofit=1, Bool_t isLS = 0){


  const TString effHistName = "hEffVsPt";
  Float_t branchingRatio = 0.666;
  Float_t pidEff = (0.954*0.954);
  Float_t pidEffTPC = (0.9999*0.9999);
  Float_t nEvents[4] = {3.069377e+06, 3.055345e+06, 3.059987e+06, 3.073845e+06}; //aod049
  Color_t color[]={kRed+1, kOrange+1, kGreen+2, kBlue+1, kBlack};
  Int_t marker[]={20, 21, 28, 22, 23};
  
  //TPC_TOF matching eff
  TCanvas * c1 = new TCanvas(Form("corrected_%i",0),"corrected spectra",600,500);
  TCanvas * ceff = new TCanvas("eff","efficiency correction",600,500);

  //K*+antK* matching eff
  for (Int_t ic = 0; ic<4;ic++){

    const TString effFileName(Form("./efficiency_RsnOut_TofMatch_centBin0%i.root",ic));
    const TString spectraHistName(Form("hRawYieldVsPt_%i",ic));
    
    TFile *fraw = TFile::Open(spectraFileName.Data());
    if (!fraw) return;
    TFile *feff = TFile::Open(effFileName.Data());
    if (!feff) return;
    TH1F * hraw = (TH1F*) fraw->Get(spectraHistName.Data());
    if (!hraw) return;
    TH1F * heff = (TH1F*) feff->Get(effHistName.Data());
    if (!heff) return;
      
    hraw->SetBinContent(1,0.0);
    //hraw->SetBinContent(2,0.0);
    Printf("raw has %i bins",hraw->GetXaxis()->GetNbins());
    Printf("eff has %i bins",heff->GetXaxis()->GetNbins());
    hraw->Scale(0.5);//TPC analysis is (K*+antiK*)/2
    hraw->Scale(1./nEvents[ic]);
    hraw->SetTitle("Raw yields normalized by the number of accepted events (EM bg)");
    hraw->GetYaxis()->SetTitle("1/N_{evt}* dN/dp_{t} (|y|<0.5)");
    heff->Scale(pidEffTPC);//scale for pi pid eff for 2 independent daughters   
    heff->Scale(pidEff);//scale for pi pid eff for 2 independent daughters   
    hraw->Divide(heff);
    hraw->Scale(1./(branchingRatio));

    hraw->SetTitle("Corrected p_{t} spectrum (EM bg)");
    hraw->GetYaxis()->SetTitle("1/N_{evt}* dN/dydp_{t} * 1/#epsilon * 1/B.R. (|y|<0.5)");
    hraw->GetYaxis()->SetRangeUser(5.e-4,1.2e1);
    hraw->GetXaxis()->SetRangeUser(0.0,6.);
    hraw->SetLineColor(color[ic]);
    hraw->SetMarkerColor(color[ic]);
    hraw->SetMarkerStyle(marker[ic]);
    hraw->SetMarkerSize(0.7);
    hraw->SetLineWidth(2);
    heff->GetXaxis()->SetRangeUser(0.0,6.);
    heff->GetYaxis()->SetRangeUser(0.0,1.);
    heff->SetLineColor(color[ic]);
    heff->SetMarkerColor(color[ic]);
    heff->SetMarkerStyle(marker[ic]);
    heff->SetLineWidth(2);
    
    c1->cd();
    gPad->SetLogy();
    if (ic==0)  hraw->Draw();
      else hraw->Draw("same");
    ceff->cd();
    if (ic==0)  heff->Draw();
    else heff->Draw("same");
    TFile * fout = new TFile(Form("correctedYields_%s%s_cent%i.root", isLS? "LS":"EM",isRoofit? "_roofit" : "", ic),"recreate");
    // if (isRoofit) {
    //   hraw->SetLineColor(kGreen+2);
    //   hraw->SetMarkerColor(kGreen+2);
    // }
    heff->Write();
    hraw->Write();
    fout->Close();
  }
  Printf("DONE");
  return;
}
