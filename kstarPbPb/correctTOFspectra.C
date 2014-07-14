void correctTOFspectra(TString spectraFileName="/Users/bellini/alice/resonances/myKstar/pwglf_train_out/compare_20mar13/normYields/normYields_tree_best_poly2_aod049_01mar13.root"){
 
  Color_t color[]={kRed, kOrange, kGreen+2, kBlue, kBlack};
  Int_t marker[]={20, 21, 28, 22, 23};
  
  
  const TString effHistName = "hEffVsPt";
  Float_t branchingRatio = 0.666;
  
  //TPC_TOF matching eff
  TCanvas * c1 = new TCanvas(Form("corrected_%i",0),"corrected spectra",600,500);
  TCanvas * ceff = new TCanvas("eff","efficiency correction",600,500);
  
  TFile *fraw = TFile::Open(spectraFileName.Data());
  if (!fraw) return;
  
  TString newName;
  if (spectraFileName.Contains("normYields")) newName = spectraFileName.ReplaceAll("normYields","corrYields");
  else newName = spectraFileName.ReplaceAll("norm","corr");
  
  newName.ReplaceAll(".root","_v05lug13.root");
  TFile * fout = new TFile(newName.Data(),"recreate");
 
  //K*+antK* matching eff
  for (Int_t ic = 0; ic<4;ic++){
    TString effFileName(Form("eff/efficiency_RsnOut_Tof20sigma_centBin0%i.root",ic));
    //TString effFileName(Form("eff/efficiency_RsnOut_Tof2s_KTrd_centBin0%i.root",ic));
    TString spectraHistName(Form("hRawYieldVsPt_%i",ic));
    
    TFile *feff = TFile::Open(effFileName.Data());
    if (!feff) return;
    TH1F * hraw = (TH1F*) fraw->Get(spectraHistName.Data());
    if (!hraw) return;
    TH1F * heff = (TH1F*) feff->Get(effHistName.Data());
    if (!heff) return;
    
    hraw->Scale(1./0.666); //BR
    hraw->Scale(1./2.);//to account for particle and antiparticle
    
    //heff->Scale(0.912);//correction for PID tail
    heff->Scale(0.95*0.95);//correction for DATA vs MC matching eff
    // Printf("raw has %i bins",hraw->GetXaxis()->GetNbins());
    // Printf("eff has %i bins",heff->GetXaxis()->GetNbins());
 

    //correct spectrum
    hraw->SetBinContent(1,0.0);
    //hraw->SetBinContent(2,0.0);
    hraw->Divide(heff);
    
    hraw->SetTitle(Form("corrected, cent bin %i",ic));
    hraw->GetYaxis()->SetTitle("1/N_{evt}* d^{2}N/dydp_{t} * 1/#epsilon * 1/B.R. (|y|<0.5)");
    hraw->GetYaxis()->SetRangeUser(1.e-4,1.5e1);
    hraw->GetXaxis()->SetRangeUser(0.0,10.);
    hraw->SetLineColor(color[ic]);
    hraw->SetMarkerColor(color[ic]);
    hraw->SetMarkerStyle(marker[ic]);
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
    
    // if (isRoofit) {
    //   hraw->SetLineColor(kGreen+2);
    //   hraw->SetMarkerColor(kGreen+2);
    // }
    fout->cd();
    heff->Write();
    hraw->Write();
   }
  fout->cd();
  c1->Write();
  ceff->Write();
  fout->Close();
  return;
}
