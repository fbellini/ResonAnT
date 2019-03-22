TH1D * getStatUncertPercentile(TString filename = "best_fit_poly2.root", TString histname="hCorrected_0", Color_t mycolor = kBlack, Int_t mymarker = 20, TString opt = "h")
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  TCanvas * c1 = new TCanvas("c1","c1", 800, 600);
  TFile * fin = TFile::Open(filename.Data());
  
  TH1F * hraw = (TH1F*) fin->Get(histname.Data());
  if (!hraw) {Printf("Error: check histo"); return 0;}

  Double_t pt[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00, 12.0, 15.0};
  Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;   

  TH1D * hrawUnc = new TH1D(Form("hRelUncert_%s",histname.Data()),"Relative stat. uncertainty", npt, pt);
  for (Int_t i=0;i<npt;i++){
    Double_t value = hraw->GetBinContent(i+1);
    Double_t stat = hraw->GetBinError(i+1);
    if (value>0) {
      hrawUnc->SetBinContent(i+1, stat*100.0/value);
      Printf("stat uncert %%: %4.2f",stat*100.0/value);
    }
  }
  hrawUnc->SetLineColor(mycolor);
  hrawUnc->SetMarkerColor(mycolor);
  hrawUnc->SetMarkerStyle(mymarker);
  hrawUnc->GetYaxis()->SetTitle("relative statistical uncertainty (%)");
  hrawUnc->GetYaxis()->SetRangeUser(0, 20.);
  
  c1->cd();
  hrawUnc->Draw(opt.Data());
  return hrawUnc;
}
