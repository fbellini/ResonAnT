void Get80to100fromSub(TString infile = "prelim_kstar_pPb_smoothSys.root", Bool_t scale2plot = 0)
{
  TCanvas * ccorr = new TCanvas("corr","corrected spectra",600,700);
  ccorr->cd();
  gPad->SetLogy();
  
  gROOT->LoadMacro("$ASD/AddPaveText.C");
  TLegend * corrleg = new TLegend(0.4,0.55,0.99,0.8,"V0A Multiplicity Classes (Pb side)");
  corrleg->SetBorderSize(0);
  corrleg->SetFillColor(kWhite);
  corrleg->SetFillStyle(0);
  
  TFile * fin = TFile::Open(infile.Data());
  TH1D * h100_stat = (TH1D*) fin->Get("hKstar_100");
  TH1D * h100_sys = (TH1D*) fin->Get("hKstar_100_sys");
  TH1D * hsub_stat = (TH1D*) h100_stat->Clone("sub80100");
  TH1D * hm_stat[5];
  for (Int_t ic = 0; ic<5 ; ic++) {
    hm_stat[ic] = (TH1D*) fin->Get(Form("hKstar_%i",ic));
    //TH1D * hm_sys = (TH1D*) fin->Get(Form("hKstar_%i_sys", ic));
    hm_stat[ic]->SetTitle(Form("%i-%i", 20*ic, 20*(ic+1)));
  }
  
  for (int ipt = 1; ipt< h100_stat->GetXaxis()->GetNbins()+1; ipt++) {
    Double_t sub = h100_stat->GetBinContent(ipt);
    Double_t suberr = h100_stat->GetBinError(ipt);  
    sub *= 5;
    suberr *=5;
    for (Int_t ic = 0; ic<4 ; ic++) {
      sub -=  hm_stat[ic]->GetBinContent(ipt);
      suberr -=  hm_stat[ic]->GetBinError(ipt);
    }
    hsub_stat->SetBinContent(ipt, sub);
    hsub_stat->SetBinError(ipt, suberr);
  }
  
  hsub_stat->SetLineColor(kRed);
  hsub_stat->SetMarkerColor(kRed);

  TCanvas * ccorr = new TCanvas("corr","corrected spectra",600,700);
  ccorr->cd();
  gPad->SetLogy();
  
  gROOT->LoadMacro("$ASD/AddPaveText.C");
  TLegend * corrleg = new TLegend(0.4,0.55,0.99,0.8,"V0A Multiplicity Classes (Pb side)");
  corrleg->SetBorderSize(0);
  corrleg->SetFillColor(kWhite);
  corrleg->SetFillStyle(0);
  corrleg->AddEntry(hsub_stat, "80-100%% from subtraction","lp");
  
  for (Int_t ic = 0; ic<5 ; ic++) {
    corrleg->AddEntry(hm_stat[ic], Form("%i-%i%%", 20*ic, 20*(ic+1)),"lp");
    ccorr->cd();
    if (ic==0) hm_stat[ic]->Draw();
    else
      hm_stat[ic]->Draw("same");
  }
  hsub_stat->Draw("same");
  
  ccorr->cd();
  corrleg->Draw();
  AddPaveText_KStar_pPb("tr");
  

  return;
}
