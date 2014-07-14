void GetInvariantSpectra(TString infile = "prelim_kstar_pPb_smoothSys.root", Bool_t scale2plot = 0)
{
  TFile * fin = TFile::Open(infile.Data());
  
  TH1D * h100_stat = (TH1D*) fin->Get("hKstar_100");
  TH1D * h100_sys = (TH1D*) fin->Get("hKstar_100_sys");
  
  TH1D * h100_stat_inv =  (TH1D*) DivideBy2piPt(h100_stat);
  TH1D * h100_sys_inv =  (TH1D*) DivideBy2piPt(h100_sys);
  TH1D * hm_stat_inv[5];
  TH1D * hm_sys_inv[5];
  

  TCanvas * ccorr = new TCanvas("corr","corrected spectra",600,700);
  ccorr->cd();
  gPad->SetLogy();

  gROOT->LoadMacro("$ASD/AddPaveText.C");
  TLegend * corrleg = new TLegend(0.4,0.55,0.99,0.8,"V0A Multiplicity Classes (Pb side)");
  corrleg->SetBorderSize(0);
  corrleg->SetFillColor(kWhite);
  corrleg->SetFillStyle(0);
  
  for (Int_t ic = 0; ic<5 ; ic++) {
    TH1D * hm_stat = (TH1D*) fin->Get(Form("hKstar_%i",ic));
    TH1D * hm_sys = (TH1D*) fin->Get(Form("hKstar_%i_sys", ic));
    
    hm_stat_inv[ic] = (TH1D*) DivideBy2piPt(hm_stat); 
    hm_sys_inv[ic] = (TH1D*) DivideBy2piPt(hm_sys); 
    hm_stat_inv[ic]->SetTitle(Form("%i-%i", 20*ic, 20*(ic+1)));
    hm_stat_inv[ic]->GetYaxis()->SetTitle("1/(2#pi#it{p}_{T}) d^{2}#it{N}/(d#it{p}_{T}d#it{y}) (GeV/#it{c})^{-2}");
    hm_sys_inv[ic]->GetYaxis()->SetTitle("1/(2#pi#it{p}_{T}) d^{2}#it{N}/(d#it{p}_{T}d#it{y}) (GeV/#it{c})^{-2}");

    hm_sys_inv[ic]->GetYaxis()->SetRangeUser(1e-8, 1.);
    if (scale2plot) {
      Float_t scaleFactor = TMath::Power(2.0, (4-ic));
      hm_stat_inv[ic]->Scale(scaleFactor);
      hm_sys_inv[ic]->Scale(scaleFactor);
      hm_stat_inv[ic]->SetTitle(Form("%i-%i%% x%2.0f", 20*ic, 20*(ic+1), scaleFactor));
      hm_sys_inv[ic]->GetYaxis()->SetRangeUser(1e-8, 10.); 
    }

    corrleg->AddEntry(hm_stat_inv[ic], Form("%i-%i%% x%2.0f", 20*ic, 20*(ic+1), scaleFactor),"lp");
   
    ccorr->cd();
    if (ic==0)  hm_sys_inv[ic]->Draw("E2");
    else hm_sys_inv[ic]->Draw("E2same");
    hm_stat_inv[ic]->Draw("same");
  }
  
  ccorr->cd();
  corrleg->Draw();
  AddPaveText_KStar_pPb("tr");
  return;
}

TH1D * DivideBy2piPt( TH1D * h) 
{
  if (!h) return 0x0;
  TH1D * dummy = (TH1D *) h->Clone(Form("%s_inv",h->GetName()));
  for (int ipt = 1; ipt< h->GetXaxis()->GetNbins()+1; ipt++) {
    Double_t binCenter = h->GetBinCenter(ipt);
    Double_t factor = 2.0*TMath::Pi()*binCenter;
    dummy->SetBinContent(ipt, h->GetBinContent(ipt)/factor);
    dummy->SetBinError(ipt, h->GetBinError(ipt)/factor);
  }
  return dummy;
}
