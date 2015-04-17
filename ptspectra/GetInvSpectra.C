void GetInvariantSpectra(TString infile = "prelim_kstar_pPb_smoothSys.root", Bool_t scale2plot = 0)
{
  TFile * fin = TFile::Open(infile.Data());
  
  TH1D * h100_stat = (TH1D*) fin->Get("hKstar_100");
  TH1D * h100_sys = (TH1D*) fin->Get("hKstar_100_sys");
  
  TH1D * h100_stat_inv =  (TH1D*) DivideBy2piPt(h100_stat);
  TH1D * h100_sys_inv =  (TH1D*) DivideBy2piPt(h100_sys);
  TH1D * hm_stat_inv[5];
  TH1D * hm_sys_inv[5];
  
  h100_stat_inv->GetYaxis()->SetTitle("1/N_{evt}1/(2#pi#it{p}_{T}) d^{2}#it{N}/(d#it{p}_{T}d#it{y}) [(GeV/#it{c})^{-2}]");
  h100_sys_inv->GetYaxis()->SetTitle("1/N_{evt}1/(2#pi#it{p}_{T}) d^{2}#it{N}/(d#it{p}_{T}d#it{y}) [(GeV/#it{c})^{-2}]");
  h100_sys_inv->GetYaxis()->SetRangeUser(1e-8, 10.); 

  TLegend * corrleg2 = new TLegend(0.5,0.72,0.89,0.8);
  corrleg2->SetBorderSize(0);
  corrleg2->SetFillColor(kWhite);
  corrleg2->SetFillStyle(0);
  corrleg2->AddEntry(h100_stat_inv, Form("ALICE NSD"),"p");

  TCanvas * ccorr2 = new TCanvas("corr2","min bias corrected spectrum",800,600);
  ccorr2->cd();
  gPad->SetLogy();
  h100_sys_inv->Draw("E2");
  h100_stat_inv->Draw("same");
  corrleg2->Draw();
  AddPaveText_KStar_pPb("tr");

  //gROOT->LoadMacro("$ASD/AddPaveText.C");
  TLegend * corrleg = new TLegend(0.4,0.55,0.99,0.8,"V0A Multiplicity Event Classes (Pb side)");
  corrleg->SetBorderSize(0);
  corrleg->SetFillColor(kWhite);
  corrleg->SetFillStyle(0);
  
  TCanvas * ccorr = new TCanvas("corr","corrected spectra",800,600);
  ccorr->cd();
  gPad->SetLogy();
  //  Float_t scaleFactor = TMath::Power(2.0, (4-ic));
  Float_t scaleFactor[5] = {2., 1., 1./2.,1./4., 1/8.};

  for (Int_t ic = 0; ic<5 ; ic++) {
    TH1D * hm_stat = (TH1D*) fin->Get(Form("hKstar_%i",ic));
    TH1D * hm_sys = (TH1D*) fin->Get(Form("hKstar_%i_sys", ic));
    
    hm_stat_inv[ic] = (TH1D*) DivideBy2piPt(hm_stat); 
    hm_sys_inv[ic] = (TH1D*) DivideBy2piPt(hm_sys); 
    hm_stat_inv[ic]->SetTitle(Form("%i-%i", 20*ic, 20*(ic+1)));
    hm_stat_inv[ic]->GetYaxis()->SetTitle("1/N_{evt}1/(2#pi#it{p}_{T}) d^{2}#it{N}/(d#it{p}_{T}d#it{y}) [(GeV/#it{c})^{-2}]");
    hm_sys_inv[ic]->GetYaxis()->SetTitle("1/N_{evt}1/(2#pi#it{p}_{T}) d^{2}#it{N}/(d#it{p}_{T}d#it{y}) [(GeV/#it{c})^{-2}]");

    hm_sys_inv[ic]->GetYaxis()->SetRangeUser(1e-8, 1.);
    if (scale2plot) {
      hm_stat_inv[ic]->Scale(scaleFactor[ic]);
      hm_sys_inv[ic]->Scale(scaleFactor[ic]);
      hm_stat_inv[ic]->SetTitle(Form("%i-%i%% x%3.1f", 20*ic, 20*(ic+1), scaleFactor[ic]));
      hm_sys_inv[ic]->GetYaxis()->SetRangeUser(1e-8, 10.); 
    }
   
    ccorr->cd();
    if (ic==0)  hm_sys_inv[ic]->Draw("E2");
    else hm_sys_inv[ic]->Draw("E2same");
    hm_stat_inv[ic]->Draw("same");
  }
  
   corrleg->AddEntry(hm_stat_inv[0], Form("0-20%% x2"),"p");
   corrleg->AddEntry(hm_stat_inv[1], Form("20-40%% x1"),"p");
   corrleg->AddEntry(hm_stat_inv[2], Form("40-60%% x1/2"),"p");
   corrleg->AddEntry(hm_stat_inv[3], Form("60-80%% x1/4"),"p");
   corrleg->AddEntry(hm_stat_inv[4], Form("80-100%% x1/8"),"p");
  ccorr->cd();
  corrleg->Draw();
  AddPaveText_KStar_pPb("bl");
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

TPaveText * AddPaveText_KStar_pPb(TString position = "bl"){
  TString text="ALICE, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV";
  TString textks = "#frac{1}{2}(K*^{0}+#bar{K*^{0}}), -0.5 < #it{y} < 0";
  //TString textcent="V0A Multiplicity Event Classes (Pb side)";
  TPaveText * pave;
  if (position.Contains("tl")) pave = new TPaveText(0.12,0.75,0.45,0.89,"NDC");
  else
    if (position.Contains("bl")) pave = new TPaveText(0.22,0.22,0.55,0.36,"NDC");
    else 
      if (position.Contains("tr")) pave = new TPaveText(0.50,0.80,0.89,0.89,"NDC");
      else
	if (position.Contains("br")) pave = new TPaveText(0.55,0.12,0.89,0.26,"NDC");
  pave->SetTextColor(kBlack);
  pave->SetTextAlign(12);
  pave->SetTextFont(42);
  //pave->SetTextSize(18);
  pave->SetBorderSize(0);
  pave->SetFillColor(kWhite);
  pave->SetFillStyle(0);
  pave->InsertText(text.Data());
  pave->InsertText(textks.Data());
  //pave->InsertText(textcent.Data());
  pave->Draw();
  return pave;
}
