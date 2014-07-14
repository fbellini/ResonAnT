void preparePlotToShow_pA(TString det="TOF", TString date="13set13")
{
  //various
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *spectra = new TCanvas(Form("%sspectra",det.Data()),Form("K*^{0} yields - %s analysis",det.Data()), 600,700);
  TLegend * l = new TLegend(0.7,0.7,0.89,0.89);
  l->SetFillColor(kWhite);
  l->SetBorderSize(0);
  l->SetTextFont(42);
  
  for (Int_t j=0;j<5;j++){
    
    TFile * fin = TFile::Open(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/systematics/finalWsyst_24sett13_%i.root",det.Data(), j)); 
    
    TString statname = Form("h%sCorrected_%i",det.Data(), j);
    TString systname = Form("h%sCorrected_%i_syst",det.Data(),j);
    
    //read Kstar histos
    TH1D * hstat = (TH1D*) fin->Get(statname.Data());
    hstat->SetTitle(Form("%i-%i%% ", j*20,(j+1)*20));
    hstat->SetMarkerStyle(0);
    hstat->SetFillStyle(0);
    hstat->SetLineWidth(1);
    
    TH1D * hsyst = (TH1D*) fin->Get(systname.Data());
    //    hsyst->SetTitle(Form("%i-%i%% - syst. uncert.", j*20,(j+1)*20));
    hsyst->SetMarkerStyle(0);
    hsyst->SetFillStyle(0);
    hsyst->SetLineWidth(1);
    
    hsyst->GetYaxis()->SetRangeUser(1e-5,2.);
    spectra->cd();
    if (j>0) hsyst->Draw("E2 same");
    else hsyst->Draw("E2");
    hstat->Draw("same");
    l->AddEntry(hstat, Form("%s %i-%i%% ", det.Data(), j*20,(j+1)*20), "lpf");
  }  
  gROOT->LoadMacro("$ASD/AddPaveText.C");
  spectra->cd();
  gPad->SetLogy();
  l->Draw("same");
  AddPaveText_pPb_cent();
  AddPaveTextErrors();
  
  spectra->SaveAs(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/FINAL/FINAL_%s_%s.root",det.Data(),det.Data(),date.Data()));
  spectra->SaveAs(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/FINAL/FINAL_%s_%s.png",det.Data(),det.Data(),date.Data()));
  spectra->SaveAs(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/FINAL/FINAL_%s_%s.C",det.Data(),det.Data(),date.Data()));

  return;
}


void preparePlotToShow_pA100(TString det="TOF", TString date="13set13")
{
  //various
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *spectra = new TCanvas(Form("%sspectra",det.Data()),Form("K*^{0} yields - %s analysis",det.Data()), 600,700);
  TLegend * l = new TLegend(0.7,0.7,0.89,0.89);
  l->SetFillColor(kWhite);
  l->SetBorderSize(0);
  l->SetTextFont(42);
  
  for (Int_t j=0;j<1;j++){
    
    TFile * fin = TFile::Open(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s100/ana2s/systematics/systUncert/finalWsyst_16set13_%i.root",det.Data(), j)); 
    
    TString statname = Form("correzione%s",det.Data());
    TString systname = Form("correzione%s0_syst",det.Data());
    
    //read Kstar histos
    TH1D * hstat = (TH1D*) fin->Get(statname.Data());
    hstat->SetTitle(Form("%i-%i%% ", j*20,(j+1)*20));
    hstat->SetMarkerStyle(20);
    hstat->SetMarkerSize(0.8);
    hstat->SetFillStyle(0);
    hstat->SetLineWidth(1);
    
    TH1D * hsyst = (TH1D*) fin->Get(systname.Data());
    //    hsyst->SetTitle(Form("%i-%i%% - syst. uncert.", j*20,(j+1)*20));
    hsyst->SetMarkerStyle(0);
    hsyst->SetFillStyle(0);
    hsyst->SetLineWidth(1);
    
    hsyst->GetYaxis()->SetRangeUser(1e-5,2.);
    hsyst->GetYaxis()->SetTitle("1/N_{evt}* d^{2}N/dydp_{t} (GeV/c)^{-1}");
    hsyst->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    spectra->cd();
    if (j>0) hsyst->Draw("E2 same");
    else hsyst->Draw("E2");
    hstat->Draw("same");
    l->AddEntry(hstat, Form("%s %i-%i%% ", det.Data(), j*20,(j+1)*20), "lpf");
  }  
  gROOT->LoadMacro("$ASD/AddPaveText.C");
  spectra->cd();
  gPad->SetLogy();
  //l->Draw("same");
  AddPaveText_pPb0to100("tr");
  AddPaveTextErrors();
  
  spectra->SaveAs(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s100/ana2s/FINAL_%s_%s.root",det.Data(),det.Data(),date.Data()));
  spectra->SaveAs(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s100/ana2s/FINAL_%s_%s.png",det.Data(),det.Data(),date.Data()));
  spectra->SaveAs(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s100/ana2s/FINAL_%s_%s.C",det.Data(),det.Data(),date.Data()));

  return;
}
