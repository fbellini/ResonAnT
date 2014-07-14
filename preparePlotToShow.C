void preparePlotToShow_pA(TString det="TOF")
{
  //various
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *spectra = new TCanvas(Form("%sspectra",det.Data()),Form("K*^{0} yields - %s analysis",det.Data()), 600,700);
  TLegend * l = new TLegend(0.5,0.75,0.89,0.89);
  l->SetFillColor(kWhite);
  l->SetBorderSize(0);
  l->SetTextFont(42);
  
  for (Int_t j=0;j<5;j++){
    
    TFile * fin = TFile::Open(Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/systematics/finalWsyst_13set13_%i.root",det.Data(), c)); 
    
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
    
    hSyst->GetYaxis()->SetRangeUser(1e-6,2.);
    spectra->cd();
    hsyst->Draw("E2");
    hstat->Draw("same");
    l->AddEntry(hstat, Form("%i-%i%% ", j*20,(j+1)*20), "lpf");
  }  
  spectra->cd();
  l->Draw("same");
  return;
}
