void compareSyst2Stat(TString file = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TPC_ANA/ana2s/systematics/systUncert/syst1_yields_allCents_ResBg.root", Bool_t isTOF=0, Bool_t save=0)
{
  TFile * fin = TFile::Open(file.Data());
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas * c1 = new TCanvas ("c1","c1", 1000,700);
  c1->Divide(3,2);

  TH1D* hsys[5];   TH1D* hstat[5]; 
  for (Int_t ic = 0;ic<5; ic++){
    hsys[ic]=(TH1D*) fin->Get(Form("hSystVsPtPercentageOfCentral_%i",ic));
    hstat[ic]=(TH1D*) fin->Get(Form("hStatVsPtPercentageOfCentral_%i",ic));
    hsys[ic]->SetTitle(Form("syst. uncert.(%%), cent.bin %i",ic));
    hstat[ic]->SetTitle(Form("stat. uncert.(%%), cent.bin %i",ic));
    if (isTOF) hsys[ic]->GetXaxis()->SetRangeUser(1.0,8.0);
    else hsys[ic]->GetXaxis()->SetRangeUser(0.5,8.0);
    hsys[ic]->GetYaxis()->SetRangeUser(0.0,100.);
    hsys[ic]->SetLineWidth(3);
    hstat[ic]->SetFillStyle(3001);
    hstat[ic]->SetFillColor(kGray);
    if (isTOF)
      hstat[ic]->GetXaxis()->SetRangeUser(1.0,8.0);
    else
      hstat[ic]->GetXaxis()->SetRangeUser(0.5,8.0);
    hstat[ic]->GetYaxis()->SetRangeUser(0.0,100.);
    hstat[ic]->GetYaxis()->SetTitle("uncertainty (% of central value)");
    c1->cd(ic+1);
    hstat[ic]->Draw();
    hsys[ic]->Draw("same");
    TLegend* leg = ((TLegend*)gPad->BuildLegend(0.12,0.76,0.86,0.89));
    leg->SetBorderSize(0);
    leg->SetFillColor(kWhite);
  }
  if (save) {
    c1->SaveAs(file.ReplaceAll(".root","_compareStat.png"));
    c1->SaveAs(file.ReplaceAll(".png",".C"));
  }
  return;
}
