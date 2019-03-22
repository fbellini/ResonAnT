void getEffVsMassPt()
{
  TFile * fin = TFile::Open("train1220.root");
  TList * lks = (TList*) fin->Get("RsnOut_tpc2s_tof3sveto");//quality2011, tof2s, 1818
  Int_t icut = 2424;//1818;
  TCanvas * c1 = new TCanvas("c1","c1",900,600);
  c1->Divide(2,2);
  TCanvas * c2 = new TCanvas("c2","c2",900,900);
 
  gROOT->LoadMacro("$ASD/SetGraphicStyle.C"); SetGraphicStyle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  
  THnSparseT<TArrayF> * hks_true = (THnSparseT<TArrayF> *) lks->FindObject(Form("TOFKStarPbPbMC_%i_kstar_TruesPM",icut));  
  TH2D * hks_true_MvsPt = (TH2D*)  hks_true->Projection(1,0);
  hks_true_MvsPt->SetTitle("True reconstructed K*");
  hks_true_MvsPt->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
  hks_true_MvsPt->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  hks_true_MvsPt->GetXaxis()->SetRangeUser(0.7,1.1);
  hks_true_MvsPt->GetXaxis()->SetLabelSize(0.05);
  hks_true_MvsPt->GetXaxis()->SetTitleSize(0.05);
  hks_true_MvsPt->GetYaxis()->SetLabelSize(0.05);
  hks_true_MvsPt->GetYaxis()->SetTitleSize(0.05);

  c1->cd(1);
  hks_true_MvsPt->Draw("colz");
  gPad->SetLogz();
  
  THnSparseT<TArrayF> * hks_gen = (THnSparseT<TArrayF> *) lks->FindObject(Form("TOFKStarPbPbMC_%i_Ks_Mother",icut));
  TH2D * hks_gen_MvsPt = (TH2D*)  hks_gen->Projection(1,0);
  hks_gen_MvsPt->SetTitle("Generated K*");
  hks_gen_MvsPt->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
  hks_gen_MvsPt->GetXaxis()->SetRangeUser(0.6,1.2);
  hks_gen_MvsPt->GetXaxis()->SetLabelSize(0.05);
  hks_gen_MvsPt->GetXaxis()->SetTitleSize(0.05);
  hks_gen_MvsPt->GetYaxis()->SetLabelSize(0.05);
  hks_gen_MvsPt->GetYaxis()->SetTitleSize(0.05);
  hks_gen_MvsPt->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  c1->cd(2);
  hks_gen_MvsPt->Draw("colz");
  gPad->SetLogz();
  
  TH2D * hks_eff = (TH2D*)  hks_true_MvsPt->Clone("eff_tof2s");  
  hks_eff->Divide(hks_gen_MvsPt);
  hks_eff->SetTitle("True reconstructed/generated");
  hks_eff->GetXaxis()->SetRangeUser(0.6,1.2);
  hks_eff->GetXaxis()->SetLabelSize(0.05);
  hks_eff->GetXaxis()->SetTitleSize(0.05);
  hks_eff->GetYaxis()->SetLabelSize(0.05);
  hks_eff->GetYaxis()->SetTitleSize(0.05);

  c1->cd(3);
  hks_eff->Draw("colz");
  gPad->SetLogz();

  //eff vs mass
  Double_t pt[] = {0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 12.0, 14.0, 16.0};
  for (Int_t i=0;i<16; i++){
    
    Float_t ptmin = pt[i];
    Float_t ptmax = pt[i+1];
    Int_t binmin = hks_true->GetAxis(1)->FindBin(ptmin);
    Int_t binmax = hks_true->GetAxis(1)->FindBin(ptmax);
    
    hks_true->GetAxis(1)->SetRange(binmin,binmax);
    hks_gen->GetAxis(1)->SetRange(binmin,binmax);
    
    TH1D * hks_true_M = (TH1D*)  hks_true->Projection(0);
    hks_true_M->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    hks_true_M->GetXaxis()->SetRangeUser(0.6,1.2);
    
    TH1D * hks_gen_M = (TH1D*)  hks_gen->Projection(0);
    hks_gen_M->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    hks_gen_M->GetXaxis()->SetRangeUser(0.6,1.2);
    
    TH1D * hks_eff_M = (TH1D*) hks_true_M->Clone(Form("effVsM_%i",i));
    hks_eff_M->SetTitle(Form("%2.1f-%2.1f GeV/c",ptmin,ptmax));
    hks_eff_M->Divide(hks_gen_M);
    hks_eff_M->GetYaxis()->SetRangeUser(1e-3,1.);
    hks_eff_M->GetYaxis()->SetTitle("efficiency");
    hks_eff_M->GetXaxis()->SetLabelSize(0.05);
    hks_eff_M->GetXaxis()->SetTitleSize(0.05);
    hks_eff_M->GetYaxis()->SetLabelSize(0.05);
    hks_eff_M->GetYaxis()->SetTitleSize(0.05);

    if (i<6) hks_eff_M->SetLineColor(kBlue+4-i);
    else hks_eff_M->SetLineColor(kCyan+3-i);
    hks_eff_M->SetLineWidth(2);

    c1->cd(4);
    if (i==0) hks_eff_M->Draw();
    else hks_eff_M->Draw("same");
    c2->cd();
    if (i==0) hks_eff_M->Draw();
    else hks_eff_M->Draw("same");
    //hks_gen_M->Draw();
  }
  c1->cd(4);
  //gPad->SetLogy();
  gPad->BuildLegend(0.75,0.62,0.99, 0.99)->SetFillColor(kWhite);
  
  return;
}
