void displayKsVsAntiKsAfterBgSub(TString filename = "_kstar_EMnorm1.30-1.50_cut1717_train215-216.root", Bool_t isLS = 1, TString pngsuffix = "tpc2s")
{
  gROOT->LoadMacro("/Users/bellini/alice/macro/SetGraphicStyle.C");
  //SetGraphicStyle(0,1,0);
  TGaxis::SetMaxDigits(3);
  Double_t pt[] = { 0.0, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };

  TFile * fks = TFile::Open(filename.Data());
  if (!fks) return;

  TFile * faks = TFile::Open(filename.ReplaceAll("_kstar_","_antikstar_"));
  if (!faks) return;
  
  TCanvas * canv[5];
  for (Int_t c=0;c<5;c++){
    
    canv[c] = new TCanvas(Form("c%i",c),Form("c%i",c),800,1000);
    canv[c]->Divide(3,4);
    
    for (Int_t p=1;p<13;p++){
      TH1D * hks = (TH1D*) fks->Get(Form("sub_norm_%s_ptBin%02i_centBin%02i",(isLS?"Like":"Mixing"),p,c));
      TString titleks(Form("%s",hks->GetTitle())); 
      // titleks.ReplaceAll("S+res.Bg (EM):","");
      // titleks.ReplaceAll("central","");
      hks->SetTitle(Form("%2.1f<p_{T}<%2.1f GeV/c (%i-%i%%)",pt[p],pt[p+1], c*20, (c+1)*20));
      hks->GetXaxis()->SetTitleSize(0.05);
      hks->GetXaxis()->SetLabelSize(0.06);
      hks->GetYaxis()->SetTitleSize(0.05);
      hks->GetYaxis()->SetLabelSize(0.06);
      hks->SetLineWidth(1);
      hks->SetMarkerColor(kRed);
      hks->SetLineColor(kRed);
      hks->SetMarkerStyle(20);
      hks->SetFillColor(kRed-10);
      hks->SetFillStyle(0);
      // hks->SetLineColor(kBlack);
      // hks->SetMarkerColor(kBlack);
      // hks->SetMarkerStyle(0);
      // hks->SetLineWidth(1);
      //HistoMakeUp(hks, 0);
      hks->Rebin(2);
      
      TH1D * haks = (TH1D*) faks->Get(Form("sub_norm_%s_ptBin%02i_centBin%02i",(isLS?"Like":"Mixing"), p,c));
      TString titleaks(Form("%s",haks->GetTitle())); 
      // titleaks.ReplaceAll("S+res.Bg (EM):","");
      // titleaks.ReplaceAll("central","");
      // haks->SetLineColor(kAzure+7);
      // haks->SetMarkerColor(kAzure+7);
      // haks->SetMarkerStyle(0);
      // haks->SetLineWidth(1);
      haks->SetTitle(Form("%2.1f<p_{T}<%2.1f GeV/c (%i-%i%%)",pt[p],pt[p+1], c*20, (c+1)*20));
      haks->GetXaxis()->SetTitleSize(0.05);
      haks->GetXaxis()->SetLabelSize(0.06);
      haks->GetYaxis()->SetTitleSize(0.05);
      haks->GetYaxis()->SetLabelSize(0.06);
      haks->SetLineWidth(1);
      haks->SetMarkerColor(kBlue+1);
      haks->SetLineColor(kBlue+1);
      haks->SetMarkerStyle(25);
      haks->SetFillColor(kBlue-10);
      haks->SetFillStyle(0);
      haks->Rebin(2);
      
      TLegend * cleg = new TLegend(0.6,0.7,0.89,0.89);
      cleg->AddEntry(hks,"K*#rightarrowK^{+}#pi^{-}","lp");
      cleg->AddEntry(haks,"#bar{K*}#rightarrowK^{-}#pi^{+}","lp");
      cleg->SetFillColor(kWhite);
      cleg->SetBorderSize(0);
      
      canv[c]->cd(p);
      //      hks->GetYaxis()->SetRangeUser(-1e3,1e4);
      hks->Draw(); haks->Draw("same");
      cleg->Draw();
    }
    canv[c]->SaveAs(Form("KsVsAntiKs_After%sBgSub_%s_c%i.png",(isLS?"Like":"Mixing"),pngsuffix.Data(),c));
  }
  return;
}
