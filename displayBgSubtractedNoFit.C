Color_t color1[4][5]={ kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-5, //tpc
		       kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2,//tof
		       kOrange+7, kPink+6, kGreen+1, kAzure+1, kBlue+4,//combined
		       kPink+4, kYellow+2, kSpring-6, kAzure+7, kBlue-3};//other 
Color_t color2[4][5]={ kTeal-5, kSpring+4, kAzure+9, kCyan+3, kBlue-5, //tpc
		       kRed-7, kOrange-7, kViolet-2, kPink-4, kViolet-9, //tof
		       kOrange-3, kPink+7, kCyan-3, kAzure+3, kBlue,//combined
		       kViolet-2, kPink+7, kGreen+1, kCyan-5, kBlue+2};//other

Double_t ptB[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00, 12.0, 15.0};

void displayBgSubtractedNoFit(TString filename, Int_t minpt1, Int_t minpt2, Bool_t isTOF, TString bg ="Like")
{
  TGaxis::SetMaxDigits(4);
  TFile * fin = TFile::Open(filename.Data());
  if (!fin) return;
  TH1D * h[10];
  TCanvas * canv[5];
  gStyle->SetOptTitle(0);
  for (int cc=0;cc<5;cc++){
    h[cc*2] = (TH1D*) fin->Get(Form("sub_norm_Like_ptBin0%i_centBin0%i", minpt1, cc));
    h[cc*2+1] = (TH1D*) fin->Get(Form("sub_norm_Like_ptBin0%i_centBin0%i", minpt2, cc));
    h[cc*2]->GetXaxis()->SetTitle("#it{M}_{K#pi} (GeV/#it{c^{2}})");
    if (minpt1<=10) 
      h[cc*2]->SetTitle(Form("%3.2f<#it{p}_{T}<%3.2f (GeV/#it{c}) (%i-%i%%)", 0.5*(minpt1-1),0.5*(minpt1), cc*20, (cc+1)*20));
    
    h[cc*2]->GetXaxis()->SetRangeUser(0.7,1.2);
    h[cc*2]->GetYaxis()->SetTitle("Counts / (0.01 GeV/#it{c^{2}})");
    h[cc*2]->Scale(2.); 
    h[cc*2]->SetLineColor(color1[isTOF][cc]);
    h[cc*2]->SetMarkerColor(color1[isTOF][cc]);
    
    h[cc*2+1]->GetXaxis()->SetTitle("#it{M}_{K#pi} (GeV/#it{c^{2}})");
    h[cc*2+1]->GetXaxis()->SetRangeUser(0.7,1.2);
    h[cc*2+1]->GetYaxis()->SetTitle("Counts / (0.01 GeV/#it{c^{2}})");
    if (minpt2<=10) 
      h[cc*2+1]->SetTitle(Form("%3.2f<#it{p}_{T}<%3.2f GeV/#it{c} (%i-%i%%)", 0.5*(minpt2-1),0.5*(minpt2), cc*20, (cc+1)*20));
    h[cc*2+1]->Scale(2.); 
    h[cc*2+1]->SetLineColor(color2[isTOF][cc]);
    h[cc*2+1]->SetMarkerColor(color2[isTOF][cc]);
 
    canv[cc] = new TCanvas(Form("cent%i",cc),Form("cent%i",cc), 600, 500);
    canv[cc]->cd();
    gPad->SetLogy();
    h[cc*2]->GetYaxis()->SetRangeUser(100.,1e6);
    h[cc*2+1]->GetYaxis()->SetRangeUser(100.,1e6);
    h[cc*2]->Draw();
    h[cc*2+1]->Draw("same");
    TLegend *leg;
    TPaveText * pave;
    if (cc>2) {
      pave = new TPaveText(0.55,0.83,0.89,0.89,"NDC");
      leg = (TLegend*) canv[cc]->BuildLegend(0.4,0.68,0.89,0.83);
    } else {
      pave = new TPaveText(0.55,0.12,0.89,0.18,"NDC");
      leg = (TLegend*) canv[cc]->BuildLegend(0.4,0.18,0.89,0.33);
    }    
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);
    pave->SetBorderSize(0);
    pave->SetFillColor(kWhite);
    pave->InsertText("Like-sign bckgd subtracted");
    pave->Draw();
  }
}

void show_distribs(TString filename, Int_t cc = 0, Int_t minpt1, Bool_t is0100 = 0, Bool_t isSavePng = 1)
{
  TGaxis::SetMaxDigits(3);
  TFile * fin = TFile::Open(filename.Data());
  if (!fin) return;
  TH1D * h[3];
  TCanvas * canv;
  gStyle->SetOptTitle(0);
  if (minpt1<=10) {
    h[0] = (TH1D*) fin->Get(Form("Signal_ptBin%02i_centBin%02i", minpt1, cc));
    h[1] = (TH1D*) fin->Get(Form("norm_Mixing_ptBin%02i_centBin%02i", minpt1, cc));
    h[2] = (TH1D*) fin->Get(Form("norm_Like_ptBin%02i_centBin%02i", minpt1, cc));
  } else {
    h[0] = (TH1D*) fin->Get(Form("Signal_ptBin%02i_centBin%02i", minpt1, cc));
    h[1] = (TH1D*) fin->Get(Form("norm_Mixing_ptBin%02i_centBin%02i", minpt1, cc));
    h[2] = (TH1D*) fin->Get(Form("norm_Like_ptBin%02i_centBin%02i", minpt1, cc));
  }
  
  h[0]->SetTitle("Unlike-Sign Pairs");
  h[1]->SetTitle("Mixed Event Background");
  h[2]->SetTitle("Like-Sign Background ");

  for (Int_t j=0;j<3;j++){
    h[j]->GetXaxis()->SetTitle("#it{M}_{K#pi} (GeV/#it{c}^{2})");
    h[j]->GetXaxis()->SetTitleOffset(1.30);
    h[j]->GetXaxis()->SetRangeUser(0.6,1.5);
    h[j]->GetYaxis()->SetTitle("Counts / (0.01 GeV/#it{c}^{2})");
    h[j]->GetYaxis()->SetTitleOffset(1.4);
  }
  h[0]->SetLineColor(kBlack);
  h[0]->SetLineWidth(2);
  h[0]->SetMarkerColor(kBlack);
  h[0]->SetMarkerStyle(20);
  h[0]->GetYaxis()->SetTitleOffset(1.3);
  h[0]->GetYaxis()->SetTitleSize(0.05);
  h[0]->GetYaxis()->SetNdivisions(510,kTRUE);
  h[0]->GetXaxis()->SetTitleOffset(1.3);
  h[0]->GetXaxis()->SetTitleSize(0.05);


  h[1]->SetLineColor(kRed);
  h[1]->SetMarkerColor(kRed);
  h[1]->SetMarkerSize(0.7);
  h[1]->SetMarkerStyle(25);

  h[2]->SetLineColor(kBlue+1);
  h[2]->SetMarkerColor(kBlue+1);
  h[2]->SetMarkerSize(0.7);
  h[2]->SetMarkerStyle(24);

  canv = new TCanvas(Form("cent%i",cc),Form("cent%i",cc), 800, 800);
  canv->cd();
  //gPad->SetLogy();
  // h[cc*2]->GetYaxis()->SetRangeUser(100.,5e4);
  // h[cc*2+1]->GetYaxis()->SetRangeUser(100.,5e4);
  h[0]->Draw();
  h[1]->Draw("same");
  h[2]->Draw("same");
  
  TLegend *leg;
  TPaveText * pave;
  //if (cc>2) {
  leg = (TLegend*) canv->BuildLegend(0.6,0.68,0.89,0.89);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  
  pave = new TPaveText(0.45,0.2,0.89,0.4,"NDC");
  pave->SetBorderSize(0);
  pave->SetFillColor(kWhite);
  pave->SetTextAlign(12);

  Float_t ptoffset = 0.5;//=0.5 for TPC and TOF stdalone
  if (is0100) {
    pave->InsertText("p-Pb  #sqrt{#it{s}_{NN}} = 5.02 TeV");
    pave->InsertText("K*^{0}#rightarrow K#pi, 0< #it{y} <0.5");  
    pave->InsertText("V0A Multiplicity Class 0-100%");
    pave->InsertText(Form("%3.1f#leq#it{p}_{T}<%3.1f GeV/#it{c}", ptB[minpt1], ptB[minpt1+1]));
    // if (minpt1==0) pave->InsertText(Form("0.0<p_{T}<0.3 GeV/#it{c} (0-100%%)"));
    // if (minpt1==1) pave->InsertText(Form("0.3<p_{T}<0.5 GeV/#it{c} (0-100%%)"));
    // if (minpt1>1 && minpt1<=10)
    //   pave->InsertText(Form("%3.2f<p_{T}<%3.2f GeV/#it{c} (0-100%%)", ptoffset+0.5*(minpt1-2), ptoffset+0.5*(minpt1-1)));
    // if (minpt1==11) pave->InsertText(Form("5.0<p_{T}<6.0 GeV/#it{c} (0-100%%)"));
    // if (minpt1==12) pave->InsertText(Form("6.0<p_{T}<7.0 GeV/#it{c} (0-100%%)"));
    // if (minpt1==13) pave->InsertText(Form("7.0<p_{T}<8.0 GeV/#it{c} (0-100%%)"));
    // if (minpt1==14) pave->InsertText(Form("8.0<p_{T}<10.0 GeV/#it{c} (0-100%%)"));
  } else {
    // pave->InsertText(Form("%3.2f<p_{T}<%3.2f GeV/#it{c} (%i-%i%%)", 0.5*(minpt1-1),0.5*(minpt1), cc*20, (cc+1)*100));
    pave->InsertText(Form("%3.1f<#it{p}_{T}<%3.1f GeV/#it{c} (%i-%i%%)", ptB[minpt1], ptB[minpt1+1], cc*20, (cc+1)*100));
  }
  pave->SetFillStyle(0);
  pave->Draw();
  if (isSavePng) canv->Print(Form("exampleDistrTpcTofVeto_pt%i.png",minpt1));
}


void show_bgSubLSvsEM_noFit(TString filename, Int_t minpt1, Int_t iPid, Bool_t is0100 = 0, Bool_t isSavePng=1)
{
  TGaxis::SetMaxDigits(4);
  TFile * fin = TFile::Open(filename.Data());
  if (!fin) return;
  TH1D * h[15];
  TCanvas * canv[5];
  gStyle->SetOptTitle(0);
  
  for (int cc=0;cc<5;cc++){
    if (is0100 && cc>0) continue;    
    if (minpt1<=10) {
      h[cc*2] = (TH1D*) fin->Get(Form("sub_norm_Like_ptBin0%i_centBin0%i", minpt1, cc));
      h[cc*2+1] = (TH1D*) fin->Get(Form("sub_norm_Mixing_ptBin0%i_centBin0%i", minpt1, cc));
    } else {
      h[cc*2] = (TH1D*) fin->Get(Form("sub_norm_Like_ptBin%i_centBin0%i", minpt1, cc));
      h[cc*2+1] = (TH1D*) fin->Get(Form("sub_norm_Mixing_ptBin%i_centBin0%i", minpt1, cc));
    }
    h[cc*2]->SetTitle("Like-sign");
    h[cc*2]->GetXaxis()->SetTitle("#it{M}_{K#pi} (GeV/#it{c^{2}})");
    h[cc*2]->GetXaxis()->SetTitleOffset(1.30);
      
    // h[cc*2]->GetXaxis()->SetRangeUser(0.7,1.2);
    h[cc*2]->GetXaxis()->SetRangeUser(0.6,1.5);
    h[cc*2]->GetYaxis()->SetTitle("Counts / (0.01 GeV/#it{c^{2}})");
    h[cc*2]->GetYaxis()->SetTitleOffset(1.4);
    // h[cc*2]->Scale(2.); 
    // h[cc*2]->SetLineColor(color1[iPid][cc]);
    // h[cc*2]->SetMarkerColor(color1[iPid][cc]);
    h[cc*2]->SetLineColor(kBlack);
    h[cc*2]->SetMarkerColor(kBlack);
    h[cc*2+1]->SetMarkerStyle(20);
    
    h[cc*2+1]->GetXaxis()->SetTitle("#it{M}_{K#pi} (GeV/#it{c^{2}})");
    //h[cc*2+1]->GetXaxis()->SetRangeUser(0.7,1.2);
    h[cc*2+1]->GetXaxis()->SetRangeUser(0.6,1.5);
    h[cc*2+1]->GetXaxis()->SetTitleOffset(1.3);
    h[cc*2+1]->GetYaxis()->SetTitle("Counts / (0.01 GeV/#it{c^{2}})");
    h[cc*2+1]->GetYaxis()->SetTitleOffset(1.4);
    h[cc*2+1]->SetTitle("Event-mixing");
    //h[cc*2+1]->Scale(2.); 
    h[cc*2+1]->SetLineColor(color1[iPid][cc]);
    h[cc*2+1]->SetMarkerColor(color1[iPid][cc]);
    h[cc*2+1]->SetMarkerStyle(24);

    canv[cc] = new TCanvas(Form("cent%i",cc),Form("cent%i",cc), 600, 500);
    canv[cc]->cd();
    //gPad->SetLogy();
    // h[cc*2]->GetYaxis()->SetRangeUser(100.,5e4);
    // h[cc*2+1]->GetYaxis()->SetRangeUser(100.,5e4);
    h[cc*2+1]->Draw();
    h[cc*2]->Draw("same");
    TLegend *leg;
    TPaveText * pave;
    //if (cc>2) {
    pave = new TPaveText(0.45,0.83,0.89,0.89,"NDC");
    leg = (TLegend*) canv[cc]->BuildLegend(0.6,0.68,0.87,0.83);
    // } else {
    //   pave = new TPaveText(0.55,0.12,0.89,0.18,"NDC");
    //   leg = (TLegend*) canv[cc]->BuildLegend(0.4,0.18,0.89,0.33);
    // }    
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);
    pave->SetBorderSize(0);
    pave->SetFillColor(kWhite);
    Float_t ptoffset = 0.5;//=0.5 for TPC and TOF stdalone
    if (is0100) {
      // if (minpt1==0) pave->InsertText(Form("0.0<p_{T}<0.3 GeV/#it{c} (0-100%%)"));
      // if (minpt1==1) pave->InsertText(Form("0.3<p_{T}<0.5 GeV/#it{c} (0-100%%)"));
      // if (minpt1>1 && minpt1<=10)
      // 	pave->InsertText(Form("%3.2f<p_{T}<%3.2f GeV/#it{c} (0-100%%)", ptoffset+0.5*(minpt1-2), ptoffset+0.5*(minpt1-1)));
      // if (minpt1==11) pave->InsertText(Form("5.0<p_{T}<6.0 GeV/#it{c} (0-100%%)"));
      // if (minpt1==12) pave->InsertText(Form("6.0<p_{T}<7.0 GeV/#it{c} (0-100%%)"));
      // if (minpt1==13) pave->InsertText(Form("7.0<p_{T}<8.0 GeV/#it{c} (0-100%%)"));
      // if (minpt1==14) pave->InsertText(Form("8.0<p_{T}<10.0 GeV/#it{c} (0-100%%)"));
      pave->InsertText(Form("%3.1f<p_{T}<%3.1f GeV/#it{c} (0-100%%)", ptB[minpt1], ptB[minpt1+1]));
    } else {
      //pave->InsertText(Form("%3.2f<p_{T}<%3.2f GeV/#it{c} (%i-%i%%)", 0.5*(minpt1-1),0.5*(minpt1), cc*20, (cc+1)*100));
      pave->InsertText(Form("%3.1f<p_{T}<%3.1f GeV/#it{c} (%i-%i%%)", ptB[minpt1], ptB[minpt1+1], cc*20, (cc+1)*100));
    }
    pave->Draw();
    if (isSavePng) canv[cc]->Print(Form("example_bgSub_LSvsEM_pt%02i_tpctofveto.png",minpt1));
  }
}
