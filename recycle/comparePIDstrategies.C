void comparePIDstrategies()
{
  for (Int_t i=0;i<22;i++){
    comparePIDstrategies(i, "Like");
    comparePIDstrategies(i, "Mixing");
  }
}

void comparePIDstrategies(Int_t pT, TString bg = "Like")
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptTitle(0);
  Color_t col[4] = {kRed+1, kMagenta+2, kBlue+2, kGreen+2};
  TFile *_file0 = TFile::Open("_sum_EMnorm-1.00--1.00_2323_tpc2s_tof4sveto.root");
  TFile *_file1 = TFile::Open("_sum_EMnorm-1.00--1.00_2323_tpc3s_tof4sveto.root");
  TFile *_file2 = TFile::Open("_sum_EMnorm-1.00--1.00_2424_tpc2s_tof3sveto.root");
  TFile *_file3 = TFile::Open("_sum_EMnorm-1.00--1.00_2424_tpc3s_tof3sveto.root");

  TH1D * h0 = (TH1D*) _file0->Get(Form("sub_norm_%s_ptBin%02i_centBin00", bg.Data(), pT));
  h0->SetLineColor(col[0]);
  h0->SetMarkerColor(col[0]);
  h0->SetMarkerStyle(20);
  TH1D * h1 = (TH1D*) _file1->Get(Form("sub_norm_%s_ptBin%02i_centBin00", bg.Data(), pT));
  h1->SetLineColor(col[1]);
  h1->SetMarkerColor(col[1]);
  h1->SetMarkerStyle(27);
  TH1D * h2 = (TH1D*) _file2->Get(Form("sub_norm_%s_ptBin%02i_centBin00", bg.Data(), pT));
  h2->SetLineColor(col[2]);
  h2->SetMarkerColor(col[2]);
  h2->SetMarkerStyle(21);
  TH1D * h3 = (TH1D*) _file3->Get(Form("sub_norm_%s_ptBin%02i_centBin00", bg.Data(), pT));
  h3->SetLineColor(col[3]);
  h3->SetMarkerColor(col[3]);
  h3->SetMarkerStyle(28);

  //binning B
  Double_t pt[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00, 12.0, 15.0};
  
  //binning C
  //Double_t pt[] = {0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 12.0, 14.0, 16.0};

  TLegend * leg = new TLegend(0.5,0.65,0.89, 0.89, Form("%3.1f< p_{T} <%3.1f [0,100]%%, %s bg.", pt[pT], pt[pT+1], bg.Data()));
  leg->SetFillColor(kWhite); 
  leg->SetBorderSize(0);
  leg->AddEntry(h0,"tpc2s_tof4sveto","lp");
  leg->AddEntry(h1,"tpc3s_tof4sveto","lp");
  leg->AddEntry(h2,"tpc2s_tof3sveto","lp");
  leg->AddEntry(h3,"tpc3s_tof3sveto","lp");
  TCanvas * c1 = new TCanvas("c1","c1", 800,600);
  c1->cd();
  h1->Draw("");
  h0->Draw("same");
  h2->Draw("same");
  h3->Draw("same");
  leg->Draw();
  
  c1->Print(Form("compare_signal_PID_pt%i_%s.png", pT, bg.Data()));
}


void comparePIDstrategiesCorr()
{
  gStyle->SetOptTitle(0);
  TString path("/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/0to100_binB/systPID");
  TFile * _file[8];   TH1D * hcorr[8]; 
  Color_t color[8] = {kRed+1, kOrange, kMagenta+2, kPink+2, kBlue+2, kAzure+8, kGreen+2, kSpring+9};
  Int_t mark[8] = {20, 24, 27, 33, 21, 25, 28, 34};
  
  _file[0] = TFile::Open(Form("%s/tpc2s_tof3sveto/fitEM_norm1_BWpoly2_fixedW/CORRECTED_br_best_fit_poly2.root", path.Data()));
  _file[1] = TFile::Open(Form("%s/tpc2s_tof3sveto/fitLS_norm0_BWpoly2_fixedW/CORRECTED_br_best_fit_poly2.root", path.Data()));
  _file[2] = TFile::Open(Form("%s/tpc2s_tof4sveto/fitEM_norm1_BWpoly2_fixedW/CORRECTED_br_best_fit_poly2.root", path.Data()));
  _file[3] = TFile::Open(Form("%s/tpc2s_tof4sveto/fitLS_norm0_BWpoly2_fixedW/CORRECTED_br_best_fit_poly2.root", path.Data()));
  _file[4] = TFile::Open(Form("%s/tpc3s_tof3sveto/fitEM_norm1_BWpoly2_fixedW/CORRECTED_br_best_fit_poly2.root", path.Data()));
  _file[5] = TFile::Open(Form("%s/tpc3s_tof3sveto/fitLS_norm0_BWpoly2_fixedW/CORRECTED_br_best_fit_poly2.root", path.Data()));
  _file[6] = TFile::Open(Form("%s/tpc3s_tof4sveto/fitEM_norm1_BWpoly2_fixedW/CORRECTED_br_best_fit_poly2.root", path.Data()));
  _file[7] = TFile::Open(Form("%s/tpc3s_tof4sveto/fitLS_norm0_BWpoly2_fixedW/CORRECTED_br_best_fit_poly2.root", path.Data()));
  
  for (Int_t j=0;j<8;j++){//loop on files
    Bool_t isLS = (j%2);
    hcorr[j] = (TH1D*) _file[j]->Get("hCorrected_0");
    hcorr[j]->SetLineColor(color[j]);
    hcorr[j]->SetMarkerColor(color[j]);
    hcorr[j]->SetMarkerStyle(mark[j]);
  }

  hcorr[0]->SetTitle("EM_tpc2s_tof3sveto");
  hcorr[1]->SetTitle("LS_tpc2s_tof3sveto");
  hcorr[2]->SetTitle("EM_tpc2s_tof4sveto");
  hcorr[3]->SetTitle("LS_tpc2s_tof4sveto");
  hcorr[4]->SetTitle("EM_tpc3s_tof3sveto");
  hcorr[5]->SetTitle("LS_tpc3s_tof3sveto");
  hcorr[6]->SetTitle("EM_tpc3s_tof4sveto");
  hcorr[7]->SetTitle("LS_tpc3s_tof4sveto");

  TH1D * ls2em[4]; 
  TCanvas *cls2em = new TCanvas("cls2em","cls2em",1000,700);
  for (Int_t j=0;j<4;j++){//loop on pairs
    ls2em[j] = (TH1D*) hcorr[j*2+1]->Clone(Form("ls2em_%i",j)); 
    ls2em[j]->Divide(ls2em[j], hcorr[j*2], 1, 1,"sumw2");
    ls2em[j]->SetLineColor(color[j*2]);
    ls2em[j]->SetMarkerColor(color[j*2]);
    ls2em[j]->SetMarkerStyle(mark[j*2]);
    ls2em[j]->SetTitle(Form("%s / %s", hcorr[j*2+1]->GetTitle(),hcorr[j*2]->GetTitle()));
    ls2em[j]->GetYaxis()->SetRangeUser(0.5,1.5);
    ls2em[j]->GetYaxis()->SetTitle("ratio");
    cls2em->cd();
    if(j==0) ls2em[j]->Draw();
    else 
      ls2em[j]->Draw("same");
  }
  cls2em->cd();
  TLegend* leg2ls = (TLegend*) gPad->BuildLegend();
  leg2ls->SetFillColor(kWhite);
  

  TH1D * r2pid[8];
  TCanvas *cr2pid = new TCanvas("cr2pid","cr2pid",1000,700);
  for (Int_t j=0;j<8;j++){//loop on pairs
    r2pid[j] = (TH1D*) hcorr[j]->Clone(Form("pid_%i",j)); 
    r2pid[j]->Divide(r2pid[j], hcorr[0], 1,1,"sumw2");
    r2pid[j]->SetLineColor(color[j]);
    r2pid[j]->SetMarkerColor(color[j]);
    r2pid[j]->SetMarkerStyle(mark[j]);
    r2pid[j]->SetTitle(Form("%s / %s", hcorr[j]->GetTitle(),hcorr[0]->GetTitle()));
    r2pid[j]->GetYaxis()->SetRangeUser(0.5,1.5);
    r2pid[j]->GetYaxis()->SetTitle("ratio");   
    cr2pid->cd();
    if(j==0) r2pid[j]->Draw();
    else r2pid[j]->Draw("same");
  }
   cr2pid->cd();
   TLegend*leg2pid = (TLegend*) gPad->BuildLegend();
   leg2pid->SetFillColor(kWhite);
  
}
