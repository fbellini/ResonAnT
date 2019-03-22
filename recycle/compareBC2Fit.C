 void compareBC2Fit(Int_t icent = 0)
{
  Color_t color[2][6]={kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack, //tpc
		       kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof
  Int_t marker[2][6]={21, 22, 23, 34, 33, 20, //tpc
		      25, 26, 32, 28, 27, 24}; //tof
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  Bool_t isTOF = 0;
  TFile *_file0 = TFile::Open(Form("best_fit_poly2_cent%i.root",icent));
  //  gROOT->LoadMacro("GetPlotRatio.C");
  gROOT->LoadMacro("AddPaveText.C");
  
  TH1D * bincount = (TH1D*) _file0->Get("bincount");
  bincount->SetLineColor(kBlack);
  bincount->SetMarkerColor(kBlack);
  bincount->SetMarkerStyle(28);
  bincount->GetXaxis()->SetRangeUser(0.5,8.0);
  bincount->GetYaxis()->SetTitle("raw counts (a.u.)");
  bincount->SetTitle(Form("bin counting"));
  TH1D * raw = (TH1D*) _file0->Get("raw");
  raw->SetLineColor(color[isTOF][icent]);
  raw->SetMarkerColor(color[isTOF][icent]);
  raw->SetMarkerStyle(marker[isTOF][icent]);
  raw->GetXaxis()->SetRangeUser(0.5,8.0);
  raw->GetYaxis()->SetTitle("raw counts (a.u.)");
  raw->GetYaxis()->SetTitleSize(0.06);
  raw->GetYaxis()->SetTitleOffset(0.8);
  raw->GetYaxis()->SetLabelSize(0.06);
  raw->GetYaxis()->SetTitleSize(0.06);   
  raw->SetTitle(Form("Breit Wigner fit"));


  TH1D * rsw2 = (TH1D*) bincount->Clone("ratio_sumw2");
  TH1D * rbin = (TH1D*) bincount->Clone("ratio_binomial");
  rbin->SetTitle("bin counting / Breit-Wigner fit (stat. err., binomial)");
  rsw2->Divide(rsw2, raw,1.,1.); // sumw2 errors
  rbin->Divide(rbin, raw,1.,1.,"B"); //binomial errors

  //cosmetics
  rbin->SetLineColor(kBlue+2);
  rbin->SetMarkerColor(kBlue+2);
  rbin->SetFillStyle(3001);
  rbin->SetFillColor(kBlue-5);
  rbin->SetMarkerStyle(marker[isTOF][icent]);
  rbin->GetYaxis()->SetRangeUser(0.8,1.2);
  rbin->GetXaxis()->SetRangeUser(0.5,8.0);
  rbin->GetYaxis()->SetLabelSize(0.05);
  rbin->GetYaxis()->SetTitle("ratio");
  rbin->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  rbin->GetYaxis()->SetRangeUser(0.68,1.32);
  rbin->GetYaxis()->SetTitleSize(0.05);
  rbin->GetYaxis()->SetTitleOffset(0.8);
  rbin->GetXaxis()->SetLabelSize(0.05);
  rbin->GetXaxis()->SetTitleSize(0.05);   

  TCanvas * cr = new TCanvas("cr","compare",600, 600);
  cr->Divide(1,2,0.0002,0.0002); 
  TPad * pad1 = (TPad*) cr->GetPad(1);
  pad1->SetFillColor(kWhite);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(0);
  pad1->SetMargin(0.1,0.05,0.005,0.05);
  // pad1->SetTicky(2);
  //pad1->SetGridx();
  // pad1->SetGridy();
  pad1->SetLogy(1);
  pad1->cd();
  raw->Draw();
  bincount->Draw("same");
  TLegend * leg = (TLegend*) pad1->BuildLegend(0.6,0.65,0.90,0.90, Form("%s analysis, %i-%i%%", isTOF?"TOF":"TPC", icent*20, (icent+1)*20));
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  
  TPad * pad2 = (TPad*) cr->GetPad(2);
  pad2->SetFillColor(kWhite);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(0);
  pad2->SetMargin(0.1,0.05,0.15,0.005);
  //   pad2->SetGridx();
  // pad2->SetTicky(2);
  pad2->SetGridy();
  //   pad1->Draw();
  //   pad2->Draw();
  pad2->cd();
  rbin->Draw("E2");
  AddPaveText(rbin->GetTitle());
  
  return;
  
}
