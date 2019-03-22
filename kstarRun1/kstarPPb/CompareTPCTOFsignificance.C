//created on 29/08/2013 for pA analysis

void CompareTPCTOFsignificance(TString tofpath = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s/fit/",
			       TString tpcpath = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TPC_ANA/ana2s/fit/",
			       TString tofname = "RAW_best_fit_poly2.root",
			       TString tpcname = "RAW_best_fit_poly2.root")
{
  Color_t color[2][6]={ kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack, //tpc
			kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kGray+1}; //tof
  Int_t marker[2][6]={21, 22, 23, 34, 33, //tpc
		      25, 26, 32, 28, 27}; //tof
  
  //TPC_TOF matching eff
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  TString tofallcent_file = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/ALLCENT/significance_kstar_approval_tof2s_centBin00.root";
  TString tpcallcent_file = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/ALLCENT/significance_kstar_approval_tpc2s_centBin00.root";
  TString hname = "hSignificanceVsPt_";
  
  TFile * intof;
  TFile * intpc;
  TH1D * Stof[6];
  TH1D * Stpc[6];
  
  intof = TFile::Open(tofallcent_file.Data());
  if (!intof) return;
  Stof[5] = (TH1D*) intof->Get("hSignif");
  Stof[5]->SetLineColor(color[1][5]);
  Stof[5]->SetLineWidth(2);
  Stof[5]->SetTitle("TOF: S/#sqrt{S+B} (3#sigma), 0-100%");
  Stof[5]->GetYaxis()->SetRangeUser(0.,150.);
  Stof[5]->GetYaxis()->SetTitle("3#sigma significance");
  Stof[5]->GetYaxis()->SetTitleOffset(1.2);
  Stof[5]->GetXaxis()->SetRangeUser(0.5, 10.0);
  Stof[5]->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  intpc = TFile::Open(tpcallcent_file.Data());
  Stpc[5] = (TH1D*) intpc->Get("hSignif");
  Stpc[5]->SetLineColor(color[0][5]);
  Stpc[5]->SetLineWidth(2);
  Stpc[5]->SetTitle("TPC: S/#sqrt{S+B} (3#sigma), 0-100%");
  Stpc[5]->GetYaxis()->SetRangeUser(0.,150.);
  Stpc[5]->GetXaxis()->SetRangeUser(0.5, 10.0);
  Stpc[5]->GetYaxis()->SetTitle("3#sigma significance");
  Stpc[5]->GetYaxis()->SetTitleOffset(1.2);
  Stpc[5]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  
  TCanvas * ctof = new TCanvas("ctof","tof",600, 600);
  ctof->cd();
  Stof[5]->Draw();
  
  TCanvas * ctpc = new TCanvas("ctpc","tpc",600, 600);
  ctpc->cd();
  Stpc[5]->Draw();

  for (Int_t ic = 0; ic<5;ic++) {
    TString opentof = tofname;
    opentof.ReplaceAll("centBin00.root", Form("centBin0%i.root",ic));
    opentof.Prepend(tofpath.Data());
    TFile * intoftmp = TFile::Open(opentof.Data());
    if (intoftmp) {
      Stof[ic] = (TH1D*) intoftmp->Get(Form("%s%i",hname.Data(),ic));
      Stof[ic]->SetLineColor(color[1][ic]);
      Stof[ic]->SetLineWidth(2);
      Stof[ic]->SetTitle(Form("S/#sqrt{S+B} (3#sigma), %i-%i%%",20*ic,20*(1+ic)));
      Stof[ic]->GetYaxis()->SetRangeUser(0.,150.);
      Stof[ic]->GetXaxis()->SetRangeUser(0.5, 8.);
 
      ctof->cd();
      Stof[ic]->Draw("same");
    }
    TString opentpc = tpcname;
    opentpc.ReplaceAll("centBin00.root", Form("centBin0%i.root",ic));
    opentpc.Prepend(tpcpath.Data());
    TFile * intpctmp = TFile::Open(opentpc.Data());
    if (intpctmp) {
      Stpc[ic] = (TH1D*) intpctmp->Get(Form("%s%i",hname.Data(),ic));
      Stpc[ic]->SetLineColor(color[0][ic]);
      Stpc[ic]->SetLineWidth(2);
      Stpc[ic]->SetTitle(Form("S/#sqrt{S+B} (3#sigma), %i-%i%%", 20*ic,20*(1+ic)));
      Stpc[ic]->GetYaxis()->SetRangeUser(0.,150.);
      Stpc[ic]->GetXaxis()->SetRangeUser(0.5, 8.);
      ctpc->cd();
      Stpc[ic]->Draw("same");
    }
  }

  TLine * line100 = new TLine(.5, 100.,8.,100.);
  line100->SetLineStyle(7);
  line100->SetLineWidth(1);
  line100->SetLineColor(kBlack);

  TLine * line50 = new TLine(.5, 50.,8.,50.);
  line50->SetLineStyle(7);
  line50->SetLineWidth(1);
  line50->SetLineColor(kBlack);

  ctof->cd();
  tofleg = (TLegend*) gPad->BuildLegend(0.45,0.65,0.88,0.88, "TOF PID analysis");
  tofleg->SetBorderSize(0);
  tofleg->SetFillColor(kWhite);
  line100->Draw("same");
  line50->Draw("same");
  ctof->SaveAs("TOFsignificance.png");

  ctpc->cd();
  tpcleg = (TLegend*) gPad->BuildLegend(0.45,0.65,0.88,0.88, "TPC PID analysis");
  tpcleg->SetBorderSize(0);
  tpcleg->SetFillColor(kWhite);
  line100->Draw("same");
  line50->Draw("same");
  ctpc->SaveAs("TPCsignificance.png");


  return;
}
