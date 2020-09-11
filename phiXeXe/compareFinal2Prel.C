#include "/Users/fbellini/alice/macros/cosmetics/SetStyle.C"
#include "/Users/fbellini/alice/macros/cosmetics/MakeUp.C"

void CompareResolution(int cent = 0);
void CompareEfficiency(int cent = 0); 
TString Centrality(int cent);
Color_t ColorCentrality(int cent = 0);
TH1F * GetPlotRatio(TH1F * num, TH1F *den, TString errorType);
TH1D * GetPlotRatioD(TH1D * num, TH1D *den, TString errorType);

void compareFinal2Prel(int cent = 0){
  TString nOld("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/preliminaryQM18/Preliminary_Spectra_phi_XeXe544TeV.root");
  TString nNew("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0503ec/phiC3_tpc2sPtDep_tof3sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070/CORRECTED_br_fitResult.root");

  TFile * fOld = TFile::Open(nOld.Data(),"read");
  TFile * fNew = TFile::Open(nNew.Data(),"read");

  TH1D * hOld = (TH1D*) fOld->Get(Form("hCorrected_%i",cent));
  TH1D * hNew = (TH1D*) fNew->Get(Form("hCorrected_%i",cent));

  TH1D * hratio = (TH1D*) GetPlotRatioD(hNew, hOld, "num");
  Beautify(hratio, ColorCentrality(cent), 1, 2, 20, 1.3);

  TLegend * leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  myLegendSetUp(leg, 0.04);
  leg->AddEntry(hratio, Centrality(cent), "lp");

  hratio->GetYaxis()->SetMaxDigits(2);
  hratio->GetYaxis()->SetRangeUser(0.,2.);
  hratio->GetXaxis()->SetRangeUser(0.,10.);

  TLine * at1 = new TLine(0, 1., 10., 1.);
  at1->SetLineStyle(7);
  at1->SetLineWidth(1);

  TLine * at120 = new TLine(0., 1.2, 10., 1.2);
  at120->SetLineStyle(2);
  at120->SetLineWidth(1);
  
  TLine * at080 = new TLine(0., 0.8, 10., 0.8);
  at080->SetLineStyle(2);
  at080->SetLineWidth(1);

  SetStyle();
  TCanvas * cs = new TCanvas("cs", "cs", 800, 600);
  cs->cd();
  gPad->SetMargin(0.2, 0.07, 0.2, 0.07);
  hratio->Draw();
  leg->Draw();
  at1->Draw("same");
  at120->Draw("same");
  at080->Draw("same");
  cs->Print(Form("new2prel_cent%i.png",cent));

}

TString Centrality(int cent){
    if (cent>4) return "Overflow";
    int centb[] = {0, 10, 30, 60, 90};
    return Form("%i-%i %%", centb[cent],centb[cent+1]); 
}

Color_t ColorCentrality(int cent){
//    Color_t col[] = {kRed, kOrange, kGreen+1, kBlue, kBlack};
    Color_t col[] = {kOrange, kSpring+5, kTeal+5, kBlue+1, kMagenta+3};
    return col[cent];
}

TH1F * GetPlotRatio(TH1F * num, TH1F *den, TString errorType)
{  
  if (!num || !den) return NULL;
  TH1F * ratio = (TH1F*) num->Clone("latest");
  ratio->SetTitle(Form("%s/%s",num->GetTitle(), den->GetTitle()));
  ratio->SetMarkerStyle(20);
  if (errorType.Contains("num")) {
    ratio->Reset();
    for (Int_t j=1; j<num->GetNbinsX()+1;j++){
      if (den->GetBinContent(j)>0 || den->GetBinContent(j)<0) {
	    ratio->SetBinContent(j, num->GetBinContent(j)/den->GetBinContent(j));
	    ratio->SetBinError(j, num->GetBinError(j)/num->GetBinContent(j));
      } else {
	    ratio->SetBinContent(j,0);
      	ratio->SetBinError(j,0);
      }
    }
  } else {
    ratio->Divide(num, den, 1., 1., errorType.Data());
  }
    ratio->GetYaxis()->SetTitle("latest / preliminary");
    return ratio;
}

TH1D * GetPlotRatioD(TH1D * num, TH1D *den, TString errorType)
{  
  if (!num || !den) return NULL;
  TH1D * ratio = (TH1D*) num->Clone("latest");
  ratio->SetTitle(Form("%s/%s",num->GetTitle(), den->GetTitle()));
  ratio->SetMarkerStyle(20);
  if (errorType.Contains("num")) {
    ratio->Reset();
    for (Int_t j=1; j<num->GetNbinsX()+1;j++){
      if (den->GetBinContent(j)>0 || den->GetBinContent(j)<0) {
	    ratio->SetBinContent(j, num->GetBinContent(j)/den->GetBinContent(j));
	    ratio->SetBinError(j, num->GetBinError(j)/num->GetBinContent(j));
      } else {
	    ratio->SetBinContent(j,0);
      	ratio->SetBinError(j,0);
      }
    }
  } else {
    ratio->Divide(num, den, 1., 1., errorType.Data());
  }
    ratio->GetYaxis()->SetTitle("latest / preliminary");
    return ratio;
}

void CompareResolution(int cent){

    TString nOld("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/simulation/res_A3_tpc2sPtDep_tof2sveto5smism.root");
    TString nNew("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0423/simulation/res_A3_tpc2sPtDep_tof3sveto5smism.root");

    TFile * fOld = TFile::Open(nOld.Data(),"read");
    TFile * fNew = TFile::Open(nNew.Data(),"read");

    TH1F * hGaus_old = (TH1F *) fOld->Get(Form("hResVsPt%i_res1", cent));
    TH1F * hRMS_old = (TH1F *) fOld->Get(Form("hResVsPtRMS%i_r1", cent));
    TH1F * hMCV_old = (TH1F *) fOld->Get(Form("hResVsPtVMC%i_r1", cent));

    TH1F * hGaus_new = (TH1F *) fNew->Get(Form("hResVsPt%i_res1", cent));
    TH1F * hRMS_new = (TH1F *) fNew->Get(Form("hResVsPtRMS%i_r1", cent));
    TH1F * hMCV_new = (TH1F *) fNew->Get(Form("hResVsPtVMC%i_r1", cent));

    TH1F * hGaus_ratio = (TH1F *) GetPlotRatio(hGaus_new, hGaus_old, "sum2");
    TH1F * hMCV_ratio = (TH1F *) GetPlotRatio(hMCV_new, hMCV_old, "sum2");
    TH1F * hRMS_ratio = (TH1F *) GetPlotRatio(hRMS_new, hRMS_old, "sum2");

    Beautify(hGaus_ratio, ColorCentrality(cent)+2, 1, 2, 21, 1.3);
    Beautify(hMCV_ratio, ColorCentrality(cent)-2, 1, 2, 24, 1.3);
    Beautify(hRMS_ratio, ColorCentrality(cent), 1, 2, 20, 1.3);

    TLegend * leg = new TLegend(0.6, 0.7, 0.9, 0.9, Centrality(cent));
    myLegendSetUp(leg, 0.04);
    leg->AddEntry(hGaus_ratio, "Gaus", "lp");
    leg->AddEntry(hMCV_ratio, "MCV", "lp");
    leg->AddEntry(hRMS_ratio, "RMS", "lp");

    hGaus_ratio->GetYaxis()->SetMaxDigits(2);
    hGaus_ratio->GetYaxis()->SetRangeUser(0.,2.);

    TLine * at1 = new TLine(0, 1., 10., 1.);
    at1->SetLineStyle(7);
    at1->SetLineWidth(1);

    TLine * at120 = new TLine(0., 1.2, 10., 1.2);
    at120->SetLineStyle(2);
    at120->SetLineWidth(1);
    
    TLine * at080 = new TLine(0., 0.8, 10., 0.8);
    at080->SetLineStyle(2);
    at080->SetLineWidth(1);


    SetStyle();
    TCanvas * cres = new TCanvas("cres", "cres", 800, 600);
    cres->cd();
    gPad->SetMargin(0.2, 0.07, 0.2, 0.07);
    hGaus_ratio->Draw();
    hMCV_ratio->Draw("same");
    hRMS_ratio->Draw("same");
    leg->Draw();
    at1->Draw("same");
    at120->Draw("same");
    at080->Draw("same");
    cres->Print(Form("new2prel_res_cent%i.png",cent));

}

void CompareEfficiency(int cent){

    TString nOld("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/simulation/eff_A3_tpc2sPtDep_tof2sveto5smism.root");
    TString nNew("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0423/simulation/eff_A3_tpc2sPtDep_tof3sveto5smism.root");

    TFile * fOld = TFile::Open(nOld.Data(),"read");
    TFile * fNew = TFile::Open(nNew.Data(),"read");

    TH1F * heff_old = (TH1F *) fOld->Get(Form("hEffVsPt%i", cent));
    TH1F * heff_new = (TH1F *) fNew->Get(Form("hEffVsPt%i", cent));
    TH1F * hratio = (TH1F *) GetPlotRatio(heff_new, heff_old, "sum2");
    Beautify(hratio, ColorCentrality(cent), 1, 2, 20, 1.3);

    TLegend * leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    myLegendSetUp(leg, 0.04);
    leg->AddEntry(hratio, Centrality(cent), "lp");

    hratio->GetYaxis()->SetMaxDigits(2);
    hratio->GetYaxis()->SetRangeUser(0.,2.);

    TLine * at1 = new TLine(0, 1., 10., 1.);
    at1->SetLineStyle(7);
    at1->SetLineWidth(1);

    TLine * at120 = new TLine(0., 1.2, 10., 1.2);
    at120->SetLineStyle(2);
    at120->SetLineWidth(1);
    
    TLine * at080 = new TLine(0., 0.8, 10., 0.8);
    at080->SetLineStyle(2);
    at080->SetLineWidth(1);

    SetStyle();
    TCanvas * ceff = new TCanvas("ceff", "ceff", 800, 600);
    ceff->cd();
    gPad->SetMargin(0.2, 0.07, 0.2, 0.07);
    hratio->Draw();
    leg->Draw();
    at1->Draw("same");
    at120->Draw("same");
    at080->Draw("same");
    ceff->Print(Form("new2prel_eff_cent%i.png",cent));
}
