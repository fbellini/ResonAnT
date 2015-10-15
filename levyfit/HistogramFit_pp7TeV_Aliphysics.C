
#include <stdio.h>
#include "TCanvas.h"
//#include "RooPlot.h"

// retrieve reference values in the database PDG
TDatabasePDG *pdg = TDatabasePDG::Instance();

enum EValue_t {
  kYield = 1,
  kYieldStat,
  kYieldSysHi,
  kYieldSysLo,
  kYieldSysData,
  kMean,
  kMeanStat,
  kMeanSysHi,
  kMeanSysLo
};

void HistogramFit_pp7TeV_Aliphysics(
			 TString infile = "Lambda-pp7TeV-Preliminary.root", //"K0s-pp7TeV-Preliminary.root",
			 TString suffix = "Lambda", //particle label
			 Double_t rangefitMin = 0.0001, //lower boundary of fit range
			 Double_t rangefitMax = 10.0,//upper boundary of fit range
			 Int_t ic = 100 //centrality or multiplicity
			 ) 
{
  gROOT->LoadMacro("$PHYS_SRC/PWGLF/SPECTRA/UTILS/SpectraUtils.C");
  gROOT->LoadMacro("$PHYS_SRC/PWGLF/SPECTRA/UTILS/YieldMean.C");
  // retrieve reference values in the database PDG
  Int_t PDG = 313;
  if (suffix.Contains("Lambda")) PDG = 3122; //K0s=310, K*=313, L=3122
  if (suffix.Contains("K0Short")) PDG = 310; //K0s=310, K*=313, L=3122
  TParticlePDG *part       = pdg->GetParticle(PDG);
  Double_t      pdgMass    = part->Mass(); // const Double_t pdgMass = 0.89594;
  Double_t      pdgWidth   = part->Width(); // 0.0487;// const Double_t pdgWidth = 0.0487;
  
  //Get Integrated YIELDS
  gStyle->SetOptStat(00001);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(1);
  gStyle->SetTextFont(42);

  //open file with histogram
  TFile *file = TFile::Open(infile.Data());
  if (!file) {Printf("Invalid file."); return;}
  TH1D *data_stat = ((TH1D*) file->Get(Form("fHistPt%sStatOnly", suffix.Data())))->Clone("data_stat");
  TH1D *data_syst = ((TH1D*) file->Get(Form("fHistPt%sStatAndSystExceptNormalization", suffix.Data())))->Clone("data_syst");
  // TH1D *data_stat = (TH1D*) file->Get("fHistPtK0ShortStatOnly")->Clone("data_stat");
  // TH1D *data_syst = (TH1D*) file->Get("fHistPtK0ShortStatAndSystExceptNormalization")->Clone("data_syst");
 
  if (!data_stat) {Printf("missing plot data_stat"); return;}
  if (!data_syst) {Printf("missing plot data_syst_uncorr"); return;}
  
  TF1 * myLevyTsallis = LevyTsallis("levy", pdgMass, 8., 0.37, 0.33);
  //fit example of stat only plot to shown T and n
  TH1D * myDummy = (TH1D*) data_stat->Clone("data_stat_copy"); 
  TFitResultPtr myDummyResult = myDummy->Fit(myLevyTsallis, "IRS","same", rangefitMin, rangefitMax);
  TPaveText *dummyResult = new TPaveText(0.5, 0.55, 0.85, 0.75,"NDC"); 
  dummyResult->SetLineWidth(1);
  dummyResult->SetBorderSize(1);
  dummyResult->SetFillColor(kWhite);
  dummyResult->SetTextAlign(12);
  //dummyResult->AddText(Form("Fit range: %4.2f-%4.2f GeV/c",rangefitMin,rangefitMax));
  dummyResult->AddText(Form("dN/dy = %6.4f #pm %6.4f", myDummyResult->Parameter(3), myDummyResult->ParError(3) ));
  dummyResult->AddText(Form("  T     = %6.4f #pm %6.4f", myDummyResult->Parameter(2), myDummyResult->ParError(2) ));
  dummyResult->AddText(Form("  n     = %6.4f #pm %6.4f", myDummyResult->Parameter(1), myDummyResult->ParError(1) ));
  dummyResult->AddText(Form("m = %6.4f", myDummyResult->Parameter(0)));

  //use spectra tools for mean pt and dN/dy
  TH1 * hresult = YieldMean(data_stat, data_syst, myLevyTsallis, rangefitMin, rangefitMax, 0.01, 0.1, "IRQ"); //IR0Q
  //  myYieldMean->Draw();
  TH1D * hSpecShiftLowSys  = YieldMean_ReturnExtremeLowHisto(data_syst);
  TH1D * hSpecShiftHighSys = YieldMean_ReturnExtremeHighHisto(data_syst);
  
  Double_t dNdy = hresult->GetBinContent(kYield);
  Double_t dNdy_stat = hresult->GetBinContent(kYieldStat);
  Double_t dNdy_sysHi = hresult->GetBinContent(kYieldSysHi);
  Double_t dNdy_sysLo = hresult->GetBinContent(kYieldSysLo);
  Double_t dNdy_sysData = hresult->GetBinContent(kYieldSysData);
  Double_t meanPt = hresult->GetBinContent(kMean);
  Double_t meanPt_stat = hresult->GetBinContent(kMeanStat);
  Double_t meanPt_sysHi = hresult->GetBinContent(kMeanSysHi);
  Double_t meanPt_sysLo = hresult->GetBinContent(kMeanSysLo);
  
  Double_t dNdy_maxSyst = TMath::Max(dNdy_sysHi,dNdy_sysLo);
  Double_t meanPt_maxSyst = TMath::Max(meanPt_sysHi,meanPt_sysLo);

  TPaveText *result = new TPaveText(0.35, 0.75, 0.88, 0.88,"NDC"); 
  result->SetLineWidth(0);
  result->SetBorderSize(0);
  result->SetFillColor(kWhite);
  result->AddText(Form("Fit range: %4.2f-%4.2f GeV/c",rangefitMin,rangefitMax));
  result->AddText(Form("dN/dy = %6.4f #pm %6.4f (stat) #pm %6.4f (sys)", dNdy, dNdy_stat, dNdy_maxSyst));
  // result->AddText(Form("T = %6.4f #pm %6.4f (stat) #pm %6.4f (sys)", ));
  // result->AddText(Form("n= %6.4f #pm %6.4f (stat) #pm %6.4f (sys)", ));
  result->AddText(Form("<p_{T}> = %6.4f #pm %6.4f (stat) #pm %6.4f (sys)",meanPt,meanPt_stat, meanPt_maxSyst));
 
  // Double_t manualdNdy= computedNdyFromData(data_stat);
  // Printf("dN/dy computed manually as integral of the data = %f", manualdNdy);
  Double_t manualStatUncert = computeStatUncertFromData(data_stat);
  Printf("stat uncertainty computed manually from data = %f --> %f%%", manualStatUncert, manualStatUncert*100./dNdy);

  Double_t piyield[3] = {16.661563, 0.002597, 0.782298};
  Double_t kyield[3] = {2.258426, 0.001456, 0.167224};
  Double_t proyield[3] = {0.934224, 0.000648, 0.062517};
  Double_t hyield[3] = {-1.,-1.,-1.};

  for (int i=0;i<3;i++){
    if (suffix.Contains("Ka")) {
      hyield[i] = kyield[i];
    } else {
      if (suffix.Contains("Pi")) {
	hyield[i] = piyield[i];
      } else 
	if (suffix.Contains("Pro")) {
	  hyield[i] = proyield[i];
	}
    }
  }
  Double_t Ks2K = 2.0*dNdy / hyield[0];
  Double_t Ks2K_stat = Ks2K * TMath::Sqrt( (2.0*dNdy_stat/dNdy)**2 + (hyield[1]/hyield[0])**2 );
  Double_t Ks_systot = TMath::Sqrt((dNdy_maxSyst/dNdy)**2 + 0.03**2 + 0.025**2);   //relative uncert. uncorr wrt K 
  Double_t K_systot = TMath::Sqrt((hyield[2]/hyield[0])**2 - 0.03**2); //relative uncert. uncorr wrt K*    
  Double_t Ks2K_sys = Ks2K * TMath::Sqrt(Ks_systot**2 + K_systot**2);
  
  Printf("ratio = %8.6f +/- %8.6f +/- %8.6f", Ks2K, Ks2K_stat, Ks2K_sys);
  
//ratio of data points to fitted function
  // TH1D *fit2hist = (TH1D*) file->Get(histogram.Data())->Clone("ratio");
  // fit2hist->SetTitle("fit/hist");
  // fit2hist->Reset("ICES");
  // fit2hist->SetLineColor(kBlack);
  // fit2hist->SetMarkerColor(kBlack);
  // fit2hist->SetMarkerStyle(1);
  // fit2hist->SetLineWidth(2);
  // fit2hist->GetYaxis()->SetRangeUser(0.7,1.3);
  
  // TLine * l11 = new TLine(rangefitMin, 1.1, rangefitMax, 1.1);
  // l11->SetLineWidth(1); l11->SetLineStyle(7); l11->SetLineColor(kRed);
  // TLine * l09 = new TLine(rangefitMin, 0.9, rangefitMax, 0.9); 
  // l09->SetLineWidth(1); l09->SetLineStyle(7); l09->SetLineColor(kRed);

  //draw
  TCanvas *c1 = new TCanvas(Form("c_%i",ic), "Fit",700,600);
  TCanvas *c2 = new TCanvas(Form("cr_%i",ic), "Ratio fit to histogram",700,600);
  data_stat->GetXaxis()->SetTitle(" p_{T} (GeV/c)");
  data_stat->GetYaxis()->SetTitle(" 1/N_{evt}*d^{2}N/dydp_{T} (GeV/c)^{-1}");   	
  c1->cd();
  gPad->SetLogy();
  data_stat->SetLineColor(kBlack);
  data_stat->SetMarkerColor(kBlack);
  data_stat->SetMarkerStyle(20);
  data_stat->SetMarkerSize(0.8);
  data_stat->GetXaxis()->SetRangeUser(0.0,15.0);
  data_stat->GetYaxis()->SetRangeUser(5e-6,2.);
  data_stat->Draw();
  result->Draw("same");
  c2->cd();
  myDummy->SetLineColor(kBlack);
  myDummy->SetMarkerColor(kBlack);
  myDummy->SetMarkerStyle(20);
  myDummy->SetMarkerSize(0.8);
  myDummy->GetXaxis()->SetRangeUser(0.0,15.0);
  myDummy->GetYaxis()->SetRangeUser(5e-6,2.);
  
  myDummy->Draw(); 
  myDummyResult->Draw("same");
  dummyResult->Draw("same");
  result->Draw("same");
  // hresult->Draw();
  return;
}

