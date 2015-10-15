
#include <stdio.h>
#include "TCanvas.h"
//#include "RooPlot.h"

// retrieve reference values in the database PDG
Int_t PDG = 313;
TDatabasePDG *pdg     = TDatabasePDG::Instance();
TParticlePDG *part    = pdg->GetParticle(PDG);
Double_t      pdgMass    = part->Mass(); // const Double_t pdgMass = 0.89594;
Double_t      pdgWidth   = part->Width(); // 0.0487;// const Double_t pdgWidth = 0.0487;

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

void HistogramFit_pp7TeV (TString infile = "K0s-pp7TeV-Preliminary.root", 
			 TString suffix = "K0s",
			 //"finalWsyst_01apr14_0.root",
			 Double_t rangefitMin = 0.0001,
			 Double_t rangefitMax = 15.0,
			 Int_t ic = 100
			 ) 
{
  //Get Integrated YIELDS
  gStyle->SetOptStat(00001);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(1);
  gStyle->SetTextFont(42);

  //open file with histogram
  TString histname = Form("fHistPt%s",suffix.Data());
  TFile *file = TFile::Open(infile.Data());
  if (!file) {Printf("Invalid file."); return;}
  TH1D *data_stat = (TH1D*) file->Get("fHistPtK0ShortStatOnly")->Clone("data_stat");
  TH1D *data_syst = (TH1D*) file->Get("fHistPtK0ShortStatAndSystExceptNormalization")->Clone("data_syst");
 
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
  dummyResult->AddText(Form("m(K*) = %6.4f", myDummyResult->Parameter(0)));

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

// gROOT->LoadMacro("/Users/bellini/alisoft/aliroot/master/src/PWGLF/SPECTRA/UTILS/SpectraUtils.C");
/*****************************************************************/
/* LEVY-TSALLIS */
/*****************************************************************/
Double_t
LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2.);
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) / n / C;
  Double_t part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 *
LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{
  
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 15., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e6);
  return fLevyTsallis;
}


// gROOT->LoadMacro("/Users/bellini/alisoft/aliroot/master/src/PWGLF/SPECTRA/UTILS/YieldMean.C");
TH1 * YieldMean(TH1 *hstat, TH1 *hsys, TF1 *f = NULL, Double_t min = 0., Double_t max = 10., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q",TString logfilename="log.root")
{
  /* set many iterations when fitting the data so we don't
     stop minimization with MAX_CALLS */
  TVirtualFitter::SetMaxIterations(1000000);

  /* create output histo */
  Double_t integral, mean;
  TH1 *hout = new TH1D("hout", "", 8, 0, 8);
  TH1 *hlo, *hhi;
  
  /* create histo with stat+sys errors */
  TH1 *htot = (TH1 *)hstat->Clone(Form("%sfittedwith%s",hstat->GetName(),f->GetName()));
  for (Int_t ibin = 0; ibin < htot->GetNbinsX(); ibin++) {
    htot->SetBinError(ibin + 1, TMath::Sqrt(hsys->GetBinError(ibin + 1) * hsys->GetBinError(ibin + 1) + hstat->GetBinError(ibin + 1) * hstat->GetBinError(ibin + 1)));
  }
  
  //fb
  htot->GetYaxis()->SetRangeUser(1.e-6, 2.);
  /*
   *   measure the central value 
   */
  
  do Int_t fitres = htot->Fit(f, opt);
  while (fitres != 0);
  TFile* filewithfits=TFile::Open(logfilename.Data(),"UPDATE");
  htot->Write();
  filewithfits->Close();		
  delete filewithfits;	 
	
  cout<<" Fit sys+stat for " <<f->GetName()<<endl;		
  cout<<"NDF="<<f->GetNDF()<<" Chi^2="<<f->GetChisquare()<<" Chi^2/NDF="<<f->GetChisquare()/f->GetNDF()<<endl;

  hlo = YieldMean_LowExtrapolationHisto(htot, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(htot, f, max, hiprecision);
  YieldMean_IntegralMean(htot, hlo, hhi, integral, mean,kTRUE);
  hout->SetBinContent(kYield, integral);
  hout->SetBinContent(kMean, mean);

  /*
   * STATISTICS
   */
  
  TCanvas *cCanvasStat = new TCanvas("cCanvasStat");
  cCanvasStat->Divide(2, 1);
  
  /*
   * measure statistical error
   */

  /* fit with stat error */
  do Int_t fitres = hstat->Fit(f, opt);
  while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hstat, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hstat, f, max, hiprecision);
  
  /* random generation with integration (coarse) */
  TH1 *hIntegral_tmp = new TH1F("hIntegral_tmp", "", 1000, 0.75 * integral, 1.25 * integral);
  TH1 *hMean_tmp = new TH1F("hMean_tmp", "", 1000, 0.75 * mean, 1.25 * mean);
  for (Int_t irnd = 0; irnd < 100; irnd++) {
    /* get random histogram */
    TH1 *hrnd = YieldMean_ReturnRandom(hstat);
    /* fit */
    TH1 *hrndlo = YieldMean_ReturnCoherentRandom(hlo);
    TH1 *hrndhi = YieldMean_ReturnCoherentRandom(hhi);
    /* integrate */
    YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean);
    hIntegral_tmp->Fill(integral);
    hMean_tmp->Fill(mean);
    delete hrnd;
    delete hrndlo;
    delete hrndhi;
  }
  /* random generation with integration (fine) */
  TH1 *hIntegral = new TH1F("hIntegral", "", 100, 
                            hIntegral_tmp->GetMean() - 10. * hIntegral_tmp->GetRMS(),
                            hIntegral_tmp->GetMean() + 10. * hIntegral_tmp->GetRMS());
  TH1 *hMean = new TH1F("hMean", "", 100,
                        hMean_tmp->GetMean() - 10. * hMean_tmp->GetRMS(),
                        hMean_tmp->GetMean() + 10. * hMean_tmp->GetRMS());
  for (Int_t irnd = 0; irnd < 1000; irnd++) {
    /* get random histogram */
    TH1 *hrnd = YieldMean_ReturnRandom(hstat);
    /* fit */
    TH1 *hrndlo = YieldMean_ReturnCoherentRandom(hlo);
    TH1 *hrndhi = YieldMean_ReturnCoherentRandom(hhi);
    /* integrate */
    YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean);
    hIntegral->Fill(integral);
    hMean->Fill(mean);
    delete hrnd;
    delete hrndlo;
    delete hrndhi;
  }
  TF1 *gaus = (TF1 *)gROOT->GetFunction("gaus");
  
  cCanvasStat->cd(1);
  hIntegral->Fit(gaus, "q");
  integral = hout->GetBinContent(kYield) * gaus->GetParameter(2) / gaus->GetParameter(1);
  hout->SetBinContent(kYieldStat, integral);
  
  cCanvasStat->cd(2);
  hMean->Fit(gaus, "q");
  mean = hout->GetBinContent(kMean) * gaus->GetParameter(2) / gaus->GetParameter(1);
  hout->SetBinContent(kMeanStat, mean);
  
  /*
   * SYSTEMATICS
   */

  TCanvas *cCanvasSys = new TCanvas("cCanvasYieldSys");
  cCanvasSys->Divide(2, 1);
  cCanvasSys->cd(1)->DrawFrame(min, 1.e-6, max, 3.e+1);
  hsys->SetMarkerStyle(20);
  hsys->SetMarkerColor(1);
  hsys->SetMarkerSize(1);
  hsys->Draw("same");
  cCanvasSys->cd(2)->DrawFrame(min, 1.e-6, max, 3.e+1);
  hsys->Draw("same");
  
  /*
   * systematic error high
   */

  TH1 *hhigh = YieldMean_ReturnExtremeHighHisto(hsys);
  do Int_t fitres = hhigh->Fit(f, opt);
  while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hhigh, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hhigh, f, max, hiprecision);
  YieldMean_IntegralMean(hhigh, hlo, hhi, integral, mean);
  integral = TMath::Abs(integral - hout->GetBinContent(kYield));
  hout->SetBinContent(kYieldSysHi, integral);

  cCanvasSys->cd(1);
  f->SetLineColor(2);
  f->DrawCopy("same");
  
  /*
   * systematic error hard
   */

  TH1 *hhard = YieldMean_ReturnExtremeHardHisto(hsys);
  do Int_t fitres = hhard->Fit(f, opt);
  while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hhard, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hhard, f, max, hiprecision);
  YieldMean_IntegralMean(hhard, hlo, hhi, integral, mean);
  mean = TMath::Abs(mean - hout->GetBinContent(kMean));
  hout->SetBinContent(kMeanSysHi, mean);

  cCanvasSys->cd(2);
  f->SetLineColor(2);
  f->DrawCopy("same");
  
  /*
   * systematic error low
   */

  TH1 *hlow = YieldMean_ReturnExtremeLowHisto(hsys);
  do Int_t fitres = hlow->Fit(f, opt);
  while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hlow, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hlow, f, max, hiprecision);
  YieldMean_IntegralMean(hlow, hlo, hhi, integral, mean);
  integral = TMath::Abs(integral - hout->GetBinContent(kYield));
  hout->SetBinContent(kYieldSysLo, integral);

  cCanvasSys->cd(1);
  f->SetLineColor(4);
  f->DrawCopy("same");

  /*
   * systematic error soft
   */

  TH1 *hsoft = YieldMean_ReturnExtremeSoftHisto(hsys);
  do Int_t fitres = hsoft->Fit(f, opt);
  while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hsoft, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hsoft, f, max, hiprecision);
  YieldMean_IntegralMean(hsoft, hlo, hhi, integral, mean);
  mean = TMath::Abs(mean - hout->GetBinContent(kMean));
  hout->SetBinContent(kMeanSysLo, mean);

  cCanvasSys->cd(2);
  f->SetLineColor(4);
  f->DrawCopy("same");

  return hout;
}

TH1 *
YieldMean_LowExtrapolationHisto(TH1 *h, TF1 *f, Double_t min, Double_t binwidth = 0.01)
{
  /* find lowest edge in histo */
  Int_t binlo;
  Double_t lo;
  for (Int_t ibin = 1; ibin < h->GetNbinsX() + 1; ibin++) {
    if (h->GetBinContent(ibin) != 0.) {
      binlo = ibin;
      lo = h->GetBinLowEdge(ibin);
      break;
    }
  }
  
  Int_t nbins = (lo - min) / binwidth;
  TH1 *hlo = new TH1F("hlo", "", nbins, min, lo);
  
  /* integrate function in histogram bins */
  Double_t cont, err, width;
  for (Int_t ibin = 0; ibin < hlo->GetNbinsX(); ibin++) {
    width = hlo->GetBinWidth(ibin + 1);
    cont = f->Integral(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), (Double_t *)0, 1.e-6);
    err = f->IntegralError(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), (Double_t *)0, (Double_t *)0, 1.e-6);
    hlo->SetBinContent(ibin + 1, cont / width);
    hlo->SetBinError(ibin + 1, err / width);
  }

  return hlo;
}

TH1 * YieldMean_HighExtrapolationHisto(TH1 *h, TF1 *f, Double_t max, Double_t binwidth = 0.1)
{
  /* find highest edge in histo */
  Int_t binhi;
  Double_t hi;
  for (Int_t ibin = h->GetNbinsX(); ibin > 0; ibin--) {
    if (h->GetBinContent(ibin) != 0.) {
      binhi = ibin + 1;
      hi = h->GetBinLowEdge(ibin + 1);
      break;
    }
  }
  
  Int_t nbins = (max - hi) / binwidth;
  TH1 *hhi = new TH1F("hhi", "", nbins, hi, max);
  
  /* integrate function in histogram bins */
  Double_t cont, err, width;
  for (Int_t ibin = 0; ibin < hhi->GetNbinsX(); ibin++) {
    width = hhi->GetBinWidth(ibin + 1);
    cont = f->Integral(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), (Double_t *)0, 1.e-6);
    err = f->IntegralError(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), (Double_t *)0, (Double_t *)0, 1.e-6);
    hhi->SetBinContent(ibin + 1, cont / width);
    hhi->SetBinError(ibin + 1, err / width);
  }

  return hhi;
}

TH1 * YieldMean_ReturnRandom(TH1 *hin)
{
  TH1 *hout = (TH1 *)hin->Clone("hout");
  hout->Reset();
  Double_t cont, err;
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    cont = hin->GetBinContent(ibin + 1);
    err = hin->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, gRandom->Gaus(cont, err));
    hout->SetBinError(ibin + 1, err);
  }
  return hout;
}

TH1 * YieldMean_ReturnCoherentRandom(TH1 *hin)
{
  TH1 *hout = (TH1 *)hin->Clone("hout");
  hout->Reset();
  Double_t cont, err, cohe;
  cohe = gRandom->Gaus(0., 1.);
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    cont = hin->GetBinContent(ibin + 1);
    err = hin->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, cont + cohe * err);
    hout->SetBinError(ibin + 1, err);
  }
  return hout;
}

TH1 * YieldMean_ReturnExtremeHighHisto(TH1 *hin)
{
  TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehigh", hin->GetName()));
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    Double_t val = hin->GetBinContent(ibin + 1);
    Double_t err = hin->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, val + err);
  }
  return hout;
}

TH1 * YieldMean_ReturnExtremeLowHisto(TH1 *hin)
{
  TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremelow", hin->GetName()));
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    Double_t val = hin->GetBinContent(ibin + 1);
    Double_t err = hin->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, val - err);
  }
  return hout;
}

TH1 * YieldMean_ReturnExtremeSoftHisto(TH1 *hin)
{
  return YieldMean_ReturnExtremeHisto(hin, -1.);
}

TH1 * YieldMean_ReturnExtremeHardHisto(TH1 *hin)
{
  return YieldMean_ReturnExtremeHisto(hin, 1.);
}

TH1 * YieldMean_ReturnExtremeHisto(TH1 *hin, Float_t sign = 1.)
{
  Double_t ptlow, pthigh;
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    ptlow = hin->GetBinLowEdge(ibin + 1);
    break;
  }
  for (Int_t ibin = hin->GetNbinsX(); ibin >= 0; ibin--) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    pthigh = hin->GetBinLowEdge(ibin + 2);
    break;
  }

  Double_t mean = hin->GetMean();
  Double_t maxdiff = 0.;
  TH1 *hmax = NULL;
  for (Int_t inode = 0; inode < hin->GetNbinsX(); inode++) {

    Double_t ptnode = hin->GetBinCenter(inode + 1);
    TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehard", hin->GetName()));
    
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      Double_t val = hin->GetBinContent(ibin + 1);
      Double_t err = hin->GetBinError(ibin + 1);
      Double_t cen = hin->GetBinCenter(ibin + 1);
      if (cen < ptnode)
        err *= -1. + (cen - ptlow) / (ptnode - ptlow);
      else
        err *= (cen - ptnode) / (pthigh - ptnode);

      hout->SetBinContent(ibin + 1, val + sign * err);
    }

    Double_t diff = TMath::Abs(mean - hout->GetMean());
    if (diff > maxdiff) {
      //      printf("found max at %f\n", ptnode);
      if (hmax) delete hmax;
      hmax = (TH1 *)hout->Clone("hmax");
      maxdiff = diff;
    }
    delete hout;
  }
  return hmax;
}

YieldMean_IntegralMean(TH1 *hdata, TH1 *hlo, TH1 *hhi, Double_t &integral, Double_t &mean,Bool_t printinfo=kFALSE)
{
  
  /*
   * compute integrals
   */

  Double_t cont, err, width, cent;
  Double_t I = 0., IX = 0., Ierr = 0., IXerr = 0., Ilerr = 0., IXlerr = 0.;
  Double_t M = 0., Merr = 0., Mlerr = 0., C;
  Double_t dataonly=0.0;

  /* integrate the data */
  for (Int_t ibin = 0; ibin < hdata->GetNbinsX(); ibin++) {
    cent = hdata->GetBinCenter(ibin + 1);
    width = hdata->GetBinWidth(ibin + 1);
    cont = width * hdata->GetBinContent(ibin + 1);
    err = width * hdata->GetBinError(ibin + 1);
    if (err <= 0.) continue;
    I += cont;
    IX += cont * cent;
  }
  dataonly=I;	
  /* integrate low */
  for (Int_t ibin = 0; ibin < hlo->GetNbinsX(); ibin++) {
    cent = hlo->GetBinCenter(ibin + 1);
    width = hlo->GetBinWidth(ibin + 1);
    cont = width * hlo->GetBinContent(ibin + 1);
    err = width * hlo->GetBinError(ibin + 1);
    if (err <= 0.) continue;
    I += cont;
    IX += cont * cent;
  }
  /* integrate high */
  for (Int_t ibin = 0; ibin < hhi->GetNbinsX(); ibin++) {
    cent = hhi->GetBinCenter(ibin + 1);
    width = hhi->GetBinWidth(ibin + 1);
    cont = width * hhi->GetBinContent(ibin + 1);
    err = width * hhi->GetBinError(ibin + 1);
    if (err <= 0.) continue;
    I += cont;
    IX += cont * cent;
  }

  /* set values */
  integral = I;
  mean = IX / I;
  if(printinfo)	
  	cout<<"data only = "<<dataonly<<" total = "<<I<<" ratio= "<<dataonly/I<<endl; 	
}

Double_t computeStatUncertFromData(TH1D * hstat)
{
  if (!hstat) return -1.0;
  Double_t stat2 = 0.0;
  for (Int_t i=1;i<hstat->GetXaxis()->GetNbins()+1;i++){
    stat2 += TMath::Power(hstat->GetBinError(i)*hstat->GetBinWidth(i), 2.0);
  }
  return TMath::Sqrt(stat2);
}



/*

  
  //*********************************************************** FIT **********************************************************
  Double_t Levy(Double_t *pt, Double_t *par)
  {
    Double_t ldNdy  = par[0];
    Double_t lTemp  = par[1];
    Double_t lPower = par[2];
    Double_t lMass  = par[3];
	     
    Double_t lBigCoef = ((lPower - 1) * (lPower - 2)) / (lPower * lTemp * (lPower * lTemp + lMass * (lPower - 2)));
    Double_t lInPower = 1 + (TMath::Sqrt(pt[0] * pt[0] + lMass * lMass) - lMass) / (lPower * lTemp);
	     
    return ldNdy * pt[0] * lBigCoef * TMath::Power(lInPower, (-1) * lPower);
  }
  TF1 *fitLevy = new TF1("fitLevy",Levy,0.0,11,4);
  fitLevy->SetLineWidth(2);
  fitLevy->SetLineStyle(7);
  fitLevy->SetLineColor(kBlue+2);
  
  Double_t    mass      = pdgMass;
  Double_t    initYield = 0.0001;
  Double_t    initTemp  = 0.3;
  Double_t    initPower = 2.0;
	
  // first try without starting values for the parameters
  // This defaults to 1 for each param.
  // this results in an ok fit for the polynomial function
  // however the non-linear part (lorenzian) does not
  // respond well.
  fitLevy->SetParName(0,"Yield");
  //fitLevy->SetParLimits(0, 0.0, 10000.0);
  fitLevy->SetParameter(0, initYield);
    
  fitLevy->SetParName(1,"Temp");
  //fitLevy->SetParLimits(1, 0.0, 10000.0);
  fitLevy->SetParameter(1, initTemp);
	
  fitLevy->SetParName(2,"n");
  fitLevy->SetParLimits(2, 0.0, 10000.0);
  fitLevy->SetParameter(2, initPower);
	
  fitLevy->SetParName(3,"Mass");
  //fitLevy->SetParLimits(3, pdgMass*0.5, 1.5*pdgMass);
  fitLevy->FixParameter(3, mass);
	
  if (correctForBR) hist->Scale(1./0.666); //commented if using BR-corrected spectra
  hist->Fit("fitLevy","VM","VM",rangefitMin,rangefitMax);
	
  TF1 *dlevy = new TF1("dlevy",Levy,0.0,11,4);

  //Double_t Integral(Double_t a, Double_t b, const Double_t* params = 0, Double_t epsilon = 1e-12)
  Double_t IntegralHist = hist->Integral(ibinmin,ibinmax,"width"); //extremes are included in integral
  cout << "The integral of the histogram between bin " << ibinmin << " and bin " << ibinmax << " is: " << IntegralHist << endl;
  cout << "The integral of the histogram between x =" << hist->GetBinLowEdge(ibinmin) << " and x =" << hist->GetBinLowEdge(ibinmax) << " is: " << IntegralHist << endl;
  Double_t IntegralFit = fitLevy->Integral(rangefitMin,rangefitMax);
  Double_t IntegralFitExtrapolatedL = fitLevy->Integral(0.0,0.5);
  Double_t IntegralFitExtrapolatedR = fitLevy->Integral(8.0,10.0);
  cout << "The integral of the fit function between " << rangefitMin << " and " << rangefitMax << " is: " << IntegralFit << endl;
  cout << "The integral of the fit function between 0.0 and 0.5 is: " << IntegralFitExtrapolatedL << endl;
  cout << "The integral of the fit function between 8.0 and 10.0 is: " << IntegralFitExtrapolatedR << endl;

  Double_t dNdyDataFit = IntegralHist + IntegralFitExtrapolatedL +IntegralFitExtrapolatedR;

  for (Int_t i= ibinmin; i<ibinmax+1; i++){
    Float_t ratio = 0.0;
    if (hist->GetBinContent(i)==0.0) {
      Printf("oooops");
      continue;
    }
    ratio = (fitLevy->Eval(hist->GetBinCenter(i)) ) / hist->GetBinContent(i);
    fit2hist->SetBinContent(i, ratio);
  }
  
  Double_t levyPar[4], levyErr[3];
  fitLevy->GetParameters(levyPar);
  levyErr[0] = fitLevy->GetParError(0);
  levyErr[1] = fitLevy->GetParError(1);
  levyErr[2] = fitLevy->GetParError(2);

  TPaveText *result = new TPaveText(0.45, 0.65, 0.88, 0.88,"NDC"); 
  result->SetLineWidth(0);
  result->SetBorderSize(0);
  result->SetFillColor(kWhite);
  result->AddText(Form("Fit range: %4.2f-%4.2f GeV/c",rangefitMin,rangefitMax));
  result->AddText(Form("dN/dy = %10.4f #pm %10.4f", levyPar[0],levyErr[0]));
  result->AddText(Form("T = %10.4f #pm %10.4f", levyPar[1],levyErr[1]));
  result->AddText(Form("n = %10.4f #pm %10.4f", levyPar[2],levyErr[2]));
  result->AddText(Form("#chi^{2}/NDF = %10.4f / %i", fitLevy->GetChisquare(),fitLevy->GetNDF()));
  result->AddText(Form("<pt> = %10.6f ", dlevy->Moment(1,0.0,8.0, fitLevy->GetParameters())));

  Printf("Fit range |  dN/dy(fit) | dN/dy (data+fit) | T (GeV) |     n    |  #chi^{2}/NDF | <p_{T}> (GeV/c)");
  Printf("%4.2f-%4.2f  &  %10.4f \pm %10.4f  & %10.4f &  %10.4f \pm %10.4f  &  %10.4f \pm %10.4f  &  %10.4f/%i & %10.6f", 
	 rangefitMin, rangefitMax, 
	 levyPar[0],levyErr[0],
	 dNdyDataFit,
	 levyPar[1],levyErr[1],
	 levyPar[2],levyErr[2],
	 fitLevy->GetChisquare(),fitLevy->GetNDF(),
	 dlevy->Moment(1,0.0,8.0, fitLevy->GetParameters()));
  
  c1->cd();
  fitLevy->Draw("same");
  fitLevy->Print();
  result->Draw();
  
  c2->cd();
  fit2hist->SetMarkerStyle(20);  //  c1->cd(2);
  fit2hist->Draw("histp");
  l11->Draw("same");
  l09->Draw("same");       
  
  //save
  TString nomeOutput = Form("cent%i_levy%3.1f-%3.1f", ic,rangefitMin,rangefitMax);
  TString OutputFit = Form("%s_%s_inNOBR.root", (inputfile.Contains("TOF")? "TOF":"TPC") ,nomeOutput.Data());
  TFile *FitOut = TFile::Open(OutputFit.Data(),"RECREATE");
  FitOut->cd();
  hist->Write();
  fitLevy->Write();
  if (fit2hist) fit2hist->Write();
  c1->Write();

  c1->SaveAs(Form("FIT_%s_%s.png", (inputfile.Contains("TOF")? "TOF":"TPC") ,nomeOutput.Data()));
  c2->SaveAs(Form("RATIOFit2Hist_%s_%s.png", (inputfile.Contains("TOF")? "TOF":"TPC") ,nomeOutput.Data()));



  return;
}
*/

