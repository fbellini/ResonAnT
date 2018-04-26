#if !defined (__CINT__) || defined (__CLING__)
#include "/Users/fbellini/alice/macros/ResonAnT/phiXeXe/functions.h"
#include "/Users/fbellini/alice/macros/ResonAnT/phiXeXe/myYieldMean.C"
#include "/Users/fbellini/alice/macros/SetStyle.C"
#include "/Users/fbellini/alice/macros/Beautify.C"
#endif

#include <stdio.h>
#include "TCanvas.h"

// retrieve reference values in the database PDG
TDatabasePDG *pdg = TDatabasePDG::Instance();

void FitSpectrum(Int_t centrality = -1, Double_t rangefitMin = 0.05,  Double_t rangefitMax = 5.0);
void FitSpectrum(TH1D *data_stat = 0x0,
		 TH1D *data_syst = 0x0,
		 TString particle = "phi",
		 TString system = "Xe-Xe",
		 TString function = "bgbw",
		 Double_t rangefitMin = 0.05, //lower boundary of fit range
		 Double_t rangefitMax = 5.0,//upper boundary of fit range
		 Bool_t refitCentral = 1,
		 Bool_t useOfficialMacro = 0
		 );


void FitSpectrum(Int_t centrality,
		 Double_t rangefitMin, //lower boundary of fit range
		 Double_t rangefitMax//upper boundary of fit range
		 )
{
  TFile * fin[3];
  fin[0] = TFile::Open("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiC3_tpc2sPtDep_tof2sveto5smism/blastWaveFit/finalWsyst_smooth1_25apr18_0.root");
  fin[1] = TFile::Open("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiC3_tpc2sPtDep_tof2sveto5smism/blastWaveFit/finalWsyst_smooth1_25apr18_1.root");
  fin[2] = TFile::Open("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiC3_tpc2sPtDep_tof2sveto5smism/blastWaveFit/finalWsyst_smooth1_25apr18_2.root");
  if (!fin[0] || !fin[1] || !fin[2]) return;
  
  TH1D * hstat[3];
  TH1D * hsys[3];
  Int_t centEdges[4] = {0, 30, 60, 90};
  
  for (int ic =0; ic<3; ic++){
    if (centrality>=0 && ic!=centrality) continue;
    hstat[ic] = (TH1D *) fin[ic]->Get(Form("hCorrected_%i", ic));
    if (!hstat[ic]) return;
    hsys[ic] = (TH1D *)  fin[ic]->Get(Form("hCorrected_%i%i_syst", ic, ic)); 
    FitSpectrum(hstat[ic], hsys[ic], "phi", Form("XeXe_%i%i", centEdges[ic], centEdges[ic+1]), "bgbw", rangefitMin, rangefitMax, 1, 0);
  }
      //hack tmp for assigning sys uncert to 10%
    //hsys[ic] = (TH1D *) hstat[ic]->Clone(); //FIXME
    //for (int j = 1; j < hsys[ic]->GetNbinsX()+1; j++){
    //  hsys[ic]->SetBinError(j, 0.02*hsys[ic]->GetBinContent(j));  
  return;
}

void FitSpectrum(TH1D *data_stat, TH1D *data_syst, TString particle, TString system, TString function, Double_t rangefitMin, Double_t rangefitMax, Bool_t refitCentral, Bool_t useOfficialMacro) 
{

  //SetStyle
  SetStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetTextFont(42);

  // retrieve particle mass in the database PDG
  Int_t PDG = -1;
  if (particle.Contains("pion")) PDG = 211;
  if (particle.Contains("kaon")) PDG = 321;
  if (particle.Contains("proton")) PDG = 2212;
  if (particle.Contains("Lambda")) PDG = 3122; 
  if (particle.Contains("K0s")) PDG = 310; 
  if (particle.Contains("phi")) PDG = 333; 
  if (particle.Contains("Kstar")) PDG = 313;
  if (PDG<0) {Printf("Particle with PDG %i not supported. Nothing done.", PDG); return;}
 
  TParticlePDG *part       = pdg->GetParticle(PDG);
  Double_t      pdgMass    = part->Mass(); // const Double_t pdgMass = 0.89594;
  Double_t      pdgWidth   = part->Width(); // 0.0487;// const Double_t pdgWidth = 0.0487;
  //check input histogram
  if (!data_stat) {Printf("missing plot data_stat"); return;}
  if (!data_syst) {Printf("missing plot data_syst_uncorr"); return;}
  
  //prepare results objects for central value refitted
  //TPaveText *paveFitParams = new TPaveText(0.5, 0.45, 0.88, 0.75,"NDC"); 
  TPaveText *paveFitParams = new TPaveText(0.02, 0.45, 0.98, 0.82,"NDC"); 
  paveFitParams->SetTextSize(0.045);
  paveFitParams->SetLineWidth(1);
  paveFitParams->SetBorderSize(1);
  paveFitParams->SetFillColor(kWhite);
  paveFitParams->SetTextAlign(12);
  paveFitParams->AddText(Form("Fit %s (1st iteration)", function.Data()));
  paveFitParams->AddText(Form("Range: %4.2f-%4.2f GeV/c", rangefitMin,rangefitMax));
  
  //prepare results objects for output of YieldMean central value fit
  //  TPaveText *paveFitParamsUff = new TPaveText(0.5, 0.45, 0.88, 0.75,"NDC"); 
  TPaveText *paveFitParamsUff = new TPaveText(0.02, 0.02, 0.98, 0.42,"NDC"); 
  paveFitParamsUff->SetTextSize(0.045);
  paveFitParamsUff->SetLineWidth(1);
  paveFitParamsUff->SetBorderSize(1);
  paveFitParamsUff->SetFillColor(kSpring+1);
  paveFitParamsUff->SetTextAlign(12);
  paveFitParamsUff->AddText(Form("Fit %s (YieldsMean central fit)", function.Data()));
  paveFitParamsUff->AddText(Form("Range: %4.2f-%4.2f GeV/c", rangefitMin,rangefitMax));

  TH1D * myDummy2 = (TH1D*) data_stat->Clone("data_stat_copy2"); 

  /* create histo with stat+sys errors */
  TH1D *data_tot = (TH1D *)data_stat->Clone("data_tot");
  data_tot->SetTitle(particle.Data());
  for (Int_t ibin = 0; ibin < data_tot->GetNbinsX(); ibin++) {
    data_stat->SetBinError(ibin + 1, TMath::Sqrt(data_syst->GetBinError(ibin + 1) * data_syst->GetBinError(ibin + 1) + data_stat->GetBinError(ibin + 1) * data_stat->GetBinError(ibin + 1)));
  }

  //Define functions
  Double_t temp = 0.17; 
  Double_t norm = 1.0;
  Double_t npar = 5.0;
  Double_t beta_max = 0.9;
  TF1 * myLevyTsallis = LevyTsallis("levy", pdgMass, npar, temp, norm);
  TF1 * myBoltzmann = Boltzmann("boltz", pdgMass, temp, norm);
  TF1 * myBGBlastWave = BGBlastWave("bgbw", pdgMass, beta_max, temp, npar, norm);
  
  Double_t extrapMin = 0; Double_t extrapMax = 50.;
  
  Printf("***************** FIT %s ******************", particle.Data());
  Printf("Function: %s", function.Data());
  Printf("Fit Range: %3.1f-%3.1f GeV/c", rangefitMin,rangefitMax );
  Printf("Extrap Range: %3.1f-%3.1f GeV/c", extrapMin,extrapMax );
  Printf("*******************************************");

  //Perform the fit
  //use spectra tools for mean pt and dN/dy
  TH1D * hresult = 0x0;

  //Set fitter options
  TVirtualFitter::SetMaxIterations(6000); 
  TFitResultPtr myDummyResult = 0x0;
  Int_t trials = 0;
  Double_t fitChi2ndf[2] = {0.0, -1.e3}; //initialization values for chi2 and NDF
  Double_t fitParMain[5] = {0.0, 0.0, 0.0, 0.0, 0.0}; // values for parameters
  Double_t fitParErrMain[5] = {0.0, 0.0, 0.0, 0.0, 0.0}; //values for parameter errors
  
  if (function.Contains("levy")) {
    if (refitCentral) {
      do {
	myDummyResult = data_tot->Fit(myLevyTsallis, "IRSq","same", rangefitMin, rangefitMax);
	Printf("My refit central trial: %d  --- chi2/ndf = %6.4f", trials++, myLevyTsallis->GetChisquare()/myLevyTsallis->GetNDF());
	if(trials > 10) {
	  Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
	break;
	}
      }  while (myDummyResult->Status() != 0);    
      paveFitParams->AddText(Form(" norm = %8.4f #pm %8.4f", myDummyResult->Parameter(3), myDummyResult->ParError(3) ));
      paveFitParams->AddText(Form("  T   = %6.4f #pm %6.4f", myDummyResult->Parameter(2), myDummyResult->ParError(2) ));
      paveFitParams->AddText(Form("  n   = %6.4f #pm %6.4f", myDummyResult->Parameter(1), myDummyResult->ParError(1) ));
      paveFitParams->AddText(Form("  m   = %6.4f", myDummyResult->Parameter(0)));
      paveFitParams->AddText(Form("#chi^{2}/ndf = %6.4f", myLevyTsallis->GetChisquare()/myLevyTsallis->GetNDF())); 
    }
    //get yield and mean
    if (useOfficialMacro) {
      hresult =  (TH1D*) YieldMean(myDummy2, data_syst, myLevyTsallis, extrapMin, extrapMax, 0.01, 0.1, "IRSq", Form("fitlog_%s_%s.root",particle.Data(), function.Data()), rangefitMin, rangefitMax); 
    } else {
      hresult =  (TH1D*) YieldMean(myDummy2, data_syst, myLevyTsallis, extrapMin, extrapMax, 0.01, 0.1, "IRSq", Form("fitlog_%s_%s.root",particle.Data(), function.Data()), rangefitMin, rangefitMax, fitChi2ndf, fitParMain, fitParErrMain); 
      paveFitParamsUff->AddText(Form(" norm = %8.4f #pm %8.4f", fitParMain[3], fitParErrMain[3]));
      paveFitParamsUff->AddText(Form("  T   = %6.4f #pm %6.4f", fitParMain[2], fitParErrMain[2]));
      paveFitParamsUff->AddText(Form("  n   = %6.4f #pm %6.4f", fitParMain[1], fitParErrMain[1]));
      paveFitParamsUff->AddText(Form("  m   = %6.4f (fixed)", fitParMain[0]));
      paveFitParamsUff->AddText(Form("#chi^{2}/ndf = %6.4f", fitChi2ndf[0]/fitChi2ndf[1])); 
    } 
  }


  else if (function.Contains("boltz")) {
    if (refitCentral) {
      myDummyResult = data_tot->Fit(myBoltzmann, "IRSq","same", rangefitMin, rangefitMax);
      paveFitParams->AddText(Form(" norm  = %8.4f #pm %8.4f", myDummyResult->Parameter(2), myDummyResult->ParError(2) ));
      paveFitParams->AddText(Form("   T    = %6.4f #pm %6.4f", myDummyResult->Parameter(1), myDummyResult->ParError(1) ));
      paveFitParams->AddText(Form("   m    = %6.4f", myDummyResult->Parameter(0)));
      paveFitParams->AddText(Form("#chi^{2}/ndf = %6.4f", myBoltzmann->GetChisquare()/myBoltzmann->GetNDF())); 
   }
    //get yield and mean
    if (useOfficialMacro) {
      hresult =  (TH1D*) YieldMean(myDummy2, data_syst, myBoltzmann, extrapMin, extrapMax, 0.01, 0.1, "IRSq",Form("fitlog_%s_%s.root",particle.Data(), function.Data()), rangefitMin, rangefitMax); 
    } else {
      hresult =  (TH1D*) YieldMean(myDummy2, data_syst, myBoltzmann, extrapMin, extrapMax, 0.01, 0.1, "IRSq",Form("fitlog_%s_%s.root",particle.Data(), function.Data()), rangefitMin, rangefitMax,fitChi2ndf, fitParMain, fitParErrMain); 
      paveFitParamsUff->AddText(Form(" norm = %8.4f #pm %8.4f", fitParMain[2], fitParErrMain[2]));
      paveFitParamsUff->AddText(Form("  T   = %6.4f #pm %6.4f", fitParMain[1], fitParErrMain[1]));
      paveFitParamsUff->AddText(Form("  m   = %6.4f (fixed)", fitParMain[0]));
      paveFitParamsUff->AddText(Form("#chi^{2}/ndf = %6.4f", fitChi2ndf[0]/fitChi2ndf[1])); 
    } 
  }
  
  else if (function.Contains("bgbw")) {
    if (refitCentral) {
      do {
	myDummyResult = data_tot->Fit(myBGBlastWave, "IRSq","same", rangefitMin, rangefitMax);
	Printf("Trial: %d  --- chi2/ndf = %6.4f", trials++, myBGBlastWave->GetChisquare()/myBGBlastWave->GetNDF());
	if(trials > 10) {
	  Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
	  break;
	}
      }  while (myDummyResult->Status() != 0);
      paveFitParams->AddText(Form(" norm   = %8.4f #pm %8.4f", myDummyResult->Parameter(4), myDummyResult->ParError(4) ));
      paveFitParams->AddText(Form("  n     = %6.4f #pm %6.4f", myDummyResult->Parameter(3), myDummyResult->ParError(3) ));
      paveFitParams->AddText(Form("  T     = %6.4f #pm %6.4f", myDummyResult->Parameter(2), myDummyResult->ParError(2) ));
      paveFitParams->AddText(Form("  #beta_{max} = %6.4f #pm %6.4f", myDummyResult->Parameter(1), myDummyResult->ParError(1) ));
      paveFitParams->AddText(Form("  m     = %6.4f", myDummyResult->Parameter(0)));
      paveFitParams->AddText(Form("#chi^{2}/ndf = %6.4f", myBGBlastWave->GetChisquare()/myBGBlastWave->GetNDF())); 
    }

    //get yield and mean
    if (useOfficialMacro) {
      hresult = (TH1D*) YieldMean(myDummy2, data_syst, myBGBlastWave, extrapMin, extrapMax, 0.01, 0.1, "IRSq", Form("fitlog_%s_%s.root",particle.Data(), function.Data()), rangefitMin, rangefitMax); 
    } else {
      hresult = (TH1D*) YieldMean(myDummy2, data_syst, myBGBlastWave, extrapMin, extrapMax, 0.01, 0.1, "IRSq", Form("fitlog_%s_%s.root",particle.Data(), function.Data()), rangefitMin, rangefitMax, fitChi2ndf, fitParMain, fitParErrMain); 
      paveFitParamsUff->AddText(Form(" norm = %8.4f #pm %8.4f", fitParMain[4], fitParErrMain[4]));
      paveFitParamsUff->AddText(Form("  n     = %6.4f #pm %6.4f", fitParMain[3], fitParErrMain[3]));
      paveFitParamsUff->AddText(Form("  T   = %6.4f #pm %6.4f", fitParMain[2], fitParErrMain[2]));
      paveFitParamsUff->AddText(Form("  #beta_{max} = %6.4f #pm %6.4f", fitParMain[1], fitParErrMain[1]));
      paveFitParamsUff->AddText(Form("  m   = %6.4f (fixed)", fitParMain[0]));
      paveFitParamsUff->AddText(Form("#chi^{2}/ndf = %6.4f", fitChi2ndf[0]/fitChi2ndf[1]));   
    }
  }

  // TH1D * hSpecShiftLowSys  = YieldMean_ReturnExtremeLowHisto(data_syst);
  // TH1D * hSpecShiftHighSys = YieldMean_ReturnExtremeHighHisto(data_syst);
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
  
  //TPaveText *paveYieldMean = new TPaveText(0.35, 0.75, 0.88, 0.88,"NDC"); 
  TPaveText *paveYieldMean = new TPaveText(0.02, 0.83, 0.98, 0.97,"NDC"); 
  paveYieldMean->SetTextSize(0.035);
  paveYieldMean->SetTextFont(62);
  paveYieldMean->SetLineWidth(0);
  paveYieldMean->SetBorderSize(0);
  paveYieldMean->SetFillColor(kWhite);
  //paveYieldMean->AddText(Form("Fit: %s, range %4.2f-%4.2f GeV/c", function.Data(), rangefitMin,rangefitMax));
  paveYieldMean->AddText(Form("dN/dy (h+fit) = %8.4f #pm %8.4f (stat) #pm %8.4f (sys)", dNdy, dNdy_stat, dNdy_maxSyst));
  paveYieldMean->AddText(Form("<p_{T}> (h+fit) = %6.4f #pm %6.4f (stat) #pm %6.4f (sys)", meanPt, meanPt_stat, meanPt_maxSyst));
  
  //ratio of data points to fitted function
  // TH1D *data2fit = (TH1D*) data_tot->Clone("ratio");
  // data2fit->Divide(myLevyTsallis, 1.0);
  // data2fit->SetTitle("data/fit");
  // data2fit->SetLineColor(kBlack);
  // data2fit->SetMarkerColor(kBlack);
  // data2fit->SetMarkerStyle(1);
  // data2fit->SetLineWidth(2);
  // data2fit->GetYaxis()->SetRangeUser(0.7,1.3);
  
  //draw
  TString titleX = "#it{p}_{T} (GeV/#it{c})";
  TString titleY = "1/#it{N}_{evt}*d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]";
  
  HistoMakeUp(data_tot, kBlack, 20, titleX.Data(), titleY.Data());
  HistoMakeUp(myDummy2, kBlack, 20, titleX.Data(), titleY.Data());
  //HistoMakeUp(data2fit, kBlack, 20, titleX.Data(), "data/fit");
  data_stat->GetXaxis()->SetTitle(" p_{T} (GeV/c)");
  data_stat->GetYaxis()->SetTitle(" 1/N_{evt}*d^{2}N/dydp_{T} (GeV/c)^{-1}");   	

  TString namefit = Form("fit_%s_%s_%s_%04.2f_%04.2f", particle.Data(), system.Data(), function.Data(), rangefitMin, rangefitMax);
  TCanvas *c1 = new TCanvas(namefit.Data(), "Fit",1000,600);

  SetStyle();

  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLogy();
  data_tot->Draw("same");
  myBGBlastWave->Draw("same");
  c1->cd(2);
  paveYieldMean->Draw("same");
  if (refitCentral) paveFitParams->Draw("same");
  if (!useOfficialMacro) paveFitParamsUff->Draw("same");
  c1->Print(Form("%s%s%s.eps",namefit.Data(), (useOfficialMacro? "_YM":""), (refitCentral?"":"woR")));
  
  // TCanvas *c2 = new TCanvas(namefit.Prepend("ratio_"), "Ratio fit to histogram",700,600);
  // c2->cd();
  // data2fit->Draw("E2");
  // c2->SaveAs(Form("%s.png",namefit.Data()));

  return;
}

