#if !defined (__CINT__) || defined (__CLING__)
#include "/Users/fbellini/alice/macros/ResonAnT/phiXeXe/functions.h"
#include "/Users/fbellini/alice/macros/ResonAnT/phiXeXe/myYieldMean.C"
#include "/Users/fbellini/alice/macros/cosmetics/SetStyle.C"
#include "/Users/fbellini/alice/macros/cosmetics/Beautify.C"
#endif

#include <stdio.h>
#include "TCanvas.h"

// retrieve reference values in the database PDG
TDatabasePDG *pdg = TDatabasePDG::Instance();
void PrintParams(TF1 * fitFunc = 0, TPaveText * paveFitParams = 0, TFitResultPtr myDummyResult = 0);
void PrintParams(TF1 * fitFunc = 0, TPaveText * paveFitParams = 0, Double_t * par = 0, Double_t * parErr = 0);
Float_t StatUncertFromData(TH1D * hStat = 0);
Float_t SystUncertFromData(TH1D * hSyst = 0);

void FitSpectrum(Int_t centrality = -1, TString function = "bgbw", Double_t rangefitMin = 0.5, Double_t rangefitMax = 10.0, TString date = "30oct20");

TF1 * FitSpectrum(TH1D *data_stat = 0x0,
	     TH1D *data_syst = 0x0,
	     TString particle = "phi",
	     TString system = "Xe-Xe",
	     TString function = "bgbw",
	     Double_t rangefitMin = 0.5, //lower boundary of fit range
	     Double_t rangefitMax = 10.0,//upper boundary of fit range
	     Bool_t refitCentral =0,
	     Bool_t useOfficialMacro = 0,
	     TH1F * fitResults = 0);

void FitSpectrum(Int_t centrality = -1,
     TString function,
		 Double_t rangefitMin, //lower boundary of fit range
		 Double_t rangefitMax,//upper boundary of fit range
		 TString date //date of spectra
		 )
{
//Fix to errors in tolerance:  
//Error in <GSLError>: Error 18 in qags.c at 548 : cannot reach tolerance because of roundoff error
//solution by L. Moneta
//https://root-forum.cern.ch/t/tolerance-problem-in-integration-can-i-solve-it-with-functor/28675/4 
  ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6); 

  const int ncentbins = 5;
  Int_t centEdges[ncentbins+1] = {0, 10, 30, 50, 70, 90};
  // fin[0] = TFile::Open("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiC3_tpc2sPtDep_tof2sveto5smism/blastWaveFit/finalWsyst_smooth1_25apr18_0.root");
  // fin[1] = TFile::Open("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiC3_tpc2sPtDep_tof2sveto5smism/blastWaveFit/finalWsyst_smooth1_25apr18_1.root");
  // fin[2] = TFile::Open("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiC3_tpc2sPtDep_tof2sveto5smism/blastWaveFit/finalWsyst_smooth1_25apr18_2.root");
  TH1F * fitResults[ncentbins];  
  TFile * fin[ncentbins];
  for (Int_t i = 0; i<ncentbins; i++){
      fitResults[i] = new TH1F(Form("fitResults%i",i),"fitResults", 7, 0., 7.);
      fitResults[i]->GetXaxis()->SetBinLabel(1, "dN/dy");
      fitResults[i]->GetXaxis()->SetBinLabel(2,"stat");
      fitResults[i]->GetXaxis()->SetBinLabel(3,"syst");
      fitResults[i]->GetXaxis()->SetBinLabel(4, "extrap");
      fitResults[i]->GetXaxis()->SetBinLabel(5, "<p_{T}>");
      fitResults[i]->GetXaxis()->SetBinLabel(6,"mpt stat");
      fitResults[i]->GetXaxis()->SetBinLabel(7,"mpt syst");
      fin[i] = TFile::Open(Form("finalWsyst_smooth1_%s_%i.root",date.Data(), i));
      if (!fin[i]) return;
  }
  
  TFile * fileres = new TFile(Form("FITSPECTRUM_%s_%3.1f-%3.1f_%s.root", function.Data(),rangefitMin,rangefitMax, date.Data()), "recreate");

  TH1D * hstat[ncentbins];
  TH1D * hsys[ncentbins];
  
  for (int ic =0; ic<ncentbins; ic++){
    if (centrality>=0 && ic!=centrality) continue;
    hstat[ic] = (TH1D *) fin[ic]->Get(Form("hCorrected_%i", ic));
    if (!hstat[ic]) return;
    hsys[ic] = (TH1D *)  fin[ic]->Get(Form("hCorrected_%i%i_syst", ic, ic)); 
    if (ic==4) {
      rangefitMin = 0.7;
      rangefitMax = 7.0;
      Printf("::::: Centrality bin 70-90: fit range restricted to 0.7-7.0 GeV/c");
    }
    TF1 * fitFcn = (TF1*) FitSpectrum(hstat[ic], hsys[ic], "phi", Form("XeXe_%i%i", centEdges[ic], centEdges[ic+1]), function.Data(), rangefitMin, rangefitMax, 1, 0, fitResults[ic]);
    fitFcn->SetName(Form("%s%i", function.Data(), ic));
    fileres->cd();
    fitResults[ic]->Write();
    fitFcn->Write();
    hstat[ic]->Write(Form("hPhiXeXe_cent%i_stat",ic));
    hsys[ic]->Write(Form("hPhiXeXe_cent%i_syst",ic));
  }

  Printf("::::: RESULTS dN/dy :::::");
  for (int ic =0; ic<ncentbins; ic++){
    Printf("dN/dy cent %i: %8.6f %8.6f %8.6f (%3.2f)", ic,
	          fitResults[ic]->GetBinContent(1), fitResults[ic]->GetBinContent(2), fitResults[ic]->GetBinContent(3), fitResults[ic]->GetBinContent(4));
  }
  
  Printf("::::: RESULTS mean pT :::::");
  for (int ic =0; ic<ncentbins; ic++){
    Printf("<pT> cent %i: %8.6f %8.6f %8.6f", ic,
	          fitResults[ic]->GetBinContent(5), fitResults[ic]->GetBinContent(6), fitResults[ic]->GetBinContent(7) );
  }
  fileres->Close();
  return;
}

TF1 * FitSpectrum(TH1D *data_stat, TH1D *data_syst, TString particle, TString system, TString function, Double_t rangefitMin, Double_t rangefitMax, Bool_t refitCentral, Bool_t useOfficialMacro, TH1F * fitResults) 
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
  if (PDG<0) {Printf("Particle with PDG %i not supported. Nothing done.", PDG); return 0;}
 
  TParticlePDG *part       = pdg->GetParticle(PDG);
  Double_t      pdgMass    = part->Mass(); // const Double_t pdgMass = 0.89594;
  Double_t      pdgWidth   = part->Width(); // 0.0487;// const Double_t pdgWidth = 0.0487;
  //check input histogram
  if (!data_stat) {Printf("missing plot data_stat"); return 0;}
  if (!data_syst) {Printf("missing plot data_syst_uncorr"); return 0;}
  
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
  TF1 * myMTexp = MTExpdNdptTimesPtFunc("mT-exp", norm, temp, pdgMass);
  TF1 * myBoseEinstein = BoseEinstein("Bose-Einstein", pdgMass, temp, norm);
  TF1 * myFermiDirac = FermiDirac("Fermi-Dirac", pdgMass, temp, norm);
  TF1 * myBylinkin = Bylinkin("Bylinkin", pdgMass);
  
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

  TF1 * fitFunc = 0; Int_t nFcnPar = 3;
  if (function.Contains("levy")) fitFunc = myLevyTsallis;
  if (function.Contains("boltz")) fitFunc = myBoltzmann;
  if (function.Contains("mtexp")) fitFunc = myMTexp;
  if (function.Contains("fermi")) fitFunc = myFermiDirac;
  if (function.Contains("bose")) fitFunc = myBoseEinstein;
  if (function.Contains("bylinkin")) fitFunc = myBylinkin;
  if (function.Contains("bgbw")) { fitFunc = myBGBlastWave; nFcnPar = 4; }

  if (refitCentral) {
    do {
      myDummyResult = data_tot->Fit(fitFunc, "IRSq","same", rangefitMin, rangefitMax);
      Printf("My refit central trial: %d  --- chi2/ndf = %6.4f", trials++, fitFunc->GetChisquare()/fitFunc->GetNDF());
      if(trials > 10) {
  	Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
  	break;
      }
    }  while (myDummyResult->Status() != 0);    
    PrintParams(fitFunc, paveFitParams, myDummyResult);
  }
  
  hresult = (TH1D*) YieldMean(myDummy2, data_syst, fitFunc, extrapMin, extrapMax, 0.001, 0.1, "IRSq", Form("fitlog_%s_%s.root",particle.Data(), function.Data()), rangefitMin, rangefitMax, fitChi2ndf, fitParMain, fitParErrMain); 
  PrintParams(fitFunc, paveFitParamsUff, fitParMain, fitParErrMain);  
  
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
  
  TPaveText *paveYieldMean = new TPaveText(0.02, 0.83, 0.98, 0.97,"NDC"); 
  paveYieldMean->SetTextSize(0.035);
  paveYieldMean->SetTextFont(62);
  paveYieldMean->SetLineWidth(0);
  paveYieldMean->SetBorderSize(0);
  paveYieldMean->SetFillColor(kWhite);
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
  TString titleY = "1/#it{N}_{evt}*d^{2}#it{N}/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}";
  
  HistoMakeUp(data_tot, kBlack, 20, titleX.Data(), titleY.Data());
  HistoMakeUp(myDummy2, kBlack, 20, titleX.Data(), titleY.Data());
  //HistoMakeUp(data2fit, kBlack, 20, titleX.Data(), "data/fit");
  data_stat->GetXaxis()->SetTitle(titleX.Data());
  data_stat->GetYaxis()->SetTitle(titleY.Data());   	

  TString namefit = Form("fit_%s_%s_%s_%04.2f_%04.2f", particle.Data(), system.Data(), function.Data(), rangefitMin, rangefitMax);
  TCanvas *c1 = new TCanvas(namefit.Data(), "Fit",1000,600);
  SetStyle();
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLogy();
  data_tot->Draw("same");
  fitFunc->Draw("same");
  c1->cd(2);
  paveYieldMean->Draw("same");
  if (refitCentral) paveFitParams->Draw("same");
  if (!useOfficialMacro) paveFitParamsUff->Draw("same");
  c1->Print(Form("%s%s%s.eps",namefit.Data(), (useOfficialMacro? "_YM":""), (refitCentral?"":"woR")));
  Printf("::::: Statistical uncertainty from data = %e", StatUncertFromData(data_stat));
  Printf("::::: Systematic uncertainty from data = %e", SystUncertFromData(data_syst));
  Int_t ibinStart = data_stat->GetXaxis()->FindBin(0.5);
  Int_t ibinStop = data_stat->GetXaxis()->FindBin(10.0);
  Double_t histIntegralErr = 0;
  Double_t histIntegral = data_stat->IntegralAndError(ibinStart, ibinStop, histIntegralErr,"width");
  
  Printf("::::: Statistical uncertainty from integral of data = %e +/- %e", histIntegral, histIntegralErr);

  if (fitResults){
    fitResults->SetBinContent(1, dNdy);
    fitResults->SetBinContent(2, dNdy_stat);
    fitResults->SetBinContent(3, dNdy_maxSyst);
    fitResults->SetBinContent(4, (1.0 - histIntegral/dNdy));
    fitResults->SetBinContent(5, meanPt);
    fitResults->SetBinContent(6, meanPt_stat);
    fitResults->SetBinContent(7, meanPt_maxSyst);
  }
  
  TCanvas *c2 = new TCanvas("r", "r",700,600);
  c2->cd();
  fitResults->Draw();
  // data2fit->Draw("E2");
  // c2->SaveAs(Form("%s.png",namefit.Data()));
  return fitFunc;
}

void PrintParams(TF1 * fitFunc, TPaveText * paveFitParams, TFitResultPtr myDummyResult)
{
  //if (!fitFunc || !paveFitParams || !myDummyResult) return;
  for (Int_t i = 0; i<fitFunc->GetNpar(); i++){
    paveFitParams->AddText(Form(" %s = %8.4f #pm %8.4f", fitFunc->GetParName(i), myDummyResult->Parameter(i), myDummyResult->ParError(i)));
  }
  return;
}

void PrintParams(TF1 * fitFunc, TPaveText * paveFitParams, Double_t * par = 0, Double_t * parErr = 0)
{
  //if (!fitFunc || !paveFitParams || !myDummyResult) return;
  for (Int_t i = 0; i<fitFunc->GetNpar(); i++){
    paveFitParams->AddText(Form(" %s = %8.4f #pm %8.4f", fitFunc->GetParName(i), par[i], parErr[i]));
  }
  return;
}

Float_t StatUncertFromData(TH1D * hStat)
{
  if (!hStat) return 0.0;
  Float_t statUnc = 0.0;
  for (int i = 1; i<hStat->GetNbinsX()+1; i++) {
    if (hStat->GetBinContent(i)>0)
      statUnc += TMath::Power(hStat->GetBinError(i)*hStat->GetXaxis()->GetBinWidth(i),2.0);
  }
  return TMath::Sqrt(statUnc);
}

Float_t SystUncertFromData(TH1D * hSyst)
{
  if (!hSyst) return 0.0;
  Float_t systUnc = 0.0;
  for (int i = 1; i<hSyst->GetNbinsX()+1; i++) {
    if (hSyst->GetBinContent(i)>0)
      systUnc += hSyst->GetBinError(i);
  }
  return systUnc;
}
