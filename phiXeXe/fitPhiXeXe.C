#include "/Users/fbellini/alice/macros/MakeUp.C"
#include "/Users/fbellini/alice/macros/ResonAnT/fit/FitHistogram.C"

void SetStyle();
void DecodeFitStatus(Int_t fitstatus = 0);

Double_t poly1(Double_t *x, Double_t *par);
Double_t poly2(Double_t *x, Double_t *par);
Double_t poly3(Double_t *x, Double_t *par);
Double_t Voigt( Double_t *x, Double_t * par);
Double_t Breit( Double_t *x, Double_t * par);
Double_t BreitB( Double_t *x, Double_t * par);

TF1 * GetVOIGTpoly0(Double_t fitMin = 0.998, Double_t fitMax = 1.068);
TF1 * GetVOIGTpoly1(Double_t fitMin = 0.998, Double_t fitMax = 1.068); //new implementation, Int_t nsig = 5.E4, Int_t nbkg = 1.e6);
TF1 * GetVOIGTpoly2(Double_t fitMin = 0.998, Double_t fitMax = 1.068);
TF1 * GetBREITpoly1(Double_t fitMin = 0.998, Double_t fitMax = 1.068);
TF1 * GetBREITpoly2(Double_t fitMin = 0.998, Double_t fitMax = 1.068);
TF1 * GetRELBWpoly1(Double_t fitMin = 0.998, Double_t fitMax = 1.068);
TF1 * GetRELBWpoly2(Double_t fitMin = 0.998, Double_t fitMax = 1.068);
TH1D * GetHistoWithoutPeak(TH1D * h = NULL, Double_t peakMin = 1.010, Double_t peakMax = 1.030);
void SetLimitsFromBgOnlyFit(TF1 * fitBg = NULL, TF1 * fitFcn = NULL);
Float_t GetResolutionFromFilePt(TFile * fileRes = NULL, Int_t centbin = 0, Int_t ptbin = 0, Int_t rangeOpt = 0, Bool_t useRMS = 0, Bool_t returnErr = 0);

//TO DO: add customised fit ranges for each pt bin
//TO DO: finx new implementation of summed functions

Int_t fitPhiXeXe(TString inName = "sub_C3.root",
		 TString imgFormat = "png",
		 TString bgType = "Mixing", //alternative "Mixing"
		 Int_t selCentBin = 0, 
		 Int_t selPtBin = -1,
		 Double_t fitMin = 0.994,
		 Double_t fitMax = 1.050,
		 Double_t desiredIMbinWidth = 0.002,
		 TString fcnSignal = "VOIGT", 
		 TString fcnBg = "poly1",
		 Double_t minMass = 1.010,
		 Double_t maxMass = 1.030,
		 Double_t minWidth = 0.001,
		 Double_t maxWidth = 0.01,
		 Double_t minRes = -1.0,
		 Double_t maxRes = -1.0,
		 Bool_t useRmsRes = 0)
{

  SetStyle();
  Color_t color[]={kRed+1, kSpring+2, kBlue+1, kBlack, kAzure+10};
  Color_t marker[]={20, 21, 33, 34, 28}; 
  
  //-----------------------
  //open input file
  //-----------------------
  TFile * fin = TFile::Open(inName.Data(),"read");
  if (!inName || !fin || !fin->IsOpen()) {
    Printf("Invalid input file or impossible to open: %s", inName.Data());
    return 1;
  }

  //-----------------------
  //open resolution file
  //-----------------------
  TString resolFileName = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/sim/ana0221mc/res_C3_tpc2s_tof3sveto.root";
  TFile * finRes = TFile::Open(resolFileName.Data(),"read");
  if (!resolFileName || !finRes || !finRes->IsOpen()) {
    Printf("Invalid input file or impossible to open: %s", resolFileName.Data());
    return 1;
  }
  
  //-----------------------
  //retrieve axis and binning
  //-----------------------
  TAxis *ptbins = (TAxis*)fin->Get("ptbins");
  const Int_t nPtBins = ptbins->GetNbins();
  TAxis *centbins = (TAxis*)fin->Get("centbins");
  const Int_t nCentBins = centbins->GetNbins();

  //-------------------------------
  //create output file and objects
  //-------------------------------
  TString folderName = Form("fit_%s_%s%s", bgType.Data(), fcnSignal.Data(), fcnBg.Data());
  if (fcnSignal.Contains("VOIGT") && (minRes<=0.0) && (maxRes<=0.0)) folderName.Append("_FixRes");
  gSystem->Exec(Form("mkdir %s", folderName.Data()));
  gSystem->Exec(Form("mkdir %s/fit_r%4.3f-%4.3f", folderName.Data(), fitMin, fitMax));
  TString foutName = Form("%s/fit_r%4.3f-%4.3f/result_c%i.root", folderName.Data(), fitMin, fitMax, selCentBin);
  TFile *fout = TFile::Open(foutName.Data(), "RECREATE");
  TString centLabel = Form("%2.0f - %2.0f", centbins->GetBinLowEdge(selCentBin+1), centbins->GetBinUpEdge(selCentBin+1));
  
  // create the histograms which will contain all the raw counts and the corrections
  TH1D  *hmass  = new TH1D("mass" , "mass; #it{p}_{T} (GeV/#it{c}); #it{M} (GeV/#it{c}^{2})" , nPtBins, ptbins->GetXbins()->GetArray());
  TH1D  *hgamma = new TH1D("gamma", "width; #it{p}_{T} (GeV/#it{c}); #Gamma (GeV/#it{c}^{2})" , nPtBins, ptbins->GetXbins()->GetArray());
  TH1D  *hsigma = new TH1D("sigma", "resolution; #it{p}_{T} (GeV/#it{c}); #sigma_{Voigt} (GeV/#it{c}^{2})", nPtBins, ptbins->GetXbins()->GetArray());
  TH1D  *hchi2  = new TH1D("chi2oNDF" , "chi2/ndf; #it{p}_{T} (GeV/#it{c}); #chi^{2}/ndf" , nPtBins, ptbins->GetXbins()->GetArray());
  TH1D  *hrawIntegral = new TH1D("rawIntegral",Form("raw yields, %s; #it{p}_{T} (GeV/#it{c}); d#it{N}_{raw}/(d#it{y}d#it{p}_{T} (GeV/#it{c})^{-1})", centLabel.Data()), nPtBins, ptbins->GetXbins()->GetArray());
  TH1D  *hrawBC = new TH1D("rawBC",Form("raw yields (BC), %s; #it{p}_{T} (GeV/#it{c}); d#it{N}_{raw}/(d#it{y}d#it{p}_{T} (GeV/#it{c})^{-1})", centLabel.Data()), nPtBins, ptbins->GetXbins()->GetArray());
  TH1D  *hBgIntegral = new TH1D("bgIntegral", Form("Res. bg., %s; #it{p}_{T} (GeV/#it{c}); d#it{B}_{raw}/(d#it{y}d#it{p}_{T} (GeV/#it{c})^{-1})", centLabel.Data()) , nPtBins, ptbins->GetXbins()->GetArray());
  TH1D * hSoverB  = new TH1D("hSoverB",Form("S/B, %s", centLabel.Data()),  nPtBins, ptbins->GetXbins()->GetArray());
  TH1D * hSignif  = new TH1D("hSignif",Form("S/#sqrt{(S+B)} in #pm3#sigma, %s", centLabel.Data()),  nPtBins, ptbins->GetXbins()->GetArray());

  Beautify(hmass, color[selCentBin], 1, 2, marker[selCentBin], 1.0);
  Beautify(hgamma, color[selCentBin], 1, 2, marker[selCentBin], 1.0);
  Beautify(hsigma, color[selCentBin], 1, 2, marker[selCentBin], 1.0);
  Beautify(hchi2, color[selCentBin], 1, 2, marker[selCentBin], 1.0);
  Beautify(hrawIntegral, color[selCentBin], 1, 2, marker[selCentBin], 1.0);
  Beautify(hrawBC, color[selCentBin], 1, 2, marker[selCentBin], 1.0);
  Beautify(hBgIntegral, color[selCentBin], 1, 2, marker[selCentBin], 1.0);
  Beautify(hSoverB, color[selCentBin], 1, 2, marker[selCentBin], 1.0);
  Beautify(hSignif, color[selCentBin], 1, 2, marker[selCentBin], 1.0);


  //-----------------------
  //Define fitting function and parameters
  //-----------------------   
  Int_t nTotPars = 0, nBgPars = 0, nSigPars = 0;
  fcnSignal.ToUpper();
  fcnBg.ToUpper();

  TF1 *fitFcn = 0x0;
  if (fcnSignal.Contains("VOIGT")) {
    nSigPars = 4;
    if (fcnBg.Contains("POLY1")) {
      nBgPars = 2;
      fitFcn = (TF1*) GetVOIGTpoly1(fitMin, fitMax);
    } else if (fcnBg.Contains("POLY2")) {
      nBgPars = 3;
      fitFcn = (TF1*) GetVOIGTpoly2(fitMin, fitMax);
    }
  }
  
  else if (fcnSignal.Contains("RELBW")) {
    nSigPars = 4;
    if (fcnBg.Contains("POLY1")) {
      nBgPars = 2;
      fitFcn = (TF1*) GetRELBWpoly1(fitMin, fitMax);
    } else if (fcnBg.Contains("POLY2")) {
      fitFcn = (TF1*) GetRELBWpoly2(fitMin, fitMax);
      nBgPars = 3;
    }
  }
  
  else if (fcnSignal.Contains("BREIT")) {
    nSigPars = 3;
    if (fcnBg.Contains("POLY1")) {
      nBgPars = 2;
      fitFcn = (TF1*) GetBREITpoly1(fitMin, fitMax);	    
    } else if (fcnBg.Contains("POLY2")) {
      nBgPars = 3;
      fitFcn = (TF1*) GetBREITpoly2(fitMin, fitMax);	    
    }
  }
  
  nTotPars = nBgPars + nSigPars;
  Printf(":::: Fitting function Npar = %i, requested Npar' = %i", fitFcn->GetNpar(), nTotPars);
  fitFcn->SetNpx(500);
  fitFcn->SetLineWidth(3);
  fitFcn->SetLineColor(kRed);

  if (!fitFcn) {
    Printf(":::: ERROR: could not retrieve fitting function. ");
    return 3;
  }

  //-----------------------
  //retrieve values in the PDG
  //-----------------------
  TDatabasePDG *pdg     = TDatabasePDG::Instance();
  TParticlePDG *part    = pdg->GetParticle(333);
  Double_t      pdgMass    = part->Mass(); 
  Double_t      pdgWidth   = part->Width();
  Double_t      massRes = 0.0018;
  Double_t      massResErr = 0.00001;

  //-----------------------
  //define ranges for bin counting and significance 
  //-----------------------
  Double_t      nsigmaPeak = 5.0;
  Double_t      peakMin = pdgMass - nsigmaPeak * pdgWidth / 2.35;
  Double_t      peakMax = pdgMass + nsigmaPeak * pdgWidth / 2.35;
  Double_t      nsigma4Signif = 3.0;
  Double_t      peakMinSignif = pdgMass - nsigma4Signif * pdgWidth / 2.35;
  Double_t      peakMaxSignif = pdgMass + nsigma4Signif * pdgWidth / 2.35;

  //Fix parameters if requested
  if ((minMass<=0.0) || (maxMass<=0.0))
    fitFcn->FixParameter(1, pdgMass);
  else 
    fitFcn->SetParLimits(1, minMass, maxMass);

  if ((minWidth<=0.0) || (maxWidth<=0.0))
    fitFcn->FixParameter(2, pdgWidth);
  else 
    fitFcn->SetParLimits(2, minWidth, maxWidth);

  Double_t integrationTolerance = 1.0e-4;
  
  for (Int_t ibin = 2; ibin<nPtBins; ibin++){     
    if (selPtBin>0 && ibin!=selPtBin) continue;
 
    //-----------------------
    //retrieve desired input
    //-----------------------
    TString histName = Form("sub_norm_%s_ptBin%02i_centBin%02i", bgType.Data(), ibin, selCentBin);
    TH1D * histo = (TH1D*) fin->Get(histName.Data());
    if (!histo) {
      Printf("Cannot find desired histogram: %s >> Doing nothing.", histName.Data());
      return 2;
    }
  
    histo->GetYaxis()->SetTitle(Form("Counts / (%4.3f GeV/#it{c}^{2})", desiredIMbinWidth));
    histo->GetXaxis()->SetTitle("#it{M}_{KK} (GeV/#it{c}^{2})");
  
    //-----------------------
    //Rebin if requested
    //-----------------------
    // Double_t invMassBinWidth = ((TAxis*) histo->GetXaxis())->GetBinWidth(10);
    // Float_t rebinFactor = desiredIMbinWidth / invMassBinWidth;
    // if (rebinFactor>2.0 && rebinFactor<3.)
    //histo->Rebin(2);
    // else if (rebinFactor>3.0 && rebinFactor<4.) histo->Rebin(3);
    //histo->Scale(1.0 / histo->GetBinWidth(2));
  
    //-----------------------
    //get initial estimate of background from
    //Fit of background without signal 
    //-----------------------
    TF1* fitBg = 0x0;
    if (fcnBg.Contains("POLY1"))
      fitBg = new TF1(Form("%s_back",fcnBg.Data()), poly1, fitMin, fitMax, nBgPars);
    else if (fcnBg.Contains("POLY2"))
      fitBg = new TF1(Form("%s_back",fcnBg.Data()), poly2, fitMin, fitMax, nBgPars);

    if (!fitBg) {
      Printf(":::: Invalid background function. >> Doing nothing.");
      return 3;
    }
  
    TH1D * histoBgOnly = (TH1D*) GetHistoWithoutPeak(histo, peakMin, peakMax);
    histoBgOnly->Fit(fitBg,"RQN");
    fitBg->SetNpx(500);
    fitBg->SetLineWidth(3);
    fitBg->SetLineColor(kBlue+1);
    fitBg->SetLineStyle(7);
      
    //SetLimitsFromBgOnlyFit(fitBg, fitFcn);

    //------------------------------
    // Set voigtian resolution
    //------------------------------
    if (fcnSignal.Contains("VOIGT")) {
      if ((minRes<=0.0) && (maxRes<=0.0)) {
	massRes = GetResolutionFromFilePt(finRes, selCentBin, ibin, 0, useRmsRes, 0);
	massResErr = GetResolutionFromFilePt(finRes, selCentBin, ibin, 0, useRmsRes, 1);
	fitFcn->SetParLimits(3, massRes-massResErr, massRes+massResErr);
      } else {
	fitFcn->SetParLimits(3, minRes, maxRes);
      }
    }
    
    //------------------------------
    //FIT
    //------------------------------
    //fit - use Minuit2 if available
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    histo->Fit(fitFcn,"SEMBR","EP");
  
    TFitResultPtr result = histo->Fit(fitFcn,"SEMBR","EP"); //add V=verbose or Q=quiet
    DecodeFitStatus((Int_t) result);
    result->Print("V");
  
    Double_t par[2][7]; // 7 is the max number of parameters foreseen for the moment
    Double_t parsig[4]; // 4 is the max number of parameters foreseen for the signal for the moment
    Double_t parbg[3]; // 3 is the max number of parameters foreseen for the signal for the moment
    for (Int_t j = 0; j< nTotPars; j++) {
      par[0][j] = fitFcn->GetParameter(j); 
      par[1][j] = fitFcn->GetParError(j); 
      Printf(":::: Parameter %i -- %s = %6.4f +/- %6.4f", j, fitFcn->GetParName(j), par[0][j], par[1][j]);
    }
    
    //-----------
    //get signal and bg functions:
    //-----------
    TF1 *signalFcn;
    if (fcnSignal.Contains("RELBW")) signalFcn = new TF1(fcnSignal.Data(), BreitB, fitMin, fitMax, nSigPars);
    else if (fcnSignal.Contains("BREIT")) signalFcn = new TF1(fcnSignal.Data(), Breit, fitMin, fitMax, nSigPars); 
    else if (fcnSignal.Contains("VOIGT")) signalFcn = new TF1(fcnSignal.Data(), Voigt, fitMin, fitMax, nSigPars);
    signalFcn->SetLineColor(kRed+1);
    signalFcn->SetLineStyle(2);
    signalFcn->SetLineWidth(2);
    signalFcn->SetNpx(500);
 
    //copy parameters and errors from fit in auxiliary functions for bg and signal
    for (Int_t j=0; j<nTotPars; j++) {
      if (j<nSigPars){
	signalFcn->SetParameter(j, par[0][j]);
	signalFcn->SetParError(j, par[1][j]);
	parsig[j] = par[0][j];
      } else {
	fitBg->SetParameter(j-nSigPars, par[0][j]);
	fitBg->SetParError(j-nSigPars, par[1][j]);
	parbg[j-nSigPars] = par[0][j];
      }
    }     

     //-----------
    //get covariance matrix to calculate errors on integrals
    //-----------
    TMatrixDSym cov = result->GetCovarianceMatrix();
    result->Print("V");

    TVirtualFitter *f = TVirtualFitter::GetFitter();  
    Double_t bgCovElements[9]; //bg = 3x3 (poly2) matrix or 2x2 (poly1)
    Double_t sigCovElements[16]; // 3x3 matrix (breit) or 4x4 matrix (voigt)
    Int_t id = 0;//reset counter
    for (Int_t k = nSigPars; k<(nBgPars+nSigPars); k++) {
      for (Int_t h = nSigPars; h<(nBgPars+nSigPars); h++) {
    	bgCovElements[id] = f->GetCovarianceMatrixElement(k,h);
    	//Printf("Bg Cov matrix element %i,%i = %e --> id #%i = %e",k,h, f->GetCovarianceMatrixElement(k,h), id,bgCovElements[id] ); 
    	id++;
      }
    }
    id = 0;//reset counter
    for (Int_t k = 0; k<nSigPars; k++) {
      for (Int_t h = 0; h<nSigPars; h++) {
    	sigCovElements[id] = f->GetCovarianceMatrixElement(k,h);
    	//Printf("Sig Cov matrix element %i,%i = %e --> id #%i = %e",k,h, f->GetCovarianceMatrixElement(k,h), id, sigCovElements[id] ); 
    	id++;
      }
    }

    //------------------------------
    //Get integrals
    //------------------------------
    Double_t integralTotal[2] = {0.0, 0.0};
    Double_t integralSignal[2] = {0.0, 0.0};
    Double_t integralBg[2] = {0.0, 0.0};
    Double_t binCountTotal[2] = {0.0, 0.0};
    Double_t binCountSignal[2] = {0.0, 0.0};
    Double_t bgIntegral4BC[2] = {0.0, 0.0};
    
     integralTotal[0] = fitFcn->Integral(fitMin, fitMax, integrationTolerance)/desiredIMbinWidth;
    integralTotal[1] = fitFcn->IntegralError(fitMin, fitMax, result->GetParams(), cov.GetMatrixArray(), integrationTolerance)/desiredIMbinWidth;
    //func->IntegralError(x1,x2,r->GetParams(), cov->GetMatrixArray()->GetSub(4, 5, 4, 5, ""));  

    integralSignal[0] = signalFcn->Integral(fitMin, fitMax, integrationTolerance)/desiredIMbinWidth;
    integralSignal[1] = signalFcn->IntegralError(fitMin, fitMax, parsig, sigCovElements, integrationTolerance)/desiredIMbinWidth;

    integralBg[0] = fitBg->Integral(fitMin, fitMax, integrationTolerance)/desiredIMbinWidth;
    integralBg[1] = fitBg->IntegralError(fitMin, fitMax, parbg, bgCovElements, integrationTolerance)/desiredIMbinWidth;

    binCountTotal[0] = histo->IntegralAndError(histo->GetXaxis()->FindBin(peakMin), histo->GetXaxis()->FindBin(peakMax), binCountTotal[1]);
    
    bgIntegral4BC[0] = fitBg->Integral(peakMin, peakMax, integrationTolerance);
    bgIntegral4BC[1] = fitBg->IntegralError(peakMin, peakMax, parbg, bgCovElements, integrationTolerance);

    binCountSignal[0] = binCountTotal[0] - bgIntegral4BC[0];
    binCountSignal[1] = TMath::Sqrt(TMath::Power(binCountTotal[1], 2.0) + TMath::Power(bgIntegral4BC[1], 2.0));


    Double_t signal3sigma = signalFcn->Integral(peakMinSignif, peakMaxSignif, integrationTolerance);
    Double_t bg3sigma = fitBg->Integral(peakMinSignif, peakMaxSignif, integrationTolerance);
    Double_t significance3sigma = signal3sigma / TMath::Sqrt(signal3sigma + bg3sigma);
    Double_t SoverB = signal3sigma / bg3sigma;
    Double_t chi2NDF = fitFcn->GetChisquare() / fitFcn->GetNDF();

    //-----------------------
    // fill output histogram
    //-----------------------
    Float_t ptBinWidth = ptbins->GetBinWidth(ibin+1);
    hmass->SetBinContent(ibin+1, par[0][1]);  hmass->SetBinError(ibin+1, par[1][1]); 
    hgamma->SetBinContent(ibin+1,par[0][2]);  hgamma->SetBinError(ibin+1,par[1][2]); 
    hsigma->SetBinContent(ibin+1,par[0][3]);  hsigma->SetBinError(ibin+1,par[1][3]);
    hchi2->SetBinContent(ibin+1,chi2NDF);  hchi2->SetBinError(ibin+1,0.0); 
    hrawIntegral->SetBinContent(ibin+1,integralSignal[0]/ptBinWidth); hrawIntegral->SetBinError(ibin+1,integralSignal[1]/ptBinWidth); 
    hBgIntegral->SetBinContent(ibin+1,integralBg[0]/ptBinWidth);  hBgIntegral->SetBinError(ibin+1,integralBg[1]/ptBinWidth); 
    hrawBC->SetBinContent(ibin+1,binCountSignal[0]/ptBinWidth);  hrawBC->SetBinError(ibin+1,binCountSignal[1]/ptBinWidth); 

    hSoverB->SetBinContent(ibin+1,SoverB);  hSoverB->SetBinError(ibin+1,0.0);
    hSignif->SetBinContent(ibin+1,significance3sigma);  hSignif->SetBinError(ibin+1,0.0);
  
    //-----------------------
    // all numbers to pave
    //-----------------------
    TPaveText * pr = new TPaveText(0.2, 0.1, 0.9, 0.9, "NDC");
    BeautifyPave(pr, 0.045);
    pr->AddText(Form("Fit %s + %s", fcnSignal.Data(), fcnBg.Data()));
    pr->AddText(Form("Range: %4.3f - %4.3f GeV/#it{c}^{2}", fitMin, fitMax));
    pr->AddText(Form("#chi^{2}/NDF = %5.2f", chi2NDF));
    pr->AddText(Form("Raw (fit) = %.4e +/- %.4e", integralSignal[0], integralSignal[1]));
    pr->AddText(Form("Tot (BC) = %.4e +/- %.4e", binCountTotal[0], binCountTotal[1]));
    pr->AddText(Form("Raw (BC) = %.4e +/- %.4e", binCountSignal[0], binCountSignal[1]));
    pr->AddText(Form("Bg (fit) = %.4e +/- %.4e", integralBg[0], integralBg[1]));
  
    for (Int_t j = 0; j< nTotPars; j++) {
      pr->AddText(Form("%s = %6.4f +/- %6.4f", fitFcn->GetParName(j), par[0][j], par[1][j]));
    }
  
    //-----------------------
    //set cosmetics
    //-----------------------
    Int_t marker_Type[2] = {20, 24}; //LSB, MEB
    Beautify(histo, kBlack, 1, 1, 20, 1.0);
    histo->GetXaxis()->SetRangeUser(0.990, 1.1);
    histo->GetYaxis()->SetRangeUser( histo->GetMinimum()*1.0, histo->GetMaximum()*1.20);
 
    //-----------------------
    //show
    //-----------------------
    TCanvas *cFit  = new TCanvas("cFit", "FIT", 0, 0, 1000, 600);
    cFit->Divide(2,1);
    cFit->cd(1);
    histo->Draw();
    fitBg->Draw("same");
    //signalFcn->Draw("same");
    cFit->cd(2);
    pr->Draw();
    TString printFitName = Form("%s/fit_r%4.3f-%4.3f/fit_c%i_pt%i.%s", folderName.Data(), fitMin, fitMax, selCentBin, ibin, imgFormat.Data()); 
    cFit->Print(printFitName.Data());
  }

  fout->cd();
  hmass->Write(); 
  hgamma->Write();
  hsigma->Write();
  hchi2->Write();
  hrawIntegral->Write();
  hrawBC->Write();
  hBgIntegral->Write();
  hSoverB->Write();
  hSignif->Write();
  
  return 0;  
}


void SetStyle()
{
  //general - set style
  // initial setup
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTextFont(42);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(1);
  gStyle->SetLabelOffset(0.005,"yx");
  gStyle->SetLabelSize(0.06,"xyz");
  gStyle->SetTitleSize(0.07,"xyz");
  gStyle->SetTitleOffset(1.1,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetEndErrorSize(0); //sets in #of pixels the lenght of the tick at the end of the error bar
  gStyle->SetTitleAlign(33);
  gStyle->SetTitleX(.95);
  gStyle->SetTitleY(.90);
  TGaxis::SetMaxDigits(3);
  return;
}


//-----------------------------------------
//Signal peak: Voigtian peak function
//-----------------------------------------
Double_t Voigt( Double_t *x, Double_t * par)
{
  return par[0] * TMath::Voigt(x[0]-par[1], par[2], par[3]);
}


//-----------------------------------------
//Signal peak: Breit Wigner Peak function
//-----------------------------------------
Double_t Breit( Double_t *x, Double_t * par)
{
  return par[0] * TMath::BreitWigner(x[0], par[1], par[2]);
}


//-----------------------------------------
//Signal peak: Rel. Breit Wigner x Boltzmann 
//-----------------------------------------
Double_t BreitB( Double_t *x, Double_t * par)
{
  //Parameters: constant = par[0], mass = par[1], width = par[2]
  //transverse momentum is par[3], to be fixed at runtime!
  const Double_t temp = 0.156; //GeV
  const Double_t mpi2 = 0.1396*0.1396;
  const Double_t mka2 = 0.4937*0.4937;
  Double_t mkpi2 = mpi2+mka2; //GeV
  Double_t arg = 0.0, arg2 = 0.0, arg3 = 0.0, arg4 = 0.0, gamma = 0.0;
  Double_t boltz = 0.0;
  
  arg = TMath::Sqrt(x[0] * x[0] + par[3] * par[3]);
  boltz = TMath::Exp(- arg / temp) / arg;
  arg2 = (x[0]*x[0] - par[1]*par[1]) * (x[0]*x[0] - par[1]*par[1]);
  arg3 = TMath::Power(x[0]*x[0] - mkpi2, 2.0) - 4.0*mpi2*mka2;
  arg4 = TMath::Power(par[1]*par[1] - mkpi2, 2.0) - 4.0*mpi2*mka2;
  gamma = par[2] * TMath::Power(par[1]/x[0], 4.0) * TMath::Power(arg3/arg4 , 1.5);
  
  return par[0] * x[0] * boltz * (x[0] * par[1] * gamma)/(arg2 + par[1]*par[1]*gamma*gamma);
}


//-----------------------------------------
//Polynomial background functions
//-----------------------------------------
Double_t poly1(Double_t *x, Double_t *par)
{
  return par[0] + par[1]*x[0];
}

//-----------------------------------------
Double_t poly2(Double_t *x, Double_t *par)
{
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

//-----------------------------------------
Double_t poly3(Double_t *x, Double_t *par)
{
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0];
}


//-----------------------------------------
// Sum of background and peak function
//-----------------------------------------
Double_t BREITpoly1(Double_t *x, Double_t *par) {
  return Breit(x, par) + poly1(x, &par[3]);
}

//-----------------------------------------
Double_t BREITpoly2(Double_t *x, Double_t *par) {
  return Breit(x, par) + poly2(x, &par[3]);
}

//-----------------------------------------
Double_t RELBWpoly1(Double_t *x, Double_t *par) {
  return BreitB(x, par) + poly1(x, &par[4]);
}

//-----------------------------------------
Double_t RELBWpoly2(Double_t *x, Double_t *par) {
  return BreitB(x, par) + poly2(x, &par[4]);
}

//-----------------------------------------
Double_t VOIGTpoly0(Double_t *x, Double_t *par) {
  return Voigt(x, par) + par[4];
}

//-----------------------------------------
Double_t VOIGTpoly1(Double_t *x, Double_t *par) {
  return Voigt(x, par) + poly1(x, &par[4]);
}

//-----------------------------------------
Double_t VOIGTpoly2(Double_t *x, Double_t *par) {
  return Voigt(x, par) + poly2(x, &par[4]);
}

//-----------------------------------------
TF1 * GetVOIGTpoly0(Double_t fitMin, Double_t fitMax)
{
  TF1* fitFcn = new TF1("VOIGTpoly1", VOIGTpoly0, fitMin, fitMax, 4);
  fitFcn->SetParNames("Norm","Mass","Width", "Resolution"); 
  return fitFcn;
}

//-----------------------------------------
TF1 * GetVOIGTpoly1(Double_t fitMin, Double_t fitMax) //, Int_t nsig, Int_t nbkg)
{
   // Int_t NEvents = nsig+nbkg;
   // Int_t NBins   = 1e3;
   // double signal_mean = 3;
   // // const Double_t doubleKaonMass = 2.0*0.493677;
   // TF1 *f_signal = new TF1("Voigt", Voigt, 0., 5., 4);
   // f_signal->SetParName(0,"norm");
   // f_signal->SetParName(1,"Mass");
   // f_signal->SetParName(2,"Width");
   // f_signal->SetParName(3,"Resolution");
   
   // TF1 *f_bg     = new TF1("Poly1", poly1, 0., 5., 2);
   // f_bg->SetParName(0,"p0");
   // f_bg->SetParName(1,"p1");

   // // CONSTRUCTION OF THE TF1NORMSUM OBJECT ........................................
   // TF1NormSum *fnorm_voigt_poly1 = new TF1NormSum(f_signal,f_bg,nsig,nbkg);
   // TF1 * f_sum = new TF1("VOIGTpoly1", *fnorm_voigt_poly1, 0., 5., fnorm_voigt_poly1->GetNpar());
   // f_sum->Draw();
   
   // for (int i = 0; i < f_sum->GetNpar(); ++i){
   //   Printf(":::: Pre - TF1NormSum parameter %i: %s = %e", i, f_sum->GetParName(i), f_sum->GetParameter(i));
   // }

   // //f_sum->SetParameters(((Double_t *)fnorm_voigt_poly1->GetParameters()));
   // f_sum->SetParName(1,"NBackground");
   // f_sum->SetParName(0,"NSignal");
   // for (int i = 2; i < f_sum->GetNpar(); ++i){
   //   f_sum->SetParName(i, fnorm_voigt_poly1->GetParName(i) );
   //   Printf(":::: Pre - TF1NormSum parameter %i: %s = %e", i,  f_sum->GetParName(i), f_sum->GetParameter(i));  
   // }
   //return f_sum;

   //old implementation
   TF1* fitFcn = new TF1("VOIGTpoly1", VOIGTpoly1, fitMin, fitMax, 6);
   fitFcn->SetParNames("Norm","Mass","Width", "Resolution", "p0","p1"); 
   return fitFcn;
}

//-----------------------------------------
TF1 * GetVOIGTpoly2(Double_t fitMin, Double_t fitMax)
{
  TF1* fitFcn = new TF1("VOIGTpoly2", VOIGTpoly2, fitMin, fitMax, 7);
  fitFcn->SetParNames("Norm","Mass","Width", "Resolution", "p0","p1","p2"); 
  return fitFcn;
}

//-----------------------------------------
TF1 * GetBREITpoly1(Double_t fitMin, Double_t fitMax)
{
  TF1* fitFcn = new TF1("BREITpoly2",  BREITpoly1, fitMin, fitMax, 5);
  fitFcn->SetParNames("Norm","Mass","Width", "p0","p1"); 
  return fitFcn;
}

//-----------------------------------------
TF1 * GetBREITpoly2(Double_t fitMin, Double_t fitMax)
{
  TF1* fitFcn = new TF1("BREITpoly2", BREITpoly2, fitMin, fitMax, 6);
  fitFcn->SetParNames("Norm","Mass","Width", "p0","p1","p2"); 
  return fitFcn;
}

//-----------------------------------------
TF1 * GetRELBWpoly1(Double_t fitMin, Double_t fitMax)
{
  TF1* fitFcn = new TF1("BREITpoly2", RELBWpoly1, fitMin, fitMax, 6);
  fitFcn->SetParNames("Norm","Mass","Width", "Pt", "p0","p1"); 
  return fitFcn;
}

//-----------------------------------------
TF1 * GetRELBWpoly2(Double_t fitMin, Double_t fitMax)
{
  TF1* fitFcn = new TF1("BREITpoly2",  RELBWpoly2, fitMin, fitMax, 7);
  fitFcn->SetParNames("Norm","Mass","Width", "Pt", "p0","p1", "p2"); 
  return fitFcn;
}


//-----------------------------------------
TH1D * GetHistoWithoutPeak(TH1D * h, Double_t peakMin, Double_t peakMax)
{
  if (!h) return 0x0; 
  //create copy of histogram h with peak removed
  TH1D* a = (TH1D*) h->Clone(Form("%s_nopeak", h->GetName()));
  
  for (Int_t j = h->GetXaxis()->FindBin(1.000001*peakMin); j<=h->GetXaxis()->FindBin(0.999999*peakMax);j++){
    a->SetBinContent(j,0.);
    a->SetBinError(j,0.);
  }
  return a;
}

void SetLimitsFromBgOnlyFit(TF1 * fitBg, TF1 * fitFcn)
{
  if(!fitBg || !fitFcn) return;
  
  Double_t parBg[4][2];
  Double_t nBgPars = fitBg->GetNpar();
  Double_t nTotPars = fitFcn->GetNpar();
  Double_t nSigPars =  nTotPars - nBgPars;
  
  for (Int_t j = 0; j < nBgPars; j++) {
    parBg[j][0] = fitBg->GetParameter(j); 
    parBg[j][1] = fitBg->GetParError(j); 
    Printf("+++ BG Parameter %i = %6.4f +/- %6.4f", j, parBg[j][0], parBg[j][1]);

    //set limits in signal+bg fit function
    fitFcn->SetParLimits(nSigPars+j, parBg[j][0]*0.5, parBg[j][0]*1.5);
  }
  
  for (Int_t j = 0; j < nTotPars; j++) {
    Double_t limLow, limUp;
    fitFcn->GetParLimits(j, limLow, limUp);
    Printf(":::: Set limits for parameter %i -- %s = %6.4f +/- %6.4f", j, fitFcn->GetParName(j), limLow, limUp);
  } 
  return;
}


//----------------------------------------------------
void DecodeFitStatus(Int_t fitstatus)
{
  // fitStatus =  migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult.
  
  if (fitstatus!=0) {
    if (fitstatus/1000 > 0) Printf("--- Improve fit result returned %i", fitstatus/1000);
    else if (fitstatus/100  > 0) Printf("--- HESSE returned %i", fitstatus/100);
    else if (fitstatus/10   > 0) Printf("--- MINOS returned %i", fitstatus/10);
    else if (fitstatus%10   > 0) Printf("--- MIGRAD returned %i", fitstatus%10);
  }
  
  return;
}

//---------------------------------------------------
Float_t GetResolutionFromFilePt(TFile * resFile, Int_t centbin, Int_t ptbin, Int_t rangeOpt, Bool_t useRMS, Bool_t returnErr)
{
  Float_t returnVal = -1.0;
  if (!resFile) returnVal = -1.0;
  TString histName;
  //FIXME: hack until you have the centrality dependent resolution
  centbin = 0;
  if (useRMS) histName = Form("hResVsPtRMS%i_r%i", centbin, rangeOpt);
  else histName = Form("hResVsPt%i_res%i", centbin, rangeOpt);
  TH1D * hRes = (TH1D *) resFile->Get(histName.Data());
  if (returnErr) returnVal = hRes->GetBinError(ptbin);
    else returnVal = hRes->GetBinContent(ptbin);
  
  return returnVal;
}


// TVirtualFitter *f = TVirtualFitter::GetFitter();  
// Double_t bgCovElements[9];
// Double_t sigCovElements[16];
// Int_t id=0;//reset counter
// for (Int_t k=nSigPars;k<(nBgPars+nSigPars);k++) {
//   for (Int_t h=nSigPars;h<(nBgPars+nSigPars);h++) {
//     bgCovElements[id]=f->GetCovarianceMatrixElement(k,h);
//     Printf("Bg Cov matrix element %i,%i = %e --> id #%i = %e",k,h, f->GetCovarianceMatrixElement(k,h), id,bgCovElements[id] ); 
//     id++;
//   }
// }
// id=0;//reset counter
// for (Int_t k=0;k<nSigPars;k++) {
//   for (Int_t h=0;h<nSigPars;h++) {
//     sigCovElements[id]=f->GetCovarianceMatrixElement(k,h);
//     Printf("Sig Cov matrix element %i,%i = %e --> id #%i = %e",k,h, f->GetCovarianceMatrixElement(k,h), id, sigCovElements[id] ); 
//     id++;
//   }
// }
