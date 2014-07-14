//Example for fitting signal/background. 
// This example can be executed with:
// root > .x FitSB.C  (using the CINT interpreter)
// root > .x FitSB.C+ (using the native complier via ACLIC)
   
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TVirtualFitter.h"
#include <assert.h>
#include "TFitResult.h"
//#include "myFitResult.C"
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/////////////////////////////////////////////////////////////
//
// AliRsnMassFitter for the fit of invariant mass distribution
// of resonances
//
// Author: fbellini@cern.ch
/////////////////////////////////////////////////////////////

#ifndef ALIRSNMASSFITTER_H
#define ALIRSNMASSFITTER_H

#include <TNamed.h>
#include <TString.h>

class TF1;
class TNtuple;
class TFile;
class TList;
class TH1F;
class TVirtualPad;

class AliRsnMassFitter : public TNamed {

 public:
  AliRsnMassFitter();
  AliRsnMassFitter(const TH1F* histoToFit, Double_t minvalue, Double_t maxvalue, Int_t rebin=1,Int_t fittypeb=0,Int_t fittypes=0);
  virtual ~AliRsnMassFitter();

  AliRsnMassFitter(const AliRsnMassFitter &mfit);
  AliRsnMassFitter& operator=(const AliRsnMassFitter &mfit);
  //BG functions
  Double_t poly1(Double_t *x, Double_t *par) { return par[0] + par[1]*x[0];} //Polynomial background functions
  Double_t poly2(Double_t *x, Double_t *par) { return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];} //Polynomial background functions
  Double_t poly3(Double_t *x, Double_t *par) { return par[0] + par[1]*x[0] + par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0];} //Polynomial background functions

  //Peak functions
  Double_t lorentzianPeak(Double_t *x, Double_t *par);  
  Double_t Breit( Double_t *x, Double_t * par);
  Double_t BreitB( Double_t *x, Double_t * par);

  // Sum of background and peak function
  Double_t BWpoly1(Double_t *x, Double_t *par) { return Breit(x, par) + poly1(x, &par[3]);}
  Double_t BWpoly2(Double_t *x, Double_t *par) { return Breit(x, par) + poly2(x, &par[3]);}
  Double_t BWpoly3(Double_t *x, Double_t *par) { return Breit(x, par) + poly3(x, &par[3]);}
  Double_t BWPSpoly2(Double_t *x, Double_t *par) { return BreitB(x, par) + poly2(x, &par[4]);}

  //setters
  void     SetHisto(const TH1F *histoToFit);
  void     SetRangeFit(Double_t fitXmin, Double_t fitXmax){fminMass=fitXmin; fmaxMass=fitXmax; CheckRangeFit();}
  void     SetMinRangeFit(Double_t fitXmin){fminMass=fitXmin;printf("CheckRangeFit after SetMaxRangeFit is also set\n");}
  void     SetMaxRangeFit(Double_t fitXmax){fmaxMass=fitXmax;printf("CheckRangeFit after SetMinRangeFit is also set\n");}
  void     SetBinN(Int_t newbinN){fNbin=newbinN;}
  void     SetType(Int_t fittypeb, Int_t fittypes);
  void     SetReflectionSigmaFactor(Int_t constant) {ffactor=constant;}
  void     SetInitialMean(Double_t mean) {fMass=mean;} // change the default value of the mean
  void     SetInitialGaussianSigma(Double_t sigma) {fSigmaSgn=sigma;} // change the default value of the sigma
  void     SetInitialWidth(Double_t gamma) {fWidthSgn=gamma;} // change the default value of the sigma
  void     SetSideBands(Bool_t onlysidebands=kTRUE) {fSideBands=onlysidebands;} // consider only side bands
  void     SetFixParam(Bool_t *fixpar){fFixPar=fixpar;}
  void     SetDefaultFixParam();
  Bool_t   SetFixThisParam(Int_t thispar,Bool_t fixpar);
  void     SetFixGaussianMean(Double_t mean=0.89595,Bool_t fixpar=kTRUE){SetInitialGaussianMean(mean); SetFixThisParam(fNFinalPars-2,fixpar);}
  void     SetFixGaussianSigma(Double_t sigma=0.005, Bool_t fixpar=kTRUE){SetInitialGaussianSigma(sigma); SetFixThisParam(fNFinalPars-1,fixpar);}
  void     SetFixWidth(Double_t gamma=0.0487, Bool_t fixpar=kTRUE){SetInitialWidth(gamma); SetFixThisParam(fNFinalPars-1,fixpar);}

  //getters
  TH1F*    GetHistoClone() const; //return the histogram
  void     GetRangeFit(Double_t &minvalue, Double_t &maxvalue) const {minvalue=fminMass; maxvalue=fmaxMass;}
  Double_t GetMinRangeFit()const {return fminMass;}
  Double_t GetMaxRangeFit()const {return fmaxMass;}
  Int_t    GetBinN()       const {return fNbin;}
  void     GetFitPars(Float_t* pars) const;
  Int_t    GetNFinalPars() const {return fNFinalPars;}
  void     GetTypeOfFit(Bool_t &background, Int_t &typeb) const {background = fWithBkg; typeb = ftypeOfFit4Bkg;}
  Int_t    GetReflectionSigmaFactor() const {return ffactor;} 
  Double_t GetMean() const {return fMass;}
  Double_t GetMeanUncertainty() const {return fMassErr;}
  Double_t GetSigma()const {return fSigmaSgn;}
  Double_t GetSigmaUncertainty()const { return fSigmaSgnErr;}
  Double_t GetWidth()const {return fWidthSgn;}
  Double_t GetWidthUncertainty()const { return fWidthSgnErr;}
  Double_t GetRawYield()const {return fRawYield;}
  Double_t GetRawYieldError()const {return fRawYieldErr;}
  Double_t GetChiSquare() const;
  Double_t GetReducedChiSquare() const;
  void     GetSideBandsBounds(Int_t& lb, Int_t& hb) const;
  Bool_t*  GetFixParam()const {return fFixPar;}
  Bool_t   GetFixThisParam(Int_t thispar)const;
  TVirtualPad* GetPad(Double_t nsigma=3,Int_t writeFitInfo=1)const;

  void     PrintParTitles() const;

  void     InitNtuParam(TString ntuname="ntupar"); // initialize TNtuple to store the parameters
  void     FillNtuParam(); //Fill the TNtuple with the current parameters
  TNtuple* GetNtuParam() const {return fntuParam;} // return the TNtuple
  TNtuple* NtuParamOneShot(TString ntuname="ntupar"); // the three functions above all together
  void     WriteHisto(TString path="./") const; // write the histogram
  void     WriteNtuple(TString path="./") const; // write the TNtuple
  void     WriteCanvas(TString userIDstring="",TString path="./",Double_t nsigma=3,Int_t writeFitInfo=1,Bool_t draw=kFALSE) const; //write the canvas in a root file
  void     DrawHere(TVirtualPad* pd,Double_t nsigma=3,Int_t writeFitInfo=1) const;
  void     DrawFit(Double_t nsigma=3) const;
  void     Reset();

  void     IntS(Float_t *valuewitherror) const;    // integral of signal given my the fit with error
  Double_t IntTot() const {return fhistoInvMass->Integral("width");}  // return total integral of the histogram
  void     Signal(Double_t nOfSigma,Double_t &signal,Double_t &errsignal) const; // signal in nsigma with error 
  void     Signal(Double_t min,Double_t max,Double_t &signal,Double_t &errsignal) const; // signal in (min, max) with error 
  void     Background(Double_t nOfSigma,Double_t &background,Double_t &errbackground) const; // backgournd in nsigma with error 
  void     Background(Double_t min,Double_t max,Double_t &background,Double_t &errbackground) const; // backgournd in (min, max) with error 
  void     Significance(Double_t nOfSigma,Double_t &significance,Double_t &errsignificance) const; // significance in nsigma with error 
  void     Significance(Double_t min,Double_t max,Double_t &significance,Double_t &errsignificance) const; // significance in (min, max) with error 

  Double_t FitFunction4MassDistr (Double_t* x, Double_t* par);
  Double_t FitFunction4Sgn (Double_t* x, Double_t* par);
  Double_t FitFunction4Bkg (Double_t* x, Double_t* par);
  Bool_t   MassFitter(Bool_t draw=kTRUE);
  Bool_t   RefitWithBkgOnly(Bool_t draw=kTRUE);
  void     RebinMass(Int_t bingroup=1);
  TF1*     GetBackgroundFullRangeFunc(){return fhistoInvMass->GetFunction("funcbkgFullRange");}
  TF1*     GetBackgroundRecalcFunc(){return fhistoInvMass->GetFunction("funcbkgRecalc");}
  TF1*     GetMassFunc(){return fhistoInvMass->GetFunction("funcmass");}

 private:

  void     PlotFit(TVirtualPad* pd,Double_t nsigma=3,Int_t writeFitInfo=1)const;

  void     ComputeParSize();
  void     ComputeNFinalPars();
  Bool_t   SideBandsBounds();
  Bool_t   CheckRangeFit();
  void     AddFunctionsToHisto();

  TH1F*     fhistoInvMass;     // histogram to fit
  Double_t  fminMass;          // lower mass limit
  Double_t  fmaxMass;          // upper mass limit
  Int_t     fminBinMass;       // bin corresponding to fminMass
  Int_t     fmaxBinMass;       // bin corresponding to fmaxMass
  Int_t     fNbin;             // number of bins
  Int_t     fParsSize;         // size of fFitPars array
  Int_t     fNFinalPars;       // number of parameters of the final function
  Float_t*  fFitPars;          //[fParsSize] array of fit parameters
  Bool_t    fWithBkg;          // signal+background (kTRUE) or signal only (kFALSE)
  Int_t     ftypeOfFit4Bkg;    // 0 = exponential; 1 = linear; 2 = pol2
  Int_t     ftypeOfFit4Sgn;    // 0 = gaus; 1 = gaus+gaus broadened
  Int_t     ffactor;           // number to multiply to the sigma of the signal to obtain the reflected gaussian
  TNtuple*  fntuParam;         // contains fit parameters
  Double_t  fMass;             // signal gaussian mean value
  Double_t  fMassErr;          // err signal gaussian mean value
  Double_t  fWidthSgn;         // signal breit wigner width
  Double_t  fWidthSgnErr;      // err signal breit wigner width
  Double_t  fSigmaSgn;         // signal gaussian sigma
  Double_t  fSigmaSgnErr;      // err signal gaussian sigma
  Double_t  fRawYield;         // signal gaussian integral
  Double_t  fRawYieldErr;      // err on signal gaussian integral
  Bool_t    fSideBands;        // kTRUE = only side bands considered
  Bool_t*   fFixPar;           //[fNFinalPars] for each par if kTRUE it is fixed in fit
  Int_t     fSideBandl;        // left side band limit (bin number)
  Int_t     fSideBandr;        // right side band limit (bin number)
  Int_t     fcounter;          // internal counter
  TList*    fContourGraph;     // TList of TGraph containing contour plots

  ClassDef(AliRsnMassFitter, 0); // class for invariant mass fit
};

#endif

inline Double_t AliRsnMassFitter::lorentzianPeak(Double_t *x, Double_t *par) 
{
  // Lorenzian Peak function
  return (0.5*par[0]*par[1]/TMath::Pi()) / 
    TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2]) 
   + .25*par[1]*par[1]);
}

inline Double_t AliRsnMassFitter::Breit( Double_t *x, Double_t * par)
{
  //Signal peak: Breit Wigner Peak function
  return par[0] * TMath::BreitWigner(x[0], par[1], par[2]);
}

inline Double_t AliRsnMassFitter::BreitB( Double_t *x, Double_t * par)
{
  //Relativistic Breit Wigner x Boltzmann factor
  //Parameters: constant = par[0], mass = par[1], width = par[2]
  //transverse momentum is par[3], to be fixed at runtime!
  const Double_t temp=0.160;
  const Double_t mpi2 = 0.1396*0.1396;
  const Double_t mka2 = 0.4937*0.4937;
  Double_t mkpi2 = mpi2+mka2; //GeV
  Double_t arg = 0.0, arg2 = 0.0, arg3=0.0, arg4=0.0, gamma = 0.0;
  Double_t boltz =0.0;
  arg = TMath::Sqrt(x[0] * x[0] + par[3] * par[3]);
  boltz = TMath::Exp(- arg / temp) / arg;
  arg2 = (x[0]*x[0] - par[1]*par[1]) * (x[0]*x[0] - par[1]*par[1]);
  arg3 = TMath::Power(x[0]*x[0] - mkpi2, 2.0) - 4.0*mpi2*mka2;
  arg4 = TMath::Power(par[1]*par[1] - mkpi2, 2.0) - 4.0*mpi2*mka2;
  gamma = par[2] * TMath::Power(par[1]/x[0], 4.0) * TMath::Power(arg3/arg4 , 1.5);
  return par[0] * x[0] * boltz * (x[0] * par[1] * gamma)/(arg2 + par[1]*par[1]*gamma*gamma);
}



/****************************************************/
/*                      FIT                         */
/****************************************************/
Int_t FitSB(TH1D * histo, 
	   Float_t *pt,
	   Int_t *cent,
	   TString fitfun = "BW+poly2",
	   Double_t fitXmin = 0.7, 
	   Double_t fitXmax = 1.1, 
	   Bool_t fixWidth = kTRUE, 
	   Bool_t fixMass = kFALSE,
	   Float_t binPt = 1.0,
	   TString gopt = "CPA",
	   TString opt = "MHIWRBC")
{
  //graphical options
  //TGaxis::SetMaxDigits(3);
  TString macroDir = "$HOME/alice/macro/kstar/fit";
  //gROOT->LoadMacro(Form("%s/myFitResult.C+",macroDir.Data()));
  TVirtualFitter::SetDefaultFitter("Minuit");
  

  Bool_t drawSignalFit=kFALSE, addPavePt=kFALSE, addPaveCent=kFALSE;
  if (gopt.Contains("B")) drawSignalFit = kTRUE;
  if (gopt.Contains("P")) addPavePt = kTRUE;
  if (gopt.Contains("C")) addPaveCent = kTRUE;
  
  TString ptLabel = Form("%3.2f #leq p_{T} < %3.2f GeV/#it{c}", pt[0],pt[1]);
  TString centLabel = Form("[%i-%i%%]", cent[0],cent[1]);
  if (gopt.Contains("A")) centLabel.Prepend("p-Pb "); 
  
  TCanvas *c1 = new TCanvas("c1","Fitting Demo",10,10,700,500);
  c1->SetFillColor(kWhite);
  c1->SetFrameFillColor(kWhite);
  c1->cd(1);
  //c1->SetGrid(0);
   
  if (!histo) {Printf("Error: invalid input histogram"); return -1;}
  //TH1F *histo = new TH1F("histo", "Breit Wigner Peak on polynomial Background",60,0,3);
  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(0.6);
  histo->SetStats(0);
  histo->GetYaxis()->SetTitle("dN/dM_{K#pi} ");
  histo->GetXaxis()->SetTitle("M_{K#pi} (GeV/#it{c}^{2})");
  histo->Scale(1.0 / histo->GetBinWidth(2));

  // create a TF1 
  TF1 *fitFcn; // = new TF1("fitFcn", fitFunction, fitXmin, fitXmax, 6);
  Int_t nBgPars = 0, nSigPars = 0;

  if (fitfun.Contains("BW+poly1")) {
    fitFcn = new TF1("fitFcn", BWpoly1, fitXmin, fitXmax, 5);
    nBgPars=2; nSigPars=3;
    fitFcn->SetParNames("constant","mass","width","p0","p1");
  }
  if (fitfun.Contains("BW+poly2")) {
    fitFcn = new TF1("fitFcn", BWpoly2, fitXmin, fitXmax, 6);
    nBgPars=3; nSigPars=3;
    fitFcn->SetParNames("constant","mass","width","p0","p1","p2");
  }
  if (fitfun.Contains("BW+poly3")) {
    fitFcn = new TF1("fitFcn", BWpoly3, fitXmin, fitXmax, 7);
    nBgPars=4; nSigPars=3;
    fitFcn->SetParNames("constant","mass","width","p0","p1","p2","p3");
  } 
  if (fitfun.Contains("BWPS+poly2")) {
    fitFcn = new TF1("fitFcn", BWPSpoly2, fitXmin, fitXmax, 7);
    nBgPars=3; nSigPars=4;
    fitFcn->SetParNames("constant","mass","width", "pT", "p0","p1","p2");
    fitFcn->FixParameter(3, binPt);
  }
  
  fitFcn->SetNpx(500);
  fitFcn->SetLineWidth(2);
  fitFcn->SetLineColor(kBlue+2);

  //-----------
  // mass 
  //-----------
   const Double_t pdgM = 0.89595;//GeV
   const Double_t pdgW = 0.0487;//GeV
   Int_t         nsigmaPeak = 5.0;  
   Double_t      peakLowLim = pdgM - nsigmaPeak * pdgW / 2.35;
   Double_t      peakUpLim  = pdgM + nsigmaPeak * pdgW / 2.35;
   
   fitFcn->SetParameter(1, pdgM);   // peak
   if (fixMass)
     fitFcn->FixParameter(1, pdgM);
   else 
     fitFcn->SetParLimits(1, peakLowLim, peakUpLim);
   
  //-----------
  // width 
  //-----------
  Double_t minW = 0.9*pdgW;
  Double_t maxW = 1.5*pdgW;
  fitFcn->SetParameter(2, pdgW); // width
  if (fixWidth)
    fitFcn->FixParameter(2, pdgW);
  else 
    fitFcn->SetParLimits(2, minW, maxW);
     
  //----------
  //FIT
  //----------
  TFitResultPtr result = histo->Fit("fitFcn","SEMBRQ","ep"); //add V=verbose or Q=quiet 
  Int_t fitstatus = result;
  // if (fitstatus>0) {
  //   // fitStatus =  migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult.
  //   Printf("*******************************\nFIT FAILED: minimization problem\n***************************");
  //   return fitstatus;
  // } else {
  //   if (fitstatus<0) {
  //     Printf("*******************************\nFIT FAILED: wrong function?\n*******************************");
  //     return fitstatus;
  //   }
  // }
  TMatrixDSym cov = result->GetCovarianceMatrix();
  result->Print("V");

  TVirtualFitter *f = TVirtualFitter::GetFitter();  
  Double_t bgCovElements[9];
  Double_t sigCovElements[16];
  Int_t id=0;//reset counter
  for (Int_t k=nSigPars;k<(nBgPars+nSigPars);k++) {
    for (Int_t h=nSigPars;h<(nBgPars+nSigPars);h++) {
      bgCovElements[id]=f->GetCovarianceMatrixElement(k,h);
      Printf("Bg Cov matrix element %i,%i = %e --> id #%i = %e",k,h, f->GetCovarianceMatrixElement(k,h), id,bgCovElements[id] ); 
      id++;
    }
  }
  id=0;//reset counter
  for (Int_t k=0;k<nSigPars;k++) {
    for (Int_t h=0;h<nSigPars;h++) {
      sigCovElements[id]=f->GetCovarianceMatrixElement(k,h);
      Printf("Sig Cov matrix element %i,%i = %e --> id #%i = %e",k,h, f->GetCovarianceMatrixElement(k,h), id, sigCovElements[id] ); 
      id++;
    }
  }
  // Double_t * covmarray = result->GetCovarianceMatrix()->GetMatrixArray();
  // for (Int_t k=0;k<(nBgPars+nSigPars)*(nBgPars+nSigPars);k++) {
  //   Printf("Cov matrix array element %i = %f",k, covmarray[k]); 
  // }
  
  Double_t par[7]; // 7 is the max number of parameters foreseen for the moment
  Double_t * parErr;//[7]; // 7 is the max number of parameters foreseen for the moment
  fitFcn->GetParameters(par); 
  parErr=fitFcn->GetParErrors(); 
 
  // writes the fit results into the par array
  Double_t fFcn[7][2]; //total fit function parameters & error
  for (Int_t j=0;j<7;j++) {
    fFcn[j][0]=result->Parameter(j);
    fFcn[j][1]=result->ParError(j);
  }

  //-----------
  //draw fitted functions:
  //-----------
   TF1 *backFcn;
   if (fitfun.Contains("poly1")) backFcn = new TF1("backFcn", poly1, fitXmin, fitXmax, nBgPars);//,"myFitFcn", "Bg");
   if (fitfun.Contains("poly2")) backFcn = new TF1("backFcn", poly2, fitXmin, fitXmax, nBgPars);//, "myFitFcn", "Bg");
   if (fitfun.Contains("poly3")) backFcn = new TF1("backFcn", poly3, fitXmin, fitXmax, nBgPars);//, "myFitFcn", "Bg");
   backFcn->SetLineColor(kRed);
   backFcn->SetLineStyle(2);
   backFcn->SetLineWidth(1);

   TF1 *signalFcn;
   if (fitfun.Contains("BWPS"))  signalFcn = new TF1("signalFcn", BreitB, fitXmin, fitXmax, nSigPars);
   else signalFcn = new TF1("signalFcn", Breit, fitXmin, fitXmax, nSigPars); 
   signalFcn->SetLineColor(kBlue);
   signalFcn->SetLineStyle(2);
   signalFcn->SetLineWidth(1);
   signalFcn->SetNpx(500);
 
   //copy parameters and errors from fit in auxiliary functions for bg and signal
   signalFcn->SetParameters(par);
   signalFcn->SetParErrors(parErr);
   if (drawSignalFit) signalFcn->Draw("same");
   if (fitfun.Contains("BWPS"))  {
     backFcn->SetParameters(&par[nSigPars]);
     backFcn->SetParErrors(&parErr[nSigPars]);
   } else {
     backFcn->SetParameters(&par[nSigPars]);
     backFcn->SetParErrors(&parErr[nSigPars]);
     // Printf("par bg #3: %f - %f", par[3], parErr[3] );
     // Printf("par set #0: %f - %f", backFcn->GetParameter(0), backFcn->GetParError(0) );
   }
   backFcn->Draw("same"); 
     
   //----------
   //Store the fit results with errors
   //----------
   Float_t min2G = pdgM - 2* pdgW;
   Float_t max2G = pdgM + 2* pdgW;
   Int_t ibinmin2G=histo->GetXaxis()->FindBin(min2G);
   Int_t ibinmax2G=histo->GetXaxis()->FindBin(max2G);
   Int_t ibinmin=histo->GetXaxis()->FindBin(fitXmin);
   Int_t ibinmax=histo->GetXaxis()->FindBin(fitXmax);

   Double_t fMass   [2]; //mass parameter & error
   fMass[0]=result->Parameter(1); fMass[1]=result->ParError(1);
  
   Double_t fGamma  [2]; //width parameter & error
   fGamma[0]=result->Parameter(2); fGamma[1]=result->ParError(2);
   //Double_t fSigma  [2]; //resolution parameter & error
  
  Double_t fChi2   [2]; //Chi2 and dof
  fChi2[0]=fitFcn->GetChisquare(); 
  fChi2[1]=fitFcn->GetNDF();

  Double_t fFcnInt [2] ; //integral of total fit function in fit range & error
  fFcnInt[0] = fitFcn->Integral(fitXmin, fitXmax);
  fFcnInt[1] = fitFcn->IntegralError(fitXmin, fitXmax);

  Double_t fBgInt  [2] = {0.0,0.0}; //integral of background fit function in fit range & error
  fBgInt[0] = backFcn->Integral(fitXmin, fitXmax);
  fBgInt[1] = backFcn->IntegralError(fitXmin, fitXmax, &par[nSigPars], bgCovElements);

  Double_t fRaw[2] = {0.0,0.0}; //integral of signal fit function in fit range & error
  fRaw[0] = fFcnInt[0] - fBgInt[0];
  fRaw[1] = fFcnInt[1] + fBgInt[1]; //simple sum because the two fits are not independent
  // fRaw[0] = signalFcn->Integral(fitXmin, fitXmax);
  // fRaw[1] = signalFcn->IntegralError(fitXmin, fitXmax, &par[0], sigCovElements);

  Double_t fHistInt[2]= {0.0,0.0}; //intqegral of histogram in fit range & error
  fHistInt[0]=histo->IntegralAndError(ibinmin, ibinmax, fHistInt[1], "width");

  Double_t fHistInt2Gamma[2]= {0.0,0.0}; //integral of histogram in +/-2Gamma range & error
  fHistInt2Gamma[0]=histo->IntegralAndError(ibinmin2G, ibinmax2G, fHistInt2Gamma[1], "width");

  Double_t fBgInt2Gamma[2] = {0.0,0.0}; //integral of background fit  function in +/-2Gamma range & error
  fBgInt2Gamma[0]=backFcn->Integral(min2G, max2G);
  fBgInt2Gamma[1]=backFcn->IntegralError(min2G, max2G, &par[nSigPars], bgCovElements);

  Double_t fRaw2Gamma[2] = {0.0,0.0}; //integral of signal fit  function in +/-2Gamma range & error
  fRaw2Gamma[0]=fHistInt2Gamma[0]-fBgInt2Gamma[0];//signalFcn->Integral(min2G, max2G);
  fRaw2Gamma[1]=TMath::Sqrt(fHistInt2Gamma[1]*fHistInt2Gamma[1]+fBgInt2Gamma[1]*fBgInt2Gamma[1]);//signalFcn->IntegralError(min2G, max2G, &par[0], sigCovElements);

  Double_t fFcnInt2Gamma[2]= {0.0,0.0}; //integral of total fit  function in +/-2Gamma range & error
  fFcnInt2Gamma[0]=fitFcn->Integral(min2G, max2G);
  fFcnInt2Gamma[1]=fitFcn->IntegralError(min2G, max2G);

  Double_t fFcnIntTail[2]= {0.0,0.0}; //integral of total fit function outside +/-2Gamma range & error
  fFcnIntTail[0]=fitFcn->Integral(0.0, min2G)+fitFcn->Integral(max2G, 2.0);
  fFcnIntTail[1]=fitFcn->IntegralError(0.0, min2G)+fitFcn->IntegralError(max2G, 2.0);

  Double_t fBgIntTail[2] = {0.0,0.0}; //integral of background fit function outside +/-2Gamma range & error
  fBgIntTail[0]=backFcn->Integral(0.0, min2G)+backFcn->Integral(max2G, 2.0);
  fBgIntTail[1]=backFcn->IntegralError(0.0, min2G,&par[nSigPars], bgCovElements)+backFcn->IntegralError(max2G, 2.0, &par[nSigPars], bgCovElements);

  Double_t fRawIntTail[2] = {0.0,0.0}; //integral of signal fit function outside +/-2Gamma range & error
  fRawIntTail[0] = fFcnIntTail[0]-fBgIntTail[0];
  fRawIntTail[1] = fFcnIntTail[1]+fBgIntTail[1];
  
  //----------
  //Print the fit results with errors
  //----------
  opt.ToUpper();
  if (opt.Contains("M")) Printf("Mass          = %7.5f +- %7.5f\n (rel. err. %6.3f%%)", fMass  [0], fMass  [1], fMass  [1]*100./fMass  [0] );
  if (opt.Contains("W")) Printf("Width         = %7.5f +- %7.5f\n (rel. err. %6.3f%%)", fGamma [0], fGamma [1], fGamma [1]*100./fGamma [0]);
  //if (opt.Contains("S")) printf("Sigma         = %7.5f +- %7.5f\n", fSigma [0], fSigma [1]);
  if (opt.Contains("H")) Printf("Hist integral = %10.1f +- %10.2f (rel. err. %6.3f%%)", fHistInt[0], fHistInt[1], fHistInt[1]*100./fHistInt[0] );
  if (opt.Contains("I")) Printf("Func integral = %10.1f +- %10.2f (rel. err. %6.3f%%)", fFcnInt [0], fFcnInt [1], fFcnInt[1]*100./fFcnInt[0] );
  if (opt.Contains("B")) Printf("BG   integral = %10.1f +- %10.2f (rel. err. %6.3f%%)", fBgInt  [0], fBgInt  [1], fBgInt  [1]*100./fBgInt  [0]);
  if (opt.Contains("R")) Printf("Raw  counts   = %10.1f +- %10.2f (rel. err. %6.3f%%)", fRaw    [0], fRaw    [1], fRaw    [1]*100./fRaw    [0]);   
  if (opt.Contains("T")) {
    Printf("Total func tails integral = %10.1f +- %10.2f (rel. err. %6.3f%%)", fFcnIntTail[0], fFcnIntTail[1],fFcnIntTail[1]*100./fFcnIntTail[0]);
    Printf(" Raw in tails = %10.1f +- %10.2f (rel. err. %6.3f%%)", fRawIntTail[0], fRawIntTail[1], fFcnIntTail[1]*100./fFcnIntTail[0]);
    Printf("Bg func tails integral = %10.1f +- %10.2f (rel. err. %6.3f%%)", fBgIntTail[0], fBgIntTail[1], fBgIntTail[1]*100./fBgIntTail[0]); 
  }
  if (opt.Contains("2G")){
    Printf("Histogram integral in +/-2Width = %10.1f +- %10.2f (rel. err. %6.3f%%)", fHistInt2Gamma[0], fHistInt2Gamma[1],  fHistInt2Gamma[1]*100./fHistInt2Gamma[0]);
    Printf("Total func integral in +/-2Width = %10.1f +- %10.2f (rel. err. %6.3f%%)", fFcnInt2Gamma[0], fFcnInt2Gamma[1], fFcnInt2Gamma[1]*100./fFcnInt2Gamma[0] );
    Printf("Signal func integral in +/-2Width = %10.1f +- %10.2f (rel. err. %6.3f%%)", fRaw2Gamma[0], fRaw2Gamma[1], fRaw2Gamma[1]*100./fRaw2Gamma[0] );
    Printf("Bg func integral in in +/-2Width = %10.1f +- %10.2f (rel. err. %6.3f%%)", fBgInt2Gamma[0], fBgInt2Gamma[1],fBgInt2Gamma[1]*100./fBgInt2Gamma[0]);
  }
  if (opt.Contains("C")) Printf("Chi2/NDF      = %4.2f/%4.2f \n", fChi2[0],fChi2[1]);
  if (opt.Contains("F")) {
    for (Int_t i = 0; i < 7; i++) {
      Printf("Fcn param %d   = %15.5f +- %15.5f (rel. err. %6.3f%%)", i, fFcn[i][0], fFcn[i][1], fFcn[i][1]*100./fFcn[i][0]);
    }
  }
  //----------
  //Print the fit results with errors in pave
  //----------
  Int_t ninfo=opt.Sizeof()-1;
  TPaveText *pave = new TPaveText(0.61, 0.61-0.04*ninfo, 0.89, 0.61, "NDC");
  if (opt.Contains("M")) pave->AddText(Form("Mass     = %7.5f #pm %7.5f", fMass  [0], fMass  [1]));
  if (opt.Contains("W")) pave->AddText(Form("Width    = %7.5f #pm %7.5f", fGamma [0], fGamma [1]));
  //   if (opt.Contains("S")) pave->AddText(Form("Sigma    = %7.5f #pm %7.5f", fSigma [0], fSigma [1]));   
  if (opt.Contains("H")) pave->AddText(Form("Hist int = %7.1f #pm %7.3f", fHistInt[0], fHistInt[1]));
  if (opt.Contains("I")) pave->AddText(Form("Func int = %7.1f #pm %7.3f", fFcnInt [0], fFcnInt [1]));
  if (opt.Contains("B")) pave->AddText(Form("BG   int = %7.1f #pm %7.3f", fBgInt  [0], fBgInt  [1]));
  if (opt.Contains("R")) pave->AddText(Form("Raw  cts = %7.1f #pm %7.3f", fRaw    [0], fRaw    [1]));
  if (opt.Contains("T")) {
    pave->AddText(Form("Tot. func tails int = %10.1f +- %10.2f", fFcnIntTail[0], fFcnIntTail[1]));
    pave->AddText(Form("S func tails int = %10.1f +- %10.2f", fRawIntTail[0], fRawIntTail[1]));
    pave->AddText(Form("Bg func tails int = %10.1f +- %10.2f", fBgIntTail[0], fBgIntTail[1]));
  }   
  if (opt.Contains("2G")) {
    pave->AddText(Form("Hist int (#pm2#Gamma) = %10.1f +- %10.2f",fHistInt2Gamma[0], fHistInt2Gamma[1]));
    pave->AddText(Form("Tot func int (#pm2#Gamma) = %10.1f +- %10.2f", fFcnInt2Gamma[0], fFcnInt2Gamma[1]));
    pave->AddText(Form("S func int (#pm2#Gamma) = %10.1f +- %10.2f", fRaw2Gamma[0], fRaw2Gamma[1]));
    pave->AddText(Form("Bg func integral (#pm2#Gamma) = %10.1f +- %10.2f", fBgInt2Gamma[0], fBgInt2Gamma[1]));  
  }
  if (opt.Contains("C")) pave->AddText(Form("Chi2/NDF = %4.2f/%4.2f", fChi2[0],fChi2[1]));
  if (opt.Contains("F")) {
    for (Int_t i = 0; i < 7; i++) {
      pave->AddText(Form("Fcn param %d   = %15.5f #pm %15.5f", i, fFcn[i][0], fFcn[i][1]));
    }
  }
  // pave->SetFillColor(kWhite);
  // pave->SetLineColor(kWhite);
  pave->SetBorderSize(0);
  pave->Draw();

  // draw the legend
  TPaveText * pave2 = new TPaveText(0.6,0.75,0.89,0.89,"NDC");
  if (addPaveCent) pave2->AddText(centLabel.Data());
  if (addPavePt) pave2->AddText(ptLabel.Data());
  //gStyle->SetOptTitle(0);
  pave2->SetFillColor(kWhite);
  pave2->SetLineColor(kWhite);
  pave2->SetBorderSize(0);
  pave2->Draw();
  TLegend *fitlegend=new TLegend(0.6,0.60,0.89,0.75);
  fitlegend->SetTextFont(42);
  fitlegend->SetFillColor(kWhite);
  fitlegend->SetLineColor(kWhite);
  fitlegend->SetTextSize(0.035);
  fitlegend->AddEntry(histo,"Data","lpe");
  fitlegend->AddEntry(backFcn,"Residual background","l");
  if (drawSignalFit) fitlegend->AddEntry(signalFcn,"Signal fit","l");
  fitlegend->AddEntry(fitFcn,"Global Fit","l");
  fitlegend->Draw();
   
  return 0;
}


