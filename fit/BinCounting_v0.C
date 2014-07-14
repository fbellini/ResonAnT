#include "Riostream.h"
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TArray.h>
#include <TMinuit.h>
#include <TFitResult.h>
#include "CPolyFit.C"
#include "myFitFcn.C"
#include "myFitResult.C"

myFitResult* BinCounting
(
   TH1D       *hist,
   myFitFcn   *fit,
   Double_t    viewMin  = 0.70,
   Double_t    viewMax  = 1.10,
   Double_t    excludeMin  = 0.896,
   Double_t    excludeMax  = 0.896,
   Double_t    BCpeakRangeMin  = 0.70,
   Double_t    BCpeakRangeMax  = 1.20,
   Double_t    peakMin  = 0.880,
   Double_t    peakMax  = 0.940,
   Double_t    pt = 1.0
 ) 
{

  //check arguments consistency
  if ( (viewMin>=excludeMin) || (viewMax<=excludeMax) || (viewMin>=BCpeakRangeMin) || (viewMax<=BCpeakRangeMax) ||
       (excludeMin>BCpeakRangeMin) || (excludeMax<BCpeakRangeMax) || (viewMin>=viewMax) || 
       (excludeMin>=excludeMax) || (BCpeakRangeMin>=BCpeakRangeMax) ){
    Printf("###### ERROR: INVALID LIMITS FOR BIN COUNTING");
    return 0;
  }
  
  Printf("###############################");
  Printf("############################### \n BIN COUNTING SETTINGS:\n############################### \n");
  Printf("Side bands for bg fits: %6.4f-%6.4f && %6.4f-%6.4f \n S+B histo integral range: %6.4f-%6.4f", 
	 viewMin, excludeMin, excludeMax, viewMax, BCpeakRangeMin, BCpeakRangeMax);
  
  // prepare output
  myFitResult *result = new myFitResult(Form("%s_result", hist->GetName()), Form("%s bin counting",fit->Functions()));
  
  //preliminary computation of bg
  TF1  *fbgtmp = fit->ComputeBgF1(hist, viewMin, viewMax, peakMin, peakMax);  
  Int_t  nbg   = fit->GetBgNumPar();
  // create function for sum + bg
  //result->fSum = new TF1("fsum", fit, &myFitFcn::Sum, viewMin, viewMax, nsig + nbg, "myFitFcn", "Sum");
  //result->fBg = new TF1("fbg", fit, &myFitFcn::Sum, viewMin, viewMax, nbg, "myFitFcn", "Bg"); 
  
  fit->EnableBgFitExclusionRange(excludeMin, excludeMax);
  TF1 * bgFitter = new TF1("my_fbg", fit, &myFitFcn::Bg , viewMin, viewMax, nbg , "myFitFcn", "Bg" );
  
  // initialize background from preliminary computation for the bg only
  for (Int_t i = 0; i < nbg; i++) bgFitter->SetParameter(i, fbgtmp->GetParameter(i));
  
  // initialize signal integral
  // Int_t    iview1 = hist->GetXaxis()->FindBin(viewMin);
  // Int_t    iview2 = hist->GetXaxis()->FindBin(viewMax);
  Int_t    ipeak1 = hist->GetXaxis()->FindBin(BCpeakRangeMin);
  Int_t    ipeak2 = hist->GetXaxis()->FindBin(BCpeakRangeMax);
  
  
  //fit only with bg function
  hist->Fit(bgFitter, "RQ0");
  gMinuit->Command("SET STRATEGY 2");
  TFitResultPtr myR = hist->Fit(bgFitter, "ERQS");   
  
  result->fBg = new TF1("fbg", fit, &myFitFcn::Sum, viewMin, viewMax, nbg, "myFitFcn", "Bg"); 
  // copy parameters
  result->CopyParams(bgFitter);
  // copy BG
  for (Int_t i = 0; i < nbg; i++) {
    result->fBg->SetParameter(i, myR->Parameter(i));
  }

  Double_t bgIntPeak = result->fBg->Integral(0.8315, 0.9605);
  Double_t bgIntPeakErr = result->fBg->IntegralError(BCpeakRangeMin, BCpeakRangeMax);
  
  result->Mass   [0]  = 0.896;//result->fSum->GetParameter(1);
  result->Mass   [1]  = 0.0;//result->fSum->GetParError(1);
  result->Gamma[0] = 0.0505; //result->fSum->GetParameter(2);
  result->Gamma[1] = 0.0; //result->fSum->GetParError(2);
  result->Sigma[0] = 0.003; //result->fSum->GetParameter(3);
  result->Sigma[1] = 0.0; //result->fSum->GetParError(3);
  result->Pt   [0]    = pt;
  result->Pt   [1]    = 0.0; //set to 0.0 - not used for the moment
  
  result->Chi2   [0]  = bgFitter->GetChisquare();
  result->Chi2   [1]  = bgFitter->GetNDF();

  //integral of the extrapolated fitted function under the peak
  result->BgInt  [0]  = bgIntPeak; //result->fBg->Integral(BCpeakRangeMin, BCpeakRangeMax);
  result->BgInt  [1]  = bgIntPeakErr ;//result->fBg->IntegralError(BCpeakRangeMin, BCpeakRangeMax);
  Printf("=============================\n I(bg function) peak %e\n", result->BgInt  [0]);

  //integral of fitted bg function in the fit range + extrapolation in excluded range
  result->FcnInt [0]  = bgFitter->Integral(viewMin, viewMax); 
  result->FcnInt [1]  = bgFitter->IntegralError(viewMin, viewMax);
  
  //integral of the histo under the peak
  result->HistInt[0]  = hist->IntegralAndError(ipeak1, ipeak2, result->HistInt[1]);
  result->HistInt[0] *= hist->GetBinWidth(2);
  result->HistInt[1] *= hist->GetBinWidth(2);
  
  //integral of histo under the peak - integralof function ==> NB: to be corrected for the tails to get the raw yields!!!
  result->Raw    [0]  = result->HistInt[0] - result->BgInt[0];
  result->Raw    [1]  = TMath::Sqrt(result->HistInt[1]*result->HistInt[1] + result->BgInt[1]*result->BgInt[1]);
  result->Print("RIBHF");
  
  return result;
}
