#include "Riostream.h"
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TArray.h>
#include <TMinuit.h>

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
  
  Printf("############################### \n BIN COUNTING SETTINGS:\n----------------------------------\n");
  Printf("Side bands for bg fits: %6.4f-%6.4f && %6.4f-%6.4f \n S+B histo integral range: %6.4f-%6.4f \n ###############################", 
	 viewMin, excludeMin, excludeMax, viewMax, BCpeakRangeMin, BCpeakRangeMax);
  
  // prepare output
  myFitResult *result = new myFitResult(Form("%s_result", hist->GetName()), Form("%s bin counting",fit->Functions()));
   
   // preliminary computation:
   // estimate background and subtract
   TF1  *fbgtmp = fit->ComputeBgF1(hist, viewMin, viewMax, peakMin, peakMax);
   
   // create function objects
   //Int_t  nsig  = fit->GetSignalNumPar();
   Int_t  nbg   = fit->GetBgNumPar();
   
   fit->EnableBgFitExclusionRange(excludeMin, excludeMax);
   result->fBg  = new TF1("_fbg", fit, &myFitFcn::BgBC , viewMin, viewMax, nbg       , "myFitFcn", "BgBC" );
   // initialize background from preliminary computation
   for (Int_t i = 0; i < nbg; i++) result->fBg->SetParameter(i, fbgtmp->GetParameter(i));
   
   // initialize signal integral
   // Int_t    iview1 = hist->GetXaxis()->FindBin(viewMin);
   // Int_t    iview2 = hist->GetXaxis()->FindBin(viewMax);
   Int_t    ipeak1 = hist->GetXaxis()->FindBin(BCpeakRangeMin);
   Int_t    ipeak2 = hist->GetXaxis()->FindBin(BCpeakRangeMax);
   
   // draw & fit
   hist->Fit(result->fBg, "RQ0");
   gMinuit->Command("SET STRATEGY 2");
   hist->Fit(result->fBg, "ERQ0");
   
   //create function with extrapolation in exclusion range region
   //   fit->DisableBgFitExclusionRange();
   TF1 * bgFitter = new TF1("cp_fbg", fit, &myFitFcn::Bg , viewMin, viewMax, nbg, "myFitFcn","Bg" );
   for (Int_t i = 0; i < nbg; i++) bgFitter->SetParameter(i, result->fBg->GetParameter(i));
   bgFitter->SetLineColor(kGray+1);
   bgFitter->SetLineStyle(7);
   hist->GetListOfFunctions()->Add(bgFitter);
   
   //store 2 separate functions for visualization
   TF1 *fleft = new TF1("leftSide_fbg", fit, &myFitFcn::Bg , viewMin, excludeMin, nbg, "myFitFcn","Bg" );
   fleft->SetParameters(bgFitter->GetParameters());
   fleft->SetLineColor(kMagenta);
   fleft->SetLineStyle(1);
   hist->GetListOfFunctions()->Add(fleft);
   //gROOT->GetListOfFunctions()->Remove(fleft);

   TF1 *fright = new TF1("rightSide_fbg", fit, &myFitFcn::Bg , excludeMax, viewMax, nbg, "myFitFcn","Bg" );
   fright->SetParameters(bgFitter->GetParameters());   
   fright->SetLineColor(kMagenta);
   fright->SetLineStyle(1);   
   hist->GetListOfFunctions()->Add(fright);
   //gROOT->GetListOfFunctions()->Remove(fright);

   //Save function attached to histo
   //hist->GetListOfFunctions()->Add(bgFitter);   
   //Printf("I(full) = %e, I(peak) = %e",bgFitter->Integral(viewMin, viewMax), bgFitter->Integral(BCpeakRangeMin, BCpeakRangeMax));

   // copy parameters
   result->CopyParams(result->fBg);
   
   result->Mass   [0]  = 0.896;//result->fSum->GetParameter(1);
   result->Mass   [1]  = 0.0;//result->fSum->GetParError(1);
   result->Gamma[0] = 0.0505; //result->fSum->GetParameter(2);
   result->Gamma[1] = 0.0; //result->fSum->GetParError(2);
   result->Sigma[0] = 0.003; //result->fSum->GetParameter(3);
   result->Sigma[1] = 0.0; //result->fSum->GetParError(3);
   result->Pt   [0]    = pt;
   result->Pt   [1]    = 0.0; //set to 0.0 - not used for the moment
   result->Chi2   [0]  = result->fBg->GetChisquare();
   result->Chi2   [1]  = result->fBg->GetNDF();
   //integral of bg function in the full range --> to be compared with B from Like-sign
   result->FcnInt [0]  = result->fBg->Integral     (viewMin, viewMax);
   result->FcnInt [1]  = result->fBg->IntegralError(viewMin, viewMax);
   //integral of bg function in the peak range 
   result->BgInt  [0]  = bgFitter->Integral(BCpeakRangeMin, BCpeakRangeMax);
   result->BgInt  [1]  = bgFitter->IntegralError(BCpeakRangeMin, BCpeakRangeMax);
   //integral of histogram in the peak range
   result->HistInt[0]  = hist->IntegralAndError(ipeak1, ipeak2, result->HistInt[1]);
   result->HistInt[0] *= hist->GetBinWidth(2);
   result->HistInt[1] *= hist->GetBinWidth(2);
   //raw signal in the peak range --> to be corrected for the tails!
   result->Raw    [0]  = result->HistInt[0]-result->BgInt[0];
   result->Raw    [1]  = TMath::Sqrt(result->HistInt[1]*result->HistInt[1] + result->BgInt[1]*result->BgInt[1]);
   
   //   result->Print("RMGSIBH");
   
   return result;
}
