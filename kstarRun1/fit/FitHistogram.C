#include "Riostream.h"
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TArray.h>
#include <TMinuit.h>
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "/Users/fbellini/alice/macros/ResonAnT/fit/CPolyFit.C"
#include "/Users/fbellini/alice/macros/ResonAnT/fit/myFitFcn.C"
#include "/Users/fbellini/alice/macros/ResonAnT/fit/myFitResult.C"

myFitResult* FitHistogram(TH1D       *hist,
			  myFitFcn   *fit,
			  Double_t    viewMin  = 0.995,
			  Double_t    viewMax  = 1.06,
			  Double_t    peakMin  = 1.005,
			  Double_t    peakMax  = 1.035,
			  Double_t    mass     = 1.019455,
			  Double_t    gamma    = 0.00426,
			  Double_t    sigma    = 0.001,
			  Bool_t      fixGamma = kFALSE,
			  Double_t    minGamma = 0.5,
			  Double_t    maxGamma = 1.5,
			  Bool_t      fixSigma = kFALSE,
			  Double_t    minSigma = 0.5,
			  Double_t    maxSigma = 1.5,
			  Double_t    pt       = 0.01,
			  TString     opt      = "",
			  Double_t    binCountMin  = 0.995,
			  Double_t    binCountMax  = 0.150)
{
   // prepare output
   myFitResult *result = new myFitResult(Form("%s_result", hist->GetName()), fit->Functions());
   
   // preliminary computation:
   // estimate background and subtract
   TF1  *fbgtmp = fit->ComputeBgF1(hist, viewMin, viewMax, peakMin, peakMax);
   TH1D *hsub   = (TH1D*)hist->Clone("sub");
   hsub->Add(fbgtmp, -1.0);
   
   // create function objects
   Int_t  nsig  = fit->GetSignalNumPar();
   Int_t  nbg   = fit->GetBgNumPar();
   result->fSum = new TF1("fsum", fit, &myFitFcn::Sum, viewMin, viewMax, nsig + nbg, "myFitFcn", "Sum");
   result->fBg  = new TF1("_fbg", fit, &myFitFcn::Bg , viewMin, viewMax, nbg       , "myFitFcn", "Bg" );
   
   // initialize background from preliminary computation
   for (Int_t i = 0; i < nbg; i++) result->fSum->SetParameter(nsig + i, fbgtmp->GetParameter(i));
   
   //define extremes for integral
   Int_t    iview1 = hsub->GetXaxis()->FindBin(viewMin);
   Int_t    iview2 = hsub->GetXaxis()->FindBin(viewMax);
   
   //define extremes for bin counting
   if ((binCountMin<0.0) && (binCountMax<0.0)) {
     binCountMin = mass - 2*gamma;
     binCountMax = mass + 2*gamma;
   }
   Int_t    ibincount1 = hsub->GetXaxis()->FindBin(binCountMin);
   Int_t    ibincount2 = hsub->GetXaxis()->FindBin(binCountMax);
   
   // initialize signal integral
   Double_t sgInt  = hsub->Integral(iview1, iview2) * hsub->GetBinWidth(2);
   result->fSum->SetParameter(0, sgInt);
   // initialize peak position (common)
   // --> we don't expect shifts extremely large
   result->fSum->SetParameter(1, mass);
   result->fSum->SetParLimits(1, peakMin, peakMax);
   
   // initialize peak width
   myFitFcn::ESignal sigType = (myFitFcn::ESignal) fit->GetSignalType();
   switch (sigType) {
   case myFitFcn::kVoigtian:
     if (fixGamma) result->fSum->FixParameter(2, gamma); 
     else {
       result->fSum->SetParameter(2, gamma);
       result->fSum->SetParLimits(2, minGamma * gamma, maxGamma * gamma);
     }
     if (fixSigma) result->fSum->FixParameter(3, sigma); 
     else {
       result->fSum->SetParameter(3, sigma);
       result->fSum->SetParLimits(3, sigma * minSigma, sigma * maxSigma);
     }
     break;
   case myFitFcn::kBreitWigner:
     if (fixGamma) result->fSum->FixParameter(2, gamma); 
     else {
       result->fSum->SetParameter(2, gamma);
       result->fSum->SetParLimits(2, minGamma * gamma, maxGamma * gamma);
     }
     break;
   case myFitFcn::kRelBreitWigner:
     if (fixGamma) result->fSum->FixParameter(2, gamma); 
     else {
       result->fSum->SetParameter(2, gamma);
       result->fSum->SetParLimits(2, minGamma * gamma, maxGamma * gamma);
     }
     break;
   case myFitFcn::kRelBreitWignerBoltzmann:
     result->fSum->FixParameter(3, pt); 
     if (fixGamma) result->fSum->FixParameter(2, gamma); 
     else {
       result->fSum->SetParameter(2, gamma);
       result->fSum->SetParLimits(2, minGamma * gamma, maxGamma * gamma);
     }
     break;
   case myFitFcn::kGaus:
     result->fSum->SetParameter(2, sigma);
     result->fSum->SetParLimits(3, sigma * minSigma, sigma * maxSigma);
     break;
   default:
     ::Error("FitHistogram:::: FitSingleMix", "Wrong function type");
     return 0x0;
   }
   
   // draw & fit
   gSystem->Load("libMathCore.so");
   ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(1000);
   ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000);
   TCanvas * ctest = new TCanvas("ctest", "ctest", 800, 600);
   ctest->cd(); hist->Draw(); ctest->Print("test.png");
   hist->Fit(result->fSum, "ERQ0M");
   gMinuit->Command("SET STRATEGY 2");
   hist->Fit(result->fSum, "ERQ0M");

   // copy parameters
   result->CopyParams(result->fSum);
   result->Chi2   [0]  = result->fSum->GetChisquare();
   result->Chi2   [1]  = result->fSum->GetNDF();
   result->Raw    [0]  = result->fSum->GetParameter(0);
   result->Raw    [1]  = result->fSum->GetParError(0);
   result->Mass   [0]  = result->fSum->GetParameter(1);
   result->Mass   [1]  = result->fSum->GetParError(1);
   
   result->Pt   [0]    = pt; //result->fSum->GetParameter(3);
   result->Pt   [1]    = 0.0; //result->fSum->GetParError(3); //set to 0.0 - not used for the moment
   
   switch (sigType) {
      case myFitFcn::kVoigtian:
         result->Gamma[0] = result->fSum->GetParameter(2);
         result->Gamma[1] = result->fSum->GetParError(2);
         result->Sigma[0] = result->fSum->GetParameter(3);
         result->Sigma[1] = result->fSum->GetParError(3);
         break;
      case myFitFcn::kBreitWigner:
         result->Gamma[0] = result->fSum->GetParameter(2);
         result->Gamma[1] = result->fSum->GetParError(2);
         break;
      case myFitFcn::kRelBreitWigner:
         result->Gamma[0] = result->fSum->GetParameter(2);
         result->Gamma[1] = result->fSum->GetParError(2);
         break;
      case myFitFcn::kRelBreitWignerBoltzmann:
         result->Gamma[0] = result->fSum->GetParameter(2);
         result->Gamma[1] = result->fSum->GetParError(2);
         break;
      case myFitFcn::kGaus:
         result->Sigma[0] = result->fSum->GetParameter(2);
         result->Sigma[1] = result->fSum->GetParError(2);
         break;
      default:
         ::Error("FitHistogram:::: FitSingleMix", "Wrong function type");
         return 0x0;
   }
   Printf("FitHistogram:::: sigma = %.6f  minSigma = %.6f maxSigma =%.6f",sigma,minSigma*sigma, maxSigma*sigma);

   //scale by bin width
   result->HistInt[0] *= hist->GetBinWidth(2);
   result->HistInt[1] *= hist->GetBinWidth(2);   

   Printf("FitHistogram:::: getting ready to integrate fitted function");
   //ROOT::Math::GSLIntegrator (1.E-9, 1E-6, 500);

   if (opt.Contains("BC")){ 
     //bin counting version!
     result->HistInt[0]  = hist->IntegralAndError(ibincount1, ibincount2, result->HistInt[1]);
     result->FcnInt [0]  = result->fSum->Integral(binCountMin,binCountMax);
     result->FcnInt [1]  = result->fSum->IntegralError(binCountMin,binCountMax);
     result->TailFcnInt [0] = result->fSum->Integral(mass-10*gamma, binCountMin) + result->fSum->Integral(binCountMax, mass+10*gamma);
     result->TailFcnInt [1] = result->fSum->IntegralError(mass-10*gamma, binCountMin) + result->fSum->IntegralError(binCountMax, mass+10*gamma);

     TF1 * bgFitter = new TF1("cp_fbg", fit, &myFitFcn::Bg , viewMin, viewMax, nbg, "myFitFcn","Bg" );
     for (Int_t i = nsig; i < nsig + nbg; i++){ 
       bgFitter->SetParameter(i-nsig, result->fSum->GetParameter(i));
     }
     Printf("FitHistogram:::: Integral = %e +/- %e", bgFitter->Integral(binCountMin,binCountMax), bgFitter->IntegralError(binCountMin,binCountMax));
     result->BgInt  [0]  = bgFitter->Integral(binCountMin,binCountMax);
     result->BgInt  [1]  = bgFitter->IntegralError(binCountMin,binCountMax);
     result->TailBgInt [0] = bgFitter->Integral(mass-10*gamma, binCountMin) + bgFitter->Integral(binCountMax, mass+10*gamma);;
     result->TailBgInt [1] = 0.0;  //bgFitter->IntegralError(mass-10*gamma, binCountMin) + bgFitter->IntegralError(binCountMax, mass+10*gamma);;
     result->TailRaw [0] = result->TailFcnInt [0] - result->TailBgInt [0];
     result->TailRaw [1] =  result->TailRaw [0]* (result->TailFcnInt [1]/result->TailFcnInt [0]);

   } else {
     
     //pure fit version!
     result->FcnInt [0]  =  result->fSum->Integral     (viewMin, viewMax);
     result->FcnInt [1]  = result->fSum->IntegralError(viewMin, viewMax);
     result->BgInt  [0]  = result->FcnInt[0] - result->Raw[0];
     result->BgInt  [1]  = TMath::Sqrt(result->FcnInt[1]*result->FcnInt[1] + result->Raw[1]*result->Raw[1]);
     result->HistInt[0]  = hist->IntegralAndError(iview1, iview2, result->HistInt[1]);
     result->TailFcnInt [0] = 0.0;
     result->TailFcnInt [1] = 0.0;
     result->TailBgInt  [0] = 0.0;
     result->TailBgInt  [1] = 0.0;
     result->TailRaw    [0] = 0.0;
     result->TailRaw    [1] = 0.0;
   }
   result->Print("RMGSIBH");
   
   // copy BG
   for (Int_t i = 0; i < nbg; i++) {
     result->fBg->SetParameter(i, result->fSum->GetParameter(i + nsig));
   }
   
   // return 
   return result;
}
