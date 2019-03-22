
#include <TF1.h>
#include <TNamed.h>
#include <TString.h>
#include <TPaveText.h>
#include "myFitOutput.h"

ClassImp(myFitOutput);

myFitOutput::myFitOutput(const char *name, const char *title) : 
   TNamed(name, title)
{
  //
  // Dummy constructor
  //
   Int_t i, j;
   for (i = 0; i < 2; i++) {
      fRaw    [i] = 0.0; 
      fMass   [i] = 0.0;
      fGamma  [i] = 0.0;
      fSigma  [i] = 0.0;
      fFcnInt [i] = 0.0;
      fBgInt  [i] = 0.0;
      fChi2   [i] = 0.0;
      fHistInt[i] = 0.0;
      fHistInt2Gamma[i] = 0.0;
      fFcnInt2Gamma[i] = 0.0;
      fRaw2Gamma[i] = 0.0;
      fBgInt2Gamma[i] = 0.0;
      fFcnIntTail[i] = 0.0;
      fRawIntTail[i] = 0.0;
      fBgIntTail[i] = 0.0;
      for (j = 0; j < 14; j++) {
         fFcn[j][i] = 0.0;
      }
   }
}

void myFitOutput::CopyFcnParams(TF1 *fcn)
{
//
// Copy parameter values from function
//
   Int_t i, npar = fcn->GetNpar();
   
   for (i = 0; i < npar; i++) {
      fFcn[i][0] = fcn->GetParameter(i);
      fFcn[i][1] = fcn->GetParError(i);
   }
}

void myFitOutput::SetFcnParams(TF1 *fcn)
{
//
// Copy parameter values from function
//

   Int_t i, npar = fcn->GetNpar();
   if (npar > 14) npar = 14;
   
   for (i = 0; i < npar; i++) {
      fcn->SetParameter(i, fFcn[i][0]);
   }
}

void myFitOutput::Print(Option_t *options) const
{
//
// Print values
//
   TString opt(options);
   opt.ToUpper();
   if (opt.Contains("M")) printf("Mass          = %7.5f +- %7.5f\n", fMass  [0], fMass  [1]);
   if (opt.Contains("W")) printf("Width         = %7.5f +- %7.5f\n", fGamma [0], fGamma [1]);
   if (opt.Contains("S")) printf("Sigma         = %7.5f +- %7.5f\n", fSigma [0], fSigma [1]);
   if (opt.Contains("B")) printf("Hist integral = %10.1f +- %10.2f\n", fHistInt[0], fHistInt[1]);
   if (opt.Contains("I")) printf("Func integral = %10.1f +- %10.2f\n", fFcnInt [0], fFcnInt [1]);
   if (opt.Contains("B")) printf("BG   integral = %10.1f +- %10.2f\n", fBgInt  [0], fBgInt  [1]);
   if (opt.Contains("R")) printf("Raw  counts   = %10.1f +- %10.2f\n", fRaw    [0], fRaw    [1]);
   
   if (opt.Contains("T")) printf("Total func tails integral = %10.1f +- %10.2f\n Signal func tails integral = %10.1f +- %10.2f\n Bg func tails integral = %10.1f +- %10.2f\n", fFcnIntTail[0], fFcnIntTail[1],fRawIntTail[0], fRawIntTail[1],fBgIntTail[0], fBgIntTail[1]);  
   
   if (opt.Contains("2G")) printf("Histogram integral in +/-2Width = %10.1f +- %10.2f\n Total func integral in +/-2Width = %10.1f +- %10.2f\n Signal func integral in +/-2Width = %10.1f +- %10.2f\n Bg func integral in in +/-2Width = %10.1f +- %10.2f\n", fHistInt2Gamma[0], fHistInt2Gamma[1], fFcnInt2Gamma[0], fFcnInt2Gamma[1], fRaw2Gamma[0], fRaw2Gamma[1], fBgInt2Gamma[0], fBgInt2Gamma[1]);  
   
   if (opt.Contains("C")) printf("Chi2/NDF      = %8.6f \n", fChi2[0]/fChi2[1]);
   if (opt.Contains("F")) {
      for (Int_t i = 0; i < 10; i++) {
         printf("Fcn param %d   = %15.5f +- %15.5f\n", i, fFcn[i][0], fFcn[i][1]);
      }
   }
}

TPaveText* myFitOutput::Pave(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *options)
{
//
// Print values
//
   TPaveText *pave = new TPaveText(x1, y1, x2, y2, "NDC");
   TString opt(options);
   opt.ToUpper();
    
   if (opt.Contains("M")) pave->AddText(Form("Mass     = %7.5f #pm %7.5f", fMass  [0], fMass  [1]));
   if (opt.Contains("W")) pave->AddText(Form("Width    = %7.5f #pm %7.5f", fGamma [0], fGamma [1]));
   if (opt.Contains("S")) pave->AddText(Form("Sigma    = %7.5f #pm %7.5f", fSigma [0], fSigma [1]));
 
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

   if (opt.Contains("C")) pave->AddText(Form("Chi2/NDF = %8.6f", fChi2[0]/fChi2[1]));
   if (opt.Contains("F")) {
      for (Int_t i = 0; i < 10; i++) {
         pave->AddText(Form("Fcn param %d   = %15.5f #pm %15.5f", i, fFcn[i][0], fFcn[i][1]));
      }
   }
   
   return pave;
}

