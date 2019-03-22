//
// Class to store fit results
//

#include <TF1.h>
#include <TNamed.h>
#include <TString.h>
#include <TPaveText.h>

#ifndef MYFITRESULT_CLASS
#define MYFITRESULT_CLASS

class myFitResult : public TNamed {

public:

  myFitResult(const char *name = "", const char *title = "");
  virtual ~myFitResult() {delete fSum; delete fBg;}
  
  Double_t Raw    [2]; //raw counts
  Double_t Mass   [2]; //mass parameter
  Double_t Gamma  [2]; // width parameter (eg. form Breit Wigner fit)
  Double_t Sigma  [2]; //resolution parameter (eg.from voigtian fit)
  Double_t Fcn    [10][2]; //fit function parameters
  Double_t FcnInt [2]; //integral of fit function in [a,b]
  Double_t BgInt  [2]; //integral of background fit function in [a,b]
  Double_t Chi2   [2]; //chi2 and dof
  Double_t HistInt[2]; //integral of histogram in [a,b]
  Double_t Pt     [2]; // transverse momentum (eg. needed for Boltzman factor fit)
  Double_t TailRaw[2]; // raw counts in  (inf, a] and [b, inf)
  Double_t TailFcnInt[2]; // integral of fit function in  (inf, a] and [b, inf)
  Double_t TailBgInt[2]; // integral of background function in  (inf, a] and [b, inf)
  
  TF1     *fSum;  //!
  TF1     *fBg;   //!
   
  void       CopyParams(TF1 *fcn);
  void       SetParams (TF1 *fcn);
  void       Print(Option_t *opt = "MWSR") const;
  TPaveText* Pave(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *opt = "MWSR");
   
private:

  ClassDef(myFitResult,1)
};

ClassImp(myFitResult);

myFitResult::myFitResult(const char *name, const char *title) : 
  TNamed(name, title),
  fSum(0x0),
  fBg(0x0)
{
  //
  // Dummy constructor
  //

  Int_t i, j;
  for (i = 0; i < 2; i++) {
    Raw    [i] = 0.0; 
    Mass   [i] = 0.0;
    Gamma  [i] = 0.0;
    Sigma  [i] = 0.0;
    FcnInt [i] = 0.0;
    BgInt  [i] = 0.0;
    Chi2   [i] = 0.0;
    HistInt[i] = 0.0;
    Pt     [i] = 0.0;
    TailRaw[i] = 0.0;
    TailFcnInt[i] = 0.0;
    TailBgInt[i] = 0.0;
    for (j = 0; j < 10; j++) {
      Fcn[j][i] = 0.0;
    }
  }
}

void myFitResult::CopyParams(TF1 *fcn)
{
  //
  // Copy parameter values from function
  //

  Int_t i, npar = fcn->GetNpar();
   
  for (i = 0; i < npar; i++) {
    Fcn[i][0] = fcn->GetParameter(i);
    Fcn[i][1] = fcn->GetParError(i);
  }
}

void myFitResult::SetParams(TF1 *fcn)
{
  //
  // Copy parameter values from function
  //

  Int_t i, npar = fcn->GetNpar();
  if (npar > 10) npar = 10;
   
  for (i = 0; i < npar; i++) {
    fcn->SetParameter(i, Fcn[i][0]);
  }
}

void myFitResult::Print(Option_t *options) const
{
  //
  // Print values
  //

  TString opt(options);
  opt.ToUpper();
    
  if (opt.Contains("P")) printf("Pt (center)   = %7.5f \n", Pt  [0]);
  if (opt.Contains("M")) printf("Mass          = %7.5f +- %7.5f\n", Mass  [0], Mass  [1]);
  if (opt.Contains("G")) printf("Width         = %7.5f +- %7.5f\n", Gamma [0], Gamma [1]);
  if (opt.Contains("S")) printf("Sigma         = %7.5f +- %7.5f\n", Sigma [0], Sigma [1]);
  cout << endl;
  if (opt.Contains("B")) printf("Hist integral = %12.6f +- %12.6f\n", HistInt[0], HistInt[1]);
  if (opt.Contains("I")) printf("Func integral = %12.6f +- %12.6f\n", FcnInt [0], FcnInt [1]);
  if (opt.Contains("B")) printf("BG   integral = %12.6f +- %12.6f\n", BgInt  [0], BgInt  [1]);
  if (opt.Contains("R")) printf("Raw  counts   = %12.6f +- %12.6f\n", Raw    [0], Raw    [1]);
  if (opt.Contains("T")) {
    printf("Tails raw   = %12.6f +- %12.6f\n", TailRaw[0], TailRaw[1]);  
    printf("Tails integral fcn  = %12.6f +- %12.6f\n", TailFcnInt[0], TailFcnInt[1]);  
    printf("Tails integral bg  = %12.6f +- %12.6f\n", TailBgInt[0], TailBgInt[1]);  
  }
  if (opt.Contains("C")) printf("Chi2/NDF      = %8.6f \n"          , Chi2[0]/Chi2[1]);
  if (opt.Contains("F")) {
    for (Int_t i = 0; i < 10; i++) {
      printf("Fcn param %d   = %15.5f +- %15.5f\n", i, Fcn[i][0], Fcn[i][1]);
    }
  }
}

TPaveText* myFitResult::Pave(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *options)
{
  //
  // Print values
  //

  TPaveText *pave = new TPaveText(x1, y1, x2, y2, "NDC");

  TString opt(options);
  opt.ToUpper();
    
  if (opt.Contains("P")) pave->AddText(Form("Pt (center)   = %7.5f \n", Pt  [0]));
  if (opt.Contains("M")) pave->AddText(Form("Mass     = %7.5f #pm %7.5f", Mass  [0], Mass  [1]));
  if (opt.Contains("G")) pave->AddText(Form("Width    = %7.5f #pm %7.5f", Gamma [0], Gamma [1]));
  if (opt.Contains("S")) pave->AddText(Form("Sigma    = %7.5f #pm %7.5f", Sigma [0], Sigma [1]));
  cout << endl;
  if (opt.Contains("B")) pave->AddText(Form("Hist int = %7.1f #pm %7.3f", HistInt[0], HistInt[1]));
  if (opt.Contains("I")) pave->AddText(Form("Func int = %7.1f #pm %7.3f", FcnInt [0], FcnInt [1]));
  if (opt.Contains("B")) pave->AddText(Form("BG   int = %7.1f #pm %7.3f", BgInt  [0], BgInt  [1]));
  if (opt.Contains("R")) pave->AddText(Form("Raw  cts = %7.1f #pm %7.3f", Raw    [0], Raw    [1]));
  if (opt.Contains("T")) pave->AddText(Form("Tails  = %7.1f #pm %7.3f", TailRaw[0], TailRaw[1]));
  if (opt.Contains("C")) pave->AddText(Form("Chi2/NDF = %8.6f", Chi2[0]/Chi2[1]));
  if (opt.Contains("F")) {
    for (Int_t i = 0; i < 10; i++) {
      pave->AddText(Form("Fcn param %d   = %15.5f #pm %15.5f", i, Fcn[i][0], Fcn[i][1]));
    }
  }
   
  return pave;
}


#endif
