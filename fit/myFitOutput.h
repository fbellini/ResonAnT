#ifndef MYFITOUTPUT
#define MYFITOUTPUT

#include <TF1.h>
#include <TNamed.h>
#include <TString.h>
#include <TPaveText.h>

class myFitOutput : public TNamed {

public:

  myFitOutput(const char *name = "", const char *title = "");
  virtual ~myFitOutput();
  void CopyFcnParams(TF1 *fcn);
  void SetFcnParams(TF1 *fcn);
  void       Print(Option_t *opt = "MWSR") const;
  TPaveText* Pave(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *opt = "MWSR");

  Double_t fRaw    [2]; //integral of signal fit function in fit range & error
  Double_t fMass   [2]; //mass parameter & error
  Double_t fGamma  [2]; //width parameter & error
  Double_t fSigma  [2]; //resolution parameter & error
  Double_t fFcn    [14][2]; //total fit function parameters & error
  Double_t fFcnInt [2]; //integral of total fit function in fit range & error
  Double_t fBgInt  [2]; //integral of background fit function in fit range & error
  Double_t fChi2   [2]; //Chi2 and dof
  Double_t fHistInt[2]; //integral of histogram in fit range & error
  Double_t fHistInt2Gamma[2]; //integral of histogram in +/-2Gamma range & error
  Double_t fFcnInt2Gamma[2]; //integral of total fit  function in +/-2Gamma range & error
  Double_t fRaw2Gamma[2]; //integral of signal fit  function in +/-2Gamma range & error
  Double_t fBgInt2Gamma[2]; //integral of background fit  function in +/-2Gamma range & error
  Double_t fFcnIntTail[2]; //integral of total fit function outside +/-2Gamma range & error
  Double_t fBgIntTail[2]; //integral of background fit function outside +/-2Gamma range & error
  Double_t fRawIntTail[2]; //integral of signal fit function outside +/-2Gamma range & error
  
private:

  ClassDef(myFitOutput,1)
};

#endif
