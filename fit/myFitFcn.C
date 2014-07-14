#include <iostream>

#include <TF1.h>
#include <TH1.h>
#include <TMath.h>
#include <TError.h>
#include <TArray.h>
#include <TString.h>
#include <TClonesArray.h>

#include "CPolyFit.C"

class myFitFcn {

public:

   enum EBackground {
      kPoly1,
      kPoly2,
      kPoly3,
      kExp,
      kBackgrounds
   };
   
   enum ESignal {
      kVoigtian,
      kBreitWigner,
      kRelBreitWigner,
      kRelBreitWignerBoltzmann,
      kGaus
   };
      
   myFitFcn(double rangeA, double rangeB, ESignal sig = kVoigtian, EBackground bg = kPoly1);
   myFitFcn(double rangeA, double rangeB, Option_t *options);
   
   int     SetSignal(ESignal type);
   int     SetBg(EBackground type);
   int     SetFunctions(Option_t *option);
   void    SetIntegralRange(double a, double b);
   void    SetIntegralDivisions(int n)        {fIntRangeDiv = n;} 
   void    SetPartIntegral(Bool_t yn = kTRUE) {fPartIntegral = yn;}
           
   int     GetSignalNumPar()  {return fNSignal;}
   int     GetBgNumPar()      {return fNBg;}
   int     GetSignalType()    {return fSigType;}
   int     GetBgType()        {return fBgType;}
   double  GetIntRangeMin()   {return fIntRangeMin;}
   double  GetIntRangeMax()   {return fIntRangeMax;}
           
  double  Signal (double *x, double *param);
  double  SigNorm(double *x, double *param);
  double  Bg     (double *x, double *param);
  double  BgBC     (double *x, double *param);
  double  Sum    (double *x, double *param);
  double  SigInt (double *param);
   
   TF1*        ComputeBgF1(TH1D *hist, Double_t fitMin, Double_t fitMax, Double_t peakMin, Double_t peakMax);
   TF1*        GetBgF1() {return fBGF1;}
   const char* Functions();
  double       RelBreitWigner(double x, double *par);
  double       RelBreitWignerBoltzmann(double x, double *par);

  //New
  double       arg(double x0, double x1, double x2);
  double       PS(double m, double pT, double T);
  double       bw(double m, double m0, double Gamma);
  double       bw1(double *x, double *par);
  double       BWPS(double *x, double *par);

  void         SetPt(double pt) {fPt = pt;}
  void         EnableBgFitExclusionRange(double min, double max);
  void         DisableBgFitExclusionRange() {fEnaExclude = kTRUE;Printf("BgFitExclusion disabled");};
private:

  Bool_t       fPartIntegral;  // if true, the signal integral is computed in the specified range
  int          fNSignal;       // number of parameters for signal
  int          fNBg;           // number of parameters for background
  ESignal      fSigType;       // chosen function for signal
  EBackground  fBgType;        // chosen function for background
  
  double       fIntRangeMin;   // lower bound of background integration range
  double       fIntRangeMax;   // upper bound of background integration range
  int          fIntRangeDiv;   // number of divisions in the integral range
  double       fPt;             // pt of the selected bin
  TF1         *fBGF1;          // temp function for BG
  Bool_t       fEnaExclude;    // if true, the backgound fit range exclusion is enablesd - for bin counting
  Double_t     fExcludeMin;    //bg fit exclusion range lower bound
  Double_t     fExcludeMax;    //bg fit exclusion range upper bound
};
//
// == IMPLEMENTATION ===============================================================================
//

//__________________________________________________________________________________________________
//
// Constructor
//
myFitFcn::myFitFcn(double a, double b, ESignal sig, EBackground bg) :
   fPartIntegral(kFALSE),
   fNSignal(0),
   fNBg(0),
   fSigType(kVoigtian),
   fBgType(kBackgrounds),
   fIntRangeMin(0.0),
   fIntRangeMax(0.0),
   fIntRangeDiv(10),
   fPt(0.01),
   fBGF1(0x0),
   fEnaExclude(kFALSE),
   fExcludeMin(0.0),
   fExcludeMax(0.0)
{
  SetSignal(sig);
  SetBg(bg);
  SetIntegralRange(a, b);
}

//__________________________________________________________________________________________________
//
// Constructor
//
myFitFcn::myFitFcn(double a, double b, Option_t *option) :
   fPartIntegral(kFALSE),
   fNSignal(0),
   fNBg(0),
   fSigType(kVoigtian),
   fBgType(kBackgrounds),
   fIntRangeMin(0.0),
   fIntRangeMax(0.0),
   fIntRangeDiv(10),
   fPt(0.01),
   fBGF1(0x0),
   fEnaExclude(kFALSE),
   fExcludeMin(0.0),
   fExcludeMax(0.0)
{
   SetFunctions(option);
   SetIntegralRange(a, b);
}

//__________________________________________________________________________________________________
//
// Set signal type.
// Defines and returns also the number of parameters
//
int myFitFcn::SetSignal(ESignal type)
{
   fSigType = type;
   
   switch (fSigType) {
   case kVoigtian:                fNSignal = 4; break;
   case kBreitWigner:             fNSignal = 3; break;
   case kRelBreitWigner:          fNSignal = 3; break;
   case kRelBreitWignerBoltzmann: fNSignal = 4; break;
   case kGaus:                    fNSignal = 3; break;
   default:                       fNSignal = 0; break;
   }
   
   return fNSignal;
}

//__________________________________________________________________________________________________
//
// Set background type.
// Defines and returns also the number of parameters
//
int myFitFcn::SetBg(EBackground type)
{
   fBgType = type;
   
   switch (fBgType) {
      case kPoly1: fNBg = 2; break;
      case kPoly2: fNBg = 3; break;
      case kPoly3: fNBg = 4; break;
      case kExp:   fNBg = 2; break;
      default    : fNBg = 0; break;
   }
   
   return fNBg;
}

//__________________________________________________________________________________________________
//
// Set any function type, using a string:
// -- signal    : "VOIGT", "BW", "REL", "BOLZ" or "GAUS"
// -- background: "POLY1", "POLY2" or "POLY3", "EXP"
//
const char* myFitFcn::Functions()
{
   TString opt;
   
   switch (fSigType) {
   case kVoigtian       : opt.Append("VOIGT"); break;
   case kBreitWigner    : opt.Append("BW"); break;
   case kRelBreitWigner : opt.Append("REL"); break;
   case kRelBreitWignerBoltzmann: opt.Append("BOLZ"); break;
   case kGaus           : opt.Append("GAUS"); break;
   default              : opt.Append("<NULL>"); break;
   }
   
   opt += '+';
   
   switch (fBgType) {
      case kPoly1: opt.Append("POLY1"); break;
      case kPoly2: opt.Append("POLY2"); break;
      case kPoly3: opt.Append("POLY3"); break;
      case kExp  : opt.Append("EXP"); break;
      default    : opt.Append("<NULL>"); break;
   }
   
   return opt.Data();
}
//__________________________________________________________________________________________________
//
// Returns a string describing the functions:
// -- signal    : "VOIGT", "BW" or "GAUS"
// -- background: "POLY1", "POLY2" or "POLY3"
//
int myFitFcn::SetFunctions(Option_t *option)
{
   TString opt(option);
   opt.ToUpper();
   
   // params
   int npar = 0;
   
   // signal
   if (opt.Contains("VOIGT"))
     npar += SetSignal(kVoigtian);
   else if (opt.Contains("BW"))
     npar += SetSignal(kBreitWigner);
   else if (opt.Contains("REL"))
     npar += SetSignal(kRelBreitWigner);
   else if (opt.Contains("BOLZ"))
     npar += SetSignal(kRelBreitWignerBoltzmann);
   else if (opt.Contains("GAUS"))
     npar += SetSignal(kGaus);
   
   // background
   if (opt.Contains("POLY1"))
      npar += SetBg(kPoly1);
   else if (opt.Contains("POLY2"))
      npar += SetBg(kPoly2);
   else if (opt.Contains("POLY3"))
      npar += SetBg(kPoly3);
   else if (opt.Contains("EXP"))
     npar += SetBg(kExp);
   return npar;
}

//__________________________________________________________________________________________________
//
// Set integral range.
// Parameters can be passed in any order
//
void myFitFcn::SetIntegralRange(double a, double b)
{
   fIntRangeMax = TMath::Max(a,b);
   fIntRangeMin = TMath::Min(a,b);
}


//__________________________________________________________________________________________________
//
// Set bg exclusion ranges
//
void myFitFcn::EnableBgFitExclusionRange(double min, double max) 
{
  if (max<min) return;
  fEnaExclude=kTRUE; 
  fExcludeMin=min;
  fExcludeMax=max;
  Printf("Bg fit exclusion range set to [%5.3f, %5.3f]", fExcludeMin, fExcludeMax );
}
//__________________________________________________________________________________________________
//
// Computation of background function
//
double myFitFcn::Bg(double *m, double *par)
{
  double x   = m[0];
  double out = 0.0;

  switch (fBgType) {
      case kExp:
	out += TMath::Exp(par[1] * x);
	 out*= par[0];
	 break;
      case kPoly3: 
         out += par[3] * x * x * x;
      case kPoly2:
         out += par[2] * x * x;
      case kPoly1:
         out += par[1] * x;
         out += par[0];
	 break;
      default:
         return 0.0;
   }
   
   return out;
}

//__________________________________________________________________________________________________
//
// Computation of background function for bin counting - with exclusion range
//
double myFitFcn::BgBC(double *m, double *par)
{
  double x   = m[0];
  double out = 0.0;
  
  if (fEnaExclude && (x>fExcludeMin) && (x<fExcludeMax) ) {
    TF1::RejectPoint();
    return 0;
  } 
  
  switch (fBgType) {
      case kExp:
	out += TMath::Exp(par[1] * x);
	 out*= par[0];
	 break;
      case kPoly3: 
         out += par[3] * x * x * x;
      case kPoly2:
         out += par[2] * x * x;
      case kPoly1:
         out += par[1] * x;
         out += par[0];
	 break;
      default:
         return 0.0;
   }
   
   return out;
}
//__________________________________________________________________________________________________
//
// Computation of signal function (normalized to 1 in full range)
//
double myFitFcn::SigNorm(double *m, double *par)
{
   double x     = m[0];
   double value = 0.0;
   
   switch (fSigType) {
      case kVoigtian: 
         value = TMath::Voigt(x - par[1], par[3], par[2], 4);
         break;
      case kBreitWigner: 
         value = TMath::BreitWigner(x, par[1], par[2]);
         break;
     case kRelBreitWignerBoltzmann: 
         //       value = RelBreitWignerBoltzmann(x, par);
         value = BWPS( m, par);
         break;
      case kRelBreitWigner: 
        value = RelBreitWigner(x, par);
	break;
      case kGaus: 
         value = TMath::Gaus(x, par[1], par[2], kTRUE);
         break;
   }
   
   return value;
}



//__________________________________________________________________________________________________
//
// Computation of signal function
//
double myFitFcn::Signal(double *m, double *par)
{
   return par[0] * SigNorm(m, par);
}

//__________________________________________________________________________________________________
//
// Computation of signal integral using the Cavalieri-Simpson rule
//
double myFitFcn::SigInt(double *param)
{
   Double_t h = (fIntRangeMax - fIntRangeMin) / fIntRangeDiv;
   Double_t sum = SigNorm(&fIntRangeMin, param) + SigNorm(&fIntRangeMax, param);
   Double_t x, y;
   
   for (int i = 1; i < fIntRangeDiv - 1; i++) {
      x = fIntRangeMin + i * h;
      y = SigNorm(&x, param);
      if (i % 2)
         sum += 2.0 * y;
      else
         sum += 4.0 * y;
   }
   
   return sum * h / 3.0;
}

//__________________________________________________________________________________________________
//
// Computation of total function, used for fit.
//
double myFitFcn::Sum(double *x, double *par)
{
   // signal: 
   // function is normalized in full real range (-inf, +inf)
   // and first parameter gives the full normalization
   double sigval, signal, sigint;
   if (fPartIntegral) {
      sigval = SigNorm(x, &par[0]);
      sigint = SigInt (&par[0]);
      signal = par[0] * sigval / sigint;
   } else {
      signal = Signal(x, &par[0]);
   }
   
   // background:
   // function is normalized in the specified integration range (fIntRangeMin, fIntRangeMax)
   // and first parameter gives this integral, which can be used for estimating counts
   double bg = Bg(x, &par[fNSignal]);
      
   // final value
   return signal + bg;
}

//__________________________________________________________________________________________________
//
// Returns a TF1 object shaped like the chosen background,
// and initialized in a proper way
//
TF1 * myFitFcn::ComputeBgF1
(TH1D *hist, Double_t fitMin, Double_t fitMax, Double_t peakMin, Double_t peakMax)
{
   // reset pointer
   if (fBGF1) delete fBGF1;
   
   // prepare variables
   CPolyFit *polyFit = 0x0;
   
   // initialize needed objects
   switch (fBgType) {
      case kPoly3: 
         polyFit = new CPolyFit(3);
         fBGF1 = new TF1("fbg", "pol3", fIntRangeMin, fIntRangeMax);
         break;
      case kPoly2:
         polyFit = new CPolyFit(2);
         fBGF1 = new TF1("fbg", "pol2", fIntRangeMin, fIntRangeMax);
         break;
      case kPoly1:
         polyFit = new CPolyFit(1);
         fBGF1 = new TF1("fbg", "pol1", fIntRangeMin, fIntRangeMax);
         break;
      case kExp:
   	fBGF1 = new TF1("fbg", "[0]*exp([1]*x)", fIntRangeMin, fIntRangeMax);
	break;
      default:
         fBGF1 = 0x0;
         ::Info("ComputeBgF1", "No background in this case");
   }
   
   // if poly fit is done, compute it here
   if (polyFit) {
      polyFit->AddPointsFromHist(hist, fitMin, fitMax, peakMin, peakMax);
      TArrayD poly = polyFit->Compute();
      fBGF1->SetParameters(poly.GetArray());
   } 
   
   // return function
   return fBGF1;
}

double myFitFcn::RelBreitWignerBoltzmann(double x, double * par)
{
  //Parameters: constant = par[0], mass = par[1], width = par[2]
  //transverse momentum is par[3], to be fixed at runtime!
  //par[0] will be use dto extract thenormalization, ie. will return the integral and error of the function
  const Double_t temp=0.160;
  const Double_t mpi2 = 0.1396*0.1396;
  const Double_t mka2 = 0.4937*0.4937;
  Double_t mkpi2 = mpi2+mka2; //GeV
  Double_t arg = 0.0, arg2 = 0.0, arg3=0.0, arg4=0.0, gamma = 0.0;
  Double_t boltz =0.0;
  arg = TMath::Sqrt(x * x + par[3] * par[3]);
  boltz = TMath::Exp(- arg / temp) / arg;
  arg2 = (x*x - par[1]*par[1]) * (x*x - par[1]*par[1]);
  arg3 = TMath::Power(x*x - mkpi2, 2.0) - 4.0*mpi2*mka2;
  arg4 = TMath::Power(par[1]*par[1] - mkpi2, 2.0) - 4.0*mpi2*mka2;
  gamma = par[2] * TMath::Power(par[1]/x, 4.0) * TMath::Power(arg3/arg4 , 1.5);
  return par[0] * x * boltz * (x * par[1] * gamma)/(arg2 + par[1]*par[1]*gamma*gamma);
}

double myFitFcn::RelBreitWigner(double x, double *par)
{
  //calculates relativistic BW x PS
  const double mkpi2 = (0.4937+0.1396)*(0.4937+0.1396); //GeV
  double arg = 0.0, gamma = 0.0;
  gamma = (par[1]*par[2]/x)*TMath::Power((x*x-mkpi2)/(par[1]*par[1]-mkpi2) , 1.5);
  arg = (x*x-par[1]*par[1])*(x*x-par[1]*par[1]);
  return (x*par[1]*gamma)/(arg + par[1]*par[1]*gamma*gamma);
}


//new 
double myFitFcn::arg(double x0, double x1, double x2)
{
  return pow(x0*x0-x1*x1-x2*x2,2.0)-4.*x1*x1*x2*x2;
}

double myFitFcn::PS(double m, double pT, double T)
{
  double mT = sqrt(m*m+pT*pT);
  return m/mT*exp(-mT/T);
}

double myFitFcn::bw(double m, double m0, double Gamma)
{
  return m*m0*Gamma/(pow(m*m-m0*m0,2.0)+m0*m0*Gamma*Gamma);
}

double myFitFcn::bw1(double *x, double *par)
{
  const double MassK = 0.49368;
  const double MassPi = 0.13957;
  double Gamma = par[2]*pow(par[1]/x[0],4.0);
  Gamma *= pow(arg(x[0],MassK,MassPi)/arg(par[1],MassK,MassPi),1.5);
  //double Gamma = par[2];
  return bw(x[0],par[1],Gamma);
}

// double myFitFcn::bw2(double *x, double *par)
// {
//   const double MassK = 0.49368;
//   const double MassPi = 0.13957;
//   double Gamma = par[2]*pow(par[1]/x[0],2.0);
//   Gamma *= pow(arg(x[0],MassK,MassPi)/arg(par[1],MassK,MassPi),1.5);
//   Gamma *= pow((pow(par[1]*LambdaPi,2)+arg(par[1],MassK,MassPi))/(pow(x[0]*LambdaPi,2)+arg(x[0],MassK,MassPi)),2);
//   //double Gamma = par[2];
//   return bw(x[0],par[1],Gamma);
// }

double myFitFcn::BWPS(double *x, double *par)
{
  const double Temp = 0.154;
  //  return par[0]*bw1(x, par)*PS(x[0], par[3], Temp)*1.e6;
  return bw1(x, par)*PS(x[0], par[3], Temp)*1.e6;
}
