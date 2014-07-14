#include <TH1.h>
#include <TError.h>
#include <TArrayD.h>
#include <TMatrixD.h>

#ifndef CPOLYFIT_CLASS
#define CPOLYFIT_CLASS

class CPolyFit {

public:

   CPolyFit(Int_t n = 2) : fOrder(n), fX(n+1, n+1), fM(n+1, 1) { }
   
   void             Reset();
   void             AddPointW(Double_t x, Double_t y, Double_t weight = 1.0);
   void             AddPointE(Double_t x, Double_t y, Double_t err    = 1.0);
   void             AddPointsFromHist(TH1D *hist, Double_t fitMin, Double_t fitMax, Double_t peakMin = 0.0, Double_t peakMax = 0.0);
   TArrayD          Compute();
   static  Double_t Pow(Double_t x, Int_t n);
   Int_t            Order() {return fOrder;}
   
private:

   Int_t     fOrder;  // polynomial order
   TMatrixD  fX;      // matrix of x powers
   TMatrixD  fM;      // vector of known terms
   
};

//__________________________________________________________________________________________________
//
// Reset all matrices to zero
//
void CPolyFit::Reset()
{


   Int_t i, j;
   
   for (i = 0; i <= fOrder; i++) {
      fM(0,i) = 0.0;
      for (j = 0; j <= fOrder; j++)
         fX(i,j) = 0.0;
   }
}

//__________________________________________________________________________________________________
//
// Include the passed point, with its weight, to the computation
//
void CPolyFit::AddPointW(Double_t x, Double_t y, Double_t weight)
{

   
   Int_t i, j;
   
   for (i = 0; i <= fOrder; i++) {
      fM(i,0) += weight * Pow(x, i) * y;
      for (j = 0; j <= fOrder; j++) 
         fX(i,j) += weight * Pow(x, i+j);
   }
}

//__________________________________________________________________________________________________
//
// Include the passed point, with its error, to the computation
//
void CPolyFit::AddPointE(Double_t x, Double_t y, Double_t err)
{
   if (err == 0.0) {
      ::Error("AddPointE", "Cannot at a point without error");
      return;
   }
   
   Double_t weight = 1.0 / err / err;
   AddPointW(x, y, weight);
}

//__________________________________________________________________________________________________
//
// Add points from a histogram, from all bins in a given range
// and with the possibility to exclude a peak range
//
void CPolyFit::AddPointsFromHist
(TH1D *hist, Double_t fitMin, Double_t fitMax, Double_t peakMin, Double_t peakMax)
{
   Double_t x, y, e;
   for (Int_t i = 1; i < hist->GetNbinsX(); i++) {
      x = hist->GetBinCenter(i);
      y = hist->GetBinContent(i);
      e = hist->GetBinError(i);
      if (x < fitMin || (x >= peakMin && x <= peakMax)) continue;
      if (x > fitMax) break;
      AddPointE(x, y, e);
   }
}

//__________________________________________________________________________________________________
//
// Computes the coefficients
//
TArrayD CPolyFit::Compute()
{
   TMatrixD invX(TMatrixD::kInverted, fX);
   TMatrixD prod(invX, TMatrixD::kMult, fM);
   
   TArrayD out(fOrder + 1);
   
   Int_t i;
   for (i = 0; i <= fOrder; i++) out[i] = prod(i,0);
   
   return out;
}

//__________________________________________________________________________________________________   
//
// Raise to an integer power in an exact way
//
Double_t CPolyFit::Pow(Double_t x, Int_t n)
{
   if (n < 1) 
      return 1;
   else
      return x * Pow(x, n-1);
}

#endif
