//
// Expected input: several THnSparseF with three axes:
//  [0] = invariant mass or resolution
//  [1] = transverse momentum
//  [2] = centrality
//
// This is converted into a collection of TH1D's 
// with invariant mass distributions for each bin
// in transverse momentum, and in a given range of centrality.
//
// The function builds a single file with all the produce histograms
// which come from two ROOT files, one for the data and another for 
// the corresponding MC.
//

#include <Riostream.h>

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TMath.h>
#include <TList.h>
#include <TError.h>
#include "TNamed.h"

class projectorInvMass_Centrality_Pt
{
public:
   
  projectorInvMass_Centrality_Pt() : fInput(0x0) { }
  
  void  SetPrefix(const char *prefix) {fPrefix = prefix;}
   
  //TH2 projectons
  TH2D* ProjectVsCentrality(Float_t cmin, Float_t cmax); //x=invMass y=pt
  TH2D* ProjectVsPt(Float_t pmin, Float_t pmax); //x=invMass y=centrality
  
  //TH1 projection
  TH1D* ProjectVsPtVsCentrality( Float_t pmin, Float_t pmax, Float_t cmin, Float_t cmax);
  void  MultiProjPtCent(Int_t npt, Double_t *pt, Int_t ncent, Double_t *cent,THnSparse *input, TList *out);
   
  //A.Pulvirenti
  TH1D* SingleProj(Float_t pmin, Float_t pmax);
  void  MultiProj(Int_t npt, Double_t *pt, THnSparse *input, TList *out);
   
private:

  TString     fPrefix;   // prefix used in the name
  Float_t     fCentMin;  // centrality range -- min
  Float_t     fCentMax;  // centrality range -- max
  Float_t     fPtMin;  // pt  range -- min
  Float_t     fPtMax;  // pt range -- max
  Float_t     fInvMassMin;  // inv mass range -- min
  Float_t     fInvMassMax;  // inv mass range -- max
  
  THnSparse       *fInput;    // processed histogram
   
};

//__________________________________________________________________________________________________

TH2D* projectorInvMass_Centrality_Pt::ProjectVsCentrality(Float_t cmin, Float_t cmax){
    // catch error
  if (!fInput) {
    ::Error("singleProj", "Null argument");
    return 0x0;
   }
 
  // compute bin
  Int_t   ic1 = fInput->GetAxis(2)->FindBin(cmin);
  Int_t   ic2 = fInput->GetAxis(2)->FindBin(cmax);
  Float_t c1  = fInput->GetAxis(2)->GetBinLowEdge(ic1);
  Float_t c2  = fInput->GetAxis(2)->GetBinUpEdge (ic2);
  if (c1 <cmin) ic1++;
  if (c2 >cmax) ic2--;
 
  ::Info("projectorInvMass_Centrality_Pt", "Centrality  bin: %d --> %d = %5.2f --> %5.2f", ic1, ic2, c1, c2);
  
  // compute projection
  TH2D *out = fInput->Projection(0, 1);
  out->SetName(Form("%s_%.2f-%.2f", fInput->GetName(), cmin, cmax));
  
  // finish
  return out;
  
};

//__________________________________________________________________________________________________

TH2D* projectorInvMass_Centrality_Pt::ProjectVsPt(Float_t pmin, Float_t pmax){
    
  // catch error
  if (!fInput) {
    ::Error("singleProj", "Null argument");
    return 0x0;
  }
  // compute bin
  Int_t   ip1 = fInput->GetAxis(1)->FindBin(pmin);
  Int_t   ip2 = fInput->GetAxis(1)->FindBin(pmax);
  Float_t p1  = fInput->GetAxis(1)->GetBinLowEdge(ip1);
  Float_t p2  = fInput->GetAxis(1)->GetBinUpEdge (ip2);
  if (p1 <pmin) ip1++;
  if (p2 >pmax) ip2--;
 
  ::Info("projectorInvMass_Centrality_Pt", "Pt  bin: %d --> %d = %5.2f --> %5.2f", ip1, ip2, p1, p2);
  
  // compute projection
  TH2D *out = fInput->Projection(0,2);
  out->SetName(Form("%s_%.2f-%.2f", fInput->GetName(), pmin, pmax));
  
  // finish
  return out;
  
};

//__________________________________________________________________________________________________
  
TH1D* projectorInvMass_Centrality_Pt::ProjectVsPtVsCentrality(Float_t pmin, Float_t pmax, Float_t cmin, Float_t cmax){
  
    // catch error
  if (!fInput) {
    ::Error("singleProj", "Null argument");
    return 0x0;
  }
  // compute bin
  Int_t   ip1 = fInput->GetAxis(1)->FindBin(pmin/*+0.0001*/);
  Int_t   ip2 = fInput->GetAxis(1)->FindBin(pmax/*-0.0001*/);
  Float_t p1  = fInput->GetAxis(1)->GetBinLowEdge(ip1);
  Float_t p2  = fInput->GetAxis(1)->GetBinUpEdge (ip2);
  if (p1 <pmin) ip1++;
  if (p2 >pmax) ip2--;
 
  // compute bin
  Int_t   ic1 = fInput->GetAxis(2)->FindBin(cmin/*+0.0001*/);
  Int_t   ic2 = fInput->GetAxis(2)->FindBin(cmax/*-0.0001*/);
  Float_t c1  = fInput->GetAxis(2)->GetBinLowEdge(ic1);
  Float_t c2  = fInput->GetAxis(2)->GetBinUpEdge (ic2);
  if (c1 <cmin) ic1++;
  if (c2 >cmax) ic2--;
  
  // ::Info("projectorInvMass_Centrality_Pt", "Pt  bin: %d --> %d = %5.2f --> %5.2f", ip1, ip2, p1, p2);
  // ::Info("projectorInvMass_Centrality_Pt", "Centrality  bin: %d --> %d = %5.2f --> %5.2f", ic1, ic2, c1, c2);

  // compute projection
  fInput->GetAxis(1)->SetRange(ip1,ip2);
  fInput->GetAxis(2)->SetRange(ic1,ic2);
  TH1D *out = fInput->Projection(0,"E0");
  out->SetName(Form("%s_pt%.2f-%.2f_cent%.2f-%.2f", fInput->GetName(), pmin, pmax,cmin, cmax));
  
  // finish
  return out;
  

};


//__________________________________________________________________________________________________
//
// Single projection from the passed THnSparseF, in a given range of pt and centrality.
// The output histogram has a name composed of the prefix and the bin number.
//
TH1D* projectorInvMass_Centrality_Pt::SingleProj(Float_t pmin, Float_t pmax)
{
   // catch error
   if (!fInput) {
      ::Error("singleProj", "Null argument");
      return 0x0;
   }
   
   // compute bin
   Int_t   ip1 = fInput->GetAxis(1)->FindBin(pmin);
   Int_t   ip2 = fInput->GetAxis(1)->FindBin(pmax);
   Float_t p1  = fInput->GetAxis(1)->GetBinLowEdge(ip1);
   Float_t p2  = fInput->GetAxis(1)->GetBinUpEdge (ip2);
   if (p1 < pmin) ip1++;
   if (p2 > pmax) ip2--;
   //p1  = fInput->GetAxis(1)->GetBinLowEdge(ip1);
   //p2  = fInput->GetAxis(1)->GetBinUpEdge(ip2);
   //::Info("ProjectSparse_IM_Pt_Cent", "Tr momentum bin: %d --> %d = %5.2f --> %5.2f", ip1, ip2, p1, p2);
   //::Info("ProjectSparse_IM_Pt_Cent", "Centrality  bin: %d --> %d = %5.2f --> %5.2f", ic1, ic2, c1, c2);
   
   // compute projection
   TH1D *out = fInput->Projection(0);//Form("%s_%.2f-%.2f", fInput->GetName(), p1, p2), ip1, ip2);
   
   // finish
   return out;
}

//__________________________________________________________________________________________________
//
// This function creates a full list of projections, by calling the above one
// several times, one per each PT bin defined in the passed array.
// All projections are added into an output list
//
void projectorInvMass_Centrality_Pt::MultiProj(Int_t npt, Double_t *pt, THnSparse *input, TList *out)
{
   // catch error
   if (!input) {
     ::Error("singleProj", "Null argument");
     return;
   } else {
     fInput = input;
   }
   
   // loop on all pt bins
   for (Int_t ip = 0; ip < npt; ip++) {
      TH1D *proj = SingleProj(pt[ip], pt[ip+1]);
      proj->SetName(Form("%s_bin%02d", fPrefix.Data(), ip));
      proj->SetTitle(Form("%.1f - %.1f", pt[ip], pt[ip+1]));
      out->Add(proj);
   }
   return;
};

//----------------------------------------------------------------------------------
void  projectorInvMass_Centrality_Pt::MultiProjPtCent(Int_t npt, Double_t *pt, Int_t ncent, Double_t *cent,THnSparse *input, TList *out){
   // catch error
   if (!input) {
     ::Error("singleProj", "Null argument");
     return;
   } else {
     fInput = input;
   }
   
   // loop on all centrality and pt bins
    for (Int_t ip = 0; ip < npt; ip++) {
      for (Int_t ic = 0; ic < ncent; ic++) {	
	TH1D *proj = ProjectVsPtVsCentrality(pt[ip], pt[ip+1],cent[ic], cent[ic+1]);
	proj->SetName(Form("%s_ptBin%02d_centBin%02d", fPrefix.Data(), ip,ic));
	proj->SetTitle(Form("pt %.1f - %.1f, cent %.1f - %.1f ", pt[ip], pt[ip+1], cent[ic], cent[ic+1] ));
	out->Add(proj);
      }
      
    }
    return;
    
    
};
