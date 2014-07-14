#include <iomanip>

enum EPart {
   kPhi,
   kKStar,
   kSigmaStar,
   kSigmaStarP,
   kSigmaStarM,
   kXiStar,
   kDelta,
   kLambda,
   kALambda,
   kAllLambda
};

void doFit(const char *fname = "spectra.root")
{
   gROOT->LoadMacro("Fit.C+");
   
   TFile *fin = TFile::Open("spectra.root");
   gROOT->LoadMacro("~/Documents/resonances2011/trunk/article_v2/Macros/HistogramsCurrent.C+");
   
   TH1D *h[4];
   
   h[0] = (TH1D*)fin->Get("Phi");
   h[1] = (TH1D*)fin->Get("KStarOld");
   h[2] = Phi(1,1,kFALSE);
   h[3] = KStarOld(1,1,kFALSE);
   
   Int_t minBin[4] = {1, 1, 0, 0};
   Int_t minBin[4] = {1, 1, 0, 0};

   Int_t    minBin, maxBin;
   Double_t mass = 0.0;
   TH1D    *histStat = 0x0;
   TH1D    *histSyst = 0x0;
   TFile   *input = 0x0;
   Double_t initY, initT, initN;
   
   initY  = 0.038;
   initT  = 0.280;
   initN  = 4.0;
   minBin = 2;
   maxBin = hrawL->GetNbinsX();
   mass   = 1.019455;

   TCanvas *c = new TCanvas("c", "FIT", 0, 0, 1600, 800);
   c->Divide(2, 2, 0.001, 0.001);
   
   c->cd(1);
   gPad->SetLogy();
   hrawL->SetLineColor(kGreen + 1);
   hrawL->SetMarkerColor(kGreen + 1);
   CFitOut out_statL = Fit(gPad, hrawL, kLevy, "EMI0", minBin, maxBin, mass, initY, initT, initN);
   out_statL.PrintValues();
   
   c->cd(2);
   gPad->SetLogy();
   hrawM->SetLineColor(kRed + 1);
   hrawM->SetMarkerColor(kRed + 1);
   CFitOut out_statM = Fit(gPad, hrawM, kLevy, "EMI0", minBin, maxBin, mass, initY, initT, initN);
   out_statM.PrintValues();
   
   c->cd(3);
   gPad->SetLogy();
   hrawF->SetLineColor(kBlue + 1);
   hrawF->SetMarkerColor(kBlue + 1);
   CFitOut out_statF = Fit(gPad, hrawF, kLevy, "EMI0", minBin, maxBin, mass, initY, initT, initN);
   out_statF.PrintValues();
   
   c->cd(4);
   gPad->SetLogy();
   hold->SetLineColor(kBlue + 1);
   hold->SetMarkerColor(kBlue + 1);
   CFitOut out_statO = Fit(gPad, hold, kLevy, "EMI0", minBin, maxBin, mass, initY, initT, initN);
   out_statO.PrintValues();
   
   c->SaveAs("fit_comparison.png");
}

/*
void Fit
(
   TH1D       *histogram,
   EFunction   fcn       = kExpPt,
   const char *optFit    = "QEMI0",
   Int_t       minBin    = 1,
   Int_t       maxBin    = 10,
   Double_t    mass      = 0.0,
   Double_t    initYield = 1.0,
   Double_t    initTemp  = 0.2,
   Double_t    initPower = 1.0
)
*/
