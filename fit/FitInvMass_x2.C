enum EFit {
   kMixing,
   kLike,
   kFunction,
   kBinCounting,
   kBinCounting2G 
};

//----------------------------------------------------------------------------------
EFit ParseMode(const char *smode)
{
   TString str(smode);
   str.ToUpper();
   
   if (!str.CompareTo("FUNCTION"))
      return kFunction;
   else if (!str.CompareTo("LIKE"))
      return kLike;
   else if (!str.CompareTo("MIXING"))
      return kMixing;
   else if (!str.CompareTo("BINCOUNTING"))
     return kBinCounting;
   else if (!str.CompareTo("BINCOUNTING2G"))
     return kBinCounting;
   else
      return -1;
}

//----------------------------------------------------------------------------------
const char *fitName(EFit mode)
{
   switch (mode) {
      case kMixing  : return "Mixing";
      case kLike    : return "LikeSign";
      case kFunction: return "Function";
      case kBinCounting : return "BinCounting";
      case kBinCounting2G : return "2GBinCounting";
      default       : return "none";
   }
}

//----------------------------------------------------------------------------------
Int_t fitColor(EFit mode)
{
   switch (mode) {
      case kMixing  : return kRed;
      case kLike    : return kBlue;
      case kFunction: return kGreen+1;
      case kBinCounting: return kMagenta+1;
      case kBinCounting: return kOrange+6;
      default       : return kBlack;
   }
}

//----------------------------------------------------------------------------------
TH1D* SumSignal(TH1D *hpm, TH1D *hmp)
{
  Int_t i, nbins = hpm->GetNbinsX();
  Double_t y1, y2, e1, e2, sig, err;   
  TH1D *out = (TH1D*)hpm->Clone();  
  out->SetName("tmp");
  out->Clear();
  out->SetEntries(hpm->GetEntries() + hmp->GetEntries());   
  for (i = 1; i <= nbins; i++) {
    y1 = hpm->GetBinContent(i);
    y2 = hmp->GetBinContent(i);
    e1 = hpm->GetBinError(i);
    e2 = hmp->GetBinError(i);
    if (y1 < 1.0 || y2 < 1.0) {
      out->SetBinContent(i, 0.0);
      out->SetBinError(i, 0.0);
    } else {
      sig = (y1+y2);//not divided by 2 in this version *0.5;
      err = e2*e2 + e1*e1;
      out->SetBinContent(i, sig);
      out->SetBinError(i, TMath::Sqrt(err));
    }
  }   
  return out;
}
//----------------------------------------------------------------------------------
TH1D* BgLike(TH1D *hpp, TH1D *hmm)
{
   Int_t i, nbins = hpp->GetNbinsX();
   Double_t y1, y2, e1, e2, bg, err;
   
   TH1D *out = (TH1D*)hpp->Clone();
   
   out->SetName("tmp");
   out->Clear();
   out->SetEntries(hpp->GetEntries() + hmm->GetEntries());
   
   for (i = 1; i <= nbins; i++) {
      y1 = hpp->GetBinContent(i);
      y2 = hmm->GetBinContent(i);
      e1 = hpp->GetBinError(i);
      e2 = hmm->GetBinError(i);
      
      if (y1 < 1.0 || y2 < 1.0) {
	out->SetBinContent(i, 0.0);
	out->SetBinError(i, 0.0);
      } else {
	//multiply by two if signal is not divided by two 
	bg = 2.0 * TMath::Sqrt(y1*y2);
	err = y1*y1*e2*e2 + y2*y2*e1*e1;
	err /= y1*y2;
	out->SetBinContent(i, bg);
	out->SetBinError(i, TMath::Sqrt(err));
      }
   }
   
   return out;
}

//----------------------------------------------------------------------------------
TH1D* BgMixing(TH1D *hpm, TH1D *hmp)
{
   Int_t i, nbins = hpm->GetNbinsX();
   Double_t y1, y2, e1, e2, bg, err;   
   TH1D *out = (TH1D*)hpm->Clone();  
   out->SetName("tmp");
   out->Clear();
   out->SetEntries(hpm->GetEntries() + hmp->GetEntries());   
   for (i = 1; i <= nbins; i++) {
      y1 = hpm->GetBinContent(i);
      y2 = hmp->GetBinContent(i);
      e1 = hpm->GetBinError(i);
      e2 = hmp->GetBinError(i);
      if (y1 < 1.0 || y2 < 1.0) {
         out->SetBinContent(i, 0.0);
         out->SetBinError(i, 0.0);
      } else {
	bg = (y1+y2);//*0.5;	//multiply by 0.5 if signal is divided by two 
	err = e2*e2 + e1*e1;
	out->SetBinContent(i, bg);
	out->SetBinError(i, TMath::Sqrt(err));
      }
   }   
   return out;
}
//----------------------------------------------------------------------------------
Double_t IntegralNormalization(TH1D *hSig, TH1D *hBg, Int_t firstBin, Int_t lastBin)
{
//
//normalize by the integral over the full range
//
  Double_t norm;
  TH1D *htmpS = (TH1D*)hSig->Clone("tmpS");
  Double_t IntS=htmpS->Integral(firstBin, lastBin);
  TH1D *htmpB = (TH1D*)hBg->Clone("tmpB");
  Double_t IntB = htmpB->Integral(firstBin, lastBin);
  
  norm = IntS / IntB;
  Printf("Normalization: range %3.2f-%3.2f -> intS = %e, intB = %e, norm = %8.4f",hSig->GetXaxis()->GetBinUpEdge(firstBin),hSig->GetXaxis()->GetBinUpEdge(lastBin), IntS, IntB, norm);
  return norm;
}

//----------------------------------------------------------------------------------
Double_t BestNormalization
(TH1D *hSig, TH1D *hBg, Double_t normMin, Double_t normMax, Double_t normStep, Int_t firstBin, Int_t lastBin)
{
//
// search best factor which leaves all subtractions positive within error
//
   Double_t norm, bestNorm = normMin;
   Double_t min, bestMin  = 1E20;
   for (norm = normMax; norm >= normMin; norm -= normStep) {
      TH1D *htmp = (TH1D*)hSig->Clone("tmp");
      htmp->Add(hBg, -norm);
      min = 1E20;
      for (Int_t k = firstBin; k <= lastBin; k++) {
         Double_t y  = htmp->GetBinContent(k);
         if (htmp->GetBinError(k) < 1E-6) continue;
         y -= htmp->GetBinError(k);
         if (y < min) min = y;
      }
      if (min < 0.0)  
         continue;
      else if (min < bestMin) {
         bestMin  = min;
         bestNorm = norm;
      }
   }
   Double_t IntS=htmp->Integral(firstBin, lastBin);
   Double_t IntB = hBg->Integral(firstBin, lastBin);
   Printf("Normalization: best -> intS = %e, intB = %e, norm = %8.4f", IntS, IntB, bestNorm);
   return bestNorm;
}
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void FitInvMass_x2
(
   const char *filein = "proj/tpc2s_centBin00.root",
   EFit        mode = kLike,
   const char *func = "BW+POLY2",
   Double_t    viewMin = 0.72,
   Double_t    viewMax = 1.10,
   Bool_t      likeNorm = kTRUE,
   Bool_t      pause = kFALSE,
   Bool_t      data = kTRUE,
   Int_t       PDG = 313,
   
   Double_t    nsigmaPeak = 7.0,
   Int_t       nrebin = 1,
   Int_t       startBin = 0,
   Int_t       stopBin = 15,
   Int_t       centBinID = 0,

   Bool_t      fixGamma  = kTRUE,
   Double_t    minGamma = 0.5,
   Double_t    maxGamma = 1.5,   
   Bool_t      fixSigma  = kTRUE,
   char *      outdirsuffix = "",
   Double_t    multSigma = 1.0,
   Double_t    minSigma = 0.0,
   Double_t    maxSigma = -999.0,
   Double_t    nBCpeak = 2.0,
   Double_t    nBCexclude = 6.0
 )
{
   // retrieve reference values in the database PDG
   TDatabasePDG *pdg     = TDatabasePDG::Instance();
   TParticlePDG *part    = pdg->GetParticle(PDG);
   Double_t      mass    = part->Mass();
   Double_t      gamma   = part->Width();
   Double_t      peakMin = mass - nsigmaPeak * gamma / 2.35;
   Double_t      peakMax = mass + nsigmaPeak * gamma / 2.35;

   //check bin counting settings
   
   if ((mode==kBinCounting) && (nBCpeak>=nBCexclude)) { 
     //ensure that peak range is always inside the edges of the exclusion range, i.e:
     // excludeMin < BCpeakRangeMin < BCpeakRangeMax < excludeMax
     nBCexclude = nBCpeak+1;
     Printf("BC WARNING::: n_peak > n_exclude ---> n_exclude = n_peak+1 = %i", nBCexclude);
   }
   
   Double_t    excludeMin = mass - nBCexclude * gamma / 2.35;
   Double_t    excludeMax = mass + nBCexclude * gamma / 2.35; 
   Double_t    BCpeakRangeMin = mass - nBCpeak * gamma / 2.35; 
   Double_t    BCpeakRangeMax = mass + nBCpeak * gamma / 2.35;
   
   if (mode==kBinCounting){
     if (viewMin>=excludeMin) {
       viewMin = excludeMin - gamma / 2.35;
       Printf("BC WARNING::: x_min > exclude_min ---> x_min = exclude_min - 1sigma = %5.3f", viewMin);
     }
     if (viewMax<=excludeMax) {
       viewMax = excludeMax + gamma / 2.35;
       Printf("BC WARNING::: x_max < exclude_max ---> x_max = exclude_max + 1sigma = %5.3f", viewMax);
     }
     
     if (TMath::Abs(viewMin-excludeMin) < 2*gamma) {
       viewMin= excludeMin-2*gamma;
       Printf("BC WARNING::: |x_min-exclude_min|<2Gamma ---> x_min = exclude_min - 2gamma = %5.3f", viewMin);
     }
     if (TMath::Abs(viewMax-excludeMax) < 2*gamma) {
       viewMax = excludeMax + 2*gamma;
       Printf("BC WARNING::: |x_max-exclude_max|<2Gamma ---> x_max = exclude_max + 2gamma = %5.3f", viewMax);
     }
   }
   
   TString macroDir = "$HOME/alice/macro/kstar/fit";
  // load and compile fit function class
   gROOT->LoadMacro(Form("%s/myFitFcn.C+",macroDir.Data()));
   gROOT->LoadMacro(Form("%s/myFitResult.C+",macroDir.Data()));
   gROOT->LoadMacro(Form("%s/FitHistogram.C+",macroDir.Data()));
   gROOT->LoadMacro(Form("%s/BinCounting.C+",macroDir.Data()));

   myFitFcn *fit  = new myFitFcn(viewMin, viewMax, func);
   Int_t     nsig = fit->GetSignalNumPar();
   Int_t     nbg  = fit->GetBgNumPar();

   // open file and get bins
   TFile  *fin   = TFile::Open(filein);
   TAxis  *bins  = (TAxis*)fin->Get("ptbins");
   Int_t   nbins = bins->GetNbins();
   
   //plots for significance
   TString signif_name(filein);
   if (signif_name.Contains("proj/")) signif_name.ReplaceAll("proj/","");
   signif_name.Prepend("significance_");
   TFile *fSignificance =new TFile(signif_name.Data(),"RECREATE");
   TH1D * hSoverB=new TH1D("hSoverB","hSoverB",  nbins, bins->GetXbins()->GetArray());
   TH1D * hSoverSB=new TH1D("hSoverSB","hSoverSB",  nbins, bins->GetXbins()->GetArray());
   TH1D * hSignif=new TH1D("hSignif","hSignif",  nbins, bins->GetXbins()->GetArray());
   
   // get lists
   TList *lUnlikePM = (TList*)fin->Get(Form("%s_UnlikePM", (data ? "Data" : "MC")));
   TList *lUnlikeMP = (TList*)fin->Get(Form("%s_UnlikeMP", (data ? "Data" : "MC")));
   TList *lLikePP = (TList*)fin->Get(Form("%s_LikePP", (data ? "Data" : "MC")));
   TList *lLikeMM = (TList*)fin->Get(Form("%s_LikeMM", (data ? "Data" : "MC")));
   TList *lMixingPM = (TList*)fin->Get(Form("%s_MixingPM", (data ? "Data" : "MC")));
   TList *lMixingMP = (TList*)fin->Get(Form("%s_MixingMP", (data ? "Data" : "MC")));
   TList *lTrues  = 0x0;
   if (!data) {
     lTrues = (TList*)fin->Get("MC_Trues");
   }
   // get resolution
   TString sFunc(func);
   TH1D *hRes = (TH1D*)fin->Get("ResGauss");
   TF1  *fRes = (TF1*) fin->Get("funcResGauss");
   if (sFunc.Contains("VOIGT") && !hRes) {
     ::Error("Cannot process: required voigtian but resolution is not present");
     return;
   }
   
   
   // using the configuration object, create the histograms
   // which will contain all the raw counts and the corrections
   TH1D  *hraw   = new TH1D("raw"  , "" , nbins, bins->GetXbins()->GetArray());
   TH1D  *hbg    = new TH1D("bg"   , "" , nbins, bins->GetXbins()->GetArray());
   TH1D  *htrue  = new TH1D("true" , "" , nbins, bins->GetXbins()->GetArray());
   TH1D  *hmass  = new TH1D("mass" , "" , nbins, bins->GetXbins()->GetArray());
   TH1D  *hsigma = new TH1D("sigma", "" , nbins, bins->GetXbins()->GetArray());
   TH1D  *hgamma = new TH1D("gamma", "" , nbins, bins->GetXbins()->GetArray());
   TH1D  *hchi2  = new TH1D("chi2" , "" , nbins, bins->GetXbins()->GetArray());
   TH1D  *hhisto   = new TH1D("histo"  , "" , nbins, bins->GetXbins()->GetArray());
   TH1D  *hfunc    = new TH1D("func"   , "" , nbins, bins->GetXbins()->GetArray());
   
   // prepare output file
   TString foutName(filein);
   foutName.ReplaceAll("proj/", Form("fit%s/", outdirsuffix));
   //  foutName.ReplaceAll(".root", Form("_%s_%s_fit%.3f-%.3f_norm%.3f-%.3f", fitName(mode), func, viewMin, viewMax, normMin, normMax));
   foutName.ReplaceAll(".root", Form("_%s_%s_fit%.3f-%.3f", fitName(mode), func, viewMin, viewMax));
   
   if (!fixGamma){
     if (maxGamma<minGamma){
       minGamma=0.5; maxGamma=1.5; }
     Printf(":::: minGamma= %6.4f, maxGamma = %6.4f - Width constrained in: 0.5*PDG - 1.5*PDG",minGamma, maxGamma); 
     foutName.Append(Form("_width%.2f-%.2f", minGamma, maxGamma));
   } else { Printf(":::: Width fixed to PDG value "); }
   
   if (sFunc.Contains("VOIGT")){
     if (multSigma < 0.9999 || multSigma > 1.0001) {
       foutName.Append(Form("_msigma%.3f", multSigma));
     } else {
       if (maxSigma<minSigma){
	 Printf("minSigma= %6.4f, maxSigma = %6.4f Setting them to default values of 0.5 and 1.5 respectively. ",minSigma, maxSigma);
	 minSigma=0.5; 
	 maxSigma=1.5;
       } 
       foutName.Append(Form("lsigma%.3f-%.3f", minSigma, maxSigma));
     } 
   }
   
   gSystem->Exec(Form("mkdir -p %s", foutName.Data()));
   //gSystem->Exec(Form("rm %s/*", foutName.Data()));
   TString fitout = Form(Form("%s/fit_bin%i-%i.root", foutName.Data(), startBin, stopBin));
   TFile *fout = TFile::Open(fitout.Data(), "RECREATE");

   //create tree to save fit result
   Double_t fitSettings[4]={-1.0,-1.0, 0.0, 0.0}; //cent bin id, pt bin id, fit range min, fit range max
   Double_t fitParams[11], normfactorCopy;
   Double_t ptinf=0.0, ptsup=0.0, SoverB = 0.0, significance = 0.0;
   Int_t  ptBinID = -1;
   Double_t nhisto=0.0, nhistoerr=0.0, nfunc=0.0, nfuncerr=0.0 ;
   TTree *tree=new TTree("tree","fit parameters tree");
   //fit settings
   tree->Branch("pt_inf", &ptinf, "pt_inf/D");
   tree->Branch("pt_sup", &ptsup, "pt_sup/D");
   tree->Branch("centBin",&centBinID,"centBin/I");
   tree->Branch("ptBin",&ptBinID,"ptBin/I");
   tree->Branch("fitrange_inf", &viewMin, "fitrange_inf/D");
   tree->Branch("fitrange_sup", &viewMax, "fitrange_sup/D");
   //fit result
   tree->Branch("signalMass",&fitParams[0],"signalMass/D");
   tree->Branch("signalMassErr",&fitParams[1],"signalMassErr/D");
   tree->Branch("signalWidth",&fitParams[2],"signalWidth/D");
   tree->Branch("signalWidthErr",&fitParams[3],"signalWidthErr/D");
   tree->Branch("nSignal",&fitParams[4],"nSignal/D");
   tree->Branch("nSignalErr",&fitParams[5],"nSignalErr/D");
   tree->Branch("nBack",&fitParams[6],"nBack/D");
   tree->Branch("nBackErr",&fitParams[7],"nBackErr/D");
   tree->Branch("chi2",&fitParams[8],"chi2/D");  
   tree->Branch("signalRes",&fitParams[9],"signalRes/D");
   tree->Branch("signalResErr",&fitParams[10],"signalResErr/D");
   tree->Branch("SoverB",&SoverB,"SoverB/D");
   tree->Branch("significance",&significance,"significance/D");
   //bg normalisation info
   tree->Branch("norm_factor", &normfactorCopy, "norm_factor/D");
   //bin counting info
   tree->Branch("nHisto",&nhisto,"nSignal/D");
   tree->Branch("nHistoErr",&nhistoerr,"nSignalErr/D");
   tree->Branch("nFunc",&nfunc,"nFunc/D");
   tree->Branch("nFuncErr",&nfuncerr,"nFuncErr/D");
   

   // prepare global canvas
   gStyle->SetOptStat("");
   TCanvas *cTmp  = new TCanvas("cTmp", "TMP", 100, 400, 640, 480);
   TCanvas *cFit  = new TCanvas("cFit", "FIT", 0, 0, 700, 650);
   TPad    *pAll  = new TPad("pAll" , "", 0.001, 0.501, 0.501, 0.999);
   TPad    *pNums = new TPad("pNums", "", 0.501, 0.501, 0.999, 0.999);
   TPad    *pFit  = new TPad("pMix" , "", 0.001, 0.001, 0.999, 0.499);
   
   TCanvas * cSummary=new TCanvas("cSummary","cSummary", 800,750);
   cSummary->Divide(5,3);   
   TCanvas * cSummarySub=new TCanvas("cSummarySub","cSummarySub", 800,750);
   cSummarySub->Divide(5,3);   
   
   pAll->SetFillColor(kWhite);
   pNums->SetFillColor(kWhite);//kYellow - 10
   pFit ->SetFillColor(kWhite);//fitColor(mode) - 10
   cFit->cd(); pAll ->Draw();
   cFit->cd(); pNums->Draw();
   cFit->cd(); pFit ->Draw();
   
   // prepare output list for results
   TList *out = new TList;
   out->SetName(Form("Results_%d", fitName(mode)));
   
   // setup fitter options
   gSystem->Load("libMathCore.so");
   ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(1000);
   ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000);
   
   // fit all bins
   Int_t    nsig, nbg;
   Double_t res, dpt, pt;
   char     dummy[100];
   if (stopBin<startBin) stopBin=nbins+1;
   for (Int_t ibin = startBin; ibin < stopBin; ibin++) {
      
     // get bin width
     dpt = bins->GetBinWidth(ibin + 1);
     pt  = bins->GetBinCenter(ibin + 1);
     Printf("===================================\n===================================\n pt = %3.2f  (bin %i, dpt = %3.2f)\n===================================\n===================================", pt, ibin, dpt);
     //set pt in fit function
     fit->SetPt(pt);
     // get resolution
     if (sFunc.Contains("VOIGT")) {
       res = multSigma * (Double_t)hRes->GetBinContent(ibin+1);
       //res = fRes->Eval(pt);	 
       //Printf("multSigma= %.6f binContent = %.6f",multSigma, hRes->GetBinContent(ibin+1));
       //Printf("Looking at bin %i with pT = %.6f, resolution = %.6f",ibin+1, pt, res);
     }
     
      // get projections
      TH1D *hSignalPM = (TH1D*)lUnlikePM->At(ibin);
      TH1D *hSignalMP = (TH1D*)lUnlikeMP->At(ibin);
      TH1D *hLikePP = (TH1D*)lLikePP->At(ibin);
      TH1D *hLikeMM = (TH1D*)lLikeMM->At(ibin);
      TH1D *hMixingPM = (TH1D*)lMixingPM->At(ibin);
      TH1D *hMixingMP = (TH1D*)lMixingMP->At(ibin);
      
      // combine PM and MP signal and mixing bg, compute like sign bg
      TH1D *hSignal   = SumSignal(hSignalPM, hSignalMP);
      TH1D *hLike   = BgLike(hLikePP, hLikeMM);
      TH1D *hMixing   = BgMixing(hMixingPM, hMixingMP);
      
      
      // compute subtraction
      TH1D *hSub = (TH1D*)hSignal->Clone(Form("%s_%s", hSignal->GetName(), fitName(mode)));
      TH1D *hBG  = 0x0;
      
      //normalization intervals
      Int_t    inorm1 = hSignal->GetXaxis()->FindBin(viewMin);
      Int_t    inorm2 = hSignal->GetXaxis()->FindBin(viewMax);
      
      Int_t    inorm3 = hSignal->GetXaxis()->FindBin(1.3);//phi 1.04
      Int_t    inorm4 = hSignal->GetXaxis()->FindBin(1.5);//kstar
      Double_t norm, normBest;
      // normalize backgrounds
      if (mode == kMixing) {
	norm = IntegralNormalization(hSignal, hMixing, inorm3, inorm4); 
	normBest = BestNormalization(hSignal, hMixing, 0.010, 0.500, 0.001, inorm1, inorm2);
	hBG = (TH1D*)hMixing->Clone(Form("%s_bg", hSub->GetName()));
	hBG->Scale(norm);
      } else if ((mode == kLike) || (mode == kBinCounting) || (mode == kBinCounting2G)) {
	norm = IntegralNormalization(hSignal, hLike, inorm3, inorm4); 
	normBest = BestNormalization(hSignal, hLike, 0.010, 0.500, 0.001, inorm1, inorm2);
	hBG = (TH1D*)hLike->Clone(Form("%s_bg", hSub->GetName()));
	// if (likeNorm) 
	//   hBG->Scale(normBest);
	// else hBG->Scale(0.5);//add if you want (K*+antiK*)/2
      }      
      if (hBG) hSub->Add(hBG, -1.0);
      normfactorCopy = normBest;

      //debug 20/09/2013
      if (ibin+1<=stopBin) {
	cSummary->cd(ibin+1);
	hSignal->SetMarkerColor(kAzure-5);
	hSignal->Draw();
	hSignalPM->Draw("same");
	hSignalMP->Draw("same");
	if (mode == kLike){
	  hLikePP->Draw("same");
	  hLikeMM->Draw("same");
	} else {
	  hMixingPM->Draw("same");
	  hMixingMP->Draw("same");	
	}
	hBG->Draw("same");
	cSummarySub->cd(ibin+1);
	hSub->Draw();
      }
      //significance   
      /*
	if (mode == kLike) { 
	//	Float_t gammaRealistic=0.006; //only for phi
	Int_t    xMinSoverB = hSignal->GetXaxis()->FindBin(mass - 3. * gamma / 2.35);
	Int_t    xMaxSoverB = hSignal->GetXaxis()->FindBin(mass + 3. * gamma / 2.35);
	Double_t    intSBSoverB= hSignal->Integral(xMinSoverB, xMaxSoverB)*1.0;
	Double_t    intBGSoverB= hBG->Integral(xMinSoverB,xMaxSoverB)*1.0;
	SoverB = (intSBSoverB - intBGSoverB) / intBGSoverB;
	Double_t SoverSB = (intSBSoverB - intBGSoverB) / intSBSoverB;
	significance =  (intSBSoverB - intBGSoverB) / TMath::Sqrt(intSBSoverB);
	Printf("Significance pt %4.2f : intSB = %f  intBG = %f  S/B = %8.4f  S/(S+B) = %8.4f  S/sqrt(S+B) = %8.4f",pt, intSBSoverB, intBGSoverB, SoverB, SoverSB, significance);
	hSoverB->SetBinContent(ibin+1,SoverB);
	hSoverSB->SetBinContent(ibin+1,SoverSB);
	hSignif->SetBinContent(ibin+1,significance);
	}
      */

      // divide by bin width
      hSub->Scale(1.0/hSub->GetBinWidth(2));
      //Printf("================= > res = %.6f",res);
      hSub->Rebin(nrebin);
      
      myFitResult *result;     
      if (mode == kBinCounting) {      //BIN COUNTING
	result = BinCounting(hSub, fit, viewMin, viewMax, excludeMin, excludeMax, BCpeakRangeMin, BCpeakRangeMax, peakMin, peakMax, pt);
	result->Raw[0] /= BinCountTailCorrection(BCpeakRangeMin, BCpeakRangeMax);
      } else {         
	if (mode == kBinCounting2G) {
	  result = FitHistogram(hSub, fit, viewMin, viewMax, peakMin, peakMax, mass, gamma, res, fixGamma, minGamma, maxGamma, fixSigma, minSigma, maxSigma, pt, "2gamma");
	} else {//FIT FUNCTION 
	  result = FitHistogram(hSub, fit, viewMin, viewMax, peakMin, peakMax, mass, gamma, res, fixGamma, minGamma, maxGamma, fixSigma, minSigma, maxSigma, pt);
	}
      }
      // name fit result
      result->SetName(Form("result_%s_%02d", fitName(mode), ibin));
      Printf("°°°°°°°°°°°°°°°°°°°°°°°°°°°°°\n°°°°°°°°°°°°°°°°°°°°°°°°°°°°°");
      result->Print("MGSRIBH");
      // store values in outputs
      // setBin(hraw  , ibin + 1, result->Raw  [0] / dpt, result->Raw  [1] / dpt);
      // setBin(hbg   , ibin + 1, result->BgInt[0] / dpt, result->BgInt[1] / dpt);
      setBin(hraw  , ibin + 1, result->Raw  [0], result->Raw  [1]);
      setBin(hbg   , ibin + 1, result->BgInt[0], result->BgInt[1]);
      setBin(hmass , ibin + 1, result->Mass [0], result->Mass [1]);
      setBin(hgamma, ibin + 1, result->Gamma[0], result->Gamma[1]);
      setBin(hsigma, ibin + 1, result->Sigma[0], result->Sigma[1]);
      setBin(hchi2 , ibin + 1, result->Chi2 [0] / result->Chi2[1], 0.0);
      setBin(hhisto  , ibin + 1, result->HistInt  [0], result->HistInt[1]);
      setBin(hfunc   , ibin + 1, result->FcnInt[0], result->FcnInt[1]);
      
      fitParams[0]=result->Mass[0];   fitParams[1]=result->Mass[1];
      fitParams[2]=result->Gamma[0];   fitParams[3]=result->Gamma[1];
      // fitParams[4]=result->Raw[0]/dpt;    fitParams[5]=result->Raw[1]/dpt;
      // fitParams[6]=result->BgInt[0]/dpt;   fitParams[7]=result->BgInt[1]/dpt;
      fitParams[4]=result->Raw[0];    fitParams[5]=result->Raw[1];
      fitParams[6]=result->BgInt[0];   fitParams[7]=result->BgInt[1];
      fitParams[8] = result->Chi2 [0]/result->Chi2[1];
      fitParams[9]=result->Sigma[0];   fitParams[10]=result->Sigma[1];
      nhisto = result->HistInt  [0]; 
      nhistoerr = result->HistInt[1];
      nfunc = result->FcnInt[0];  
      nfuncerr=result->FcnInt[1];
      ptinf = bins->GetBinLowEdge(ibin+1);
      ptsup = bins->GetBinUpEdge(ibin+1);
      Printf("*************\n ibin %i\n*************", ibin);
      ptBinID = (Int_t) ibin;
      
      // save result into the list
      out->Add(result);
      
      //variables to estimate S/B for kstar
      Int_t    xMinSoverB = hSignal->GetXaxis()->FindBin(mass - 3. * gamma / 2.35);
      Int_t    xMaxSoverB = hSignal->GetXaxis()->FindBin(mass + 3. * gamma / 2.35);
      Double_t intSBSoverB= hSignal->Integral(xMinSoverB, xMaxSoverB)*1.0;
      Double_t intSraw=result->Raw[0];
      SoverB = intSraw / (intSBSoverB-intSraw);
      Double_t SoverSB = intSraw / intSBSoverB;
      significance = intSraw / TMath::Sqrt(intSBSoverB);
      Printf("Significance wrt ResBg (pt %4.2f): intSB = %f  intSraw = %f  S/B = %8.4f  S/(S+B) = %8.4f  S/sqrt(S+B) = %8.4f",pt, intSBSoverB, intSraw, SoverB, SoverSB, significance);
      hSoverB->SetBinContent(ibin+1,SoverB);
      hSoverSB->SetBinContent(ibin+1,SoverSB);
      hSignif->SetBinContent(ibin+1,significance);
      
      //fill the tree
      tree->Fill();
      
      // draw
      hSignal->GetXaxis()->SetRangeUser(0.60, 1.30);
      //hSignal->GetXaxis()->SetRangeUser(viewMin, viewMax);
      hSignal->GetXaxis()->SetTitle("#it{M}_{K#pi} (GeV/#it{c^{2}})");
      hSignal->GetYaxis()->SetTitle("Counts / (0.01 GeV/#it{c^{2}})");
      hSignal->SetMarkerStyle(20);
      hSignal->SetMarkerColor(kBlack);
      hSignal->SetLineColor(kBlack);
      if (hBG) hBG->SetLineColor(fitColor(mode));
      pAll->cd();
      hSignal->SetMinimum(0.0);
      hSignal->Draw();
      if (hBG) hBG->Draw("histsame");
      pFit->cd();
      hSub->SetTitle(fitName(mode));
      hSub->GetXaxis()->SetRangeUser(0.60, 1.30);
      //hSub->GetXaxis()->SetRangeUser(viewMin, viewMax);
      hSub->GetXaxis()->SetTitle("#it{M}_{K#pi} (GeV/#it{c^{2}})");
      hSub->GetYaxis()->SetTitle("dN/d#it{M}_{K#pi}");
      hSub->SetMarkerStyle(20);
      hSub->SetMarkerColor(kBlack);
      hSub->SetLineColor(kBlack);
      hSub->Draw();
      TPaveText *pavept = new TPaveText(0.12, 0.67, 0.40, 0.85,"NDC");
      pavept->SetLineWidth(0);
      pavept->SetLineColor(kWhite);
      pavept->SetFillStyle(0);
      pavept->AddText(Form("centrality %i-%i%%",20*centBinID,20*(centBinID+1)));
      pavept->AddText(Form("%3.2f<p_{T}<%3.2f GeV/#it{c}",bins->GetBinLowEdge(ibin+1),bins->GetBinUpEdge(ibin+1)));
      pavept->Draw("same");
      
      if (mode != kBinCounting) {
	result->fSum->SetLineColor(fitColor(mode)+2);
	result->fSum->SetNpx(100000);
	result->fSum->SetLineWidth(2.0);
	result->fSum->Draw("same");
	//result->fBg->SetLineColor(kMagenta);
	result->fBg ->SetNpx(100000);
	result->fBg ->SetLineWidth(1.0);
	result->fBg ->SetLineStyle(2);
	result->fBg ->Draw("same");
      }
      
      // draw a text with results
      pNums->cd();
      pNums->Clear();
      TPaveText *pave = result->Pave(0.05, 0.05, 0.95, 0.95, "RMGIBHC");
      pave->SetBorderSize(0);
      pave->SetTextFont(42);
      pave->SetTextAlign(21);
      pave->SetTextSize(0.06);
      pave->SetFillColor(pNums->GetFillColor());
      pave->SetLineColor(kBlack);
      pave->Draw();
      

      // store figure
      cFit->Update();
      cFit->SaveAs(Form("%s/%s_%s_bin_%02d_all.png", foutName.Data(), fitName(mode), func, ibin));
      
      // save functions and histograms
      cTmp->cd();
      cTmp->Clear();
      hSignal->Draw("PE");
      pavept->Draw("same");
      if (hBG)  hBG->Draw("histsame");
      cTmp->SaveAs(Form("%s/%s_%s_bin_%02d_distrib.png", foutName.Data(), fitName(mode), func, ibin));
      cTmp->SaveAs(Form("%s/%s_%s_bin_%02d_distrib.C", foutName.Data(), fitName(mode), func, ibin));
      cTmp->Clear();
      hSub->Draw();
      if (mode != kBinCounting) result->fSum->Draw("same");
      result->fBg->Draw("same");
      pavept->Draw("same");
      cTmp->SaveAs(Form("%s/%s_%s_bin_%02d_sub.png", foutName.Data(), fitName(mode), func, ibin));
      cTmp->SaveAs(Form("%s/%s_%s_bin_%02d_sub.C", foutName.Data(), fitName(mode), func, ibin));
      
      // pause
      if (pause) {
        cout << "Pause (q=quit, c=continue): " << endl;
	cin >> dummy;
	if (!strcmp(dummy, "q")) return;
	if (!strcmp(dummy, "c")) pause = kFALSE;
      }
   }

   // save histograms
   fout->cd();
   hraw  ->Write();
   hbg   ->Write();
   htrue ->Write();
   hmass ->Write();
   hgamma->Write();
   hsigma->Write();
   hchi2 ->Write();
   tree  ->Write();
   fout->Close();

   fSignificance->cd();
   hSoverB->Write();
   hSoverSB->Write();
   hSignif->Write();
   fSignificance->Close();
}

//----------------------------------------------------------------------------------
void setBin(TH1D *hist, Int_t bin, Double_t value, Double_t error)
{
   if (!hist) return;
   
   hist->SetBinContent(bin, value);
   hist->SetBinError  (bin, error);
}
//----------------------------------------------------------------------------------
Float_t BinCountTailCorrection(Double_t BCpeakRangeMin = 0.7, Double_t BCpeakRangeMax=1.2, Bool_t display = kFALSE)
{
  //correction for the tail left outside the nominal peak range in bin counting 
  TF1 *bw = new TF1("bw","TMath::BreitWigner(x,[0],[1])",0.0,2.0);
  bw->SetParameters(0.896,0.0487);
  Float_t intPeak = bw->Integral(BCpeakRangeMin,BCpeakRangeMax);
  Float_t intTot = bw->Integral(0.0,2.0);
  if (intTot<=0) return -1.0;
  Double_t corrfactor = intPeak/intTot;
  Printf("BinCountTailCorrection:::: intPeak(%5.3f-%5.3f) = %e, intTot(0.0-2.0) = %e \n BinCountTailCorrection:::: Correction factor = %6.4f", 
	 BCpeakRangeMin,BCpeakRangeMax, intPeak, intTot, intPeak/intTot);
  if (display) {
    bw->Draw();
    TLine * lmin = new TLine(BCpeakRangeMin, 0.0, BCpeakRangeMin, 12.0);
    TLine * lmax = new TLine(BCpeakRangeMax, 0.0, BCpeakRangeMax, 12.0);
    lmin->SetLineStyle(7);
    lmax->SetLineStyle(7);
    lmin->Draw("same");
    lmax->Draw("same");
  }
  return corrfactor;
}
