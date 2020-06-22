#include "/Users/fbellini/alice/macros/ResonAnT/fit/FitHistogram.C"
#include "/Users/fbellini/alice/macros/cosmetics/MakeUp.C"

enum EFit { kMixing = 0,
	    kLike = 1,
	    kFunction,
	    kBinCounting,
	    kEMBinCounting,
	    kNfits};

enum ENorm {  kNoNorm = 0,
	      kBestNorm,
	      kNorm106108,
	      kNorm1112,
	      kNorm115125};

enum EHistStyle {kSig = 0,
		 kLSBPP,
		 kLSBMM,
		 kLSB,
		 kLSBnorm,
		 kLSBsub,
		 kMEB,
		 kMEBnorm,
		 kMEBsub};

EFit    ParseMode(const char *smode);
void    setBin(TH1D *hist, Int_t bin, Double_t value, Double_t error);
TString fitName(EFit mode);
Int_t   fitColor(EFit mode);

TH1D * SumSignal(TH1D *hpm, TH1D *hmp);
TH1D * BgLike(TH1D *hpp, TH1D *hmm);
TH1D * BgMixing(TH1D *hpm, TH1D *hmp);
TH1D * GetNormalisedBg(TH1D * hSignal, TH1D * hBg, ENorm normMethod, Double_t * normFactors);

Double_t BestNormalization(TH1D *hSig, TH1D *hBg, Double_t normMin, Double_t normMax, Double_t normStep, Int_t firstBin, Int_t lastBin);
Double_t IntegralNormalization(TH1D *hSig, TH1D *hBg, Int_t firstBin, Int_t lastBin);


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void fitInvMassPhiXeXe(TString     filein = "sub_20180123_RsnOut.root",
		       EFit        mode = EFit::kMixing,
		       TString     func = "VOIGT+POLY1",
		       Double_t    fitMin = 0.994,
		       Double_t    fitMax = 1.050,
		       ENorm       normMethod = ENorm::kNorm106108,
		       Int_t       startBin = 3,
		       Int_t       stopBin = 10,
		       Int_t       centBinID = -1,
		       Bool_t      fixGamma  = kFALSE,
		       Bool_t      fixSigma  = kTRUE,
		       Double_t    minGamma = 0.25,
		       Double_t    maxGamma = 1.75,   
		       Double_t    minSigma = 1.0,
		       Double_t    maxSigma = 1.0,
		       Bool_t      scaleByBinWidth = 0,
		       Bool_t      pause = kFALSE,
		       Bool_t      data = kTRUE,
		       Int_t       PDG = 333,
		       Double_t    desiredIMbinWidth = 0.002) //in GeV/cˆ2
{

  // load and compile fit function class
  // TString macroDir = "$HOME/alice/macro/kstar/fit";
  // gROOT->LoadMacro(Form("%s/myFitFcn.C+",macroDir.Data()));
  // gROOT->LoadMacro(Form("%s/myFitResult.C+",macroDir.Data()));
  // gROOT->LoadMacro(Form("%s/FitHistogram.C+",macroDir.Data()));

  //Set style
  TGaxis::SetMaxDigits(3);
  gStyle->SetLineWidth(1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.07,"xyz");
  gStyle->SetTitleOffset(1.1,"y");
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetEndErrorSize(0); //sets in #of pixels the lenght of the tick at the end of the error bar
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);

  //histo make up settings
  Float_t Marker_Size = .8;
  Float_t Marker_Style[] = {    20,      24,      25,     21,     21,     25,       28,      28,     20};
  Color_t color[]        = {kBlack, kOrange, kCyan+2, kRed-1, kRed+1, kGray+2, kAzure-7, kBlue+1, kBlack};

  // retrieve reference values in the database PDG
  TDatabasePDG *pdg     = TDatabasePDG::Instance();
  TParticlePDG *part    = pdg->GetParticle(PDG);
  Double_t      mass    = part->Mass(); 
  Double_t      gamma   = part->Width();
  Double_t      massRes = 0.0018;
  Float_t       nsigmaPeak = 7.0;
  Double_t      peakMin = mass - nsigmaPeak * gamma / 2.35;
  Double_t      peakMax = mass + nsigmaPeak * gamma / 2.35;
  
  Float_t       nsigma4Signif = 3.0;
  Double_t      peakMinSignif = mass - nsigma4Signif * gamma / 2.35;
  Double_t      peakMaxSignif = mass + nsigma4Signif * gamma / 2.35;
  
  myFitFcn *fit  = 0;
  func.ToUpper();
  if (func.Contains("VOIGT+POLY1")) fit = new myFitFcn(fitMin, fitMax);
  if (func.Contains("VOIGT+POLY2")) fit = new myFitFcn(fitMin, fitMax);
  if (func.Contains("BW+POLY1")) fit = new myFitFcn(fitMin, fitMax);
  if (func.Contains("BW+POLY2")) fit = new myFitFcn(fitMin, fitMax);
  if (func.Contains("BW+POLY3")) fit = new myFitFcn(fitMin, fitMax);

  // open file and get bins
  Printf(":::: Opening file %s", filein.Data());
  TFile  *fin   = TFile::Open(filein.Data(), "read"); 
  if (!fin) return;
  
  TAxis  *bins  = (TAxis*)fin->Get("ptbins");
  Int_t   nbins = bins->GetNbins();
  TAxis  *centbins  = (TAxis*)fin->Get("centbins");
  Int_t   nCentBins = centbins->GetNbins();
  
  // get lists
  TList *lUnlikePM = (TList*)fin->Get("hUnlikePM");
  TList *lLikePP   = (TList*)fin->Get("hLikePP");
  TList *lLikeMM   = (TList*)fin->Get("hLikeMM");
  TList *lMixingPM = (TList*)fin->Get("hMixingPM");  
  
  //-------------------------------
  //create output file and objects
  //-------------------------------
  TString folderName = Form("fit%s_%s", fitName(mode).Data(), func.Data());
  gSystem->Exec(Form("mkdir %s", folderName.Data()));

  TString foutName = Form("%s/fit%4.3f-%4.3f_W%s_S%s", folderName.Data(), fitMin, fitMax, (fixGamma? "fixed" : "free"), (fixSigma? "fixed" : "free"));
  gSystem->Exec(Form("mkdir -p %s", foutName.Data()));

  TString fitout = Form("%s/fit_c%i_bin%i-%i.root", foutName.Data(), centBinID, startBin, stopBin);
  TFile *fout = TFile::Open(fitout.Data(), "RECREATE");
  
  // create the histograms which will contain all the raw counts and the corrections
  TH1D  *hmass  = new TH1D("mass" , "mass" , nbins, bins->GetXbins()->GetArray());
  TH1D  *hgamma = new TH1D("gamma", "width" , nbins, bins->GetXbins()->GetArray());
  TH1D  *hsigma = new TH1D("sigma", "resolution" , nbins, bins->GetXbins()->GetArray());
  TH1D  *hchi2  = new TH1D("chi2" , "chi2/ndf" , nbins, bins->GetXbins()->GetArray());
  TH1D  *hraw   = new TH1D("raw"  , "raw yields" , nbins, bins->GetXbins()->GetArray());
  TH1D  *hbg    = new TH1D("bg"   , "background entries" , nbins, bins->GetXbins()->GetArray());
  TH1D  *hhisto = new TH1D("histo"  , "" , nbins, bins->GetXbins()->GetArray());
  TH1D  *hBgIntegral = new TH1D("func"   , "" , nbins, bins->GetXbins()->GetArray());
  TH1D  *hbincount = new TH1D("bincount"  , "" , nbins, bins->GetXbins()->GetArray());

  // create output and histograms for significance (before bg subtraction)
  TH1D * hSoverB  = new TH1D("hSoverB","hSoverB",  nbins, bins->GetXbins()->GetArray());
  TH1D * hSoverSB = new TH1D("hSoverSB","hSoverSB",  nbins, bins->GetXbins()->GetArray());
  TH1D * hSignif  = new TH1D("hSignif","hSignif",  nbins, bins->GetXbins()->GetArray());
  
  
  //create tree to save fit result // FIX ME
  Int_t    ptBinID = -1; //index of pt bin
  Double_t fitParams[11]; //output from the fit
  Double_t normFactors[2] = {0.0, 0.0}; //bg normalization factor //integral, best
  Double_t ptLow = 0.0, ptHigh = 0.0; //extremes of pt interval
  Double_t centLow = 0.0, centHigh = 0.0; //extremes of cent interval
  Double_t SoverB = 0.0, significance = 0.0; //significance
  Double_t nhisto = 0.0, nhistoerr = 0.0; 
  Double_t nRawBC = 0.0, nRawBCErr = 0.0;
  Double_t nrawtail = 0.0, nrawtailErr = 0.0;
  Double_t nBgFcn = 0.0, nBgFcnErr = 0.0 ; //integralof fcn in BC (fit) range if bin counting (fit)

  TTree *tree = new TTree("tree","fit parameters tree");
  //fit settings
  tree->Branch("centBin",&centBinID,"centBin/I");
  tree->Branch("centLow",&centLow,"centLow/I");
  tree->Branch("centHigh",&centHigh,"centHigh/I");
  tree->Branch("ptBin",&ptBinID,"ptBin/I");
  tree->Branch("ptLow", &ptLow, "ptLow/D");
  tree->Branch("ptHigh", &ptHigh, "ptHigh/D");
  tree->Branch("fitRangeLow", &fitMin, "fitRangeLow/D");
  tree->Branch("fitRangeHigh", &fitMax, "fitRangeHigh/D");
  //fit result
  tree->Branch("Mass",&fitParams[0],"Mass/D");
  tree->Branch("MassErr",&fitParams[1],"MassErr/D");
  tree->Branch("Width",&fitParams[2],"Width/D");
  tree->Branch("WidthErr",&fitParams[3],"WidthErr/D");
  tree->Branch("nSignal",&fitParams[4],"nSignal/D");
  tree->Branch("nSignalErr",&fitParams[5],"nSignalErr/D");
  tree->Branch("nBg",&fitParams[6],"nBg/D");
  tree->Branch("nBgErr",&fitParams[7],"nBgErr/D");
  tree->Branch("chi2oNDF",&fitParams[8],"chi2oNDF/D");  
  tree->Branch("VoigtSigma",&fitParams[9],"VoigtSigma/D");
  tree->Branch("VoigtSigmaErr",&fitParams[10],"VoigtSigmaErr/D");
  tree->Branch("SoB",&SoverB,"SoB/D");
  tree->Branch("significance",&significance,"significance/D");
  //bg normalisation info
  tree->Branch("normIntegral", &normFactors[0], "normIntegral/D");
  tree->Branch("normBest", &normFactors[1], "normBest/D");
  //bin counting info
  tree->Branch("nHisto",&nhisto,"nHisto/D");
  tree->Branch("nHistoErr",&nhistoerr,"/D");
  tree->Branch("nBgFcn",&nBgFcn,"nBgFcn/D");
  tree->Branch("nBgFcnErr",&nBgFcnErr,"nBgFcnErr/D");
  tree->Branch("nRawBC",&nRawBC,"nRawBC/D");
  tree->Branch("nRawBCErr",&nRawBCErr,"nRawBCErr/D");
  
  //prepare to display
  TCanvas *cTmp  = new TCanvas("cTmp", "TMP", 800, 600);
  TCanvas *cFit  = new TCanvas("cFit", "FIT", 0, 0, 700, 650);
  TPad    *pAll  = new TPad("pAll" , "", 0.001, 0.501, 0.501, 0.999);
  TPad    *pNums = new TPad("pNums", "", 0.501, 0.501, 0.999, 0.999);
  TPad    *pFit  = new TPad("pMix" , "", 0.001, 0.001, 0.999, 0.499);
  
  pAll->SetFillColor(kWhite);
  pNums->SetFillColor(kWhite);
  pFit->SetFillColor(kWhite);
  cFit->cd(); pAll ->Draw();
  cFit->cd(); pNums->Draw();
  cFit->cd(); pFit ->Draw();
   
  // setup fitter options
  gSystem->Load("libMathCore.so");
  ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(5000);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(5000);

  //----------------------------
  // FIT PARAMETERS
  //----------------------------
  // width
  if (fixGamma) Printf(":::: FITTING Width fixed to PDG value = %6.4f GeV/c", gamma);
  else {
    if (maxGamma < minGamma){
      minGamma = 0.5;
      maxGamma = 1.5;
    } else {    
      Printf(":::: FITTING Width: MinGamma = %6.4f, maxGamma = %6.4f", minGamma*gamma, maxGamma*gamma); 
    }
  }
    
  // resolution
  //-- Voigtian resolution
  Double_t voigtRes = 0.0018;
  if (func.Contains("VOIGT")){
    if (fixSigma) {
      Printf(":::: FITTING Resolution (Voigt) FIXED");
      voigtRes = 0.0018;
    } else {
      if (maxSigma < minSigma){
	minSigma = 0.5;
	maxSigma = 1.5;
      } else {
	Printf(":::: FITTING Resolution (Voigt): minSigma factor = %6.4f, maxSigma factor = %6.4f", minSigma, maxSigma);
      }
    }
  }

  // fit
  Double_t res, dpt, pt;
  char dummy[100];
  if (stopBin<startBin) stopBin = nbins+1;
  for (Int_t ibin = startBin; ibin < stopBin; ibin++) {
    
    // get pt bin 
    dpt = bins->GetBinWidth(ibin+1);
    ptLow  = bins->GetBinLowEdge(ibin+1);
    ptHigh  = bins->GetBinUpEdge(ibin+1);
    pt = bins->GetBinCenter(ibin+1);
    ptBinID = ibin;
    //fit->SetPt(pt);
    centLow  = centbins->GetBinLowEdge(ibin+1);
    centHigh  = centbins->GetBinUpEdge(ibin+1);
    Printf("======================================================================");
    Printf(":::: FITTING %3.2f < pt < %3.2f  (bin %i, dpt = %3.2f)", ptLow, ptHigh, ibin, dpt);
    Printf("======================================================================");

    
    //--------------------------
    // get input histogram
    //--------------------------
    TH1D *hSignal = (TH1D*)lUnlikePM->FindObject(Form("hUnlikePM_ptBin%02i_centBin%02i", ibin, centBinID));
    if (hSignal) Printf(">>>> Reading histo %s", hSignal->GetName());
    Beautify(hSignal, color[EHistStyle::kSig], 1, 1, Marker_Style[EHistStyle::kSig], Marker_Size);

    TH1D *hMixing = (TH1D*)lMixingPM->FindObject(Form("hMixingPM_ptBin%02i_centBin%02i", ibin, centBinID));
    if (hMixing) Printf(">>>> Reading histo %s", hMixing->GetName());
    Beautify(hMixing, color[EHistStyle::kMEB], 1, 1, Marker_Style[EHistStyle::kMEB], Marker_Size);
    
    TH1D *hLikePP = (TH1D*)lLikePP->FindObject(Form("hLikePP_ptBin%02i_centBin%02i", ibin, centBinID));
    if (hLikePP) Printf(">>> Reading histo %s", hLikePP->GetName());
    Beautify(hLikePP, color[EHistStyle::kLSBPP], 1, 1, Marker_Style[EHistStyle::kLSBPP], Marker_Size);

    TH1D *hLikeMM = (TH1D*)lLikeMM->FindObject(Form("hLikeMM_ptBin%02i_centBin%02i", ibin, centBinID)); 
    if (hLikeMM) Printf(">>>> Reading histo %s", hLikeMM->GetName());
    Beautify(hLikeMM, color[EHistStyle::kLSBMM], 1, 1, Marker_Style[EHistStyle::kLSBMM], Marker_Size);

    TH1D *hLike = BgLike(hLikePP, hLikeMM);
    hLike->SetName("hComputedLSB");
    if (hLike) Printf(">>>> Computed histo %s", hLike->GetName());
    hLike->SetTitle(hLikePP->GetTitle());
    Beautify(hLike, color[EHistStyle::kLSB], 1, 1, Marker_Style[EHistStyle::kLSB], Marker_Size);

    //-----------------------
    //Chose and normalise background
    //-----------------------
    TH1D *hBg  = 0x0;
    if (mode == EFit::kMixing) hBg = (TH1D *) GetNormalisedBg(hSignal, hMixing, normMethod, normFactors);
    if (mode == EFit::kLike) hBg = (TH1D *) GetNormalisedBg(hSignal, hLike, normMethod, normFactors);
    //-----------------------
    //Get invariant mass binning and rebin if requested
    //-----------------------
    Double_t invMassBinWidth = ((TAxis*) hSignal->GetXaxis())->GetBinWidth(1);
    Int_t rebinFactor = desiredIMbinWidth / invMassBinWidth;
    
    if (rebinFactor>0) {
      Printf(":::: REBINNING: Minv bin width: %4.3f --> %4.3f via rebin by factor %i", invMassBinWidth, desiredIMbinWidth, rebinFactor);
      hSignal->Rebin(rebinFactor);
      hBg->Rebin(rebinFactor);
    } else {
      Printf(":::: NO rebinning performed");
    }

    //scale by bin width
    if (scaleByBinWidth) {
      hSignal->Scale(1.0/invMassBinWidth);
      hBg->Scale(1.0/invMassBinWidth);
    }
    
    //-----------------------
    // compute subtraction
    //-----------------------
    TH1D *hSub = (TH1D*) hSignal->Clone(Form("%s_%s", hSignal->GetName(), fitName(mode).Data()));    
    if (hBg) hSub->Add(hBg, -1.0);
    else Printf(":::: Warning: histogram to be subtracted cannot be found");
    Beautify(hSub, color[EHistStyle::kSig], 1, 1, Marker_Style[EHistStyle::kSig], Marker_Size);

    //-----------------------
    // Run Fit
    //-----------------------
    TString fitopt = "EMR";
    myFitResult * result = FitHistogram(hSub, 
					fit, 
					fitMin, fitMax, 
					peakMin, peakMax, 
					mass, gamma, voigtRes, 
					fixGamma, minGamma, maxGamma, 
					fixSigma, minSigma, maxSigma, 
					pt, fitopt.Data());
    
    // name fit result
    result->SetName(Form("result_%s_%02d", fitName(mode).Data(), ibin));
    Printf("°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°");
    result->Print("MGSRIBHT");
    Printf("°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°");
    
    //store values in output tree
    fitParams[0]=result->Mass[0]; fitParams[1]=result->Mass[1];
    fitParams[2]=result->Gamma[0]; fitParams[3]=result->Gamma[1];
    fitParams[4]=result->Raw[0]/dpt; fitParams[5]=result->Raw[1]/dpt;
    fitParams[6]=result->BgInt[0]/dpt; fitParams[7]=result->BgInt[1]/dpt;
    fitParams[8]=result->Chi2 [0]/result->Chi2[1];
    fitParams[9]=result->Sigma[0]; fitParams[10]=result->Sigma[1];

    //integral of histo in fit(BC) range
    nhisto    = result->HistInt[0]/dpt; //divide already by dpT ==> dN/dpT
    nhistoerr = result->HistInt[1]/dpt;
    //integral of tot fit function in fit(BC) range
    nBgFcn     = result->FcnInt[0]/dpt; //divide already by dpT ==> dN/dpT
    nBgFcnErr  = result->FcnInt[1]/dpt;
    
    //raw yields resulting form bin counting
    nRawBC    = (result->HistInt[0]-result->BgInt[0]);  //divide already by dpT ==> dN/dpT
    nRawBCErr = TMath::Sqrt(result->HistInt[1]*result->HistInt[1]+result->BgInt[1]*result->BgInt[1]);
    nRawBC /= dpt;
    nRawBCErr /= dpt;
    
    setBin(hraw  , ibin + 1, result->Raw  [0]/dpt, result->Raw  [1]/dpt);//<21set
    setBin(hbg   , ibin + 1, result->BgInt[0]/dpt, result->BgInt[1]/dpt);//<21set
    setBin(hmass , ibin + 1, result->Mass [0], result->Mass [1]);
    setBin(hgamma, ibin + 1, result->Gamma[0], result->Gamma[1]);
    setBin(hsigma, ibin + 1, result->Sigma[0], result->Sigma[1]);
    setBin(hchi2 , ibin + 1, result->Chi2 [0] / result->Chi2[1], 0.0);  
    setBin(hhisto, ibin + 1, nhisto, nhistoerr);
    setBin(hBgIntegral,  ibin + 1, nBgFcn, nBgFcnErr);
    setBin(hbincount, ibin + 1, nRawBC, nRawBCErr);
    
    
    //variables to estimate S/Bg for kstar (use signal int. before bg subtraction)
    Int_t    xMinSoverB = hSignal->GetXaxis()->FindBin(mass - 3.0 * gamma / 2.35);
    Int_t    xMaxSoverB = hSignal->GetXaxis()->FindBin(mass + 3.0 * gamma / 2.35);
    Double_t intSigBg = hSignal->Integral(xMinSoverB, xMaxSoverB);
    Double_t intSraw = result->Raw[0];
    Double_t SoverSB = intSraw / intSigBg; // S/(S+B)
    SoverB = intSraw / (intSigBg-intSraw); // S/B
    significance = intSraw / TMath::Sqrt(intSigBg); // S/sqrt(S+B)
    Printf("fitInvMassPhiXeXe:::: Significance after bg subtraction pt %4.2f:\n intSB = %f  intSraw = %f  S/B = %8.4f  S/(S+B) = %8.4f  S/sqrt(S+B) = %8.4f",pt, intSigBg, intSraw, SoverB, SoverSB, significance);
    hSoverB->SetBinContent(ibin+1,SoverB);
    hSoverSB->SetBinContent(ibin+1,SoverSB);
    hSignif->SetBinContent(ibin+1,significance);
    
    //fill the tree
    tree->Fill();

    // draw 
    hSignal->GetXaxis()->SetRangeUser(0.90, 1.20);
    hSignal->GetXaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c^{2}})");
    hSignal->GetYaxis()->SetTitle("Counts / (0.001 GeV/#it{c^{2}})");

    pAll->cd();
    hSignal->SetMinimum(0.0);
    hSignal->Draw();
    hBg->Draw("same");

    pFit->cd();
    hSub->SetTitle(fitName(mode).Data());
    hSub->GetXaxis()->SetRangeUser(fitMin, fitMax);
    hSub->GetXaxis()->SetTitle("#it{M}_{K^{+}K^{-}} (GeV/#it{c^{2}})");
    hSub->GetYaxis()->SetTitle("d#it{N}/d#it{M}_{KK}");
    hSub->Draw();

    TPaveText *pavept = new TPaveText(0.65, 0.67, 0.95, 0.85,"NDC");
    pavept->SetLineWidth(0);
    pavept->SetLineColor(kWhite);
    pavept->SetFillStyle(0);
    pavept->AddText(Form("V0M %2.0f-%2.0f%%", centbins->GetBinLowEdge(centBinID+1), centbins->GetBinLowEdge(centBinID+2)));
    pavept->AddText(Form("%3.1f#leq #it{p}_{T} < %3.1f GeV/#it{c}", bins->GetBinLowEdge(ibin+1),bins->GetBinUpEdge(ibin+1)));
    pavept->Draw("same");
      
    TPaveText *pavept2 = new TPaveText(0.6, 0.75, 0.87, 0.89,"NDC");
    pavept2->SetLineWidth(0);
    pavept2->SetLineColor(kWhite);
    pavept2->SetFillStyle(0);
    pavept2->SetTextFont(42);
    pavept2->AddText(Form("V0M %2.0f-%2.0f%%", centbins->GetBinLowEdge(centBinID+1), centbins->GetBinLowEdge(centBinID+2)));
    pavept2->AddText(Form("%3.1f #leq #it{p}_{T} < %3.1f GeV/#it{c}",bins->GetBinLowEdge(ibin+1),bins->GetBinUpEdge(ibin+1)));
      
    result->fSum->SetLineColor(fitColor(mode));
    result->fSum->SetNpx(100000);
    result->fSum->SetLineWidth(3.0);
    result->fSum->Draw("same");
    result->fBg ->SetNpx(100000);
    result->fBg ->SetLineWidth(3.0);
    result->fBg ->SetLineColor(kBlue+1);
    result->fBg ->SetLineStyle(7);
    result->fBg ->Draw("same");
          
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
    cFit->SaveAs(Form("%s/c%2.0f%2.0f_pt%02i_%s_%s_all.png", foutName.Data(), centbins->GetBinLowEdge(centBinID+1), centbins->GetBinLowEdge(centBinID+2), ibin, fitName(mode).Data(), func.Data()));
      
    // save functions and histograms
    cTmp->cd();
    cTmp->Clear();
    hSignal->Draw("PE");
    if (hBg)  hBg->Draw("same");
    cTmp->SaveAs(Form("%s/c%02.0f%02.0f_pt%02i_%s_sub.png", foutName.Data(), centbins->GetBinLowEdge(centBinID+1), centbins->GetBinLowEdge(centBinID+2), ibin,  fitName(mode).Data()));
    
    cTmp->Clear();
    gStyle->SetOptTitle(0);
    hSub->Draw();      
    result->fSum->Draw("same");
    result->fBg->Draw("same");
    pavept2->Draw("same");

    TString funLeg = "Breit-Wigner peak fit";
    if (func.Contains("VOIGT"))  funLeg = "Voigtian peak fit";
    
    TLegend * tmpleg = new TLegend(0.56, 0.6, 0.88, 0.75);
    tmpleg->SetLineWidth(0);
    tmpleg->SetLineColor(kWhite);   
    tmpleg->SetFillColor(kWhite);
    tmpleg->AddEntry(hSub, "Data (stat. uncert.)","lp");
    tmpleg->AddEntry(result->fSum, funLeg.Data(),"lp");
    tmpleg->AddEntry(result->fBg, "Residual background","lp");
    tmpleg->SetTextSize(0.04);
    tmpleg->Draw("same");
    cTmp->SaveAs(Form("%s/c%02.0f%02.0f_pt%02i_%s_%s_fit.png", foutName.Data(), centbins->GetBinLowEdge(centBinID+1), centbins->GetBinLowEdge(centBinID+2), ibin,  fitName(mode).Data(), func.Data()));
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
  tree->Write();
  hraw  ->Write();
  hbg   ->Write();
  hmass ->Write();
  hgamma->Write();
  hsigma->Write();
  hchi2 ->Write();
  hhisto->Write();
  hBgIntegral->Write();
  hbincount->Write();
  hSoverB->Write();
  hSoverSB->Write();
  hSignif->Write();
 
}

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
  else if (!str.CompareTo("EMBINCOUNTING"))
    return kEMBinCounting;
  else
    return kNfits;
}

//----------------------------------------------------------------------------------
TString fitName(EFit mode)
{
  switch (mode) {
  case kMixing  : return "MEB";
  case kLike    : return "LSB";
  case kFunction: return "Funct";
  case kBinCounting : return "BinCount";
  case kEMBinCounting : return "EMBinCount";
  default       : return "none";
  }
}

//----------------------------------------------------------------------------------
Int_t fitColor(EFit mode)
{
  switch (mode) {
  case kMixing  : return kRed;
  case kLike    : return kRed;
  case kFunction: return kMagenta+2;
  case kBinCounting: return kBlue;
  case kEMBinCounting: return kBlue;
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
      sig = (y1+y2);//*0.5;
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
   
   out->SetName("LSB");
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
	bg =  2.0*TMath::Sqrt(y1*y2);//multiply by two if signal is not divided by two bg = 2.0 * TMath::Sqrt(y1*y2);
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
      bg = (y1+y2);//*0.5;
      err = e2*e2 + e1*e1;
      out->SetBinContent(i, bg);
      out->SetBinError(i, TMath::Sqrt(err));
    }
  }   
  return out;
}

//----------------------------------------------------------------------------------
TH1D * GetNormalisedBg(TH1D * hSignal, TH1D * hBg, ENorm normMethod, Double_t * normFactors)
{
  
  if (!hBg || !hSignal) {
    Printf(":::: GetNormalisedBg missing inputs. Doing nothing.");
    return 0x0;
  }
  
  //normalization intervals
  Int_t inormLow, inormHigh;
  switch (normMethod) {
  case ENorm::kNorm106108:
    inormLow = hSignal->GetXaxis()->FindBin(1.060);
    inormHigh = hSignal->GetXaxis()->FindBin(1.080);
    break;
  case ENorm::kNorm1112 :
    inormLow = hSignal->GetXaxis()->FindBin(1.1);
    inormHigh = hSignal->GetXaxis()->FindBin(1.2);
    break;
  case ENorm::kNorm115125:
    inormLow = hSignal->GetXaxis()->FindBin(1.15);
    inormHigh = hSignal->GetXaxis()->FindBin(1.25);
    break;
  case ENorm::kBestNorm:
    inormLow = hSignal->GetXaxis()->FindBin(0.99);
    inormHigh = hSignal->GetXaxis()->FindBin(1.2);
    break;
  default :
    inormLow = hSignal->GetXaxis()->FindBin(1.05);
    inormHigh = hSignal->GetXaxis()->FindBin(1.15);
    break;
  }
  
  // get normalization factors for backgrounds
  Double_t norm, normBest;
  norm = IntegralNormalization(hSignal, hBg, inormLow, inormHigh); 
  normBest = BestNormalization(hSignal, hBg,  0.005, 0.500, 0.001, inormLow, inormHigh);
  
  TH1D * hNorm = (TH1D*) hBg->Clone(Form("norm_%s", hBg->GetName()));
  hNorm->SetTitle(Form("norm %s", hBg->GetTitle()));
  
  if (normMethod == ENorm::kBestNorm) {
    hNorm->Scale(normBest);
    Printf(":::: Normalisation of BG by best factor %e", normBest);
  } else {
    if (normMethod>ENorm::kNoNorm) {
      hBg->Scale(norm);
      Printf("fitInvMassPhiXeXe::::Normalisation of BG by factor %e", norm);
    }
  }

  normFactors[0] = norm;
  normFactors[1] = normBest;
  
  return hNorm;

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
  Printf("fitInvMassPhiXeXe::::Normalization: range %3.2f-%3.2f -> intS = %e, intB = %e, norm = %8.4f",hSig->GetXaxis()->GetBinUpEdge(firstBin),hSig->GetXaxis()->GetBinUpEdge(lastBin), IntS, IntB, norm);
  return norm;
}

//----------------------------------------------------------------------------------
Double_t BestNormalization(TH1D *hSig, TH1D *hBg, Double_t normMin, Double_t normMax, Double_t normStep, Int_t firstBin, Int_t lastBin)
{
  //
  // search best factor which leaves all subtractions positive within error
  //
  Double_t norm, bestNorm = normMin;
  Double_t min, bestMin  = 1E20;
  for (norm = normMax; norm >= normMin; norm -= normStep) {
    TH1D * htmp = (TH1D*)hSig->Clone("tmp");
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
  Double_t IntS = hSig->Integral(firstBin, lastBin);
  Double_t IntB = hBg->Integral(firstBin, lastBin);
  Printf("fitInvMassPhiXeXe:::: Normalization: best -> intS = %e, intB = %e, norm = %8.4f", IntS, IntB, bestNorm);
  return bestNorm;
}


//----------------------------------------------------------------------------------
void setBin(TH1D *hist, Int_t bin, Double_t value, Double_t error)
{
  if (!hist) return;
   
  hist->SetBinContent(bin, value);
  hist->SetBinError  (bin, error);
}
