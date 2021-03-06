//
// Splits a TH2F using the myHistSplig_sparse_IM_PT_CENT class
//
#include "/Users/fbellini/alice/macros/ResonAnT/core/projectorTH3_InvMass_Centrality_Pt.C"
#include "/Users/fbellini/alice/macros/cosmetics/MakeUp.C"
#include "TFile.h"

void GetEfficiencyFromBinnedMinv(TH1F* trueMinv=NULL, TH1F* momMinv=NULL, Float_t* effAndErr=NULL);
void GetGausSigma(TH1F * htmp = NULL, Float_t * sigmaAndErr = 0, Float_t absRange = 0.01);
void GetSigmaRMS(TH1F * htmp = NULL, Float_t * sigmaAndErr = 0, Int_t absRange = 3);
Float_t GetTruncationCorrection(Int_t absRange = 3);
Double_t Voigt( Double_t *x, Double_t * par);
TF1 * GetVOIGT(Double_t fitMin, Double_t fitMax);

int projectMC(TString nameData = "RsnOut.root",
	       TString listName = "RsnOut_default_LowBdca",
	       TString cutLabel  = "default_LowBdca",
	       TString binning  = "final",
	       Bool_t doEffOnly = 0,
	       Color_t customColor = kBlack)
{
   // initial setup
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTextFont(42);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(1);
  gStyle->SetLabelSize(0.06,"xyz");
  gStyle->SetLabelOffset(0.005,"yx");
  gStyle->SetTitleSize(0.07,"xyz");
  gStyle->SetTitleOffset(1.1,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetEndErrorSize(0); //sets in #of pixels the lenght of the tick at the end of the error bar
  gStyle->SetTitleAlign(33);
  gStyle->SetTitleX(.91);
  gStyle->SetTitleY(.91);
  TGaxis::SetMaxDigits(2);
    
  //Color_t color[]={kRed+1, kOrange+1, kSpring-5, kBlue+1, kBlack};
  //Color_t marker[]={20, 21, 24, 25, 28}; 
  Color_t color[] = {kOrange, kSpring+5, kAzure+2, kBlue+1, kMagenta+2, kBlack};
  Int_t  marker[] = {20, 21, 24, 25, 34, 33}; 
  Color_t histoColor = customColor;
  //const Int_t kNhistosData = 4;
    
  // open input file  
  TFile * fileData = TFile::Open(nameData.Data());
  if (!(fileData && fileData->IsOpen())) { 
    Printf("ERROR: cannot open file");
    return 1; 
  }
    
  TList * listData = (TList*)fileData->Get(listName.Data());
  if (! listData) return 2;
//  TH2F * hVtx = (TH2F*) listData->FindObject("PhiXeXeMC_eventVtx");
  TH1F * hCounters = (TH1F*) listData->FindObject("hEventStat");
  TH1F * hMulti = (TH1F*) listData->FindObject("hAEventsVsMulti");
  TH1F * hMultiAccepted = (TH1F*) listData->FindObject("PhiXeXeMC_eventMult");
  TH3F * hTrueIn =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_True%s", cutLabel.Data()));
  if (!hTrueIn) return 3;
  TH3F * hTrueYIn =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_TrueY%s", cutLabel.Data()));
  if (!hTrueYIn) return 4;
  TH3F * hMother =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_Mother%s", cutLabel.Data()));
  if (!hMother) return 5;
  TH3F * hMotherY =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_MotherY%s", cutLabel.Data()));
  if (!hMotherY && !doEffOnly) return 6;
  // TH3F * hPhaseSpace =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_PhaseSpace%s", cutLabel.Data()));
  // if (!hPhaseSpace) return;

  TCanvas * cev = new TCanvas("cev","events", 800,600);
  cev->cd(); gPad->SetLogy();
  Beautify(hMulti, kBlue+2, 1, 2, 24, 1.0);
  if (hMultiAccepted) Beautify(hMultiAccepted, kRed+2, 1, 2, 20, 1.0);
  hMulti->GetXaxis()->SetRangeUser(0.0, 100.0);
  hMulti->GetYaxis()->SetRangeUser(1., hMulti->GetMaximum()*2);
  hMulti->SetTitle("events; V0M percentile; events");
  hMulti->Draw("hist");
  if (hMultiAccepted) hMultiAccepted->Draw("same E0");
  cev->Print(Form("events_%s.png",nameData.Data()));
  
  //------------------------------
  //define binning
  //------------------------------
  Double_t centA[] = {0.0, 10.0, 30.0, 60.0, 90.0};   
  Double_t centB[] = {0.0, 5.0, 10.0, 30.0, 50.0, 70.0, 90.0};
  Double_t centC[] = {0.0, 10.0, 30.0, 50.0, 70.0, 90.0};
  Double_t centD[] = {0.0, 30.0, 60.0, 90.0};
  Double_t cento[] = {0.0, 10.0, 30.0, 90.0};
  Double_t   pt1[100]; 
  pt1[0] = 0.0; 
  for(int j = 1; j<100; j++){ pt1[j] = pt1[j-1]+0.1;}
  Double_t   pt2[] = {0.0, 0.3, 0.5, 0.7, 1.00, 1.50, 2.00, 2.50, 3.00, 3.5, 4.00, 4.5, 5.0, 7.0, 10.0};
  Double_t   pt3[] = {0.0, 0.3, 0.5, 0.7, 0.9, 1.10, 1.30, 1.50, 2.00, 3.00, 4.00, 5.0, 7.0, 10.0};
  Double_t   rap[] = {-0.5, 0.5};//{-0.5, -0.3, -0.1, 0.1, 0.3, 0.5};

  Double_t minbias[] = {0.0, 90.0};
  
  Int_t     npt   = 0;
  Int_t   ncent   = 0;
  Int_t    nrap   = sizeof(rap) / sizeof(rap[0]) - 1;
  TAxis *ptbins   = 0;
  TAxis *centbins = 0;
  Double_t * selectedPtBinning;
  
  if (binning.Contains("C")) {
      ncent = sizeof(centC) / sizeof(centC[0]) - 1;
      centbins = new TAxis(ncent, centC);
    } else
  if (binning.Contains("A")) {
    ncent = sizeof(centA) / sizeof(centA[0]) - 1;
    centbins = new TAxis(ncent, centA);
  }  else
    if (binning.Contains("B")) {
      ncent = sizeof(centB) / sizeof(centB[0]) - 1;
      centbins = new TAxis(ncent, centB);
    } if (binning.Contains("final")) {
      ncent = sizeof(cento) / sizeof(cento[0]) - 1;
      centbins = new TAxis(ncent, cento);
    } 

  if (binning.Contains("1")) {
    npt = sizeof(pt1) / sizeof(pt1[0]) - 1;
    ptbins = new TAxis(npt, pt1);
    selectedPtBinning = pt1;
  } else if (binning.Contains("2")) {
    npt = sizeof(pt2) / sizeof(pt2[0]) - 1;
    ptbins = new TAxis(npt, pt2);
    selectedPtBinning = pt2;
  } else if (binning.Contains("3") || binning.Contains("final")) {
    npt = sizeof(pt3) / sizeof(pt3[0]) - 1;
    ptbins = new TAxis(npt, pt3);
    selectedPtBinning = pt3;
  }

  Printf(":::::::::::: Binning scheme: %s -- ncent = %i, npt = %i", binning.Data(), ncent, npt);
  //output file
  TString efffilename = Form("eff_%s_%s.root", binning.Data(), cutLabel.Data());
  TFile * efffile = new TFile(efffilename.Data(), "recreate");

  //create projections  
  TList * lTrue = new TList(); lTrue->SetName(Form("true_%s", cutLabel.Data()));
  TList * lMother = new TList(); lMother->SetName(Form("mother_%s", cutLabel.Data()));
  TList * lTrueY = new TList(); lTrueY->SetName(Form("trueY_%s", cutLabel.Data()));
  TList * lMotherY = new TList(); lMotherY->SetName(Form("motherY_%s", cutLabel.Data()));

  projectorTH3_InvMass_Centrality_Pt projector;

  TH3F * hInput[] = {hTrueIn, hMother, hTrueYIn, hMotherY};
  if (hInput[0]) hInput[0]->SetName("True");
  if (hInput[1]) hInput[1]->SetName("Mother");
  if (hInput[2]) hInput[2]->SetName("TrueY");
  if (hInput[3]) hInput[3]->SetName("MotherY");

  //---------------------------------
  //efficiency vs pt
  //---------------------------------
  //project pt vs cent
  for (Int_t i = 0; i < 2; i++) {
    if (!hInput[i]) continue;
    projector.SetPrefix(hInput[i]->GetName());
    TList * out = (i==0? out = lTrue : out = lMother);    
    if (binning.Contains("A")) projector.MultiProjPtCent(npt, selectedPtBinning, ncent, centA,  hInput[i], out);
    else if (binning.Contains("B")) projector.MultiProjPtCent(npt, selectedPtBinning, ncent, centB,  hInput[i], out);
    else if (binning.Contains("C")) projector.MultiProjPtCent(npt, selectedPtBinning, ncent, centC,  hInput[i], out);
      //hack for final efficiency: 30-90% bins are merged to reduce statistical fluctuations, since effs are compatible
    else if (binning.Contains("final")) projector.MultiProjPtCent(npt, selectedPtBinning, ncent, cento,  hInput[i], out);
    projector.SetPrefix(Form("%smb", hInput[i]->GetName()));
    projector.MultiProjPtCent(npt, selectedPtBinning, 1, minbias,  hInput[i], out);
  }
  Printf(":::: Projected according to binning strategy %s", binning.Data());
  if (!lTrue || !lMother)  7;  

  //Print efficiency plot on canvas
  TCanvas *ctrue[6];
  TCanvas * ceffo = new TCanvas("ceffo","efficiency", 800,600);
  TCanvas * cratio = new TCanvas("cratio","efficiency ratio", 800,600);

  for (Int_t icentbin = 0; icentbin<ncent+1;icentbin++){

    TString centLabel = Form("(%2.0f-%02.0f %%)", centbins->GetBinLowEdge(icentbin+1), centbins->GetBinLowEdge(icentbin+2));   
    //Create histos for efficiency
    TH1F * hTrueCountsPM = new TH1F(Form("hTrueCounts%i",icentbin),"True Counts", npt, selectedPtBinning);
    TH1F * hMotherCounts = new TH1F(Form("hMotherCounts%i",icentbin),"hMotherCounts", npt, selectedPtBinning);
    TH1F * hEffVsPt = new TH1F(Form("hEffVsPt%i", icentbin), "; #it{p}_{T} (GeV/#it{c}); A #times #varepsilon", npt, selectedPtBinning);
    if (icentbin == ncent) {
      centLabel = "(0-90%)";
      hEffVsPt->SetName("hEffVsPtMB");
    }
    hEffVsPt->SetTitle(Form("#epsilon_{#phi}, %s %s", cutLabel.Data(), centLabel.Data()));
    hEffVsPt->GetYaxis()->SetRangeUser(0.0,1.0);
    Beautify(hTrueCountsPM, color[icentbin], 1, 2, 24, 1.0);
    Beautify(hMotherCounts, color[icentbin], 1, 2, 25, 1.0);
    Beautify(hEffVsPt, color[icentbin], 1, 2, marker[icentbin], 1.0);

    //loop over pt bins
    if (icentbin==ncent) ctrue[icentbin] = new TCanvas(Form("true_090"), Form("true_090"),1200,800);
    else ctrue[icentbin] = new TCanvas(Form("true_%i",icentbin),Form("true_%i",icentbin),1200,800);
    ctrue[icentbin]->Divide(4, 3);

    for (Int_t iptbin = 1; iptbin<npt+1; iptbin++){
      
      TH1F * htmp;
      if (icentbin == ncent) htmp = (TH1F*) lTrue->FindObject(Form("Truemb_ptBin%02i_centBin00", iptbin-1));
      else htmp = (TH1F*) lTrue->FindObject(Form("True_ptBin%02i_centBin%02i", iptbin-1, icentbin));
      Beautify(htmp, kBlack, 1, 2, 20, 1.0);
      Double_t ntrue = htmp->Integral();
      Double_t ntrueerr = TMath::Sqrt(ntrue);				       
      hTrueCountsPM->SetBinContent(iptbin, ntrue);
      hTrueCountsPM->SetBinError(iptbin, ntrueerr);
      htmp->GetXaxis()->SetRangeUser(0.98, 1.08);
      htmp->GetXaxis()->SetTitle("#it{M}_{KK} GeV/#it{c}");
      htmp->GetYaxis()->SetTitle("counts / (1 MeV/#it{c}^{2})");
      htmp->SetTitle(Form("%3.1f < #it{p}_{T} < %3.1f GeV/#it{c} %s", selectedPtBinning[iptbin-1], selectedPtBinning[iptbin], centLabel.Data()));

      TH1F * htmp2;
      if (icentbin == ncent) htmp2 = (TH1F*) lMother->FindObject(Form("Mothermb_ptBin%02i_centBin00", iptbin-1));
      else htmp2 = (TH1F*) lMother->FindObject(Form("Mother_ptBin%02i_centBin%02i", iptbin-1, icentbin));
      Beautify(htmp2, kOrange, 1, 2, 20, 1.0);
      Double_t nmum = htmp2->Integral();
      Double_t nmumerr = TMath::Sqrt(nmum);				       
      hMotherCounts->SetBinContent(iptbin, nmum);
      hMotherCounts->SetBinError(iptbin, nmumerr);
      htmp2->GetXaxis()->SetRangeUser(0.98, 1.08);
      htmp2->SetTitle(Form("%3.1f < #it{p}_{T} < %3.1f GeV/#it{c} %s", selectedPtBinning[iptbin-1], selectedPtBinning[iptbin], centLabel.Data()));
      htmp2->GetXaxis()->SetTitle("#it{M}_{KK} GeV/#it{c}");
      htmp2->GetYaxis()->SetTitle("counts / (1 MeV/#it{c}^{2})");

      ctrue[icentbin]->cd(iptbin);
      gPad->SetLogy();
      htmp2->GetYaxis()->SetRangeUser(1., htmp2->GetMaximum()*4.0);
      htmp2->Draw();
      htmp->Draw("same");
    
      Double_t ratio = (nmum > 0 ? ntrue/nmum : 0.0); 
      Double_t err = ratio * TMath::Sqrt(TMath::Power(ntrueerr/ntrue, 2.0) + TMath::Power(nmumerr/nmum, 2.0));
      if (ntrue == 0) { ratio = -1.0; err = 0.0;}
      hEffVsPt->SetBinContent(iptbin, ratio);
      hEffVsPt->SetBinError(iptbin, err);
      //Printf("Efficiency  %3.1f < pt < %3.1f = %6.5f +/- %6.5f", pt[iptbin-1], pt[iptbin], ratio, err);
    }
    TString namePrinted = Form("trueVsGen_%s_%i_%s.png", binning.Data(), icentbin, cutLabel.Data());
    ctrue[icentbin]->Print(namePrinted.Data());

    //display legend for eff
    ceffo->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    TString opt = (icentbin>0 ? "same" : "");
    hEffVsPt->GetYaxis()->SetNdivisions(509);
    hEffVsPt->GetXaxis()->SetRangeUser(0.0001, selectedPtBinning[npt]);
    hEffVsPt->Draw(opt.Data());
    
    // save into a file
    efffile->cd();
    hMotherCounts->Write();
    hTrueCountsPM->Write();
    hEffVsPt->Write();
    ctrue[icentbin]->Write();   
  }
  
  for (Int_t icentbin = 0; icentbin < ncent; icentbin++){
    TH1F * hEffMB = (TH1F*)((TH1F*) efffile->Get("hEffVsPtMB"))->Clone("hEffMB");
    TH1F * hRatio2MB = (TH1F*)((TH1F*) efffile->Get(Form("hEffVsPt%i", icentbin)))->Clone(Form("hRatio2MB%i", icentbin));
    Beautify(hRatio2MB, color[icentbin], 1, 2, marker[icentbin], 1.0);
    hRatio2MB->Divide(hRatio2MB, hEffMB, 1.0, 1.0, "B");
    efffile->cd();
    hRatio2MB->Write();
    cratio->cd();
    TString opt = (icentbin>0 ? "same" : "");
    hRatio2MB->Draw(opt.Data());
    hRatio2MB->GetYaxis()->SetNdivisions(509);
    hRatio2MB->GetYaxis()->SetRangeUser(0.5, 1.5);
  }

  gStyle->SetOptTitle(0);
  ceffo->cd();
  TLegend * legeff = (TLegend*) gPad->BuildLegend(0.6, 0.75, 0.92, 0.98);
  legeff->SetFillStyle(1001); legeff->SetFillColor(kWhite);
  efffilename.ReplaceAll(".root",".png");
  ceffo->Print(efffilename.Data());
  cratio->cd();
  legeff->Draw();
  cratio->Print(efffilename.Prepend("ratio2MB"));

  efffile->cd();
  lMother->Write("Mother");
  lTrue->Write("True");

  if (doEffOnly) return 8;
  //-------------------- -------------
  // efficiency vs y
  //---------------------------------
  //Project pt vs Y
  // projector.SetPrefix(hInput[2]->GetName());
  // projector.MultiProjPtCent(npt, selectedPtBinning, nrap, rap,  hTrueYIn, lTrueY);
  // projector.SetPrefix(hInput[3]->GetName());
  // projector.MultiProjPtCent(npt, selectedPtBinning, nrap, rap,  hMotherY, lMotherY);
  // if (!lTrueY || !lMotherY) {
  //   Printf("invalid lists");
  //   return;    
  // } else {
  //   efffile->cd();
  //   lMotherY->Write();
  //   lTrueY->Write();
  // }
  // //display
  // TCanvas * ceffy = new TCanvas("ceffy","efficiency", 800,600);

  // //get efficiency vs y and pt
  // for (Int_t irapbin = 0; irapbin<nrap;irapbin++){
    
  //   TString rapLabel = Form("(%2.1f < #it{y} < %2.1f)", rap[irapbin], rap[irapbin+1]);
  //   TH1F * hEffVsY = new TH1F(Form("hEffVsY%i", irapbin), "; #it{p}_{T} (GeV/#it{c}); A #times #varepsilon", npt, selectedPtBinning);
  //   hEffVsY->SetTitle(Form("#epsilon_{#phi}, %s %s", cutLabel.Data(), rapLabel.Data()));
  //   hEffVsY->GetYaxis()->SetRangeUser(0.0,1.0);
  //   Beautify(hEffVsY, color[irapbin], 1, 2, marker[irapbin], 1.0);
    
  //   Float_t effAndErr[2] = {0.,0.};
  //   //loop over pt bins
  //   for (Int_t iptbin = 1; iptbin<npt+1; iptbin++){
  //     TH1F *  htmp = (TH1F*) lTrueY->FindObject(Form("TrueY_ptBin%02i_centBin%02i", iptbin-1, irapbin));
  //     TH1F *  htmp2 = (TH1F*) lMotherY->FindObject(Form("MotherY_ptBin%02i_centBin%02i", iptbin-1, irapbin));
  //     GetEfficiencyFromBinnedMinv(htmp, htmp2, effAndErr);
  //     hEffVsY->SetBinContent(iptbin, effAndErr[0]);
  //     hEffVsY->SetBinError(iptbin, effAndErr[1]);
  //   }
  //   //show
  //   ceffy->cd();
  //   TString opt = (irapbin>0 ? "same" : "");
  //   hEffVsY->Draw(opt.Data());
  //   //save to file
  //   efffile->cd();
  //   hEffVsY->Write();
  // }
  // TLegend * legeffY = (TLegend*) gPad->BuildLegend(0.6, 0.75, 0.92, 0.98);
  // legeffY->SetFillStyle(1001); legeffY->SetFillColor(kWhite);
  // ceffy->Print(efffilename.ReplaceAll("ratio2MB","rapidity"));
 
  //---------------------------------
  // mass, pt and y resolution
  //---------------------------------
  //display
  TCanvas * cres = new TCanvas("cres","pt, y resolution", 800,600);
  cres->Divide(1,2);

  TCanvas * cresMethods = new TCanvas("cresMethods","resolution - methods", 900,1500);
  cresMethods->Divide(3,2);

  TCanvas * cresGausCent = new TCanvas("cresGausCent","mass resolution, Gaussian fit", 800,600);
  TCanvas * cresRMSCent = new TCanvas("cresRMSCent","mass resolution, RMS", 800,600);
  TCanvas * cresVMCCent = new TCanvas("cresVMCCent","mass resolution, Voigt fit to MC true", 800,600);

  TCanvas * cres_cent[6];
  for (Int_t icent = 0; icent<ncent; icent++) {
    cres_cent[icent] = new TCanvas(Form("cres_cent%i",icent),"resolution", 1000,600);
    cres_cent[icent]->Divide(2,2);
  }
  
  //Get input
  TH3F * hInputRes[4] = {0x0, 0x0, 0x0, 0x0};
  hInputRes[0] =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_ResCent%s", cutLabel.Data()));
  //hInputRes[1] =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_Res%s", cutLabel.Data()));
  //hInputRes[2] =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_ResPt%s", cutLabel.Data()));
  //hInputRes[3] =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_ResY%s", cutLabel.Data()));
 
  hInputRes[0]->SetName("MassResC");
  //hInputRes[1]->SetName("MassRes");
  //hInputRes[2]->SetName("PtRes");
  //hInputRes[3]->SetName("RapRes");

  //project mass resolution vs pt and Y
  TFile * fres = new TFile(Form("res_%s_%s.root", binning.Data(), cutLabel.Data()),"recreate");
  TList * lResolution = new TList();
  lResolution->SetName(Form("Resolution_%s", cutLabel.Data()));

  for (int ii = 0; ii<1; ii++) {
    if (ii == 1) {
      projector.SetPrefix(hInputRes[ii]->GetName());
      if (binning.Contains("A")) projector.MultiProjPtCent(npt, selectedPtBinning, ncent, centA, hInputRes[ii], lResolution);
      else if (binning.Contains("B")) projector.MultiProjPtCent(npt, selectedPtBinning, ncent, centB, hInputRes[ii], lResolution);
      else if (binning.Contains("C") || binning.Contains("final")) projector.MultiProjPtCent(npt, selectedPtBinning, ncent, centC, hInputRes[ii], lResolution);
      projector.SetPrefix(Form("%smb", hInputRes[ii]->GetName()));
      projector.MultiProjPtCent(npt, selectedPtBinning, 1, minbias,  hInputRes[ii], lResolution);
    } else {
      projector.SetPrefix(hInputRes[ii]->GetName());
      projector.MultiProjPtCent(npt, selectedPtBinning, nrap, rap,  hInputRes[ii], lResolution);   
    }
  }
  fres->cd();
  lResolution->Write();
  
  //--------------------------------
  //get mass resolution vs y and pt
  //--------------------------------  
  TH1F * htmp = NULL;
  Float_t sigmaAndErr[2] = {0., 0.};
  TString opt = "";
  Color_t rescol[6] = {kRed+1, kOrange, kSpring-5, kBlue+1, kMagenta+2, kBlack};
  Float_t massResolRange[4] = {0.002, 0.003, 0.004, 0.005};
  Float_t fitLowRange[4] = {1.008, 1.007, 1.006, 1.005};
  Float_t fitHighRange[4] = {1.032, 1.033, 1.034, 1.035};

  TF1* fitFcn = GetVOIGT(1.00, 1.1);// norm, mean, res, width
  Printf("A");
  fitFcn->SetParLimits(1, 1.0, 1.030);
  Printf("B");
  fitFcn->SetParLimits(2, 0.001, 0.01);//0.004266); // width fixed to the pdg value
  Printf("C");
  fitFcn->SetParLimits(3, 0.001, 0.01);//0.004266); // width fixed to the pdg value
  Printf("D");
  fitFcn->SetParameter(3, 0.004266); // width fixed to the pdg value
  if (!fitFcn) return 9;
  
  for (Int_t icent = 0; icent<ncent;icent++){
    Printf ("\n\n Centrality bin %i ----------- %f - %f", icent, centC[icent], centC[icent+1]);
    TH1F * hResVsPt[4];
    TH1F * hResVsPtRMS[4];
    TH1F * hResVsPtVMC[4];
    TH1F * hMassVsPtVMC[4];
    TH1F * hWidthVsPtVMC[4];
    TString centLabel = Form("%2.0f-%2.0f %%", centC[icent], centC[icent+1]);
    if (icent == ncent) centLabel = "0-90%";
    
    for (Int_t i = 0; i<4; i++){
      hMassVsPtVMC[i] = new TH1F(Form("hMassVsPtVMC%i_r%i", icent, i), "; #it{p}_{T} (GeV/#it{c}); M_{#phi, VMC} (GeV/#it{c}^{2})", npt, selectedPtBinning);
      hMassVsPtVMC[i]->SetTitle(Form("M_{#phi,VMC}, Fit %4.3f<#it{M}<%4.3f", fitLowRange[i], fitHighRange[i]));
      Beautify(hMassVsPtVMC[i], rescol[icent]-i, 1, 2, marker[i], 1.0);
      
      hWidthVsPtVMC[i] = new TH1F(Form("hWidthVsPtVMC%i_r%i", icent, i), "; #it{p}_{T} (GeV/#it{c}); #Gamma_{#phi, VMC} (GeV/#it{c}^{2})", npt, selectedPtBinning);
      hWidthVsPtVMC[i]->SetTitle(Form("#Gamma_{#phi,VMC}, Fit %4.3f<#it{M}<%4.3f", fitLowRange[i], fitHighRange[i]));
      Beautify(hWidthVsPtVMC[i], rescol[icent]-i, 1, 2, marker[i], 1.0);

      hResVsPtVMC[i] = new TH1F(Form("hResVsPtVMC%i_r%i", icent, i), "; #it{p}_{T} (GeV/#it{c}); #sigma_{#phi, VMC} (GeV/#it{c}^{2})", npt, selectedPtBinning);
      hResVsPtVMC[i]->SetTitle(Form("#sigma_{#phi,VMC}, Fit %4.3f<#it{M}<%4.3f", fitLowRange[i], fitHighRange[i]));
      Beautify(hResVsPtVMC[i], rescol[icent]-i, 1, 2, marker[i], 1.0);
      
      hResVsPt[i] = new TH1F(Form("hResVsPt%i_res%i", icent, i), "; #it{p}_{T} (GeV/#it{c}); #sigma_{#phi, Gaus} (GeV/#it{c}^{2})", npt, selectedPtBinning);
      hResVsPt[i]->SetTitle(Form("#sigma_{#phi,Gaus}, |#Delta#it{M}| < %4.3f", massResolRange[i]));
      Beautify(hResVsPt[i], rescol[icent]-i, 1, 2, marker[i], 1.0);
      
      hResVsPtRMS[i] = new TH1F(Form("hResVsPtRMS%i_r%i", icent, i), "; #it{p}_{T} (GeV/#it{c}); #sigma_{#phi, RMS} (GeV/#it{c}^{2})", npt, selectedPtBinning);
      hResVsPtRMS[i]->SetTitle(Form("#sigma_{#phi,RMS}, |#Delta#it{M}| < %i#sigma_{ext}", i+2));
      Beautify(hResVsPtRMS[i], rescol[icent]-i, 1, 2, marker[i], 1.3);
      
      //loop over pt bins
      for (Int_t iptbin = 1; iptbin<npt+1; iptbin++){
        if (icent == ncent) htmp = (TH1F*) lResolution->FindObject(Form("MassResCmb_ptBin%02i_centBin%02i", iptbin-1, 0));
        else htmp = (TH1F*) lResolution->FindObject(Form("MassResC_ptBin%02i_centBin%02i", iptbin-1, icent));
        Beautify(htmp, kAzure+7-iptbin, 1, 2, 1, 1.0);	

        //get resolution from gaussian fit of mrec - mgen
        GetGausSigma(htmp, sigmaAndErr, massResolRange[i]);       
        hResVsPt[i]->SetBinContent(iptbin, sigmaAndErr[0]);
        hResVsPt[i]->SetBinError(iptbin, sigmaAndErr[1]);

        //get resolution from RMS of mrec - mgen, by also correcting for the truncated tails
        GetSigmaRMS(htmp, sigmaAndErr, i+2);       
        hResVsPtRMS[i]->SetBinContent(iptbin, sigmaAndErr[0]/GetTruncationCorrection(i+2));
        hResVsPtRMS[i]->SetBinError(iptbin, sigmaAndErr[1]/GetTruncationCorrection(i+2));
        //Printf("::::: Truncation correction factor = 1./%f",GetTruncationCorrection(i+2));

        //voigtian fit of true distribution
        if (icent < ncent) htmp = (TH1F*) lTrue->FindObject(Form("True_ptBin%02i_centBin%02i", iptbin-1, icent));
        else htmp = (TH1F*) lTrue->FindObject(Form("Truemb_ptBin%02i_centBin00", iptbin-1));

        if (!htmp) {
          // Printf(":::: WARNING : no input to fit Voigt with!");
          hMassVsPtVMC[i]->SetBinContent(iptbin, 0.0);
          hWidthVsPtVMC[i]->SetBinContent(iptbin, 0.0);
          hResVsPtVMC[i]->SetBinContent(iptbin, 0.0);
        } else {
          TFitResultPtr result = htmp->Fit(fitFcn, "SBRNQ", "", fitLowRange[i],fitHighRange[i]);
          //Printf("mean = %f, \n sigma = %f, \n res = %f", fitFcn->GetParameter(1),fitFcn->GetParameter(2), fitFcn->GetParameter(3));
          hMassVsPtVMC[i]->SetBinContent(iptbin, fitFcn->GetParameter(1));
          hWidthVsPtVMC[i]->SetBinContent(iptbin, fitFcn->GetParameter(3));
          hResVsPtVMC[i]->SetBinContent(iptbin, fitFcn->GetParameter(2));
          
          hMassVsPtVMC[i]->SetBinError(iptbin, fitFcn->GetParError(1));
          hWidthVsPtVMC[i]->SetBinError(iptbin, fitFcn->GetParError(3));
          hResVsPtVMC[i]->SetBinError(iptbin, fitFcn->GetParError(2));
        }
      } //end loop on pt

      fres->cd();
      hResVsPt[i]->Write();
      hResVsPtRMS[i]->Write();
      hResVsPtVMC[i]->Write();
      hMassVsPtVMC[i]->Write();
      hWidthVsPtVMC[i]->Write();
      
      //draw
      if (i>0) opt = "same";
      cres_cent[icent]->cd(1);
      hResVsPt[i]->GetYaxis()->SetRangeUser(0.0001, 0.0045);
      hResVsPt[i]->Draw(opt.Data());

      cres_cent[icent]->cd(2);
      hResVsPtRMS[i]->GetYaxis()->SetRangeUser(0.0001, 0.0045);
      hResVsPtRMS[i]->Draw(opt.Data());

      cres_cent[icent]->cd(3);
      hResVsPtVMC[i]->GetYaxis()->SetRangeUser(0.0001, 0.0045);
      hResVsPtVMC[i]->Draw(opt.Data());

      cres_cent[icent]->cd(4);
      hMassVsPtVMC[i]->GetYaxis()->SetRangeUser(1.018, 1.022);
      hMassVsPtVMC[i]->Draw(opt.Data());
      
      opt = "";      
    } //end loop on range

    if (icent>0) opt = "same";
    cresGausCent->cd();
    hResVsPt[0]->Draw(opt.Data());
    hResVsPt[3]->Draw("same");
    cresRMSCent->cd();
    hResVsPtRMS[0]->Draw(opt.Data());
    hResVsPtRMS[3]->Draw("same");
    cresVMCCent->cd();
    hResVsPtVMC[0]->Draw(opt.Data());
    hResVsPtVMC[3]->Draw("same");

    cresMethods->cd(icent+1);
    hResVsPt[1]->Draw(opt.Data());
    hResVsPtRMS[1]->SetMarkerStyle(30);
    hResVsPtRMS[1]->Draw("same");
    hResVsPtVMC[1]->SetMarkerStyle(28);
    hResVsPtVMC[1]->Draw("same");
    gPad->BuildLegend(0.25, 0.7, 0.6, 0.9, centLabel);
    //draw legends
    cres_cent[icent]->cd(1);
    TLegend * legRes = (TLegend*) gPad->BuildLegend(0.2, 0.6, 0.62, 0.9, centLabel.Data(), "p");
    legRes->SetFillStyle(0);
    legRes->SetBorderSize(0);
    cres_cent[icent]->cd(2);
    TLegend * legRes2 = (TLegend*) gPad->BuildLegend(0.2, 0.6, 0.62, 0.9, centLabel.Data(), "p");
    legRes2->SetFillStyle(0);
    legRes2->SetBorderSize(0);
    cres_cent[icent]->cd(3);
    TLegend * legRes3 = (TLegend*) gPad->BuildLegend(0.2, 0.6, 0.62, 0.9, centLabel.Data(), "p");
    legRes3->SetFillStyle(0);
    legRes3->SetBorderSize(0);
    TString resSave = Form("%s_cent%i.pdf", fres->GetName(), icent);
    resSave.ReplaceAll(".root","");
    cres_cent[icent]->Print(resSave.Data());    
  }//end loop on centrality

  cresGausCent->cd();
  TLegend * legResG = (TLegend*) gPad->BuildLegend(0.2, 0.7, 0.7, 0.9, "", "p");
  legResG->SetFillStyle(0);
  legResG->SetBorderSize(0);
  legResG->SetNColumns(2);
  legResG->SetTextSize(0.035);
  cresGausCent->Print("centDep_GaussRes.pdf");

  cresRMSCent->cd();
  TLegend * legResR = (TLegend*) gPad->BuildLegend(0.2, 0.7, 0.7, 0.9, "", "p");
  legResR->SetFillStyle(0);
  legResR->SetBorderSize(0);
  legResR->SetNColumns(2);
  legResR->SetTextSize(0.035);
  cresRMSCent->Print("centDep_RMSRes.pdf");

  cresVMCCent->cd();
  TLegend * legResV = (TLegend*) gPad->BuildLegend(0.2, 0.7, 0.8, 0.9, "", "p");
  legResV->SetFillStyle(0);
  legResV->SetBorderSize(0);
  legResV->SetNColumns(2);
  legResV->SetTextSize(0.025);
  cresVMCCent->Print("centDep_VMCRes.pdf");
  
//loop on rapidity bins
/*
  for (Int_t irapbin = 0; irapbin<nrap;irapbin++){
    TString rapLabel = Form("(%2.1f < #it{y} < %2.1f)", rap[irapbin], rap[irapbin+1]);
    //pt and rapidity resolution
    TH1F * hPtResVsPt = new TH1F(Form("hPtResVsPt%i", irapbin), "#it{p}_{T} resolution from gaussian fit; #it{p}_{T} (GeV/#it{c}); #sigma_{pT}", npt, selectedPtBinning);
    hPtResVsPt->SetTitle(Form("#sigma_{#phi}, Gaus fit %s", rapLabel.Data()));
    Beautify(hPtResVsPt, color[irapbin], 1, 2, marker[irapbin], 1.0);
    
    TH1F * hPtResVsPtRMS = new TH1F(Form("hPtResVsPtRMS%i", irapbin), "#it{p}_{T} resolution from RMS; #it{p}_{T} (GeV/#it{c}); #sigma_{pT}", npt, selectedPtBinning);
    hPtResVsPtRMS->SetTitle(Form("#sigma_{#phi}, RMS %s", rapLabel.Data()));
    Beautify(hPtResVsPtRMS, color[irapbin]-2, 1, 2, marker[irapbin], 1.0);
    
    TH1F * hYResVsPt = new TH1F(Form("hYResVsPt%i", irapbin), "#it{y} resolution from gaussian fit; #it{p}_{T} (GeV/#it{c}); #sigma_{y}", npt, selectedPtBinning);
    hYResVsPt->SetTitle(Form("#sigma_{y}, Gaus Fit %s", rapLabel.Data()));
    Beautify(hYResVsPt, color[irapbin], 1, 2, marker[irapbin], 1.0);
    
    TH1F * hYResVsPtRMS = new TH1F(Form("hYResVsPtRMS%i", irapbin), "#it{y} resolution from RMS; #it{p}_{T} (GeV/#it{c}); #sigma_{y}", npt, selectedPtBinning);
    hYResVsPtRMS->SetTitle(Form("#sigma_{y}, RMS %s", rapLabel.Data()));
    Beautify(hYResVsPtRMS, color[irapbin]-2, 1, 2, marker[irapbin], 1.0);
    
    for (Int_t iptbin = 1; iptbin<npt+1; iptbin++){
      //Pt resolution
      htmp = (TH1F*) lResolution->FindObject(Form("PtRes_ptBin%02i_centBin%02i", iptbin-1, irapbin));
      GetGausSigma(htmp, sigmaAndErr);       
      hPtResVsPt->SetBinContent(iptbin, sigmaAndErr[0]);
      hPtResVsPt->SetBinError(iptbin, sigmaAndErr[1]);
     
      hPtResVsPtRMS->SetBinContent(iptbin, htmp->GetRMS());
      hPtResVsPtRMS->SetBinError(iptbin, htmp->GetRMSError());

      //Y resolution
      htmp = (TH1F*) lResolution->FindObject(Form("RapRes_ptBin%02i_centBin%02i", iptbin-1, irapbin));
      GetGausSigma(htmp, sigmaAndErr);       
      hYResVsPt->SetBinContent(iptbin, sigmaAndErr[0]);
      hYResVsPt->SetBinError(iptbin, sigmaAndErr[1]);

      hYResVsPtRMS->SetBinContent(iptbin, htmp->GetRMS());
      hYResVsPtRMS->SetBinError(iptbin, htmp->GetRMSError());
    }
    //save to file
    fres->cd();
    hPtResVsPt->Write();
    hPtResVsPtRMS->Write();
    hYResVsPt->Write();
    hYResVsPtRMS->Write();
    
    //show
    TString opt = ""; 
    opt = (irapbin>0 ? "same" : "");
    cres->cd(1);
    hPtResVsPt->Draw(opt.Data());
    cres->cd(2);
    hYResVsPt->Draw(opt.Data());
  }//end loop on rapidity bins
*/
  return 0;
 
}




void GetEfficiencyFromBinnedMinv(TH1F* trueMinv, TH1F* momMinv, Float_t* effAndErr)
{
  //
  //compute efficiency for a given bin by integrating the projections in Minv for that bin 
  //
  if (!trueMinv || !momMinv) return;
  Double_t ntrue = trueMinv->Integral();
  Double_t ntrueerr = TMath::Sqrt(ntrue);				       
  
  Double_t nmum = momMinv->Integral();
  Double_t nmumerr = TMath::Sqrt(nmum);				       
  
  Double_t ratio = (nmum > 0 ? ntrue/nmum : 0.0); 
  Double_t err = ratio * TMath::Sqrt(TMath::Power(ntrueerr/ntrue, 2.0) + TMath::Power(nmumerr/nmum, 2.0));
  if (ntrue == 0) { ratio = -1.0; err = 0.0;}

  effAndErr[0] = ratio;
  effAndErr[1] = err;
  return;
}

void GetGausSigma(TH1F * htmp, Float_t * sigmaAndErr, Float_t absRange)
{
  //returns sigma and its error from gaussian fit
  if (!htmp || !sigmaAndErr) return;
  if (htmp->GetEntries()<1) {
    sigmaAndErr[0] = 0; sigmaAndErr[1] = 0;
  }
  htmp->Fit("gaus","S0RQ","", -absRange, absRange);
  TF1* fitg = (TF1*) htmp->GetFunction("gaus");
  sigmaAndErr[0] = (fitg ? fitg->GetParameter(2) : 0.0);
  sigmaAndErr[1] = (fitg ? fitg->GetParError(2) : 0.0);
  return;
}

void GetSigmaRMS(TH1F * htmp, Float_t * sigmaAndErr, Int_t absRange)
{
  //returns sigma and its error from gaussian fit
  if (!htmp || !sigmaAndErr) return;
  if (htmp->GetEntries()<1) {
    sigmaAndErr[0] = 0; sigmaAndErr[1] = 0;
  }
  Float_t rmsAll = htmp->GetRMS();
  Int_t ibin1 = htmp->GetXaxis()->FindBin(-absRange*rmsAll);
  Int_t ibin2 = htmp->GetXaxis()->FindBin(absRange*rmsAll);
  htmp->GetXaxis()->SetRange(ibin1, ibin2);
  sigmaAndErr[0] = htmp->GetRMS();
  sigmaAndErr[1] = htmp->GetRMSError();
  htmp->GetXaxis()->SetRange(); 
  return;
}

Float_t GetTruncationCorrection(Int_t absRange)
{
  Float_t x;
  TF1 g1("g1", "TMath::Gaus(x,[0],[1])", -10., 10.);
  g1.SetParameter(0, 0.0);
  g1.SetParameter(1, 1.0);
  Float_t variance = g1.Variance(-absRange, absRange);
  return TMath::Sqrt(variance);
}


//-----------------------------------------
//Signal peak: Voigtian peak function
//-----------------------------------------
Double_t Voigt( Double_t *x, Double_t * par)
{
  return par[0] * TMath::Voigt(x[0]-par[1], par[2], par[3]);
}


//-----------------------------------------
TF1 * GetVOIGT(Double_t fitMin, Double_t fitMax)
{
  TF1* fitFcn = new TF1("VOIGTpoly0", Voigt, fitMin, fitMax, 4);
  fitFcn->SetParNames("Norm","Mass","Width", "Resolution"); 
  return fitFcn;
}
