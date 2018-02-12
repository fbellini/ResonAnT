//
// Splits a TH2F using the myHistSplig_sparse_IM_PT_CENT class
//
#include "/Users/fbellini/alice/macros/ResonAnT/projectorTH3_InvMass_Centrality_Pt.C"
#include "/Users/fbellini/alice/macros/MakeUp.C"
#include "TFile.h"

void GetEfficiencyFromBinnedMinv(TH1F* trueMinv=NULL, TH1F* momMinv=NULL, Float_t* effAndErr=NULL);
void GetGausSigma(TH1D * htmp = NULL, Float_t * sigmaAndErr = 0);

void projectMC(TString nameData = "phiXeXe_LHC17j7_RsnOut.root",
			TString listName = "RsnOut_default",
			TString cutLabel  = "tpc3s",
			TString binning  = "A3",
			Color_t customColor = kBlack)
{
   // initial setup
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTextFont(42);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
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
  
  Color_t color[]={kCyan+1, kAzure+2, kBlue+2, kAzure+4, kBlack};
  Color_t marker[]={21, 22, 32, 28, 20};
  Color_t histoColor = customColor;
  const Int_t kNhistosData = 4;
    
  // open input file  
  TFile * fileData = TFile::Open(nameData.Data());
  if (!(fileData && fileData->IsOpen())) { 
    Printf("ERROR: cannot open file");
    return; 
  }
    
  TList * listData = (TList*)fileData->Get(listName.Data());
  if (! listData) return;
  TString cutString("tpc3s");
  TH2F * hVtx = (TH2F*) listData->FindObject("PhiXeXeMC_eventVtx");
  TH1F * hCounters = (TH1F*) listData->FindObject("hEventStat");
  TH1F * hMulti = (TH1F*) listData->FindObject("hAEventsVsMulti");
  TH1F * hMultiAccepted = (TH1F*) listData->FindObject("PhiXeXeMC_eventMult");
  TH3F * hTrueIn =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_True_%s", cutString.Data()));
  TH3F * hTrueYIn =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_TrueY_%s", cutString.Data()));
  TH3F * hMother =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_Mother_%s", cutString.Data()));
  TH3F * hMotherY =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_MotherY_%s", cutString.Data()));
  TH3F * hPhaseSpace =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_PhaseSpace_%s", cutString.Data()));

  TCanvas * cev = new TCanvas("cev","events", 800,600);
  cev->cd(); gPad->SetLogy();
  Beautify(hMulti, kBlue+2, 1, 2, 24, 1.0);
  Beautify(hMultiAccepted, kRed+2, 1, 2, 20, 1.0);
  hMulti->GetXaxis()->SetRangeUser(0.0, 100.0);
  hMulti->GetYaxis()->SetRangeUser(1., hMulti->GetMaximum()*2);
  hMulti->SetTitle("events; V0M percentile; events");
  hMulti->Draw("hist");
  hMultiAccepted->Draw("same E0");
  cev->Print(Form("events_%s.png",nameData.Data()));
  
  //------------------------------
  //define binning
  //------------------------------
  Double_t centA[] = {0.0, 10.0, 30.0, 60.0, 90.0};   
  Double_t centB[] = {0.0, 20.0, 40.0, 60.0, 80.0};
  Double_t centC[] = {0.0, 30.0, 50.0, 70.0, 90.0};
  Double_t   pt1[] = {0.0, 0.3, 0.5, 1.00, 1.50, 2.00, 2.50, 3.00, 3.5, 4.00, 4.5, 5.0};// 6.0, 8.0, 10.0};
  Double_t   pt2[] = {0.0, 0.3, 0.5, 0.7, 1.00, 1.50, 2.00, 2.50, 3.00, 3.5, 4.00, 4.5, 5.0};// 7.0, 10.0};
  Double_t   pt3[] = {0.0, 0.3, 0.5, 0.7, 0.9, 1.10, 1.30, 1.50, 2.00, 3.00, 4.00, 5.0};//, 6.0, 8.0, 10.0};
  Double_t   rap[] = {-0.5, -0.3, -0.1, 0.1, 0.3, 0.5};

  Double_t minbias[] = {0.0, 90.0};
  
  Int_t     npt   = 0;
  Int_t   ncent   = 0;
  Int_t    nrap   = sizeof(rap) / sizeof(rap[0]) - 1;
  TAxis *ptbins   = 0;
  TAxis *centbins = 0;
  Double_t * selectedPtBinning;
  
  if (binning.Contains("1")) {
    npt = sizeof(pt1) / sizeof(pt1[0]) - 1;
    ptbins = new TAxis(npt, pt1);
    selectedPtBinning = pt1;
  } else if (binning.Contains("2")) {
    npt = sizeof(pt2) / sizeof(pt2[0]) - 1;
    ptbins = new TAxis(npt, pt2);
    selectedPtBinning = pt2;
  } else if (binning.Contains("3")) {
    npt = sizeof(pt3) / sizeof(pt3[0]) - 1;
    ptbins = new TAxis(npt, pt3);
    selectedPtBinning = pt3;
  }

  if (binning.Contains("A")) {
    ncent = sizeof(centA) / sizeof(centA[0]) - 1;
    centbins = new TAxis(ncent, centA);
  }  else
    if (binning.Contains("B")) {
      ncent = sizeof(centB) / sizeof(centB[0]) - 1;
      centbins = new TAxis(ncent, centB);
    } else
      if (binning.Contains("C")) {
	ncent = sizeof(centC) / sizeof(centC[0]) - 1;
	centbins = new TAxis(ncent, centC);
      }

  //output file
  TString efffilename = Form("eff_%s_%s", binning.Data(), nameData.Data());
  TFile * efffile = new TFile(efffilename.Data(), "recreate");

  //create projections  
  TList * lTrue = new TList(); lTrue->SetName(Form("true_%s", cutLabel.Data()));
  TList * lMother = new TList(); lMother->SetName(Form("mother_%s", cutLabel.Data()));
  TList * lTrueY = new TList(); lTrueY->SetName(Form("trueY_%s", cutLabel.Data()));
  TList * lMotherY = new TList(); lMotherY->SetName(Form("motherY_%s", cutLabel.Data()));

  projectorTH3_InvMass_Centrality_Pt projector;

  TH3F * hInput[] = {hTrueIn, hMother, hTrueYIn, hMotherY};
  hInput[0]->SetName("True");
  hInput[1]->SetName("Mother");
  hInput[2]->SetName("TrueY");
  hInput[3]->SetName("MotherY");

  //---------------------------------
  //efficiency vs pt
  //---------------------------------
  //project pt vs cent
  for (Int_t i = 0; i < 2; i++) {
    projector.SetPrefix(hInput[i]->GetName());
    TList * out = (i==0? out = lTrue : out = lMother);    
    if (binning.Contains("A")) projector.MultiProjPtCent(npt, selectedPtBinning, ncent, centA,  hInput[i], out);
    else if (binning.Contains("B")) projector.MultiProjPtCent(npt, selectedPtBinning, ncent, centB,  hInput[i], out);
    else if (binning.Contains("C")) projector.MultiProjPtCent(npt, selectedPtBinning, ncent, centC,  hInput[i], out);
    projector.SetPrefix(Form("%smb", hInput[i]->GetName()));
    projector.MultiProjPtCent(npt, selectedPtBinning, 1, minbias,  hInput[i], out);
  }
  Printf(":::: Projected according to binning strategy %s", binning.Data());
  if (!lTrue || !lMother) return;  

  //Print efficiency plot on canvas
  TCanvas *ctrue[5];
  TCanvas * ceffo = new TCanvas("ceffo","efficiency", 800,600);
  TCanvas * cratio = new TCanvas("cratio","efficiency ratio", 800,600);

  for (Int_t icentbin = 0; icentbin<ncent+1;icentbin++){

    TString centLabel = Form("(%2.0f-%02.0f%%)", centbins->GetBinLowEdge(icentbin+1), centbins->GetBinLowEdge(icentbin+2));
   
    //Create histos for efficiency
    TH1F * hTrueCountsPM = new TH1F(Form("hTrueCounts%i",icentbin),"True Counts", npt, selectedPtBinning);
    TH1F * hMotherCounts = new TH1F(Form("hMotherCounts%i",icentbin),"hMotherCounts", npt, selectedPtBinning);
    TH1F * hEffVsPt = new TH1F(Form("hEffVsPt%i", icentbin), "; #it{p}_{T} (GeV/#it{c}); A #times #varepsilon", npt, selectedPtBinning);
    if (icentbin == ncent) {
      centLabel = "( 0-90%)";
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
      Beautify(htmp, kOrange+2, 1, 2, 20, 1.0);
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
    TString namePrinted = Form("trueVsGen_%i_%s_%s", icentbin, binning.Data(), nameData.Data());
    namePrinted.ReplaceAll(".root",".png");
    ctrue[icentbin]->Print(namePrinted.Data());

    //display legend for eff
    ceffo->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    TString opt = (icentbin>0 ? "same" : "");
    hEffVsPt->GetYaxis()->SetNdivisions(509);
    hEffVsPt->GetXaxis()->SetRangeUser(0.001, selectedPtBinning[npt]);
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
  
  //---------------------------------
  // efficiency vs y
  //---------------------------------
  //Project pt vs Y
  projector.SetPrefix(hInput[2]->GetName());
  projector.MultiProjPtCent(npt, selectedPtBinning, nrap, rap,  hTrueYIn, lTrueY);
  projector.SetPrefix(hInput[3]->GetName());
  projector.MultiProjPtCent(npt, selectedPtBinning, nrap, rap,  hMotherY, lMotherY);
  if (!lTrueY || !lMotherY) {
    Printf("invalid lists");
    return;    
  } else {
    efffile->cd();
    lMotherY->Write();
    lTrueY->Write();
  }
  //display
  TCanvas * ceffy = new TCanvas("ceffy","efficiency", 800,600);

  //get efficiency vs y and pt
  for (Int_t irapbin = 0; irapbin<nrap;irapbin++){
    TString rapLabel = Form("(%2.1f < #it{y} < %2.1f)", rap[irapbin], rap[irapbin+1]);
    TH1F * hEffVsY = new TH1F(Form("hEffVsY%i", irapbin), "; #it{p}_{T} (GeV/#it{c}); A #times #varepsilon", npt, selectedPtBinning);
    hEffVsY->SetTitle(Form("#epsilon_{#phi}, %s %s", cutLabel.Data(), rapLabel.Data()));
    hEffVsY->GetYaxis()->SetRangeUser(0.0,1.0);
    Beautify(hEffVsY, color[irapbin], 1, 2, marker[irapbin], 1.0);
    
    Float_t effAndErr[2] = {0.,0.};
    //loop over pt bins
    for (Int_t iptbin = 1; iptbin<npt+1; iptbin++){
      TH1F *  htmp = (TH1F*) lTrueY->FindObject(Form("TrueY_ptBin%02i_centBin%02i", iptbin-1, irapbin));
      TH1F *  htmp2 = (TH1F*) lMotherY->FindObject(Form("MotherY_ptBin%02i_centBin%02i", iptbin-1, irapbin));
      GetEfficiencyFromBinnedMinv(htmp, htmp2, effAndErr);
      hEffVsY->SetBinContent(iptbin, effAndErr[0]);
      hEffVsY->SetBinError(iptbin, effAndErr[1]);
    }
    //show
    ceffy->cd();
    TString opt = (irapbin>0 ? "same" : "");
    hEffVsY->Draw(opt.Data());
    //save to file
    efffile->cd();
    hEffVsY->Write();
  }
  TLegend * legeffY = (TLegend*) gPad->BuildLegend(0.6, 0.75, 0.92, 0.98);
  legeffY->SetFillStyle(1001); legeffY->SetFillColor(kWhite);
  ceffy->Print(efffilename.ReplaceAll("ratio2MB","rapidity"));
 
  //---------------------------------
  // mass, pt and y resolution
  //---------------------------------
  TH3F * hInputRes[3];
  hInputRes[0] =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_Res_%s", cutString.Data()));
  hInputRes[1] =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_ResPt_%s", cutString.Data()));
  hInputRes[2] =  (TH3F*) listData->FindObject(Form("PhiXeXeMC_ResY_%s", cutString.Data()));
 
  hInputRes[0]->SetName("MassRes");
  hInputRes[1]->SetName("PtRes");
  hInputRes[2]->SetName("RapRes");

  //project mass resolution vs pt and Y
  TFile * fres = new TFile(Form("res_%s_%s", binning.Data(), nameData.Data()),"recreate");
  TList * lResolution = new TList();
  lResolution->SetName(Form("Resolution_%s", cutLabel.Data()));

  for (int ii = 0; ii<3; ii++) {
    projector.SetPrefix(hInputRes[ii]->GetName());
    projector.MultiProjPtCent(npt, selectedPtBinning, nrap, rap,  hInputRes[ii], lResolution);   
  }

  fres->cd();
  lResolution->Write();

  //display
  TCanvas * cres = new TCanvas("cres","resolution", 800,600);
  cres->Divide(2,2);

  //get resolution vs y and pt
  for (Int_t irapbin = 0; irapbin<nrap;irapbin++){
    TString rapLabel = Form("(%2.1f < #it{y} < %2.1f)", rap[irapbin], rap[irapbin+1]);
    TH1F * hResVsPt = new TH1F(Form("hResVsPt%i", irapbin), "mass resolution from gaussian fit; #it{p}_{T} (GeV/#it{c}); #sigma", npt, selectedPtBinning);
    hResVsPt->SetTitle(Form("#sigma_{#phi}, %s %s", cutLabel.Data(), rapLabel.Data()));
    Beautify(hResVsPt, color[irapbin], 1, 2, marker[irapbin], 1.0);
    
    TH1F * hResVsPtRMS = new TH1F(Form("hResVsPtRMS%i", irapbin), "mass resolution from RMS; #it{p}_{T} (GeV/#it{c}); #sigma", npt, selectedPtBinning);
    hResVsPtRMS->SetTitle(Form("#sigma_{#phi}, %s %s", cutLabel.Data(), rapLabel.Data()));
    Beautify(hResVsPtRMS, color[irapbin], 1, 2, marker[irapbin], 1.0);

    TH1F * hPtResVsPt = new TH1F(Form("hPtResVsPt%i", irapbin), "#it{p}_{T} resolution from gaussian fit; #it{p}_{T} (GeV/#it{c}); #sigma_{pT}", npt, selectedPtBinning);
    hPtResVsPt->SetTitle(Form("#sigma_{#phi}, %s %s", cutLabel.Data(), rapLabel.Data()));
    Beautify(hPtResVsPt, color[irapbin], 1, 2, marker[irapbin], 1.0);
    
    TH1F * hPtResVsPtRMS = new TH1F(Form("hPtResVsPtRMS%i", irapbin), "#it{p}_{T} resolution from RMS; #it{p}_{T} (GeV/#it{c}); #sigma_{pT}", npt, selectedPtBinning);
    hPtResVsPtRMS->SetTitle(Form("#sigma_{#phi}, %s %s", cutLabel.Data(), rapLabel.Data()));
    Beautify(hPtResVsPtRMS, color[irapbin], 1, 2, marker[irapbin], 1.0);

    TH1F * hYResVsPt = new TH1F(Form("hYResVsPt%i", irapbin), "#it{y} resolution from gaussian fit; #it{p}_{T} (GeV/#it{c}); #sigma_{y}", npt, selectedPtBinning);
    hYResVsPt->SetTitle(Form("#sigma_{y}, %s %s", cutLabel.Data(), rapLabel.Data()));
    Beautify(hYResVsPt, color[irapbin], 1, 2, marker[irapbin], 1.0);
    
    TH1F * hYResVsPtRMS = new TH1F(Form("hYResVsPtRMS%i", irapbin), "#it{y} resolution from RMS; #it{p}_{T} (GeV/#it{c}); #sigma_{y}", npt, selectedPtBinning);
    hYResVsPtRMS->SetTitle(Form("#sigma_{y}, %s %s", cutLabel.Data(), rapLabel.Data()));
    Beautify(hYResVsPtRMS, color[irapbin], 1, 2, marker[irapbin], 1.0);
    
    //loop over pt bins
    for (Int_t iptbin = 1; iptbin<npt+1; iptbin++){

      //Mass resolution
      TH1F *  htmp = (TH1F*) lResolution->FindObject(Form("MassRes_ptBin%02i_centBin%02i", iptbin-1, irapbin));
      Float_t sigmaAndErr[2] = {0., 0.};
      GetGausSigma(htmp, sigmaAndErr);       
      hResVsPt->SetBinContent(iptbin, sigmaAndErr[0]);
      hResVsPt->SetBinError(iptbin, sigmaAndErr[1]);
     
      hResVsPtRMS->SetBinContent(iptbin, htmp->GetRMS());
      hResVsPtRMS->SetBinError(iptbin, htmp->GetRMSError());
      
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
    hResVsPt->Write();
    hResVsPtRMS->Write();
    hPtResVsPt->Write();
    hPtResVsPtRMS->Write();
    hYResVsPt->Write();
    hYResVsPtRMS->Write();

    //show
    TString opt = (irapbin>0 ? "same" : "");
    cres->cd(1);
    hResVsPt->Draw(opt.Data());
    cres->cd(2);
    hPtResVsPt->Draw(opt.Data());
    cres->cd(3);
    hYResVsPt->Draw(opt.Data());
  }

  
  cres->cd(1);
  TLegend * legRes = (TLegend*) gPad->BuildLegend(0.6, 0.75, 0.92, 0.98);
  legRes->SetFillStyle(1001); legRes->SetFillColor(kWhite);
  TString resSave = fres->GetName();
  resSave.ReplaceAll(".root",".png");
  cres->Print(resSave.Data());
 
  
  //Save lists to file and print canvases
  efffile->cd();
  lMother->Write();
  lTrue->Write();
  lMotherY->Write();
  lTrueY->Write();

  return;
 
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

void GetGausSigma(TH1D * htmp, Float_t * sigmaAndErr)
{
  //returns sigma and its error from gaussian fit
  if (!htmp || !sigmaAndErr) return;
  if (htmp->GetEntries()<1) {
    sigmaAndErr[0] = 0; sigmaAndErr[1] = 0;
  }
  htmp->Fit("gaus","S0","");
  TF1* fitg = (TF1*) htmp->GetFunction("gaus");
  sigmaAndErr[0] = (fitg ? fitg->GetParameter(2) : 0.0);
  sigmaAndErr[1] = (fitg ? fitg->GetParError(2) : 0.0);
  return;
}
