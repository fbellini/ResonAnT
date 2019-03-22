/******************************************************************
         Author: fbellini@cern.ch - Created on 25.11.2016

Macro dedicated to the checks on TPC and TOF PID QA
using the output produced by the Rsn package monitor framework
e.g. $ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C

Usage:
- 

*******************************************************************/
void checkRsnPIDqa(TString filename = "/Users/fbellini/alice/resonances/kstar_pA5.02TeV/LF_pPb17-21/train1920.root",
		   TString foldername = "RsnOut_cut5",
		   Bool_t savePng = 1);
  
void MakeUpHisto(TH1* histo = NULL, TString titleY = "", Int_t marker = 20, Color_t color = kBlue+2);

Double_t TOFsignal(Double_t *x, Double_t *par)
{
  //Define function to fit TOF signal t-texp
  //as a gaussian + exponential tail
  Double_t norm = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];
  Double_t tail = par[3];
  Double_t a = par[4];
  Double_t b = par[5];
  //Double_t c = par[6];
  if (x[0] <= (tail + mean))
    return norm * TMath::Gaus(x[0], mean, sigma) + a + b * x[0]; /*+ c * x[0] * x[0]*/
  else
    return norm * TMath::Gaus(tail + mean, mean, sigma) * TMath::Exp(-tail * (x[0] - tail - mean) / (sigma * sigma)) + a + b * x[0];
    /*+ c * x[0] * x[0]*/
}

void checkRsnPIDqa(TString filename, TString foldername, Bool_t checkPIDqaTaskOutput, Bool_t savePng)
{
  //Open input file 
  TFile * fin = TFile::Open(filename.Data());
  if (!fin) return 0x0;
  
  //Access output of specific wagon
  TList * list = (TList*) fin->Get(foldername.Data());
  if (!list) return 0x0;

  //Set range for fit
  Float_t RangeFitMomMin = 0.1; //range in momentum where to check the mean and pull
  Float_t RangeFitMomMax = 2.0;
  Int_t xbinFitMin = 0;
  Int_t xbinFitMax = -1;
  Float_t RangeFitNsigmaPIDmin = -2.0; //range in Nsigma where the fit is to be performed
  Float_t RangeFitNsigmaPIDmax = 2.0;

  //Set range for visualisation
  Float_t RangeShowTPC[2] = {0.1, 2.0}; 
  Float_t RangeShowTOF[2] = {0.25, 2.0};
  
  //Histogram names for TPC PID QA --- 
  TString plotTPCpi = "cutQ_bit5.TPC_nsigmaPi_VsPtpc_pTPC_pi";
  TString plotTPCka  = "cutQ_bit5.TPC_nsigmaK_VsPtpc_pTPC_K";
  TString plotTPCpro = "cutQ_bit5.TPC_nsigmaPro_VsPtpc_pTPC_p";
  
  //Histogram names for TPC PID QA --- 
  TString plotTOFpi  = "cutQ_bit5.TOF_nsigmaPi_vsP_p_pi";
  TString plotTOFka  = "cutQ_bit5.TOF_nsigmaK_vsP_p_K";
  TString plotTOFpro = "cutQ_bit5.TOF_nsigmaPro_vsP_p_p";
  
  //--------------------------
  // TPC PID Nsigma
  // fit with simple gaussian
  //--------------------------
  //Gaussian function
  TF1 *fGaus = new TF1("f","gaus", -7.0, 7.0);

  //--- pions
  TH2F * hTPCsigmaPi = (TH2F*)list->FindObject(plotTPCpi.Data());
  hTPCsigmaPi->RebinX(2);
  hTPCsigmaPi->SetTitle("TPC Pions");
  MakeUpHisto(hTPCsigmaPi,"p_{TPC} (GeV/c)", "N#sigma_{TPC}", 1, kBlack, 2);
  hTPCsigmaPi->GetYaxis()->SetRangeUser(-5.1,5.1);
  hTPCsigmaPi->GetXaxis()->SetRangeUser(RangeShowTPC[0], RangeShowTPC[1]);
  xbinFitMin = hTPCsigmaPi->GetXaxis()->FindBin(RangeFitMomMin);
  xbinFitMax = hTPCsigmaPi->GetXaxis()->FindBin(RangeFitMomMax);
  hTPCsigmaPi->FitSlicesY(fGaus, xbinFitMin, xbinFitMax );
  TH1D * hTPCsigmaPi_mean = ((TH1D*)gDirectory->FindObject(Form("%s_1", plotTPCpi.Data())))->Clone("hNsigmaTPCpi_mean");
  TH1D * hTPCsigmaPi_pull = ((TH1D*)gDirectory->FindObject(Form("%s_2", plotTPCpi.Data())))->Clone("hNsigmaTPCpi_pull");
  MakeUpHisto(hTPCsigmaPi_mean, "", "", 1, kBlack, 2);
  MakeUpHisto(hTPCsigmaPi_pull, "", "", 1, kRed+2, 2);

  //--- kaons
  TH2F * hTPCsigmaKa = (TH2F*)list->FindObject(plotTPCka.Data());
  hTPCsigmaKa->RebinX(2);
  hTPCsigmaKa->SetTitle("TPC Kaons");
  hTPCsigmaKa->GetYaxis()->SetRangeUser(-5.1,5.1);
  hTPCsigmaKa->GetXaxis()->SetRangeUser(RangeShowTPC[0], RangeShowTPC[1]);
  hTPCsigmaKa->FitSlicesY(fGaus, xbinFitMin, xbinFitMax );
  MakeUpHisto(hTPCsigmaKa,"p_{TPC} (GeV/c)", "N#sigma_{TPC}", 1, kBlack, 2);
  TH1D * hTPCsigmaKa_mean = ((TH1D*)gDirectory->FindObject(Form("%s_1", plotTPCka.Data())))->Clone("hNsigmaTPCka_mean");
  TH1D * hTPCsigmaKa_pull = ((TH1D*)gDirectory->FindObject(Form("%s_2", plotTPCka.Data())))->Clone("hNsigmaTPCka_pull");
  MakeUpHisto(hTPCsigmaKa_mean, "", "", 1, kBlack, 2);
  MakeUpHisto(hTPCsigmaKa_pull, "", "", 1, kRed+2, 2);

  //--- protons
  TH2F * hTPCsigmaPro = (TH2F*)list->FindObject(plotTPCpro.Data());
  hTPCsigmaPro->RebinX(2);
  hTPCsigmaPro->SetTitle("TPC Protons");
  MakeUpHisto(hTPCsigmaPro,"p_{TPC} (GeV/c)", "N#sigma_{TPC}", 1, kBlack, 2);
  hTPCsigmaPro->GetYaxis()->SetRangeUser(-5.1,5.1);
  hTPCsigmaPro->GetXaxis()->SetRangeUser(RangeShowTPC[0], RangeShowTPC[1]);
  hTPCsigmaPro->FitSlicesY(fGaus, xbinFitMin, xbinFitMax );
  TH1D * hTPCsigmaPro_mean = ((TH1D*)gDirectory->FindObject(Form("%s_1", plotTPCpro.Data())))->Clone("hNsigmaTPCpro_mean");
  TH1D * hTPCsigmaPro_pull = ((TH1D*)gDirectory->FindObject(Form("%s_2", plotTPCpro.Data())))->Clone("hNsigmaTPCpro_pull");
  MakeUpHisto(hTPCsigmaPro_mean, "", "", 1, kBlack, 2);
  MakeUpHisto(hTPCsigmaPro_pull, "", "", 1, kRed+2, 2);

   //--- plot TPC
  TLine *l11=new TLine(RangeShowTPC[0],0.,RangeShowTPC[1],0.); l11->SetLineWidth(1); l11->SetLineStyle(7);
  TLine *l12=new TLine(RangeShowTPC[0],1.,RangeShowTPC[1],1.); l12->SetLineWidth(1); l12->SetLineStyle(7);

  gStyle->SetOptStat(0);
  TCanvas *cPidPerformance4 = new TCanvas("cPIDperformance4","TPC PID",1200,500);
  cPidPerformance4->Divide(3,1);
  cPidPerformance4->cd(1);
  gPad->SetLogz(); gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
  hTPCsigmaPi->DrawCopy("colz");
  hTPCsigmaPi_mean->DrawCopy("same");
  hTPCsigmaPi_pull->DrawCopy("same");
  l11->Draw("same"); l12->Draw("same");

  cPidPerformance4->cd(2);
  gPad->SetLogz(); gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
  hTPCsigmaKa->DrawCopy("colz");
  hTPCsigmaKa_mean->DrawCopy("same");
  hTPCsigmaKa_pull->DrawCopy("same");
  l11->Draw("same"); l12->Draw("same");

  cPidPerformance4->cd(3);
  gPad->SetLogz(); gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
  hTPCsigmaPro->DrawCopy("colz");
  hTPCsigmaPro_mean->DrawCopy("same");
  hTPCsigmaPro_pull->DrawCopy("same");
  l11->Draw("same"); l12->Draw("same");

  TLegend * pidLegTPC = new TLegend(0.15,0.8,0.88,0.88);
  pidLegTPC->SetBorderSize(0); pidLegTPC->SetFillStyle(1001); pidLegTPC->SetFillColor(kWhite);
  pidLegTPC->SetTextSize(0.04); pidLegTPC->SetNColumns(2);
  pidLegTPC->AddEntry(hTPCsigmaPro_mean,"Mean","lp");
  pidLegTPC->AddEntry(hTPCsigmaPro_pull,Form("#sigma, Gaus fit (%2.1f,%2.1f)",RangeFitNsigmaPIDmin,RangeFitNsigmaPIDmax),"lp");
  pidLegTPC->Draw("same");

  if (savePng) cPidPerformance4->SaveAs("RsnQA_TPC_Nsigma.png");
  
  //----------------------------------------------------
  // TOF
  // fit with signal model = gaussian + exponential tail
  //----------------------------------------------------
  //Signal model for TOF signal = gaus + exp tail
  const Int_t npars = 6;
  TF1 *fSignalModel = new TF1("fSignalModel", TOFsignal, -7.0, 7.0, npars);
  fSignalModel->SetTitle("TOF Signal");
  fSignalModel->SetParameter(0, 1.);
  fSignalModel->SetParameter(1, 0.);
  fSignalModel->SetParLimits(1, -2., 1.);
  fSignalModel->SetParameter(2, 1.);
  fSignalModel->SetParLimits(2, 0.5, 2.);
  fSignalModel->SetParameter(3, 1.);
  fSignalModel->SetParLimits(3, 0.5, 1.5);
  fSignalModel->SetParameter(4, 1.);
  fSignalModel->SetParLimits(4, 0., 1.e8);
  fSignalModel->SetParameter(5, 0.);
  fSignalModel->SetParLimits(5, -10., 10.);
  fSignalModel->SetNpx(2000);
  fSignalModel->SetParNames("Norm", "Mean", "Sigma", "Tail", "Shift", "Slope"/*, "Square"*/);
  fSignalModel->SetLineColor(kRed+1);

  //results
  TObjArray *results[3];
  for(Int_t i = 0; i < 3; i++){
    results[i] = new TObjArray(10);
  }
  TH1D * par[3][npars];
  //--- pions
  TH2F * hTOFsigmaPi = (TH2F*)list->FindObject(plotTOFpi.Data());
  hTOFsigmaPi->SetTitle("TOF Pions");
  hTOFsigmaPi->RebinX(2);
  MakeUpHisto(hTOFsigmaPi,"p (GeV/c)", "N#sigma_{TOF}", 1, kBlack, 2);
  hTOFsigmaPi->GetYaxis()->SetRangeUser(-5.1,5.1);
  hTOFsigmaPi->GetXaxis()->SetRangeUser(RangeShowTOF[0], RangeShowTOF[1]);
  fSignalModel->SetParLimits(4, 0., hTOFsigmaPi->GetMaximum()*0.5);
  fSignalModel->SetParLimits(0, 0., hTOFsigmaPi->GetMaximum()*1.2);
  hTOFsigmaPi->FitSlicesY(fSignalModel, xbinFitMin, xbinFitMax, 0, "QR", results[0] );
  for(Int_t cc = 0; cc < npars ; cc++) {
    par[0][cc] = (TH1D*)gDirectory->FindObject(Form("%s_%i", plotTOFpi.Data(), cc));
  }
  MakeUpHisto(par[0][1], "", "", 1, kBlue, 2);
  MakeUpHisto(par[0][2], "", "", 1, kMagenta+2, 2);
  
  //--- KAONS
  TH2F * hTOFsigmaKa = (TH2F*)list->FindObject(plotTOFka.Data());
  hTOFsigmaKa->SetTitle("TOF Kaons");
  hTOFsigmaKa->RebinX(2);
  MakeUpHisto(hTOFsigmaKa,"p (GeV/c)", "N#sigma_{TOF}", 1, kBlack, 2);
  hTOFsigmaKa->GetYaxis()->SetRangeUser(-5.1,5.1);
  hTOFsigmaKa->GetXaxis()->SetRangeUser(RangeShowTOF[0], RangeShowTOF[1]);
  fSignalModel->SetParLimits(4, 0., hTOFsigmaKa->GetMaximum()*0.5);
  fSignalModel->SetParLimits(0, 0., hTOFsigmaKa->GetMaximum()*1.2);
  hTOFsigmaKa->FitSlicesY(fSignalModel, xbinFitMin, xbinFitMax, 0, "QR", results[0] );
  for(Int_t cc = 0; cc < npars ; cc++) {
    par[1][cc] = (TH1D*)gDirectory->FindObject(Form("%s_%i", plotTOFka.Data(), cc));
  }
  MakeUpHisto(par[1][1], "", "", 1, kBlue, 2);
  MakeUpHisto(par[1][2], "", "", 1, kMagenta+2, 2);

  //--- protons
  TH2F * hTOFsigmaPro = (TH2F*)list->FindObject(plotTOFpro.Data());
  hTOFsigmaPro->SetTitle("TOF Protons");
  hTOFsigmaPro->RebinX(2);
  MakeUpHisto(hTOFsigmaPro,"p (GeV/c)", "N#sigma_{TOF}", 1, kBlack, 2);
  hTOFsigmaPro->GetYaxis()->SetRangeUser(-5.1,5.1);
  hTOFsigmaPro->GetXaxis()->SetRangeUser(RangeShowTOF[0], RangeShowTOF[1]);
  fSignalModel->SetParLimits(4, 0., hTOFsigmaPro->GetMaximum()*0.5);
  fSignalModel->SetParLimits(0, 0., hTOFsigmaPro->GetMaximum()*1.2);
  hTOFsigmaPro->FitSlicesY(fSignalModel, xbinFitMin, xbinFitMax, 0, "QR", results[0] );
  for(Int_t cc = 0; cc < npars ; cc++) {
    par[2][cc] = (TH1D*)gDirectory->FindObject(Form("%s_%i", plotTOFpro.Data(), cc));
  }
  MakeUpHisto(par[2][1], "", "", 1, kBlue, 2);
  MakeUpHisto(par[2][2], "", "", 1, kMagenta+2, 2);

  //--- plot TOF
  gStyle->SetOptStat(0);
  TCanvas *cPidPerformance3 = new TCanvas("cPidPerformance3","TOF PID performance",1200,500);
  cPidPerformance3->Divide(3,1);
  cPidPerformance3->cd(1);
  gPad->SetLogz(); gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
  hTOFsigmaPi->DrawCopy("colz");
  if(par[0][1]) par[0][1]->DrawCopy("same");
  if(par[0][2]) par[0][2]->DrawCopy("same");
  l11->Draw("same"); l12->Draw("same");

  cPidPerformance3->cd(2);
  gPad->SetLogz(); gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
  hTOFsigmaKa->DrawCopy("colz");
  if(par[1][1]) par[1][1]->DrawCopy("same");
  if(par[1][2]) par[1][2]->DrawCopy("same");
  l11->Draw("same"); l12->Draw("same");

  cPidPerformance3->cd(3);
  gPad->SetLogz(); gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
  hTOFsigmaPro->DrawCopy("colz");
  if(par[2][1]) par[2][1]->DrawCopy("same");
  if(par[2][2]) par[2][2]->DrawCopy("same");
  l11->Draw("same"); l12->Draw("same");

  TLegend * pidLegTOF = new TLegend(0.15,0.8,0.88,0.88);
  pidLegTOF->SetBorderSize(0); pidLegTOF->SetFillStyle(1001); pidLegTOF->SetFillColor(kWhite);
  pidLegTOF->SetTextSize(0.04); pidLegTOF->SetNColumns(2);
  pidLegTOF->AddEntry(par[0][1],"Mean","lp");
  pidLegTOF->AddEntry(par[0][2], Form("#sigma, Gaus+Tail fit (%2.1f,%2.1f)",RangeFitNsigmaPIDmin, RangeFitNsigmaPIDmax),"lp");
  pidLegTOF->Draw("same");
  
  if (savePng) cPidPerformance3->Print("RsnQA_TOF_Nsigma.png");
  return;
}

//-------------------------------------------------------------------------------------
void MakeUpHisto(TH1* histo, TString titleX, TString titleY, Int_t marker, Color_t color, Int_t lineWidth)
{
  if (!histo) return;
  histo->SetMarkerStyle(marker);
  histo->SetMarkerSize(0.7);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->SetLineWidth(lineWidth);
  histo->SetFillColor(kWhite);
  histo->SetFillStyle(0);
  histo->GetYaxis()->SetNdivisions(515);
  if (!titleX.IsNull()) histo->GetXaxis()->SetTitle(titleX.Data());
  if (!titleY.IsNull()) histo->GetYaxis()->SetTitle(titleY.Data());
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.045);
  histo->GetXaxis()->SetLabelOffset(-0.003);
  histo->GetXaxis()->SetTitleOffset(1.1);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.045);
  histo->GetYaxis()->SetLabelOffset(0.007);
  histo->GetYaxis()->SetTitleOffset(1.2);
  return;
}


/*
//--------------------------------- NSigma PID from PIDqa ----------------------------------//
TH2F * hSigmaPiT0 = 0x0;
TH2F * hSigmaKaT0 = 0x0;
TH2F * hSigmaProT0 = 0x0;

//results
TObjArray *resultsT0[3];
for(Int_t i = 0; i < 3; i++) resultsT0[i] = new TObjArray(10);
TH1D * parT0[3][npars];

//------- PIONS ------//
hSigmaPiT0=(TH2F*)tofPidListT0->FindObject("hNsigmaP_TOF_pion");
//hSigmaPiT0->SetName("hSigmaPiT0");
hSigmaPiT0->GetYaxis()->SetRangeUser(-5.,5.);
hSigmaPiT0->GetXaxis()->SetRangeUser(0.2, 5.);
//fit with simple gaussian
hSigmaPiT0->FitSlicesY(f);
TH1D * hSigmaPiT0_mean = (TH1D*)gDirectory->Get("hNsigmaP_TOF_pion_1")->Clone("hNsigmaP_TOF_pion_mean");
TH1D * hSigmaPiT0_pull = (TH1D*)gDirectory->Get("hNsigmaP_TOF_pion_2")->Clone("hNsigmaP_TOF_pion_pull");
MakeUpHisto(hSigmaPiT0_mean, "", 1, kBlack, 2);
MakeUpHisto(hSigmaPiT0_pull, "", 1, kRed, 2);

//fit with signal model
if (fitSignalModel) {
fSignalModel->SetParLimits(4, 0., hSigmaPiT0->GetMaximum()*0.5);
fSignalModel->SetParLimits(0, 0., hSigmaPiT0->GetMaximum()*1.2);
hSigmaPiT0->FitSlicesY(fSignalModel, 0, -1, 0, "QR", resultsT0[0] );
for(Int_t cc = 0; cc < npars; cc++) {
parT0[0][cc] = (TH1D*)gDirectory->Get(Form("hNsigmaP_TOF_pion_%i", cc));
}
MakeUpHisto(parT0[0][1], "", 1, kBlue, 2);
MakeUpHisto(parT0[0][2], "", 1, kMagenta+2, 2);
}

//------- KAONS ------//
hSigmaKaT0=(TH2F*)tofPidListT0->FindObject("hNsigmaP_TOF_kaon");
hSigmaKaT0->GetYaxis()->SetRangeUser(-5.,5.);
hSigmaKaT0->GetXaxis()->SetRangeUser(0.2, 5.);
//fit with simple gausssian
hSigmaKaT0->FitSlicesY(f);
TH1D * hSigmaKaT0_mean = (TH1D*)gDirectory->Get("hNsigmaP_TOF_kaon_1")->Clone("hNsigmaP_TOF_kaon_mean");
TH1D * hSigmaKaT0_pull = (TH1D*)gDirectory->Get("hNsigmaP_TOF_kaon_2")->Clone("hNsigmaP_TOF_kaon_pull");
MakeUpHisto(hSigmaKaT0_mean, "", 1, kBlack, 2);
MakeUpHisto(hSigmaKaT0_pull, "", 1, kRed, 2);

//fit with signal model
if (fitSignalModel) {
fSignalModel->SetParLimits(4, 0., hSigmaKaT0->GetMaximum()*0.5);
fSignalModel->SetParLimits(0, 0., hSigmaKaT0->GetMaximum()*1.2);
hSigmaKaT0->FitSlicesY(fSignalModel, 0, -1, 0, "QR", resultsT0[1] );
for(Int_t cc = 0; cc < npars; cc++) {
parT0[1][cc] = (TH1D*)gDirectory->Get(Form("hNsigmaP_TOF_kaon_%i", cc));
}
MakeUpHisto(parT0[1][1], "", 1, kBlue, 2);
MakeUpHisto(parT0[1][2], "", 1, kMagenta+2, 2);
}

//------- PROTONS ------//
hSigmaProT0 = (TH2F*)tofPidListT0->FindObject("hNsigmaP_TOF_proton");
hSigmaProT0->GetYaxis()->SetRangeUser(-5.,5.);
hSigmaProT0->GetXaxis()->SetRangeUser(0.2, 5.);
//fit with simple gausssian
hSigmaProT0->FitSlicesY(f);
TH1D * hSigmaProT0_mean = (TH1D*)gDirectory->Get("hNsigmaP_TOF_proton_1")->Clone("hNsigmaP_TOF_proton_mean");
TH1D * hSigmaProT0_pull = (TH1D*)gDirectory->Get("hNsigmaP_TOF_proton_2")->Clone("hNsigmaP_TOF_proton_pull");
MakeUpHisto(hSigmaProT0_mean, "", 1, kBlack, 2);
MakeUpHisto(hSigmaProT0_pull, "", 1, kRed, 2);

//fit with signal model
if (fitSignalModel) {
fSignalModel->SetParLimits(4, 0., hSigmaProT0->GetMaximum()*0.5);
fSignalModel->SetParLimits(0, 0., hSigmaProT0->GetMaximum()*1.2);
hSigmaProT0->FitSlicesY(fSignalModel, 0, -1, 0, "QR", resultsT0[2] );
for(Int_t cc = 0; cc < npars; cc++) {
parT0[2][cc] = (TH1D*)gDirectory->Get(Form("hNsigmaP_TOF_proton_%i", cc));
}
MakeUpHisto(parT0[2][1], "", 1, kBlue, 2);
MakeUpHisto(parT0[2][2], "", 1, kMagenta+2, 2);
}

//Show in canvas
TLine *l1=new TLine(0.,0.,5.,0.);
TLine *l2=new TLine(0.,1.,5.,1.);
TCanvas *cPidPerformance3T0 = new TCanvas("cPidPerformance3T0","PID performance from PIDqa - N_{#sigma}^{TOF} with StartTime",1200,500);
cPidPerformance3T0->Divide(3,1);
cPidPerformance3T0->cd(1);
gPad->SetLogz(); gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
hSigmaPiT0->DrawCopy("colz");
hSigmaPiT0_mean->DrawCopy("same");
hSigmaPiT0_pull->DrawCopy("same");
l1->Draw("same");
l2->Draw("same");

cPidPerformance3T0->cd(2);
gPad->SetLogz(); gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
hSigmaKaT0->DrawCopy("colz");
hSigmaKaT0_mean->DrawCopy("same");
hSigmaKaT0_pull->DrawCopy("same");
l1->Draw("same");
l2->Draw("same");

cPidPerformance3T0->cd(3);
gPad->SetLogz(); gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
hSigmaProT0->DrawCopy("colz");
hSigmaProT0_mean->DrawCopy("same");
hSigmaProT0_pull->DrawCopy("same");
l1->Draw("same");
l2->Draw("same");

if (fitSignalModel) {
for (Int_t jj = 0; jj<3; jj++){
cPidPerformance3T0->cd(jj+1);
if(parT0[jj][1]) parT0[jj][1]->DrawCopy("same");
if(parT0[jj][2]) parT0[jj][2]->DrawCopy("same");
}
}

cPidPerformance3T0->cd(1);
TLegend * pidLegT0 = new TLegend(0.15,0.76,0.88,0.88);
pidLegT0->SetBorderSize(0); pidLegT0->SetFillStyle(1001); pidLegT0->SetFillColor(kWhite);
pidLegT0->SetTextSize(0.04);
pidLegT0->SetNColumns(2);
pidLegT0->AddEntry(hSigmaPiT0_mean,"Mean","lp");
pidLegT0->AddEntry(hSigmaPiT0_pull,Form("#sigma, Gaus fit (%2.1f,%2.1f)",RangeFitNsigmaPIDmin,RangeFitNsigmaPIDmax),"lp");
if (fitSignalModel && parT0[0][1] && parT0[0][2]) {
pidLegT0->AddEntry(parT0[0][1],"Mean","lp");
pidLegT0->AddEntry(parT0[0][2],Form("#sigma, Gaus+Tail fit (%2.1f,%2.1f)", ModelRangeFitNsigmaPIDmin, ModelRangeFitNsigmaPIDmax),"lp");
}
pidLegT0->Draw("same");

if (savePng) cPidPerformance3T0->Print(Form("%s/%i%s_PID_sigmasStartTime.png", plotDir.Data(), runNumber, dirsuffix.Data()));
}
*/


//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
TH1 * projectCutPlot(TString filename="/Users/fbellini/alice/resonances/kstar_pA5.02TeV/LF_pPb17-21/train1920.root",
		     TString foldername="RsnOut_cut5",
		     TString plot="dcaxy",
		     TString projopt="y",
		     TString newname = "cut5")
{
  TFile * fin = TFile::Open(filename.Data());
  if (!fin) return 0x0;
  TList * list = (TList*) fin->Get(foldername.Data());
  if (!list) return 0x0;
  TH2F * hist = (TH2F*) list->FindObject(getPlotName(plot.Data(),5));
  if (!hist) return;

  TH1F * proj = 0x0;
  if (projopt.Contains("x")) proj = (TH1F*) hist->ProjectionX(Form("%s_px",newname.Data()));
  else proj = (TH1F*) hist->ProjectionY(Form("%s_py",newname.Data()));

  return proj;
}

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
TString getPlotName(TString whatToPlot="dcaxy", Int_t bit=5)
{

  TString speciename = "cutQ";
  TString suffix;
  /*
    if (what.Contains("mom")) suffix = "P_p";
    if (what.Contains("pt")) suffix = "Pt_pt"; // at: 0x7fdbd2cfaa00
    if (what.Contains("eta")) suffix = "Eta_eta"; // at: 0x7fdbd2cfadd0
    if (what.Contains("phi")) suffix = "Phi_phi"; // at: 0x7fdbd2cfb1a0
    if (what.Contains("tpcout")) suffix = "PhiOuterTPC_phiOuterTPC"; // at: 0x7fdbd2cfb570
  */
  if (what.Contains("phi")) suffix = "PhiVsPt_pt_phi"; // at: 0x7fdbd2cfb940
  if (what.Contains("dcaxy")) suffix = "DCAxyVsPt_pt_DCAxy"; // at: 0x7fdbd4800000
  if (what.Contains("dcaz")) suffix = "DCAzVsP_p_DCAz"; // at: 0x7fdbd48003f0
  if (what.Contains("itscls")) suffix = "ITSclsVsPt_pt_ITScls"; // at: 0x7fdbd48007e0
  //    if (what.Contains("tpccls")) suffix = "TPCclsVsPt_pt_TPCcls"; // at: 0x7fdbd4800bd0
  if (what.Contains("tpccls")) suffix = "TPCclsVsPtpc_pTPC_TPCcls"; // at: 0x7fdbd4800fc0
  if (what.Contains("crossedrows")) suffix = "TPCcrossedRowsVsPtpc_pTPC_TPCcrossedRows"; // at: 0x7fdbd48013b0
  //    if (what.Contains("")) suffix = "TPCcrossedRows2FclsVsPtpc_pTPC_TPCcrossedRows2Fcls"; // at: 0x7fdbd48017a0
  if (what.Contains("itschi2")) suffix = "ITSchi2VsPt_pt_ITSchi2"; // at: 0x7fdbd4801b90
  //    if (what.Contains("tpcchi2")) suffix = "TPCchi2VsPt_pt_TPCchi2"; // at: 0x7fdbd4802350
  if (what.Contains("tpcchi2")) suffix = "TPCchi2VsPtpc_pTPC_TPCchi2"; // at: 0x7fdbd4802740
  if (what.Contains("dedx")) suffix = "dEdx_VsPtpc_pTPC_sTPC"; // at: 0x7fdbd4802e60
  if (what.Contains("nstpcpi")) suffix = "TPC_nsigmaPi_VsPtpc_pTPC_pi"; // at: 0x7fdbd4803250
  if (what.Contains("nstpck")) suffix = "TPC_nsigmaK_VsPtpc_pTPC_K"; // at: 0x7fdbd4803640
  if (what.Contains("nstpcpro")) suffix = "TPC_nsigmaPro_VsPtpc_pTPC_p"; // at: 0x7fdbd4803a30
  if (what.Contains("allsigmapi")) suffix = "TOFnsigmaPi_TPCnsigmaPi_pi_pi"; // at: 0x7fdbd4803e20
  if (what.Contains("allsigmak")) suffix = "TOFnsigmaK_TPCnsigmaK_K_K"; // at: 0x7fdbd4804240
  if (what.Contains("allsigmapro")) suffix = "TOFnsigmaP_TPCnsigmaP_p_p"; // at: 0x7fdbd48046b0
  if (what.Contains("nstofpi")) suffix = "TOF_nsigmaPi_vsP_p_pi"; // at: 0x7fdbd4804b20
  if (what.Contains("nstofk")) suffix = "TOF_nsigmaK_vsP_p_K"; // at: 0x7fdbd4804f60
  if (what.Contains("nstofpro")) suffix = "TOF_nsigmaPro_vsP_p_p"; // at: 0x7fdbd48053a0
  if (what.Contains("deltapi")) suffix = "TOF_deltaPi_vsP_p_Dpi"; // at: 0x7fdbd48057e0
  if (what.Contains("deltak")) suffix = "TOF_deltaK_vsP_p_DK"; // at: 0x7fdbd4805c20
  if (what.Contains("deltapro")) suffix = "TOF_deltaPro_vsP_p_Dp"; // at: 0x7fdbd4806060

  //cutK24_2.0sigma.TOF_deltaPro_vsP_p_Dp
  TString name = Form("%s_bit%i.%s", speciename.Data(), bit, suffix.Data());
  return name;

}
