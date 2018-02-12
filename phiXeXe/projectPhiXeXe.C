//
// Splits a TH2F using the myHistSplig_sparse_IM_PT_CENT class
//
#include "/Users/fbellini/alice/macros/ResonAnT/projectorTH3_InvMass_Centrality_Pt.C"
#include "TFile.h"


void FillListForCentrality(TList * listIn, TList * listOut, Int_t icentBin);

void projectPhiXeXe( const char *nameData = "20180123_RsnOut.root",
		     TString outName  = "phi",
		     const char *listName = "RsnOut_default",
		     TString binning  = "B2",
		     Bool_t isMC = 0,
		     Bool_t saveSeparatelyByCent = 0)
{
  //graphics for display
  gStyle->SetOptStat("1111");
  gStyle->SetTextFont(42);
  
  //generalia
  const Int_t kNhistosData = 4;
  outName.Append(binning.Data());

  //------------------------------
  //define binning
  //------------------------------
  Double_t centA[] = {0.0, 30.0, 60.0, 90.0};   
  Double_t centB[]  = {0.0, 20.0, 40.0, 60.0, 80.0};
  Double_t   pt1[] = {0.5, 1.00, 1.50, 2.00, 2.50, 3.00, 3.5, 4.00, 4.5, 5.0, 6.0, 7.0, 10.00};
  Double_t   pt2[] = {0.3, 0.5,  1.00, 1.50, 2.00, 2.50, 3.00, 3.5, 4.00, 4.5, 5.0, 7.0, 10.00};
  Double_t   pt3[] = {0.3, 0.5,  0.7, 0.9, 1.10, 1.30, 1.50, 2.00, 3.00, 4.00, 5.0, 7.0, 10.00};

  Int_t   npt  = 0;
  Int_t   ncent  = 0;
  TAxis *ptbins = 0;
  TAxis *centbins = 0;

  if (binning.Contains("1")) {
    npt = sizeof(pt1) / sizeof(pt1[0]) - 1;
    ptbins = new TAxis(npt, pt1);
  } else if (binning.Contains("2")) {
    npt = sizeof(pt2) / sizeof(pt2[0]) - 1;
    ptbins = new TAxis(npt, pt2);
  } else if (binning.Contains("3")) {
    npt = sizeof(pt3) / sizeof(pt3[0]) - 1;
    ptbins = new TAxis(npt, pt3);
  }

  if (binning.Contains("A")) {
    ncent = sizeof(centA) / sizeof(centA[0]) - 1;
    centbins = new TAxis(ncent, centA);
  }  else
    if (binning.Contains("B")) {
      ncent = sizeof(centB) / sizeof(centB[0]) - 1;
      centbins = new TAxis(ncent, centB);
    }

  //------------------------------
  // open input file
  //------------------------------
  TFile *fileData = TFile::Open(nameData);
  if (!fileData || !fileData->IsOpen()) { Printf("Invalid file passed as input. Doing nothing."); return;}
  
  TList *listData = listData = (TList*)fileData->Get(listName);
  if (!listData) { Printf("Invalid list passed as input. Doing nothing."); return;}
    
  //get input histograms
  TH3F* hInput[kNhistosData] = {0,0,0,0}; 
  hInput[0] = (TH3F*)listData->FindObject("PhiXeXeData_Unlike_tpc3s");
  hInput[1] = (TH3F*)listData->FindObject("PhiXeXeData_LikePP_tpc3s");
  hInput[2] = (TH3F*)listData->FindObject("PhiXeXeData_LikeMM_tpc3s");
  hInput[3] = (TH3F*)listData->FindObject("PhiXeXeData_Mixing_tpc3s");

  if (!hInput[0] || !hInput[1] || !hInput[2] || !hInput[3]) { Printf("Invalid histogram requested as input. Doing nothing."); return;}

  // rename
  hInput[0]->SetName("hUnlikePM");
  hInput[1]->SetName("hLikePP");
  hInput[2]->SetName("hLikeMM");
  hInput[3]->SetName("hMixingPM");

  
  //------------------------------
  //output
  //------------------------------
  gSystem->Exec(Form("mkdir %s", outName.Data()));
  TFile *fileOut = new TFile(Form("%s/proj_%s", outName.Data(), nameData),"RECREATE");
  //TList *out = new TList();
  fileOut->cd();
  ptbins->Write("ptbins");
  centbins->Write("centbins");

  //output lists - define once and then clear them
  TList * lUnlikePM = new TList();  lUnlikePM->SetName("UnlikePM");
  TList * lLikePP = new TList();    lLikePP->SetName("LikePP");
  TList * lLikeMM = new TList();    lLikeMM->SetName("LikeMM");
  TList * lMixingPM  = new TList();  lMixingPM->SetName("MixingPM");

  //------------------------------
  // project histograms
  //------------------------------
  projectorTH3_InvMass_Centrality_Pt projector;
  TList * out = lUnlikePM; // support list
  
  // loop on inputs
  for (Int_t i = 0; i < kNhistosData; i++) {
    if (!hInput[i]) continue;
    projector.SetPrefix(hInput[i]->GetName());
    
    if (i == 1) out = lLikePP;
    if (i == 2) out = lLikeMM;
    if (i == 3) out = lMixingPM;
    
    if (binning.Contains("A1")) projector.MultiProjPtCent(npt, pt1, ncent, centA,  hInput[i], out);
    else if (binning.Contains("B1")) projector.MultiProjPtCent(npt, pt1, ncent, centB,  hInput[i], out);
    else if (binning.Contains("A2")) projector.MultiProjPtCent(npt, pt2, ncent, centA,  hInput[i], out);
    else if (binning.Contains("B2")) projector.MultiProjPtCent(npt, pt2, ncent, centB,  hInput[i], out);
    else if (binning.Contains("A3")) projector.MultiProjPtCent(npt, pt3, ncent, centA,  hInput[i], out);
    else if (binning.Contains("B3")) projector.MultiProjPtCent(npt, pt3, ncent, centB,  hInput[i], out);
  }
  Printf(":::: Projected according to binning strategy %s", binning.Data());

  for (Int_t ih = 0; ih<lUnlikePM->GetEntries();ih++){
    ((TH1F*)lUnlikePM->At(ih))->SetLineColor(kBlack);
    ((TH1F*)lUnlikePM->At(ih))->SetMarkerColor(kBlack);
    ((TH1F*)lUnlikePM->At(ih))->SetMarkerStyle(20);
    ((TH1F*)lUnlikePM->At(ih))->SetMarkerSize(0.8); 
    ((TH1F*)lUnlikePM->At(ih))->SetLineWidth(2);
    ((TH1F*)lUnlikePM->At(ih))->SetFillColor(kGray);
    ((TH1F*)lUnlikePM->At(ih))->SetFillStyle(3002);
    ((TH1F*)lUnlikePM->At(ih))->SetDrawOption("b");
  }
  for (Int_t ih = 0; ih<lMixingPM->GetEntries();ih++){  
    ((TH1F*)lMixingPM->At(ih))->SetLineColor(kMagenta-2);
    ((TH1F*)lMixingPM->At(ih))->SetMarkerColor(kMagenta-2);
    ((TH1F*)lMixingPM->At(ih))->SetMarkerStyle(1);
    ((TH1F*)lMixingPM->At(ih))->SetLineWidth(2);
  }
  for (Int_t ih = 0; ih<lLikePP->GetEntries();ih++){
    ((TH1F*)lLikePP->At(ih))->SetLineColor(kRed);
    ((TH1F*)lLikePP->At(ih))->SetMarkerColor(kRed);
    ((TH1F*)lLikePP->At(ih))->SetMarkerStyle(1);
    ((TH1F*)lLikePP->At(ih))->SetLineWidth(2);
  }
  for (Int_t ih = 0; ih<lLikeMM->GetEntries();ih++){
    ((TH1F*)lLikeMM->At(ih))->SetLineColor(kBlue);
    ((TH1F*)lLikeMM->At(ih))->SetMarkerColor(kBlue);
    ((TH1F*)lLikeMM->At(ih))->SetMarkerStyle(1);
    ((TH1F*)lLikeMM->At(ih))->SetLineWidth(2);
  }          
  fileOut->cd();
  lUnlikePM->Write(hInput[0]->GetName(), TObject::kSingleKey);
  lLikePP->Write(hInput[1]->GetName(), TObject::kSingleKey);
  lLikeMM->Write(hInput[2]->GetName(), TObject::kSingleKey);
  lMixingPM->Write(hInput[3]->GetName(), TObject::kSingleKey);
  fileOut->Close();
  
  Printf(":::: Output written in %s", fileOut->GetName());
  return;
   
}




 