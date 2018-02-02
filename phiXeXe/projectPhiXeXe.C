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
  
  //output
  gSystem->Exec(Form("mkdir %s", outName.Data()));
  TFile *dummy = new TFile(Form("%s/proj_%s", outName.Data(), nameData),"RECREATE");
  TList *out = new TList();
  dummy->cd();
  ptbins->Write("ptbins");
  centbins->Write("centbins");

  //output lists - define once and then clear them
  TList lUnlikePM;
  lUnlikePM.SetName(Form("%s_Unlike", (!isMC ? "Data" : "MC")));
  TList lLikePP;
  lLikePP.SetName(Form("%s_LikePP", (!isMC ? "Data" : "MC")));
  TList lLikeMM;
  lLikeMM.SetName(Form("%s_LikeMM", (!isMC ? "Data" : "MC")));
  TList lMixingPM;
  lMixingPM.SetName(Form("%s_Mixing", (!isMC ? "Data" : "MC")));
  
  // open input file
  TFile *fileData;
  TList *listData;
  fileData = TFile::Open(nameData);
  if (fileData && fileData->IsOpen()) listData = (TList*)fileData->Get(listName);
  
  //get input histograms
  TH3F* hInput[kNhistosData]= {0,0,0,0};
  if (listData) {
    hInput[ 0] = (TH3F*)listData->FindObject("PhiXeXeData_Unlike_tpc3s");
    hInput[ 1] = (TH3F*)listData->FindObject("PhiXeXeData_Mixing_tpc3s");
    hInput[ 2] = (TH3F*)listData->FindObject("PhiXeXeData_LikePP_tpc3s");
    hInput[ 3] = (TH3F*)listData->FindObject("PhiXeXeData_LikeMM_tpc3s");
    
    // rename
    hInput[ 0]->SetName("Data_UnlikePM");
    hInput[ 1]->SetName("Data_MixingPM");
    hInput[ 2]->SetName("Data_LikePP");
    hInput[ 3]->SetName("Data_LikeMM");
  }


  
  // define binning
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
    ptbins = new TAxis(ncent, centA);
  }  else
    if (binning.Contains("B")) {
      ncent = sizeof(centB) / sizeof(centB[0]) - 1;
      ptbins = new TAxis(ncent, centB);
    }
  
  // split histograms
  projectorTH3_InvMass_Centrality_Pt projector;
   
  // loop on inputs
  for (Int_t i = 0; i < kNhistosData; i++) {
    if (!hInput[i]) continue;
    // projector.SetPrefix(hInput[i]->GetName());
    // projector.MultiProj(npt, pt, hInput[i], out);
    projector.SetPrefix(hInput[i]->GetName());
    if (binning.Contains("A1")) projector.MultiProjPtCent(npt, pt1, ncent, centA,  hInput[i], out);
    else if (binning.Contains("B1")) projector.MultiProjPtCent(npt, pt1, ncent, centB,  hInput[i], out);
    else if (binning.Contains("A2")) projector.MultiProjPtCent(npt, pt2, ncent, centA,  hInput[i], out);
    else if (binning.Contains("B2")) projector.MultiProjPtCent(npt, pt2, ncent, centB,  hInput[i], out);
    else if (binning.Contains("A3")) projector.MultiProjPtCent(npt, pt3, ncent, centA,  hInput[i], out);
    else if (binning.Contains("B3")) projector.MultiProjPtCent(npt, pt3, ncent, centB,  hInput[i], out);
  }

  //save separate files for each centrality bin
  if (saveSeparatelyByCent) {
    for (Int_t icentbin=0; icentbin<ncent;icentbin++){
      lUnlikePM.Clear();
      lMixingPM.Clear();
      lLikePP.Clear();
      lLikeMM.Clear();
      FillListForCentrality(out,&lUnlikePM,icentbin);
      FillListForCentrality(out,&lMixingPM,icentbin);
      FillListForCentrality(out,&lLikePP,icentbin);
      FillListForCentrality(out,&lLikeMM,icentbin);

      //output file
      TFile *fout = TFile::Open(Form("%s/proj_%i.root", outName.Data(), icentbin), "RECREATE");
      // save into a file
      fout->cd();   
      ptbins->Write("ptbins");
      centbins->Write("centbins");
      //Printf("=========################ DEBUG %i",lUnlikePM.GetEntries());
      //make-up
      for (Int_t ih = 0; ih<lUnlikePM.GetEntries();ih++){
	((TH1F*)lUnlikePM.At(ih))->SetLineColor(kBlack);
	((TH1F*)lUnlikePM.At(ih))->SetMarkerColor(kBlack);
	((TH1F*)lUnlikePM.At(ih))->SetMarkerStyle(20);
	((TH1F*)lUnlikePM.At(ih))->SetMarkerSize(0.8); 
	((TH1F*)lUnlikePM.At(ih))->SetLineWidth(2);
	((TH1F*)lUnlikePM.At(ih))->SetFillColor(kGray);
	((TH1F*)lUnlikePM.At(ih))->SetFillStyle(3002);
	((TH1F*)lUnlikePM.At(ih))->SetDrawOption("b");
      }
      for (Int_t ih = 0; ih<lMixingPM.GetEntries();ih++){  
	((TH1F*)lMixingPM.At(ih))->SetLineColor(kMagenta-2);
	((TH1F*)lMixingPM.At(ih))->SetMarkerColor(kMagenta-2);
	((TH1F*)lMixingPM.At(ih))->SetMarkerStyle(1);
	((TH1F*)lMixingPM.At(ih))->SetLineWidth(2);
      }
      for (Int_t ih = 0; ih<lLikePP.GetEntries();ih++){
	((TH1F*)lLikePP.At(ih))->SetLineColor(kRed);
	((TH1F*)lLikePP.At(ih))->SetMarkerColor(kRed);
	((TH1F*)lLikePP.At(ih))->SetMarkerStyle(1);
	((TH1F*)lLikePP.At(ih))->SetLineWidth(2);
      }
      for (Int_t ih = 0; ih<lLikeMM.GetEntries();ih++){
	((TH1F*)lLikeMM.At(ih))->SetLineColor(kBlue);
	((TH1F*)lLikeMM.At(ih))->SetMarkerColor(kBlue);
	((TH1F*)lLikeMM.At(ih))->SetMarkerStyle(1);
	((TH1F*)lLikeMM.At(ih))->SetLineWidth(2);
      }          
      lUnlikePM.Write(hInput[0]->GetName(), TObject::kSingleKey);
      lMixingPM.Write(hInput[1]->GetName(), TObject::kSingleKey);
      lLikePP.Write(hInput[2]->GetName(), TObject::kSingleKey);
      lLikeMM.Write(hInput[3]->GetName(), TObject::kSingleKey);
      fout->Close();
    }
  }
  // write the full list on dummy file
  dummy->cd();
  out->Write();
  dummy->Close();
  return;
   
}
//---------------------------------------------
void FillListForCentrality(TList * listIn, TList * listOut, Int_t icentBin)
{
  //selects all type of histos per centrality bin and saves them in separate lists  
  TString listName = listOut->GetName();
  TString type;
  
  if (listName.Contains("Unlike")) type = "Unlike";
  if (listName.Contains("LikeMM")) type = "LikeMM";
  if (listName.Contains("LikePP")) type = "LikePP";
  if (listName.Contains("Mixing")) type = "Mixing";
  
  Int_t nhisto = listIn->GetEntries();
  for (Int_t ih=0;ih<nhisto;ih++){
    TH1F * dummy = (TH1F *)listIn->At(ih);
    TString hName = dummy->GetName();
    if (hName.Contains(type.Data())){
      if (hName.Contains(Form("centBin%02d",icentBin))){
	listOut->AddLast(dummy);
	Printf("Histo found with name: %s",hName.Data());
      }
    }
  } 
  
  Printf("============== Output");
  listOut->ls();
  return;
}






 
