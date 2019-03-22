//
// Splits a TH2F using the myHistSplig_sparse_IM_PT_CENT class
//
Bool_t kComputeEfficiency=kFALSE;
const Int_t kNhistosData=8;
TString macroDir ="$ASD/ResonAnT/"; //"/Users/bellini/alice/macro/kstar";

void projectKStar_TPCbins
(
 const char *nameData = "AnalysisResults.root",
 const char *listName = "RsnOut_Tof20sigma",
 const char *outName  = "proj/kstar",
 Char_t *icut = "77",
 Bool_t isPP= 0,
 Bool_t isMC = 0,
 Bool_t doMixLS = 0
)
{
  // initial setup
  gStyle->SetOptStat("1111");
  gStyle->SetTextFont(42);

   // load utility macro
  gROOT->LoadMacro(Form("%s/projectorInvMass_Centrality_Pt.C+g",macroDir.Data()));
  gSystem->MakeDirectory("proj");
  
  // output lists - define once and then clear them
  TList lUnlikePM;
  TList lUnlikeMP;
  // lUnlikePM.SetOwner(kTRUE);
  // lUnlikeMP.SetOwner(kTRUE);
  lUnlikePM.SetName(Form("%s_UnlikePM", (!isMC ? "Data" : "MC")));
  lUnlikeMP.SetName(Form("%s_UnlikeMP", (!isMC ? "Data" : "MC")));
  
  TList lLikePP;
  TList lLikeMM;
  // lLikePP.SetOwner(kTRUE);
  // lLikeMM.SetOwner(kTRUE);
  lLikePP.SetName(Form("%s_LikePP", (!isMC ? "Data" : "MC")));
  lLikeMM.SetName(Form("%s_LikeMM", (!isMC ? "Data" : "MC")));
  
  TList lMixingMP;
  TList lMixingPM;
  // lMixingMP.SetOwner(kTRUE);
  // lMixingPM.SetOwner(kTRUE);
  lMixingPM.SetName(Form("%s_MixingPM", (!isMC ? "Data" : "MC")));
  lMixingMP.SetName(Form("%s_MixingMP", (!isMC ? "Data" : "MC")));
  
  TList lMixingPP;
  TList lMixingMM;
  // lMixingPP.SetOwner(kTRUE);
  // lMixingMM.SetOwner(kTRUE);
  lMixingPP.SetName(Form("%s_MixingPP", (!isMC ? "Data" : "MC")));
  lMixingMM.SetName(Form("%s_MixingMM", (!isMC ? "Data" : "MC")));
  
  TList lTrues;
  lTrues.SetName("MC_Trues");
  //lTrues.SetOwner(kTRUE);
  // open input file
  TFile *fileData;// = 0x0;//, *fileMC = 0x0;
  TList *listData;// = 0x0;//, *listMC = 0x0;
  fileData = TFile::Open(nameData);
  if (fileData && fileData->IsOpen()) listData = (TList*)fileData->Get(listName);
  
  // names of computed histograms (one per settings)
  THnSparse* hInput[kNhistosData]= {0,0,0,0,0,0,0,0};
  if (listData) {  
    //Official train jobs PbPb
    hInput[ 0] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_UnlikePM", (isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
    hInput[ 1] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_UnlikeMP",(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
    hInput[ 2] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_MixingPM",(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
    hInput[ 3] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_MixingMP",(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
    hInput[ 4] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_LikePP",(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
    hInput[ 5] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_LikeMM",(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
     
    // rename
     hInput[ 0]->SetName("Data_UnlikePM");
     hInput[ 1]->SetName("Data_UnlikeMP");
     hInput[ 2]->SetName("Data_MixingPM");
     hInput[ 3]->SetName("Data_MixingMP");
     hInput[ 4]->SetName("Data_LikePP");
     hInput[ 5]->SetName("Data_LikeMM");
 
     if (doMixLS){
       hInput[ 6] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_MixingPP",(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
       hInput[ 7] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_MixingMM",(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
       hInput[ 6]->SetName("Data_MixingPP");
       hInput[ 7]->SetName("Data_MixingMM");
     }
  }
   
   // define binning in pT    
   Double_t cent[]={ 0.0, 20.0, 40.0, 60.0, 80.0};   
   Double_t pt[] = { 0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0};

   Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;   
   Int_t   ncent  = sizeof(cent) / sizeof(cent[0]) - 1;
 
   //create axis to reproduce the binning
   TAxis *ptbins = new TAxis(npt, pt);
   TAxis *centbins = new TAxis(ncent, cent);
   
   // split histograms
   projectorInvMass_Centrality_Pt projector;
   
   //output file with all histos
   TList * out = new TList();
   out->SetOwner(kTRUE);
     
   // loop on inputs
   for (Int_t i = 0; i < kNhistosData; i++) {
     if (!hInput[i]) continue;
     projector.SetPrefix(hInput[i]->GetName());
     projector.MultiProjPtCent(npt, pt, ncent, cent, hInput[i], out);
   }

   //save separate files for each centrality bin
   for (Int_t icentbin=0; icentbin<ncent;icentbin++){
     lUnlikePM.Clear();
     lUnlikeMP.Clear();
     lMixingPM.Clear();
     lMixingMP.Clear();
     lLikePP.Clear();
     lLikeMM.Clear();
     FillListForCentrality(out,&lUnlikePM,icentbin);
     FillListForCentrality(out,&lUnlikeMP,icentbin);
     FillListForCentrality(out,&lMixingPM,icentbin);
     FillListForCentrality(out,&lMixingMP,icentbin);
     FillListForCentrality(out,&lLikePP,icentbin);
     FillListForCentrality(out,&lLikeMM,icentbin);
     
     if (doMixLS){
       lMixingPP.Clear();
       lMixingMM.Clear();
       FillListForCentrality(out,&lMixingPP,icentbin);
       FillListForCentrality(out,&lMixingMM,icentbin);
     }
     //output file
     TFile *fout = TFile::Open(Form("%s_centBin%02d.root",outName,icentbin), "RECREATE");
     // save into a file
     fout->cd();   
     ptbins->Write("ptbins");
     centbins->Write("centbins");
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
     for (Int_t ih = 0; ih<lUnlikeMP.GetEntries();ih++){
       ((TH1F*)lUnlikeMP.At(ih))->SetLineColor(kYellow+1);
       ((TH1F*)lUnlikeMP.At(ih))->SetMarkerColor(kYellow+1);
       ((TH1F*)lUnlikeMP.At(ih))->SetMarkerStyle(21);
       ((TH1F*)lUnlikeMP.At(ih))->SetMarkerSize(0.8); 
       ((TH1F*)lUnlikeMP.At(ih))->SetLineWidth(2);
       ((TH1F*)lUnlikeMP.At(ih))->SetFillColor(kYellow+1);
       ((TH1F*)lUnlikeMP.At(ih))->SetFillStyle(3002);
       ((TH1F*)lUnlikeMP.At(ih))->SetDrawOption("b");
     }
     for (Int_t ih = 0; ih<lMixingPM.GetEntries();ih++){  
       ((TH1F*)lMixingPM.At(ih))->SetLineColor(kMagenta-2);
       ((TH1F*)lMixingPM.At(ih))->SetMarkerColor(kMagenta-2);
       ((TH1F*)lMixingPM.At(ih))->SetMarkerStyle(1);
       ((TH1F*)lMixingPM.At(ih))->SetLineWidth(2);
     }
     for (Int_t ih = 0; ih<lMixingMP.GetEntries();ih++){
       ((TH1F*)lMixingMP.At(ih))->SetLineColor(kBlue-2);
       ((TH1F*)lMixingMP.At(ih))->SetMarkerColor(kBlue-2);
       ((TH1F*)lMixingMP.At(ih))->SetMarkerStyle(1);
       ((TH1F*)lMixingMP.At(ih))->SetLineWidth(2);
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
     if (doMixLS){
       for (Int_t ih = 0; ih<lMixingPP.GetEntries();ih++){
	 ((TH1F*)lMixingPP.At(ih))->SetLineColor(kOrange+5);
	 ((TH1F*)lMixingPP.At(ih))->SetMarkerColor(kOrange+5);
	 ((TH1F*)lMixingPP.At(ih))->SetMarkerStyle(34);
	 ((TH1F*)lMixingPP.At(ih))->SetMarkerSize(0.8);
	 ((TH1F*)lMixingPP.At(ih))->SetLineWidth(2);
       }
       for (Int_t ih = 0; ih<lMixingMM.GetEntries();ih++){
	 ((TH1F*)lMixingMM.At(ih))->SetLineColor(kGreen+1);
	 ((TH1F*)lMixingMM.At(ih))->SetMarkerColor(kGreen+1);
	 ((TH1F*)lMixingMM.At(ih))->SetMarkerStyle(34);
	 ((TH1F*)lMixingMM.At(ih))->SetMarkerSize(0.8);
	 ((TH1F*)lMixingMM.At(ih))->SetLineWidth(2);
       }   
     }
     lUnlikePM.Write(hInput[0]->GetName(), TObject::kSingleKey);
     lUnlikeMP.Write(hInput[1]->GetName(), TObject::kSingleKey);
     lMixingPM.Write(hInput[2]->GetName(), TObject::kSingleKey);
     lMixingMP.Write(hInput[3]->GetName(), TObject::kSingleKey);
     lLikePP.Write(hInput[4]->GetName(), TObject::kSingleKey);
     lLikeMM.Write(hInput[5]->GetName(), TObject::kSingleKey);
     if (doMixLS){
       lMixingPP.Write(hInput[6]->GetName(), TObject::kSingleKey);
       lMixingMM.Write(hInput[7]->GetName(), TObject::kSingleKey);
     }
     fout->Close();
   }
   // write the full list on dummy file
   TFile *dummy = new TFile(Form("proj/proj_%s",nameData),"RECREATE");
   dummy->cd();
   ptbins->Write("ptbins");
   centbins->Write("centbins");
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
  
  if (listName.Contains("UnlikePM")) type = "UnlikePM";
  if (listName.Contains("UnlikeMP")) type = "UnlikeMP";
  if (listName.Contains("LikeMM")) type = "LikeMM";
  if (listName.Contains("LikePP")) type = "LikePP";
  if (listName.Contains("MixingPM")) type = "MixingPM";
  if (listName.Contains("MixingMP")) type = "MixingMP";
  if (listName.Contains("MixingPP")) type = "MixingPP";
  if (listName.Contains("MixingMM")) type = "MixingMM";
  
  Int_t nhisto = listIn->GetEntries();
  for (Int_t ih=0;ih<nhisto;ih++){
    TH1F * dummy = listIn->At(ih);
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

//----------------------------------------------------
Int_t definePtBins(Int_t npt = 0, Double_t *pt){

  //notice that pt must be an array of npt+1 elements
  // define several choices of binning in pT   
  Float_t pt0_default[2] = {0.0, 10.0};
  Float_t pt1[2] = {1.0, 10.0};
  Float_t pt2[3] = {0.0, 1.0, 10.0};
  Float_t pt14[15] = {0.0, 0.5, 1.00, 1.50, 2.00, 2.50, 3.00, 3.5, 4.00, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };
  Float_t pt18[19] = {0.0, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.5, 4.00, 4.5, 5.0, 6.0, 7.0, 10.00 };
  
  switch (npt) {
  case 1 :
    Printf("===== TR.MOMENTUM BINNING: \n  1.0,10.0");
    for (Int_t j=0;j<npt+1;j++){
      pt[j]=pt1[j];
    }
    break;
  case 2 :
    Printf("===== TR.MOMENTUM BINNING: \n  0.0, 1.0, 10.0");
    for (Int_t j=0;j<npt+1;j++){
      pt[j]=pt2[j];
    }
    break;
  case 14 :
    Printf("===== TR.MOMENTUM BINNING: \n  0.0, 0.5, 1.00, 1.50, 2.00, 2.50, 3.00, 3.5, 4.00, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 ");
    for (Int_t j=0;j<npt+1;j++){
      pt[j]=pt14[j];
    }
    break;
  case 18:
    Printf("===== TR.MOMENTUM BINNING: \n  0.0, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.5, 4.00, 4.5, 5.0, 6.0, 7.0, 10.00");
    for (Int_t j=0;j<npt+1;j++){
      pt[j]=pt18[j];
    }
    break;
  default:
    Printf("===== TR.MOMENTUM BINNING: \n  0.0,10.0");
    for (Int_t j=-2;j<npt;j++){
      pt[j+2]=pt0_default[j+2];
    }
    break;
  }
  return npt;
}

//----------------------------------------------------
Int_t defineCentBins(Int_t nc = 0,Double_t *cent){

  //define binning in centrality
  Float_t cent0_default[2]={ 0.0, 90.0};
  Float_t cent1[2]={ 0.0, 80.0};
  Float_t cent2[3]={ 0.0, 80.0, 100.0};
  Float_t cent5[6]={ 0.0, 20.0, 40.0, 60.0, 80.0, 90.0};   
  Float_t cent9[10]={ 0.0, 10.0 , 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0};   
  Float_t cent10[11]={ 0.0, 5.0, 10.0 , 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0};   
  
  switch (nc) {
     case 1 :
      Printf("===== CENTRALITY BINNING: \n  0.0, 80.0");
      for (Int_t j=0;j<nc+1;j++){
	cent[j]=cent1[j];
      }
      break;
    case 2 :
      Printf("===== CENTRALITY BINNING: \n  0.0, 80.0, 100.0");
      for (Int_t j=0;j<nc+1;j++){
	cent[j]=cent2[j];
      }
      break;
    case 5 :
      Printf("===== CENTRALITY BINNING: \n   0.0, 20.0, 40.0, 60.0, 80.0, 90.0");
      for (Int_t j=0;j<nc+1;j++){
	cent[j]=cent5[j];
      }
      break;
  case 9 :
      Printf("===== CENTRALITY BINNING: \n   0.0, 10.0 , 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0");
      for (Int_t j=0;j<nc+1;j++){
	cent[j]=cent9[j];
      }
      break;
    case 10 :
      Printf("===== CENTRALITY BINNING: \n   0.0, 5.0, 10.0 , 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0");
      for (Int_t j=0;j<nc+1;j++){
	cent[j]=cent10[j];
      }
      break;
  default:
    Printf("===== CENTRALITY BINNING: \n  0.0, 80.0");
    for (Int_t j=0;j<2;j++){
      cent[j]=cent0_default[j];
    }
    break;
  }
  return nc;
}
 
