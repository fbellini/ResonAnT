//
// Splits a TH2F using the myHistSplig_sparse_IM_PT_CENT class
//
Bool_t kComputeEfficiency=kFALSE;
const Int_t kNhistosData=8;
TString macroDir ="/Users/bellini/alice/macro/kstar";
Color_t color[6] = {kRed+3, kMagenta+3, kOrange, kGreen+1, kBlue+3, kCyan+3};
void projectKStar
(
 const char *nameData = "train89.root",
 const char *listName = "RsnOut_tpc3s_tof3sveto",
 const char *outName  = "tpc3s_tof3sveto",
 Char_t *icut = "2424",
 Bool_t saveProj = kFALSE,
 Bool_t isTPC=0,
 Bool_t isPP= 0,
 Bool_t isMC = 0,
 Bool_t doMixLS = 0)
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
    hInput[ 0] = (THnSparse*)listData->FindObject(Form("%sKStar%s%s_%s_kstar_UnlikePM", (isTPC?"TPC":"TOF"),(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
    hInput[ 1] = (THnSparse*)listData->FindObject(Form("%sKStar%s%s_%s_kstar_UnlikeMP", (isTPC?"TPC":"TOF"),(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
    hInput[ 2] = (THnSparse*)listData->FindObject(Form("%sKStar%s%s_%s_kstar_MixingPM", (isTPC?"TPC":"TOF"),(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
    hInput[ 3] = (THnSparse*)listData->FindObject(Form("%sKStar%s%s_%s_kstar_MixingMP", (isTPC?"TPC":"TOF"),(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
    hInput[ 4] = (THnSparse*)listData->FindObject(Form("%sKStar%s%s_%s_kstar_LikePP", (isTPC?"TPC":"TOF"),(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
    hInput[ 5] = (THnSparse*)listData->FindObject(Form("%sKStar%s%s_%s_kstar_LikeMM", (isTPC?"TPC":"TOF"),(isPP? "pp" : "PbPb"), (isMC? "MC" : "Data"), icut));
     
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
  //old ana  Double_t pt[] = { 0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.00 };
  //new ana signed Double_t pt[] = { 0.0, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.00 };
  //integrated Double_t pt[] = { 1.0, 10.0}; 
  //old tpc ana Double_t pt[] = { 0.0, 0.5, 1.0, 1.50, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0,10.0};

  //trd ana signed 
  //Double_t pt[] = {1.0, 2.0, 3.0, 4.0, 5.0, 7.0};
  
  /****************************/
  //pA analysis
  /****************************/
  Double_t cent[]={ 0.0, 100};
  //Double_t cent[]={0.0, 20.0, 40.0, 60.0, 80.0, 100.0};   
  //TOF and TPC standalone analyses - binning 0
  Double_t pt[] = { 0.0, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };
  //TOF and TPC standalone analyses+high pT
  //Double_t pt[] = { 0.0, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00, 12.0, 15., 18.0, 20.0};
  //define binning in pT  - 300MeV bins - binning A 
  //Double_t pt[] = {0.0, 0.15, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };
  //define binning in pT  - 200MeV bins - binning B
  // Double_t pt[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00, 12.0, 15.0};
  //define binning in pT  - 200MeV bins - binning C
  //Double_t pt[] = {0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 12.0, 14.0, 16.0};
  //define binning in pT  - 200MeV bins - binning D
  // Double_t pt[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00, 12.0, 15.0 };

  //*/PILE UP REJECTION CHECK
  //Double_t cent[]={ 0.0, 20.0, 40.0, 60.0, 80.0, 100.0};   
  //Double_t pt[] = { 0.0, 1.0, 5.0};
  //*/

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
     //make-up
     for (Int_t ih = 0; ih<lUnlikePM.GetEntries();ih++){
       ((TH1F*)lUnlikePM.At(ih))->SetLineColor(color[0]-ih);
       ((TH1F*)lUnlikePM.At(ih))->SetMarkerColor(color[0]-ih);
       ((TH1F*)lUnlikePM.At(ih))->SetMarkerStyle(20);
       ((TH1F*)lUnlikePM.At(ih))->SetMarkerSize(0.7); 
       ((TH1F*)lUnlikePM.At(ih))->SetLineWidth(1);
       ((TH1F*)lUnlikePM.At(ih))->SetFillColor(kGray);
       ((TH1F*)lUnlikePM.At(ih))->SetFillStyle(0);// or 3002
       ((TH1F*)lUnlikePM.At(ih))->SetDrawOption("b");
     }
     for (Int_t ih = 0; ih<lUnlikeMP.GetEntries();ih++){
       ((TH1F*)lUnlikeMP.At(ih))->SetLineColor(color[1]-ih);
       ((TH1F*)lUnlikeMP.At(ih))->SetMarkerColor(color[1]-ih);
       ((TH1F*)lUnlikeMP.At(ih))->SetMarkerStyle(24);
       ((TH1F*)lUnlikeMP.At(ih))->SetMarkerSize(0.7); 
       ((TH1F*)lUnlikeMP.At(ih))->SetLineWidth(1);
       ((TH1F*)lUnlikeMP.At(ih))->SetFillColor(color[1]-ih);
       ((TH1F*)lUnlikeMP.At(ih))->SetFillStyle(9);//or 3002
       ((TH1F*)lUnlikeMP.At(ih))->SetDrawOption("b");
     }
     for (Int_t ih = 0; ih<lMixingPM.GetEntries();ih++){  
       ((TH1F*)lMixingPM.At(ih))->SetLineColor(color[4]-ih);
       ((TH1F*)lMixingPM.At(ih))->SetMarkerColor(color[4]-ih);
       ((TH1F*)lMixingPM.At(ih))->SetMarkerStyle(1);
       ((TH1F*)lMixingPM.At(ih))->SetLineWidth(2);
     }
     for (Int_t ih = 0; ih<lMixingMP.GetEntries();ih++){
       ((TH1F*)lMixingMP.At(ih))->SetLineColor(color[5]-ih);
       ((TH1F*)lMixingMP.At(ih))->SetMarkerColor(color[5]-ih);
       ((TH1F*)lMixingMP.At(ih))->SetMarkerStyle(1);
       ((TH1F*)lMixingMP.At(ih))->SetLineWidth(2);
     }
     for (Int_t ih = 0; ih<lLikePP.GetEntries();ih++){
       ((TH1F*)lLikePP.At(ih))->SetLineColor(color[2]);
       ((TH1F*)lLikePP.At(ih))->SetMarkerColor(color[2]);
       ((TH1F*)lLikePP.At(ih))->SetMarkerStyle(1);
       ((TH1F*)lLikePP.At(ih))->SetLineWidth(2);
     }
     for (Int_t ih = 0; ih<lLikeMM.GetEntries();ih++){
       ((TH1F*)lLikeMM.At(ih))->SetLineColor(color[3]);
       ((TH1F*)lLikeMM.At(ih))->SetMarkerColor(color[3]);
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
	 ((TH1F*)lMixingMM.At(ih))->SetLineColor(kPink+1);
	 ((TH1F*)lMixingMM.At(ih))->SetMarkerColor(kPink+1);
	 ((TH1F*)lMixingMM.At(ih))->SetMarkerStyle(34);
	 ((TH1F*)lMixingMM.At(ih))->SetMarkerSize(0.8);
	 ((TH1F*)lMixingMM.At(ih))->SetLineWidth(2);
       }   
     }
     if (saveProj){
       //output file
       TFile *fout = TFile::Open(Form("%s_centBin%02d.root",outName,icentbin), "RECREATE");
       // save into a file
       fout->cd();   
       ptbins->Write("ptbins");
       centbins->Write("centbins");
       
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
   }
   // write the full list on dummy file
   TFile *dummy = new TFile(Form("proj_%s_%s.root",icut, outName, nameData),"RECREATE");
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

