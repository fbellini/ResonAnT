//
// Splits a TH2F using the myHistSplig_sparse_IM_PT_CENT class
//
Bool_t kComputeEfficiency=1;
const Int_t kNhistosData=4;
TString macroDir ="/Users/bellini/alice/macro/kstar";

void projectKStarMCtrue_old
(
   const char *nameData = "analysisAOD.root",
   TString listName = "RsnOut_TofPid",
   const char *outName  = "mcproj/eff",
   char * cutID = "77",
   Bool_t isPP=0,
   Bool_t isMC=1
)
{
   // initial setup
  gStyle->SetOptStat("1111");
  gStyle->SetTextFont(42);

   // load utility macro
  gROOT->LoadMacro(Form("%s/projectorInvMass_Centrality_Pt.C+g",macroDir.Data()));
  
  TFile *dummy = new TFile(Form("proj/projMC_%s",nameData),"RECREATE");
  TList * out = new TList();
  
  // output lists - define once and then clear them
  TList lTruesPM;
  TList lTruesMP;
  lTruesPM.SetName(Form("%s_TruesPM", (!isMC ? "Data" : "MC")));
  lTruesMP.SetName(Form("%s_TruesMP", (!isMC ? "Data" : "MC")));
  
  TList lMother;//,lAntiMother;
  lMother.SetName(Form("%s_Mother", (!isMC ? "Data" : "MC")));
  // lAntiMother.SetName(Form("%s_AntiMother", (!isMC ? "Data" : "MC")));

   // open input file
  TFile *fileData;// = 0x0;//, *fileMC = 0x0;
  TList *listData;
  TList *listDataAnti;// = 0x0;//, *listMC = 0x0;
  // TString antiListName=listName.Data();
  // antiListName.ReplaceAll("RsnOut_","RsnOut_anti_");
  fileData = TFile::Open(nameData);
  if (fileData && fileData->IsOpen()) listData = (TList*)fileData->Get(listName.Data());
  // if (fileData && fileData->IsOpen()) listDataAnti = (TList*)fileData->Get(antiListName.Data());
  
   // names of computed histograms (one per settings)
  THnSparse* hInput[kNhistosData]= {0,0,0,0};
   if (listData) {
    
     hInput[ 0] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_TruesPM", (isPP ? "pp" : "PbPb"), (isMC? "MC" : "Data"), cutID));
     hInput[ 1] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_TruesMP", (isPP ? "pp" : "PbPb"), (isMC? "MC" : "Data"), cutID));
     hInput[ 2] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_Mother", (isPP ? "pp" : "PbPb"), (isMC? "MC" : "Data"), cutID));
     // hInput[ 3] = (THnSparse*)listDataAnti->FindObject(Form("TOFKStar%s%s_%s_kstar_Mother", (isPP ? "pp" : "PbPb"), (isMC? "MC" : "Data"), cutID));

     // rename
     hInput[ 0]->SetName("TruesPM");
     hInput[ 1]->SetName("TruesMP");
     hInput[ 2]->SetName("Mother");
     // hInput[ 3]->SetName("AntiMother");

   }
   
   // define binning in pT    
   Double_t cent[]={ 0.0, 20.0, 40.0, 60.0, 80.0};   
   Double_t pt[] = {0.0, 0.5, 1.00, 1.50, 2.00, 2.50, 3.00, 3.5, 4.00, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };
    Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;
   
   //define binning in centrality
   //Double_t cent[]={ 0.0, 100.0};
   Int_t   ncent  = sizeof(cent) / sizeof(cent[0]) - 1;
 
   //create axis to reproduce the binning
   TAxis *ptbins = new TAxis(npt, pt);
   TAxis *centbins = new TAxis(ncent, cent);
   
     // split histograms
   projectorInvMass_Centrality_Pt projector;
   
   // loop on inputs
   for (Int_t i = 0; i < kNhistosData; i++) {
     if (!hInput[i]) continue;
     // projector.SetPrefix(hInput[i]->GetName());
     // projector.MultiProj(npt, pt, hInput[i], out);
     projector.SetPrefix(hInput[i]->GetName());
     projector.MultiProjPtCent(npt, pt, ncent, cent, hInput[i], out);
   }

   TCanvas *c[5];
   c[0]= new TCanvas(Form("c_0"),Form("cent_0"),1200,600);
   c[1]= new TCanvas(Form("c_1"),Form("cent_1"),1200,600);
   c[2]= new TCanvas(Form("c_2"),Form("cent_2"),1200,600);
   c[3]= new TCanvas(Form("c_3"),Form("cent_3"),1200,600);
   c[4]= new TCanvas(Form("c_4"),Form("cent_4"),1200,600);
  
   //save separate files for each centrality bin
   for (Int_t icentbin=0; icentbin<ncent;icentbin++){
     c[icentbin]->Divide(5,3);
     lTruesPM.Clear();
     lTruesMP.Clear();
     lMother.Clear();
     //     lAntiMother.Clear();
     FillListForCentrality(out,&lTruesPM,icentbin);
     FillListForCentrality(out,&lTruesMP,icentbin);
     FillListForCentrality(out,&lMother,icentbin);
     //  FillListForCentrality(out,&lAntiMother,icentbin);
     //output file
     TFile *efffile = TFile::Open(Form("efficiency_%s_centBin%02d.root",listName.Data(),icentbin), "RECREATE");
      //Create histos for efficiency
     TH1F * hTrueCountsPM = new TH1F("hTrueCountsPM","True Counts K+#pi-", npt, pt);
     TH1F * hTrueCountsMP = new TH1F("hTrueCountsMP","True Counts K-#pi+", npt, pt);
     TH1F * hMotherCounts = new TH1F("hMotherCounts","hMotherCounts", npt, pt);
     //   TH1F * hAntiMotherCounts = new TH1F("hAntiMotherCounts","hAntiMotherCounts", npt, pt);
     TH1F * hEffVsPt = new TH1F("hEffVsPt", "Efficiency for K*+#bar{K*} in PbPb@2.76 TeV; p_{t} (GeV/c); #epsilon = reco true / generated", npt, pt);
     hTrueCountsPM->Sumw2();
     hTrueCountsMP->Sumw2();
     hMotherCounts->Sumw2();
     //    hAntiMotherCounts->Sumw2();
     hEffVsPt->Sumw2();
     //make-up
     for (Int_t ih = 0; ih<lTruesPM.GetEntries();ih++){
       ((TH1F*)lTruesPM.At(ih))->SetLineColor(kMagenta+2);
       ((TH1F*)lTruesPM.At(ih))->SetMarkerColor(kMagenta+2);
       ((TH1F*)lTruesPM.At(ih))->SetMarkerStyle(20);
       ((TH1F*)lTruesPM.At(ih))->SetMarkerSize(0.5); 
       ((TH1F*)lTruesPM.At(ih))->SetLineWidth(1);
       ((TH1F*)lTruesPM.At(ih))->SetFillColor(kMagenta-8);
       ((TH1F*)lTruesPM.At(ih))->SetFillStyle(3002);
       ((TH1F*)lTruesPM.At(ih))->SetDrawOption("H");
       Int_t counts = ((TH1F*)lTruesPM.At(ih))->Integral();
       hTrueCountsPM->SetBinContent(ih+1, counts);
       hTrueCountsPM->SetBinError(ih+1, TMath::Sqrt(counts));

       ((TH1F*)lTruesMP.At(ih))->SetLineColor(kOrange-2);
       ((TH1F*)lTruesMP.At(ih))->SetMarkerColor(kOrange-2);
       ((TH1F*)lTruesMP.At(ih))->SetMarkerStyle(21);
       ((TH1F*)lTruesMP.At(ih))->SetMarkerSize(0.5); 
       ((TH1F*)lTruesMP.At(ih))->SetLineWidth(1);
       ((TH1F*)lTruesMP.At(ih))->SetFillColor(kYellow+10);
       ((TH1F*)lTruesMP.At(ih))->SetFillStyle(3002);
       //((TH1F*)lTruesMP.At(ih))->SetDrawOption("bar");
       Int_t counts = ((TH1F*)lTruesMP.At(ih))->Integral();
       hTrueCountsMP->SetBinContent(ih+1, counts);
       hTrueCountsMP->SetBinError(ih+1, TMath::Sqrt(counts));

       ((TH1F*)lMother.At(ih))->SetLineColor(kCyan-2);
       ((TH1F*)lMother.At(ih))->SetMarkerColor(kCyan-2);
       ((TH1F*)lMother.At(ih))->SetMarkerStyle(1);
       ((TH1F*)lMother.At(ih))->SetLineWidth(2);
       Int_t counts = ((TH1F*)lMother.At(ih))->Integral();
       hMotherCounts->SetBinContent(ih+1, counts);
       hMotherCounts->SetBinError(ih+1, TMath::Sqrt(counts));

       // ((TH1F*)lAntiMother.At(ih))->SetLineColor(kBlue+2);
       // ((TH1F*)lAntiMother.At(ih))->SetMarkerColor(kBlue+2);
       // ((TH1F*)lAntiMother.At(ih))->SetMarkerStyle(1);
       // ((TH1F*)lAntiMother.At(ih))->SetLineWidth(2);
       // Int_t counts = ((TH1F*)lAntiMother.At(ih))->Integral();
       // hAntiMotherCounts->SetBinContent(ih+1, counts);
       // hAntiMotherCounts->SetBinError(ih+1, TMath::Sqrt(counts));
       
       c[icentbin]->cd(ih+1);
       gPad->SetLogy();
       ((TH1F*)lMother.At(ih))->GetYaxis()->SetRangeUser(0.1,5.e5);
       ((TH1F*)lMother.At(ih))->Draw();
       //  ((TH1F*)lAntiMother.At(ih))->Draw("same");
       ((TH1F*)lTruesMP.At(ih))->Draw("same");
       ((TH1F*)lTruesPM.At(ih))->Draw("same");
     }

     //efficiency estimation
     TH1D* halltrue = (TH1D*)hTrueCountsPM->Clone();  
     halltrue->Add(hTrueCountsMP,1);
     
     TH1D* hallmother = (TH1D*)hMotherCounts->Clone();  
     // hallmother->Add(hAntiMotherCounts,1);
     
     // hEffVsPt = (TH1F*) halltrue->Clone("efficiency");
     // hEffVsPt->Divide(hallmother);
     for (Int_t ii = 1; ii<npt+1; ii++){
       Double_t alltrue = halltrue->GetBinContent(ii);
       Double_t errtrue = halltrue->GetBinError(ii);
       Double_t mothers = hMotherCounts->GetBinContent(ii);//+hAntiMotherCounts->GetBinContent(ii);
       Double_t errMothers = hMotherCounts->GetBinError(ii);
       //Double_t errAntiMothers = hAntiMotherCounts->GetBinError(ii);
       Double_t errmum = errMothers;//TMath::Sqrt(errAntiMothers*errAntiMothers+errMothers*errMothers);
       Double_t ratio;
       if (mothers>0) 
	 ratio = alltrue/mothers;
       else ratio = 0.0;
       Double_t err = ratio*TMath::Sqrt((errtrue/alltrue)*(errtrue/alltrue)+(errmum/mothers)*(errmum/mothers));
       hEffVsPt->SetBinContent(ii, ratio);
       hEffVsPt->SetBinError(ii, err);
       Printf("Efficiency  %3.1f < pt < %3.1f = %6.5f +/- %6.5f", pt[ii-1], pt[ii], ratio, err);
       }
     
     hEffVsPt->GetYaxis()->SetRangeUser(0.0,1.0);
     // save into a file
     efffile->cd();
     ptbins->Write("ptbins");
     centbins->Write("centbins");
     hMotherCounts->Write();
     //hAntiMotherCounts->Write();
     hTrueCountsMP->Write();
     hTrueCountsPM->Write();
     halltrue->SetName("AllTrue");
     halltrue->Write();
     hallmother->SetName("AllMother");
     hallmother->Write();
     hEffVsPt->Write();
     efffile->Close();


     //output file
     TFile *fout = TFile::Open(Form("%s_centBin%02d.root",outName,icentbin), "RECREATE");
     fout->cd();   
     lTruesPM.Write(hInput[0]->GetName(), TObject::kSingleKey);
     lTruesMP.Write(hInput[1]->GetName(), TObject::kSingleKey);
     lMother.Write(hInput[2]->GetName(), TObject::kSingleKey);
     // lAntiMother.Write(hInput[3]->GetName(), TObject::kSingleKey);
     fout->Close();
   }
   // write the full list on dummy file
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
  if (listName.Contains("TruesPM")) type = "TruesPM";
  if (listName.Contains("TruesMP")) type = "TruesMP";
  if (listName.Contains("Mother")) type = "Mother";
  
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
