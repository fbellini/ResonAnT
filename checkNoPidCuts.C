TString infileName = "analysisAOD.root";
TString listPrefix = "RsnOut";
TString tofListSuffix = "_noPid";
Bool_t isQsameFile = 1;
enum specie{ kAnyQ, kPi, kK, /* kPro, kAnyMatch,*/ kNtot};
TString cutNames[species::kNtot]={ "cutQuality", "cutPionTOFPbPb2010_10.0sigma", "cutKaonTOFPbPb2010_10.02sigma"/*, "cutProtonTof","cutMatchTof"*/}; //_matchNoPid
Color_t color[3] = {kRed, kBlue, kGreen+1};
Color_t fillcolor[3] = {kOrange, kAzure+6, kGreen-9};
Int_t marker[3] = {20, 21, 22};
Int_t ptRebinFactor = 20;
//-------------------------------------------
void checkNoPidCuts(TString finName, Bool_t display = kTRUE, Bool_t display=1){

  //input
  if (finName) infileName = finName.Data();
  TFile * fin = TFile::Open(infileName.Data());
  TList * ltof     = (TList*) fin->Get(Form("%s%s", listPrefix.Data(),tofListSuffix.Data()));
  TList * lquality;
  if (!isQsameFile) lquality = (TList*) fin->Get(listPrefix.Data());
  else  lquality = (TList*) fin->Get(Form("%s%s", listPrefix.Data(),tofListSuffix.Data()));

  TH1D * hPt[species::kNtot]; 
  TH1D * hP[species::kNtot]; 
  TH1D * hEta[species::kNtot]; 

  //output
  TFile * fOut = new TFile("checkNoPidCuts.root","recreate");

  /* CHECK CUTS FOR PID */
  for (Int_t icut = 0; icut<species::kNtot; icut++){

    TList * pidPerformance = new TList();
    //quality
    if (icut==species::kAnyQ){
      hP[icut]=(TH1D*)(lquality->FindObject(Form("%s.P_p",cutNames[icut].Data())) )->Clone();
      hPt[icut]=(TH1D*)(lquality->FindObject(Form("%s.Pt_pt",cutNames[icut].Data())) )->Clone();
      hEta[icut]=(TH1D*)(lquality->FindObject(Form("%s.Eta_eta",cutNames[icut].Data())) )->Clone();
    } else { 
      //tof pid
      hP[icut]=(TH2F*)(ltof->FindObject(Form("%s.P_p",cutNames[icut].Data())) )->Clone();
      hPt[icut]=(TH2F*)(ltof->FindObject(Form("%s.Pt_pt",cutNames[icut].Data())) )->Clone();
      hEta[icut]=(TH2F*)(ltof->FindObject(Form("%s.Eta_eta",cutNames[icut].Data())) )->Clone();
    }
    
    hP[icut]->GetXaxis()->SetTitle("p (GeV/c)");
    hPt[icut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hEta[icut]->GetXaxis()->SetTitle("#eta");
    
  }//END LOOP ON CUTS
    
  //compute ratios between tracks passing Quality cuts and PID cuts
  TH1D * ratioPionCut_P =  (TH1D*) hP[species::kPi]->Clone();
  TH1D * ratioPionCut_Pt = (TH1D*) hPt[species::kPi]->Clone();
  TH1D * ratioPionCut_Eta = (TH1D*) hEta[species::kPi]->Clone();
  ratioPionCut_P->Rebin(ptRebinFactor/2);
  ratioPionCut_Pt->Rebin(ptRebinFactor/2);
  hP[species::kAnyQ]->Rebin(ptRebinFactor/2);
  hPt[species::kAnyQ]->Rebin(ptRebinFactor/2);
  ratioPionCut_P->Divide( hP[species::kAnyQ] );
  ratioPionCut_Pt->Divide( hPt[species::kAnyQ] );
  ratioPionCut_Eta->Divide( hEta[species::kAnyQ] );
    
  TH1D * ratioKaonCut_P =  (TH1D*) hP[species::kK]->Clone();
  TH1D * ratioKaonCut_Pt = (TH1D*) hPt[species::kK]->Clone();
  TH1D * ratioKaonCut_Eta = (TH1D*) hEta[species::kK]->Clone();
  ratioKaonCut_P->Rebin(ptRebinFactor/2);
  ratioKaonCut_Pt->Rebin(ptRebinFactor/2);
  ratioKaonCut_P->Divide( hP[species::kAnyQ] );
  ratioKaonCut_Pt->Divide( hPt[species::kAnyQ] );
  ratioKaonCut_Eta->Divide( hEta[species::kAnyQ] );
    
  for (Int_t j=0;j<3;j++){
    ratioPionCut_P->GetYaxis()->SetTitle("PID cut /quality cut");
    ratioPionCut_Pt->GetYaxis()->SetTitle("PID cut /quality cut");
    ratioPionCut_Eta->GetYaxis()->SetTitle("PID cut /quality cut");
    ratioKaonCut_P->GetYaxis()->SetTitle("PID cut /quality cut");
    ratioKaonCut_Pt->GetYaxis()->SetTitle("PID cut /quality cut");
    ratioKaonCut_Eta->GetYaxis()->SetTitle("PID cut /quality cut");
  }

  fOut->cd();
  ratioPionCut_P->Write();
  ratioPionCut_Pt->Write();
  ratioPionCut_Eta->Write();
  ratioKaonCut_P->Write();
  ratioKaonCut_Pt->Write();
  ratioKaonCut_Eta->Write();

  fOut->Close();
  return;
  }
