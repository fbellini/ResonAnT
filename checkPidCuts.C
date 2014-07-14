TString infileName = "analysisAOD.root";
TString listPrefix = "RsnOut";
Bool_t isQsameFile = 1;
enum species{ kAnyQ, kPi, kK, /* kPro, kAnyMatch,*/ kNtot};

//pp
//TString cutNames[species::kNtot]={ "cutQuality", "cutPionTOFPbPb2010_2.0sigma", "cutKaonTOFPbPb2010_2.02sigma"/*, "cutProtonTof","cutMatchTof"*/}; //_matchNoPid
//PbPb
//TString cutNames[species::kNtot]={ "cutQuality", "cutPionTOFPbPb2010_2.0sigma", "cutKaonTOFPbPb2010_2.02sigma"/*, "cutProtonTof","cutMatchTof"*/}; //_matchNoPid
TString cutNames[species::kNtot]={"cutQuality", "cutPionTOFPbPb2010_2.0sigma", "cutKaonTOFPbPb2010_2.0sigma"/*, "cutProtonTof","cutMatchTof"*/};
Color_t color[species::kNtot] = {kRed, kBlue, kGreen+1};
Color_t fillcolor[species::kNtot] = {kOrange, kAzure+6, kGreen-9};
Int_t marker[species::kNtot] = {20, 21, 22};
Int_t ptRebinFactor = 20;
//-------------------------------------------
void checkPidCuts(TString finName, TString tofListSuffix = "_Q5Tof20sigma", Bool_t display = kTRUE){

  //input
  if (finName) infileName = finName.Data();
  TFile * fin = TFile::Open(infileName.Data());
  TList * ltof     = (TList*) fin->Get(Form("%s%s", listPrefix.Data(),tofListSuffix.Data()));
  TList * lquality;
  if (!isQsameFile) lquality = (TList*) fin->Get(listPrefix.Data());
  else  lquality = (TList*) ltof; //fin->Get(Form("%s%s", listPrefix.Data(),tofListSuffix.Data()));
  if (!lquality) return;

  //output
  TFile * fOut = new TFile(Form("checkPidCuts%s.root",tofListSuffix.Data()),"recreate");
  TList * pidPerformance = new TList();
  TList * lProjNsigmaTPC = new TList();
  lProjNsigmaTPC->SetOwner(1);
  TList * lProjNsigmaTOF = new TList();
  lProjNsigmaTOF->SetOwner(1);

  TH1D * hPt[species::kNtot]; 
  TH1D * hP[species::kNtot]; 
  TH1D * hEta[species::kNtot]; 
  TH2F * hNsigmaTPCVsPtpc[3];
  TH2F * hNsigmaTOFVsP[3];

  for (Int_t ip = 0; ip<3;ip++){ //loop on pi, k, p 
    hNsigmaTPCVsPtpc[ip] = new TH2F();
    hNsigmaTOFVsP[ip] = new TH2F();
  }	
  
  /* CHECK CUTS FOR PID */
  for (Int_t icut = 0; icut<species::kNtot; icut++){
    for (Int_t j=0;j<3;j++){
      hNsigmaTPCVsPtpc[j]->Reset();
      hNsigmaTOFVsP[j]->Reset();
    }
    //quality
    if (icut==species::kAnyQ){
      hP[icut]=(TH1D*)(lquality->FindObject(Form("%s.P_p",cutNames[icut].Data())) )->Clone();
      hPt[icut]=(TH1D*)(lquality->FindObject(Form("%s.Pt_pt",cutNames[icut].Data())) )->Clone();
      hEta[icut]=(TH1D*)(lquality->FindObject(Form("%s.Eta_eta",cutNames[icut].Data())) )->Clone();
      
      hNsigmaTPCVsPtpc[0]=(TH2F*)(lquality->FindObject(Form("%s.TPC_nsigmaPi_VsPtpc_pTPC_pi",cutNames[icut].Data())) )->Clone();
      hNsigmaTPCVsPtpc[1]=(TH2F*)(lquality->FindObject(Form("%s.TPC_nsigmaK_VsPtpc_pTPC_K",cutNames[icut].Data())) )->Clone();
      hNsigmaTPCVsPtpc[2]=(TH2F*)(lquality->FindObject(Form("%s.TPC_nsigmaPro_VsPtpc_pTPC_p",cutNames[icut].Data())) )->Clone();
      
      hNsigmaTOFVsP[0]=(TH2F*)(lquality->FindObject(Form("%s.TOF_nsigmaPi_vsP_p_pi",cutNames[icut].Data())) )->Clone();
      hNsigmaTOFVsP[1]=(TH2F*)(lquality->FindObject(Form("%s.TOF_nsigmaK_vsP_p_K",cutNames[icut].Data())) )->Clone();
      hNsigmaTOFVsP[2]=(TH2F*)(lquality->FindObject(Form("%s.TOF_nsigmaPro_vsP_p_p",cutNames[icut].Data())) )->Clone();
    } else { 
      //tof pid
      hP[icut]=(TH2F*)(ltof->FindObject(Form("%s.P_p",cutNames[icut].Data())) )->Clone();
      hPt[icut]=(TH2F*)(ltof->FindObject(Form("%s.Pt_pt",cutNames[icut].Data())) )->Clone();
      hEta[icut]=(TH2F*)(ltof->FindObject(Form("%s.Eta_eta",cutNames[icut].Data())) )->Clone();
    
      hNsigmaTPCVsPtpc[0]=(TH2F*)(ltof->FindObject(Form("%s.TPC_nsigmaPi_VsPtpc_pTPC_pi",cutNames[icut].Data())) )->Clone();
      hNsigmaTPCVsPtpc[1]=(TH2F*)(ltof->FindObject(Form("%s.TPC_nsigmaK_VsPtpc_pTPC_K",cutNames[icut].Data())) )->Clone();
      hNsigmaTPCVsPtpc[2]=(TH2F*)(ltof->FindObject(Form("%s.TPC_nsigmaPro_VsPtpc_pTPC_p",cutNames[icut].Data())) )->Clone();

      hNsigmaTOFVsP[0]=(TH2F*)(ltof->FindObject(Form("%s.TOF_nsigmaPi_vsP_p_pi",cutNames[icut].Data())) )->Clone();
      hNsigmaTOFVsP[1]=(TH2F*)(ltof->FindObject(Form("%s.TOF_nsigmaK_vsP_p_K",cutNames[icut].Data())) )->Clone();
      hNsigmaTOFVsP[2]=(TH2F*)(ltof->FindObject(Form("%s.TOF_nsigmaPro_vsP_p_p",cutNames[icut].Data())) )->Clone();
    }

    //make nsigma projections
    //note: TPCnsigma on y axis, TOFnsigma on Xaxis
    for (Int_t ip = 0; ip<3;ip++){ //loop on pi, k, p      
      //make nsigma projections for slices of Dpt = 0.2 GeV/c for nsigma_TPC
      MakeSliceProjections(hNsigmaTPCVsPtpc[ip],lProjNsigmaTPC, 1, ptRebinFactor, -1, ip, "n#sigma_{TPC}",  "p_{TPC} (GeV/c)" ); //TPCnsigma on y axis, slices of pT along x
      //make nsigma projections for slices of Dpt = 0.2 GeV/c for nsigma_TOF
      MakeSliceProjections(hNsigmaTOFVsP[ip],lProjNsigmaTOF, 1, ptRebinFactor, -1, ip, "n#sigma_{TOF}", "p (GeV/c)"); // TOFnsigma on x axis, slices of pT along y
      fOut->cd();
      lProjNsigmaTPC->Write(Form("proj_NsigmaTPC_%s_%i",cutNames[icut].Data(),ip), 1);      
      lProjNsigmaTOF->Write(Form("proj_NsigmaTOF_%s_%i",cutNames[icut].Data(),ip), 1);
    }
    hP[icut]->GetXaxis()->SetTitle("p (GeV/c)");
    hPt[icut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hEta[icut]->GetXaxis()->SetTitle("#eta");
  }//END LOOP ON CUTS
    
  //compute ratios between tracks passing Quality cuts and PID cuts
  // for (Int_t j=0;j<3;j++){
  // }
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


//-------------------------------------------------------------------
Int_t MakeSliceProjections(TH2F * histo, TList * projList, Bool_t sliceXaxis = 0, Int_t rebinFactor=1, Int_t selectedBin = -1, Int_t colorID = 0, TString titleXaxis = "projected var", TString titleYaxis = "p (GeV/c)")
{
  //make nsigma projections for slices of pT
  if (!histo) return 1;
  if (!projList) return 2;
  //define axis
  TAxis * slicedAxis;
  TAxis * projectedAxis;
  if (sliceXaxis) {
    histo->RebinX(rebinFactor);
    slicedAxis = (TAxis*) histo->GetXaxis();
    projectedAxis = (TAxis*) histo->GetYaxis();
  } else {
    histo->RebinY(rebinFactor);
    slicedAxis = (TAxis*) histo->GetYaxis();
    projectedAxis = (TAxis*) histo->GetXaxis();
  }
  TString histoName(histo->GetName());
  Int_t nslices = slicedAxis->GetNbins();
  Int_t nbins = projectedAxis->GetNbins();
  Printf("Projecting %s for %i slices the content of %i bins", histo->GetName(),nslices,nbins);
  
  TH1D * projDummy = new TH1D(); 
  
  for (Int_t islice = 1; islice < nslices; islice++){
    //select only one bin 
    if (selectedBin>0 && islice!=selectedBin) continue;
    
    //check entries per slice, if <20 do not project
    Bool_t isFilled = kFALSE;
    Int_t nentries = 0;  
    for (Int_t ibin=1;ibin<nbins;ibin++){
      if (sliceXaxis) nentries+=histo->GetBinContent(islice,ibin);
      else nentries+=histo->GetBinContent(ibin,islice);     
      if (nentries>20) isFilled = kTRUE;
    }
    //    Printf("Entries %i in slice %i", nentries, islice);
    //project slices
    if (isFilled){
      Float_t sliceLowEdge=slicedAxis->GetBinLowEdge(islice);
      Float_t sliceUpEdge=slicedAxis->GetBinUpEdge(islice);
      if (sliceXaxis) 
	projDummy = (TH1D*) histo->ProjectionY(Form("proj_%i",islice), islice, islice);
      else 
	projDummy = (TH1D*) histo->ProjectionX(Form("proj_%5.2f-%5.2f",sliceLowEdge,sliceUpEdge), islice, islice);
      //     projDummy->SetTitle(Form("proj_%i:%5.2f #leq %s < %5.2f", islice, sliceLowEdge,titleYaxis.Data(),  sliceUpEdge));
      projDummy->SetTitle(Form("%5.2f #leq %s < %5.2f", sliceLowEdge,titleYaxis.Data(),sliceUpEdge));
      projDummy->GetXaxis()->SetTitle(titleXaxis.Data());
      projDummy->SetLineColor(color[colorID]);
      projDummy->SetMarkerColor(color[colorID]);
      if (!histoName.Contains("Quality")){
	projDummy->SetFillColor(fillcolor[colorID]);
	projDummy->SetDrawOption("h");
      }
      projList->AddLast(projDummy);
    }
  }
  projList->AddLast(histo);
  return 0;
}

//-------------------------------------------------------------------
void DisplayOnCanvas(TString finName="checkPidCuts.root")
{
  
  Int_t id4tv[4]={3,6,8,13};

  TFile * fin = TFile::Open(finName.Data());
  TList * lProjNsigmaQPi = (TList*) fin->Get("proj_NsigmaTOF_cutQuality_0");
  TList * lProjNsigmaTOFPi = (TList*) fin->Get(Form("proj_NsigmaTOF_%s_0",cutNames[1].Data()));
  TList * lProjNsigmaQK = (TList*) fin->Get("proj_NsigmaTOF_cutQuality_1");
  TList * lProjNsigmaTOFK = (TList*) fin->Get(Form("proj_NsigmaTOF_%s_1",cutNames[2].Data()));
  TLegend * legPi = new TLegend(0.6,0.7,0.88,0.88);
  legPi->SetFillColor(kWhite);
  legPi->SetBorderSize(0.0);
  TLegend * legK = new TLegend(0.6,0.7,0.88,0.88);
  legK->SetFillColor(kWhite);
  legK->SetBorderSize(0.0);
    //show on canvas pion ID
  TCanvas *tv = new TCanvas("cPi","cPi",1200,500);
  tv->Divide(2,2);
  TH1F* hQPi_05 = (TH1F*) lProjNsigmaQPi->FindObject(Form("proj_%i",id4tv[0]));
  TH1F* hQPi_1 = (TH1F*) lProjNsigmaQPi->FindObject(Form("proj_%i",id4tv[1]));
  TH1F* hQPi_15 = (TH1F*) lProjNsigmaQPi->FindObject(Form("proj_%i",id4tv[2]));
  TH1F* hQPi_25 = (TH1F*) lProjNsigmaQPi->FindObject(Form("proj_%i",id4tv[3]));
  hQPi_05->GetXaxis()->SetRangeUser(-10.,10.);
  hQPi_05->GetXaxis()->SetTitle("n#sigma_{#pi}^{TOF}");
  hQPi_05->GetYaxis()->SetTitle("tracks");
  hQPi_05->GetXaxis()->SetLabelSize(0.05);
  hQPi_05->GetXaxis()->SetTitleSize(0.05);
  hQPi_05->GetYaxis()->SetLabelSize(0.05);
  hQPi_05->GetYaxis()->SetTitleSize(0.05);
  hQPi_1->GetXaxis()->SetRangeUser(-10.,10.);
  hQPi_1->GetXaxis()->SetTitle("n#sigma_{#pi}^{TOF}");
  hQPi_1->GetYaxis()->SetTitle("tracks");
  hQPi_15->GetXaxis()->SetRangeUser(-10.,10.);
  hQPi_1->GetXaxis()->SetLabelSize(0.05);
  hQPi_1->GetXaxis()->SetTitleSize(0.05);
  hQPi_1->GetYaxis()->SetLabelSize(0.05);
  hQPi_1->GetYaxis()->SetTitleSize(0.05);
  hQPi_15->GetXaxis()->SetTitle("n#sigma_{#pi}^{TOF}");
  hQPi_15->GetYaxis()->SetTitle("tracks");
  hQPi_15->GetXaxis()->SetLabelSize(0.05);
  hQPi_15->GetXaxis()->SetTitleSize(0.05);
  hQPi_15->GetYaxis()->SetLabelSize(0.05);
  hQPi_15->GetYaxis()->SetTitleSize(0.05);
  hQPi_25->GetXaxis()->SetRangeUser(-10.,10.);
  hQPi_25->GetXaxis()->SetTitle("n#sigma_{#pi}^{TOF}");
  hQPi_25->GetYaxis()->SetTitle("tracks");
  hQPi_25->GetXaxis()->SetLabelSize(0.05);
  hQPi_25->GetXaxis()->SetTitleSize(0.05);
  hQPi_25->GetYaxis()->SetLabelSize(0.05);
  hQPi_25->GetYaxis()->SetTitleSize(0.05);
  legPi->AddEntry(hQPi_05,"std quality cuts","lpf");
  tv->cd(1);   hQPi_05->Draw();
  tv->cd(2);   hQPi_1->Draw();
  tv->cd(3);   hQPi_15->Draw();
  tv->cd(4);   hQPi_25->Draw();
  TH1F* htofPi_05 = (TH1F*) lProjNsigmaTOFPi->FindObject(Form("proj_%i",id4tv[0]));
  TH1F* htofPi_1 = (TH1F*) lProjNsigmaTOFPi->FindObject(Form("proj_%i",id4tv[1]));
  TH1F* htofPi_15 = (TH1F*) lProjNsigmaTOFPi->FindObject(Form("proj_%i",id4tv[2]));
  TH1F* htofPi_25 = (TH1F*) lProjNsigmaTOFPi->FindObject(Form("proj_%i",id4tv[3]));
  legPi->AddEntry(htofPi_05,"2#sigma_{#pi}^{TOF} && 5#sigma_{#pi}^{TPC} cut","lpf");
  tv->cd(1);   gPad->SetLogy(); htofPi_05->Draw("sameH"); legPi->Draw();
  tv->cd(2);   gPad->SetLogy(); htofPi_1->Draw("sameH"); legPi->Draw();
  tv->cd(3);   gPad->SetLogy(); htofPi_15->Draw("sameH"); legPi->Draw();
  tv->cd(4);   gPad->SetLogy(); htofPi_25->Draw("sameH"); legPi->Draw();

  //show on canvas kaon ID
  TCanvas *tvK = new TCanvas("cK","cK",1200,500);
  tvK->Divide(2,2);
  TH1F* hQK_05 = (TH1F*) lProjNsigmaQK->FindObject(Form("proj_%i",id4tv[0]));
  TH1F* hQK_1 = (TH1F*) lProjNsigmaQK->FindObject(Form("proj_%i",id4tv[1]));
  TH1F* hQK_15 = (TH1F*) lProjNsigmaQK->FindObject(Form("proj_%i",id4tv[2]));
  TH1F* hQK_25 = (TH1F*) lProjNsigmaQK->FindObject(Form("proj_%i",id4tv[3]));
  hQK_05->GetXaxis()->SetRangeUser(-10.,10.);
  hQK_05->GetXaxis()->SetTitle("n#sigma_{K}^{TOF}");
  hQK_05->GetYaxis()->SetTitle("tracks");
  hQK_05->GetXaxis()->SetLabelSize(0.05);
  hQK_05->GetXaxis()->SetTitleSize(0.05);
  hQK_05->GetYaxis()->SetLabelSize(0.05);
  hQK_05->GetYaxis()->SetTitleSize(0.05);
  hQK_1->GetXaxis()->SetRangeUser(-10.,10.);
  hQK_1->GetXaxis()->SetTitle("n#sigma_{K}^{TOF}");
  hQK_1->GetYaxis()->SetTitle("tracks");
  hQK_1->GetXaxis()->SetLabelSize(0.05);
  hQK_1->GetXaxis()->SetTitleSize(0.05);
  hQK_1->GetYaxis()->SetLabelSize(0.05);
  hQK_1->GetYaxis()->SetTitleSize(0.05);
  hQK_15->GetXaxis()->SetRangeUser(-10.,10.);
  hQK_15->GetXaxis()->SetTitle("n#sigma_{K}^{TOF}");
  hQK_15->GetYaxis()->SetTitle("tracks");
  hQK_15->GetXaxis()->SetLabelSize(0.05);
  hQK_15->GetXaxis()->SetTitleSize(0.05);
  hQK_15->GetYaxis()->SetLabelSize(0.05);
  hQK_15->GetYaxis()->SetTitleSize(0.05);
  hQK_25->GetXaxis()->SetRangeUser(-10.,10.);
  hQK_25->GetXaxis()->SetTitle("n#sigma_{K}^{TOF}");
  hQK_25->GetYaxis()->SetTitle("tracks");
  hQK_25->GetXaxis()->SetLabelSize(0.05);
  hQK_25->GetXaxis()->SetTitleSize(0.05);
  hQK_25->GetYaxis()->SetLabelSize(0.05);
  hQK_25->GetYaxis()->SetTitleSize(0.05);
  legK->AddEntry(hQK_05,"std quality cuts","lpf");
  tvK->cd(1);   gPad->SetLogy(); hQK_05->Draw();
  tvK->cd(2);   gPad->SetLogy(); hQK_1->Draw();
  tvK->cd(3);   gPad->SetLogy(); hQK_15->Draw();
  tvK->cd(4);   gPad->SetLogy(); hQK_25->Draw();
      
  TH1F* htofK_05 = (TH1F*) lProjNsigmaTOFK->FindObject(Form("proj_%i",id4tv[0]));
  TH1F* htofK_1 = (TH1F*) lProjNsigmaTOFK->FindObject(Form("proj_%i",id4tv[1]));
  TH1F* htofK_15 = (TH1F*) lProjNsigmaTOFK->FindObject(Form("proj_%i",id4tv[2]));
  TH1F* htofK_25 = (TH1F*) lProjNsigmaTOFK->FindObject(Form("proj_%i",id4tv[3]));
  legK->AddEntry(htofK_05,"2#sigma_{K}^{TOF} && 5#sigma_{K}^{TPC} cut","lpf");
  tvK->cd(1);   htofK_05->Draw("sameH"); legK->Draw();
  tvK->cd(2);   htofK_1->Draw("sameH");  legK->Draw();
  tvK->cd(3);   htofK_15->Draw("sameH"); legK->Draw();
  tvK->cd(4);   htofK_25->Draw("sameH"); legK->Draw();
      
  //Display contamination from tpc pid to TOF PIONS
  TList * lProjNsigmaTPCPi_pi = (TList*) fin->Get(Form("proj_NsigmaTPC_%s_0",cutNames[1].Data()));
  TList * lProjNsigmaTPCPi_k = (TList*) fin->Get(Form("proj_NsigmaTPC_%s_1",cutNames[1].Data()));
  TList * lProjNsigmaTPCPi_p = (TList*) fin->Get(Form("proj_NsigmaTPC_%s_2",cutNames[1].Data()));
  //pion for tof, pion for tpc
  TH1F* htpcpi_Pi05 = (TH1F*) lProjNsigmaTPCPi_pi->FindObject(Form("proj_%i",id4tv[0]));
  TH1F* htpcpi_Pi1 = (TH1F*) lProjNsigmaTPCPi_pi->FindObject(Form("proj_%i",id4tv[1]));
  TH1F* htpcpi_Pi15 = (TH1F*) lProjNsigmaTPCPi_pi->FindObject(Form("proj_%i",id4tv[2]));
  TH1F* htpcpi_Pi25 = (TH1F*) lProjNsigmaTPCPi_pi->FindObject(Form("proj_%i",id4tv[3]));
  //pion for tof, kaon for tpc
  TH1F* htpck_Pi05 = (TH1F*) lProjNsigmaTPCPi_k->FindObject(Form("proj_%i",id4tv[0]));
  TH1F* htpck_Pi1 = (TH1F*) lProjNsigmaTPCPi_k->FindObject(Form("proj_%i",id4tv[1]));
  TH1F* htpck_Pi15 = (TH1F*) lProjNsigmaTPCPi_k->FindObject(Form("proj_%i",id4tv[2]));
  TH1F* htpck_Pi25 = (TH1F*) lProjNsigmaTPCPi_k->FindObject(Form("proj_%i",id4tv[3]));
    //pion for tof, proton for tpc
  TH1F* htpcp_Pi05 = (TH1F*) lProjNsigmaTPCPi_p->FindObject(Form("proj_%i",id4tv[0]));
  TH1F* htpcp_Pi1 = (TH1F*) lProjNsigmaTPCPi_p->FindObject(Form("proj_%i",id4tv[1]));
  TH1F* htpcp_Pi15 = (TH1F*) lProjNsigmaTPCPi_p->FindObject(Form("proj_%i",id4tv[2]));
  TH1F* htpcp_Pi25 = (TH1F*) lProjNsigmaTPCPi_p->FindObject(Form("proj_%i",id4tv[3]));
 //show on canvas kaon ID
  TCanvas *tvPiCont = new TCanvas("tvPiCont","tvPiCont",1200,500);
  tvPiCont->Divide(4,1);
  tvPiCont->cd(1); 
  gPad->SetLogy();
  htpcpi_Pi05->GetXaxis()->SetRangeUser(-5.,5.);
  htpcpi_Pi05->Draw();
  htpck_Pi05->Draw("same");
  htpcp_Pi05->Draw("same");

  tvPiCont->cd(2);
  gPad->SetLogy();
  htpcpi_Pi1->GetXaxis()->SetRangeUser(-5.,5.);
  htpcpi_Pi1->Draw();
  htpck_Pi1->Draw("same");
  htpcp_Pi1->Draw("same");

  tvPiCont->cd(3);
  gPad->SetLogy();
  htpcpi_Pi15->GetXaxis()->SetRangeUser(-5.,5.);
  htpcpi_Pi15->Draw();
  htpck_Pi15->Draw("same");
  htpcp_Pi15->Draw("same");

  tvPiCont->cd(4);
  gPad->SetLogy();
  htpcpi_Pi25->GetXaxis()->SetRangeUser(-5.,5.);
  htpcpi_Pi25->Draw();
  htpck_Pi25->Draw("same");
  htpcp_Pi25->Draw("same");

  //Display contamination from tpc pid to TOF KAONS
  TList * lProjNsigmaTPCK_pi = (TList*) fin->Get(Form("proj_NsigmaTPC_%s_0",cutNames[2].Data()));
  TList * lProjNsigmaTPCK_k = (TList*) fin->Get(Form("proj_NsigmaTPC_%s_1",cutNames[2].Data()));
  TList * lProjNsigmaTPCK_p = (TList*) fin->Get(Form("proj_NsigmaTPC_%s_2",cutNames[2].Data()));
  //pion for tof, pion for tpc
  TH1F* htpcpi_K05 = (TH1F*) lProjNsigmaTPCK_pi->FindObject(Form("proj_%i",id4tv[0]));
  TH1F* htpcpi_K1 = (TH1F*) lProjNsigmaTPCK_pi->FindObject(Form("proj_%i",id4tv[1]));
  TH1F* htpcpi_K15 = (TH1F*) lProjNsigmaTPCK_pi->FindObject(Form("proj_%i",id4tv[2]));
  TH1F* htpcpi_K25 = (TH1F*) lProjNsigmaTPCK_pi->FindObject(Form("proj_%i",id4tv[3]));
  //pion for tof, kaon for tpc
  TH1F* htpck_K05 = (TH1F*) lProjNsigmaTPCK_k->FindObject(Form("proj_%i",id4tv[0]));
  TH1F* htpck_K1 = (TH1F*) lProjNsigmaTPCK_k->FindObject(Form("proj_%i",id4tv[1]));
  TH1F* htpck_K15 = (TH1F*) lProjNsigmaTPCK_k->FindObject(Form("proj_%i",id4tv[2]));
  TH1F* htpck_K25 = (TH1F*) lProjNsigmaTPCK_k->FindObject(Form("proj_%i",id4tv[3]));
    //pion for tof, proton for tpc
  TH1F* htpcp_K05 = (TH1F*) lProjNsigmaTPCK_p->FindObject(Form("proj_%i",id4tv[0]));
  TH1F* htpcp_K1 = (TH1F*) lProjNsigmaTPCK_p->FindObject(Form("proj_%i",id4tv[1]));
  TH1F* htpcp_K15 = (TH1F*) lProjNsigmaTPCK_p->FindObject(Form("proj_%i",id4tv[2]));
  TH1F* htpcp_K25 = (TH1F*) lProjNsigmaTPCK_p->FindObject(Form("proj_%i",id4tv[3]));
 //show on canvas kaon ID
  TCanvas *tvKCont = new TCanvas("tvKCont","tvKCont",1200,500);
  tvKCont->Divide(4,1);
  tvKCont->cd(1); 
  gPad->SetLogy();
  htpcpi_K05->GetXaxis()->SetRangeUser(-5.,5.);
  htpcpi_K05->Draw();
  htpck_K05->Draw("same");
  htpcp_K05->Draw("same");

  tvKCont->cd(2);
  gPad->SetLogy();
  htpcpi_K1->GetXaxis()->SetRangeUser(-5.,5.);
  htpcpi_K1->Draw();
  htpck_K1->Draw("same");
  htpcp_K1->Draw("same");

  tvKCont->cd(3);
  gPad->SetLogy();
  htpcpi_K15->GetXaxis()->SetRangeUser(-5.,5.);
  htpcpi_K15->Draw();
  htpck_K15->Draw("same");
  htpcp_K15->Draw("same");

  tvKCont->cd(4);
  gPad->SetLogy();
  htpcpi_K25->GetXaxis()->SetRangeUser(-5.,5.);
  htpcpi_K25->Draw();
  htpck_K25->Draw("same");
  htpcp_K25->Draw("same");

  return;
}

