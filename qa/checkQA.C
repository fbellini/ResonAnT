/**********************************************/
/*  Macro to display and analyse the output   */
/*  of the AddMonitorOutput.C @ Rsn package   */
/*  fbellini@cern.ch,   04/12/2013            */
/**********************************************/
EColor col[9] = {kBlack, kGray+2, kYellow+3, kRed+1, kMagenta+2, kOrange-3, kBlue+1, kTeal-6, kViolet+6};

void MakeUp(TH1* h, TColor color, Int_t mark=20){
  if (!h) return;
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(mark);
  h->SetMarkerSize(0.6);
  return;
}

void checkQA(TString filename = "train1.root",  TString listsuffix = "cut2", Bool_t doTrk = 1, Bool_t  doTPCpid = 1, Bool_t doTOFpid=1)
{
  //gSystem->Load("libPWGLFresonances.so");
  // gROOT->LoadMacro("AddPaveText.C");
  // gROOT->LoadMacro("GetProfilePlot.C");
  // gROOT->LoadMacro("fitSlices.C");

  /********************************************************************/
  // Open input file with Rsn QA
  /********************************************************************/
  TFile * fin = TFile::Open(filename.Data());
  TString listname = Form("RsnOut_%s", listsuffix.Data());
  TList * list = (TList *) fin->Get(listname.Data());
  
  /********************************************************************/
  // Retrieve cut set from list name - needed for histos names
  /********************************************************************/
  Int_t cutPiCandidate = 0;  Int_t cutKaCandidate = 0;
  Float_t nsigmaPi = -1.0;  Float_t nsigmaKa = -1.0;
  Int_t aodFilterBit = 5;
  if (listname.Contains("cut2")){
    cutPiCandidate=2; cutKaCandidate=2; 
    nsigmaPi=3.0; nsigmaKa=3.0;
  }
  if (listname.Contains("comb")){
    cutPiCandidate=19; cutKaCandidate=19; 
    nsigmaPi=-1.0; nsigmaKa=-1.0;
  }
  if (listname.Contains("tpc2s")){
    cutPiCandidate=17; cutKaCandidate=17; 
    nsigmaPi=2.0; nsigmaKa=2.0;
  }
  if (listname.Contains("tof2s")){
    cutPiCandidate=18; cutKaCandidate=18; 
    nsigmaPi=2.0; nsigmaKa=2.0;
  }

  //Cut names
  TString cutQName =  Form("cutQ_bit%i",aodFilterBit);
  TString cutPiName = Form("cutPi%i_%2.1fsigma",cutPiCandidate, nsigmaPi);
  TString cutKaName = Form("cutK%i_%2.1fsigma",cutKaCandidate, nsigmaKa); 
  Char_t cutname[3][20];
  sprintf(cutname[0], "%s", cutQName.Data()); Printf("Cut 0 == %s", cutname[0]);
  sprintf(cutname[1], "%s", cutPiName.Data()); Printf("Cut 1 == %s", cutname[1]);
  sprintf(cutname[2], "%s", cutKaName.Data()); Printf("Cut 2 == %s", cutname[2]);
 
  /********************************************************************/
  // Retrieve histos from list
  /********************************************************************/

  //Events QA
  TH1F*	hEventStat          = (TH1F*) list->FindObject("hEventStat"); 
  TH1F*	hAEventsVsMulti     = (TH1F*) list->FindObject("hAEventsVsMulti");
  TH2F*	hVzVsCent           = (TH2F*) list->FindObject("hVzVsCent");
  TH2F*	hMultiVsCent        = (TH2F*) list->FindObject("MultiVsCent");
  TH1F*	hAcceptedEventVtx   = (TH1F*) list->FindObject(Form("TOFKStarPbPbData_%i%i_eventVtx",cutPiCandidate,cutKaCandidate));
  TH1F* hAcceptedEventMulti = (TH1F*) list->FindObject(Form("TOFKStarPbPbData_%i%i_eventMult",cutPiCandidate,cutKaCandidate));
  TH2F*	hPtVsMult_pt_multi  = (TH2F*) list->FindObject(Form("%s.PtVsMult_pt_mult", cutname[0]));

  //Single-track QA plots
  TH2F*	 hP_charge[3], hPt_charge[3], hEta_charge[3], hPhi_charge[3], hPhiOuterTPC_charge[3], hDCAxyVsPt[3], hDCAzVsP[3], hITSclsVsPt[3], hTPCclsVsPt[3], hTPCclsVsPtpc[3], hITSchi2VsPt[3], hTPCchi2VsPt[3], hTPCchi2VsPtpc[3];
  TH3F * hPhiVsPt_charge[3];

  TH2F * hdEdxVsPtpc[3], hTPCnsigmaPiVsPtpc[3], hTPCnsigmaKVsPtpc[3], hTPCnsigmaProVsPtpc[3];
  TH2F * hTOFnsigmaPiVsTPCnsigmaPi[3], hTOFnsigmaPiVsP[3], hTOFnsigmaKVsP[3], hTOFnsigmaProVsP[3];

  TH1D * hPt[3], hPt_pos[3], hPt_neg[3], hPt_neg2pos[3], hEta[3], hEta_pos[3], hEta_neg[3], hEta_neg2pos[3], hPhi[3],hPhi_pos[3], hPhi_neg[3], hPhi_neg2pos[3];
  //hP[3], hPhiOuterTPC[3], hDCAxy[3], hDCAz[3], hITSclsVsPt[3], hTPCclsVsPt[3], hTPCclsVsPtpc[3], hITSchi2VsPt[3], hTPCchi2VsPt[3], hTPCchi2VsPtpc[3];  TH3F * hPhiVsPt_charge[3];


 
  TCanvas *cptqa = new TCanvas("cptqa","cptqa", 800,600);
  cptqa->Divide(4,2);
  TCanvas *cetaqa = new TCanvas("cetaqa","cetaqa", 800,600);
  cetaqa->Divide(4,2);
  TCanvas *cphiqa = new TCanvas("cphiqa","cphiqa", 800,600);
  cphiqa->Divide(4,2);
  
  for (Int_t j=0;j<3;j++) {  
    Printf("Doing cut %i", j);
    TString opt = ((j>0)?"same":"");
    //Tracking QA
    if (doTrk) {
      Printf("Doing cut %i tracking", j);
      hP_charge[j] = (TH2F*) list->FindObject(Form("%s.P_p_charge", cutname[j]));
      
      /***** Pt *****/
      hPt_charge[j] = (TH2F*) list->FindObject(Form("%s.Pt_pt_charge", cutname[j]));
      if (hPt_charge[j]) {
	hPt[j] = (TH1D*) hPt_charge[j]->ProjectionX(Form("hPt%i",j));
	hPt_pos[j] = (TH1D*) hPt_charge[j]->ProjectionX(Form("hPt%i_pos",j),3,3);
	hPt_neg[j] = (TH1D*) hPt_charge[j]->ProjectionX(Form("hPt%i_neg",j),1,1);
	gROOT->LoadMacro("$ASD/GetPlotRatio.C");
	hPt_neg2pos[j]  = (TH1D*) GetPlotRatio(hPt_neg[j], hPt_pos[j], 0, "neg", "pos", "entries");// 0, 0, "", 0.0, -1.e10, "sum2");
	hPt_neg2pos[j]->SetName(Form("hPt%i_neg2pos",j));
	MakeUp(hPt[j], col[j*3], 1);
	MakeUp(hPt_neg[j], col[j*3+1], 28);
	MakeUp(hPt_pos[j], col[j*3+2], 20);
	MakeUp(hPt_neg2pos[j], col[j*3], 21);
	cptqa->cd(1); hPt[j]->Draw(opt.Data()); 
	cptqa->cd(1+j); hPt[j]->Draw(); hPt_pos[j]->Draw("same"); hPt_neg[j]->Draw("same");
	cptqa->cd(6+j); hPt_neg2pos[j]->Draw(); 
      } else {Printf("failing to get Pt plot for cut #i",j);}
      /***** Eta *****/
      hEta_charge[j] = (TH2F*) list->FindObject(Form("%s.Eta_eta_charge", cutname[j]));
      if (hEta_charge[j]) {
	hEta[j]     = (TH1F*) hEta_charge[j]->ProjectionX(Form("hEta%i",j),1,3);
	hEta_pos[j] = (TH1F*) hEta_charge[j]->ProjectionX(Form("hEta%i_pos",j),3,3);
	hEta_neg[j] = (TH1F*) hEta_charge[j]->ProjectionX(Form("hEta%i_neg",j),1,1);
	hEta_neg2pos[j]  = (TH1D*) GetPlotRatio(hEta_neg[j], hEta_pos[j], 0, "neg", "pos", "entries", 0, 0, "", 0.0, -1.e10, "sum2");
	hEta_neg2pos[j]->SetName(Form("hEta%i_neg2pos",j));
	MakeUp(hEta[j], col[j*3], 1);
	MakeUp(hEta_neg[j], col[j*3+1], 28);
	MakeUp(hEta_pos[j], col[j*3+2], 20);
	MakeUp(hEta_neg2pos[j], col[j*3], 21);
	cetaqa->cd(1); hEta[j]->Draw(opt.Data()); 
	cetaqa->cd(1+j); hEta[j]->Draw(); hEta_pos[j]->Draw("same"); hEta_neg[j]->Draw("same");
	cetaqa->cd(6+j); hEta_neg2pos[j]->Draw(); 
      } else {Printf("failing to get Eta plot for cut #i",j);}
      /***** Phi *****/
      hPhi_charge[j]  = (TH2F*) list->FindObject(Form("%s.Phi_phi_charge", cutname[j]));
      if (hPhi_charge[j]) {
	hPhi[j]     = (TH1F*) hPhi_charge[j]->ProjectionX(Form("hPhi%i",j),1,3);
	hPhi_pos[j] = (TH1F*) hPhi_charge[j]->ProjectionX(Form("hPhi%i_pos",j),3,3);
	hPhi_neg[j] = (TH1F*) hPhi_charge[j]->ProjectionX(Form("hPhi%i_neg",j),1,1);
	hPhi_neg2pos[j]  = (TH1D*) GetPlotRatio(hPhi_neg[j], hPhi_pos[j], 0, "neg", "pos", "entries", 0, 0, "", 0.0, -1.e10, "sum2");
	hPhi_neg2pos[j]->SetName(Form("hPhi%i_neg2pos",j));
	MakeUp(hPhi[j], col[j*3], 1);
	MakeUp(hPhi_neg[j], col[j*3+1], 28);
	MakeUp(hPhi_pos[j], col[j*3+2], 20);
	MakeUp(hPhi_neg2pos[j], col[j*3], 21);
	cphiqa->cd(1); hPhi[j]->Draw(opt.Data()); 
	cphiqa->cd(1+j); hPhi[j]->Draw(); hPhi_pos[j]->Draw("same"); hPhi_neg[j]->Draw("same");
	cphiqa->cd(6+j); hPhi_neg2pos[j]->Draw(); 
      }  else {Printf("failing to get Phi plot for cut #i",j);}

      hPhiOuterTPC_charge[j] = (TH2F*) list->FindObject(Form("%s.PhiOuterTPC_phiOuterTPC_charge", cutname[j]));
      hPhiVsPt_charge[j]     = (TH3F*) list->FindObject(Form("%s.PhiVsPt_pt_phi_charge", cutname[j])); 
      hDCAxyVsPt[j]          = (TH2F*) list->FindObject(Form("%s.DCAxyVsPt_pt_DCAxy", cutname[j]));
      hDCAzVsP[j]            = (TH2F*) list->FindObject(Form("%s.DCAzVsP_p_DCAz", cutname[j]));
      hITSclsVsPt[j]         = (TH2F*) list->FindObject(Form("%s.ITSclsVsPt_pt_ITScls", cutname[j]));
      hTPCclsVsPt[j]         = (TH2F*) list->FindObject(Form("%s.TPCclsVsPt_pt_TPCcls", cutname[j]));
      hTPCclsVsPtpc[j]       = (TH2F*) list->FindObject(Form("%s.TPCclsVsPtpc_pTPC_TPCcls", cutname[j]));
      hITSchi2VsPt[j]        = (TH2F*) list->FindObject(Form("%s.ITSchi2VsPt_pt_ITSchi2", cutname[j]));
      hTPCchi2VsPt[j]        = (TH2F*) list->FindObject(Form("%s.TPCchi2VsPt_pt_TPCchi2", cutname[j]));
      hTPCchi2VsPtpc[j]      = (TH2F*) list->FindObject(Form("%s.TPCchi2VsPtpc_pTPC_TPCchi2", cutname[j]));
    }
    if (doTPCpid) {
      //PID QA - TPC
      Printf("Doing cut %i TPC pid", j);
      hdEdxVsPtpc[j]         = (TH2F*) list->FindObject(Form("%s.dEdx_VsPtpc_pTPC_sTPC", cutname[j]));
      hTPCnsigmaPiVsPtpc[j]  = (TH2F*) list->FindObject(Form("%s.TPC_nsigmaPi_VsPtpc_pTPC_pi", cutname[j]));
      hTPCnsigmaKVsPtpc[j]   = (TH2F*) list->FindObject(Form("%s.TPC_nsigmaK_VsPtpc_pTPC_K", cutname[j]));
      hTPCnsigmaProVsPtpc[j] = (TH2F*) list->FindObject(Form("%s.TPC_nsigmaPro_VsPtpc_pTPC_p", cutname[j]));
    }
    if (doTOFpid) {
      //PID QA - TOF
      Printf("Doing cut %i TOF pid", j);
      hTOFnsigmaPiVsTPCnsigmaPi[j] = (TH2F*) list->FindObject(Form("%s.TOFnsigmaPi_TPCnsigmaPi_pi_pi", cutname[j]));
      hTOFnsigmaPiVsP[j]           = (TH2F*) list->FindObject(Form("%s.TOF_nsigmaPi_vsP_p_pi", cutname[j]));
      hTOFnsigmaKVsP[j]            = (TH2F*) list->FindObject(Form("%s.TOF_nsigmaK_vsP_p_K", cutname[j]));
      hTOFnsigmaProVsP[j]          = (TH2F*) list->FindObject(Form("%s.TOF_nsigmaPro_vsP_p_p", cutname[j]));
    }
  }
  //output file
  //  TFile * fout = new TFile("checkRsnQA.root","recreate");
 
  /********************************************************************/
  /************************ PT DISTRIBS *******************************/
  /********************************************************************/
  /*
  //primaries distribs - projections
  TH2D * hpt = (TH2D *) ltrd->FindObject("cutQ_bit5.Pt_pt_charge");
  TH1D * hpt_neg = (TH1D *) hpt->ProjectionX("hpt_neg", 0, 1, "eo"); 
  TH1D * hpt_pos = (TH1D *) hpt->ProjectionX("hpt_pos", 3, 4, "eo"); 
  hpt_neg->Rebin(10);
  hpt_pos->Rebin(10);
  hpt_neg->SetTitle("primary negative");
  hpt_pos->SetTitle("primary positive");
  fout->cd();
  hpt_neg->Write();
  hpt_pos->Write();
  // Printf("********************************************************************");
  // Printf("ALL positive = %e \n negative = %e \n totali = %e", hpt_pos->GetIntegral(), hpt_neg->GetIntegral(), hpt->GetIntegral());

  // trd distribs - projections
  TH2D * hptKtrd = (TH2D *) ltrd->FindObject("cutK9_-1.0sigma.Pt_pt_charge");
  TH1D * hptKtrd_neg = (TH1D *) hptKtrd->ProjectionX("hptKtrd_neg", 0, 1, "eo");
  TH1D * hptKtrd_pos = (TH1D *) hptKtrd->ProjectionX("hptKtrd_pos", 3, 4, "eo"); 
  hptKtrd_neg->SetTitle("negative, TRD");
  hptKtrd_pos->SetTitle("positive, TRD");
  hptKtrd_neg->Rebin(10);
  hptKtrd_pos->Rebin(10);

  */
  return;
}



TString ltrdn = "RsnOut_TofMatch_Trd";
TString lnotrdn = "RsnOut_TofMatch_NoTrd";
TString l2sn = "RsnOut_Tof20sigma";
TString l2strdn = "RsnOut_Tof2s_KaTrd";
TString l2snotrdn = "RsnOut_Tof2s_KaNoTrd";


void checkQATrd(TString fn = "AnalysisResults_Trd.root")
{
  //analyse output from rsn
  TFile * fin = TFile::Open(fn.Data());
  TList * ltrd = (TList *) fin->Get(ltrdn.Data());
  TList * lnotrd = (TList *) fin->Get(lnotrdn.Data());
  
  //output file
  TFile * fout = new TFile("out.root","recreate");

 
  /********************************************************************/
  /************************ PT DISTRIBS *******************************/
  /********************************************************************/

  //primaries distribs - projections
  TH2D * hpt = (TH2D *) ltrd->FindObject("cutQ_bit5.Pt_pt_charge");
  TH1D * hpt_neg = (TH1D *) hpt->ProjectionX("hpt_neg", 0, 1, "eo"); 
  TH1D * hpt_pos = (TH1D *) hpt->ProjectionX("hpt_pos", 3, 4, "eo"); 
  hpt_neg->Rebin(10);
  hpt_pos->Rebin(10);
  hpt_neg->SetTitle("primary negative");
  hpt_pos->SetTitle("primary positive");
  fout->cd();
  hpt_neg->Write();
  hpt_pos->Write();
  // Printf("********************************************************************");
  // Printf("ALL positive = %e \n negative = %e \n totali = %e", hpt_pos->GetIntegral(), hpt_neg->GetIntegral(), hpt->GetIntegral());

  // trd distribs - projections
  TH2D * hptKtrd = (TH2D *) ltrd->FindObject("cutK9_-1.0sigma.Pt_pt_charge");
  TH1D * hptKtrd_neg = (TH1D *) hptKtrd->ProjectionX("hptKtrd_neg", 0, 1, "eo");
  TH1D * hptKtrd_pos = (TH1D *) hptKtrd->ProjectionX("hptKtrd_pos", 3, 4, "eo"); 
  hptKtrd_neg->SetTitle("negative, TRD");
  hptKtrd_pos->SetTitle("positive, TRD");
  hptKtrd_neg->Rebin(10);
  hptKtrd_pos->Rebin(10);
  //  Printf("TRD positive = %e \n negative = %e \n totali = %e", hptKtrd_pos->GetIntegral(), hptKtrd_neg->GetIntegral(), hptKtrd->GetIntegral());
  TH1D * rpt_negOverPos = (TH1D *) hpt_neg->Clone("rpt_negOverPos");
  rpt_negOverPos->Divide(hpt_pos);
  fout->cd();
  hptKtrd_neg->Write();
  hptKtrd_pos->Write();
  rpt_negOverPos->Write();
  
  // trd neg/pos
  TH1D * rKtrd_pos = (TH1D *) hptKtrd_pos->Clone("rKtrd_pos");
  rKtrd_pos->Divide(hpt_pos);
  TH1D * rKtrd_neg = (TH1D *) hptKtrd_neg->Clone("rKtrd_neg");
  rKtrd_neg->Divide(hpt_neg);
  TH1D * rKtrd_negOverPos = (TH1D *) hptKtrd_neg->Clone("rKtrd_negOverPos");
  rKtrd_negOverPos->Divide(hptKtrd_pos);

  rKtrd_neg->SetTitle("negative matched/primary, TRD");
  rKtrd_pos->SetTitle("positive matched/primary, TRD");
  rKtrd_negOverPos->SetTitle("negative/positive matched, TRD");

  fout->cd();
  rKtrd_pos->Write();
  rKtrd_neg->Write();
  rKtrd_negOverPos->Write();

  // no trd distribs - projections
  TH2D * hptKnotrd = (TH2D *) lnotrd->FindObject("cutK10_-1.0sigma.Pt_pt_charge");
  TH1D * hptKnotrd_neg = (TH1D *) hptKnotrd->ProjectionX("hptKnotrd_neg", 0, 1, "eo"); 
  TH1D * hptKnotrd_pos = (TH1D *) hptKnotrd->ProjectionX("hptKnotrd_pos", 3, 4, "eo"); 
  hptKnotrd_neg->SetTitle("negative, no TRD");
  hptKnotrd_pos->SetTitle("positive, no TRD");
  hptKnotrd_neg->Rebin(10);
  hptKnotrd_pos->Rebin(10);
    //Printf("NO TRD positive = %e \n negative = %e \n totali = %e", hptKnotrd_pos->GetIntegral(), hptKnotrd_neg->GetIntegral(), hptKnotrd->GetIntegral());
  fout->cd();
  hptKnotrd_neg->Write();
  hptKnotrd_pos->Write();

  TH1D * rKnotrd_pos = (TH1D *) hptKnotrd_pos->Clone("rKnotrd_pos");
  rKnotrd_pos->Divide(hpt_pos);
  TH1D * rKnotrd_neg = (TH1D *) hptKnotrd_neg->Clone("rKnotrd_neg");
  rKnotrd_neg->Divide(hpt_neg);
  TH1D * rKnotrd_negOverPos = (TH1D *) hptKnotrd_neg->Clone("rKnotrd_negOverPos");
  rKnotrd_negOverPos->Divide(hptKnotrd_pos);
  rKnotrd_neg->SetTitle("negative matched/primary, no TRD");
  rKnotrd_pos->SetTitle("positive matched/primary, no TRD");
  rKnotrd_negOverPos->SetTitle("negative/positive matched, no TRD");

  fout->cd();
  rKnotrd_pos->Write();
  rKnotrd_neg->Write();
  rKnotrd_negOverPos->Write();

  //matched neg/pos
  TH1D * hptmatch_neg = (TH1D*) hptKnotrd_neg->Clone("hptmatch_neg");
  hptmatch_neg->Add(hptKtrd_neg,1.0);
  TH1D * hptmatch_pos = (TH1D*) hptKnotrd_pos->Clone("hptmatch_pos");
  hptmatch_pos->Add(hptKtrd_pos,1.0);  
 
  TH1D * rmatch_negOverPos = (TH1D *) hptmatch_neg->Clone("rmatch_negOverPos");
  rmatch_negOverPos->Divide(hptmatch_pos);

  hptmatch_neg->SetTitle("negative matched, TRD+no TRD");
  hptmatch_pos->SetTitle("positive matched/primary, TRD+no TRD");
  rmatch_negOverPos->SetTitle("negative/positive matched, TRD+no TRD");

  fout->cd();
  hptmatch_neg->Write();
  hptmatch_pos->Write();
  rmatch_negOverPos->Write();
  
  /********************************************************************/
  /************************ Phi DISTRIBS ******************************/
  /********************************************************************/
    //primaries distribs - projections
  TH2D * hphiTPC = (TH2D *) ltrd->FindObject("cutQ_bit5.PhiOuterTPC_phiOuterTPC_charge");
  TH1D * hphiTPC_neg = (TH1D *) hphiTPC->ProjectionX("hphiTPC_neg", 0, 1, "eo"); 
  TH1D * hphiTPC_pos = (TH1D *) hphiTPC->ProjectionX("hphiTPC_pos", 3, 4, "eo"); 
  hphiTPC_pos->Rebin(20);
  hphiTPC_neg->Rebin(20);
  TH1D * rphiTPC_negOverPos = (TH1D *) hphiTPC_neg->Clone("rphiTPC_negOverPos");
  rphiTPC_negOverPos->Divide(hphiTPC_pos);
  hphiTPC_neg->SetTitle("primary negative");
  hphiTPC_pos->SetTitle("primary positive");
  rphiTPC_negOverPos->SetTitle("primary negative/positive");
  
  fout->cd();
  hphiTPC_neg->Write();
  hphiTPC_pos->Write();
  rphiTPC_negOverPos->Write();

  // trd distribs - projections
  TH2D * hphiTPCtrd = (TH2D *) ltrd->FindObject("cutPi9_-1.0sigma.PhiOuterTPC_phiOuterTPC_charge");
  TH1D * hphiTPCtrd_neg = (TH1D *) hphiTPCtrd->ProjectionX("hphiTPCtrd_neg", 0, 1, "eo");
  TH1D * hphiTPCtrd_pos = (TH1D *) hphiTPCtrd->ProjectionX("hphiTPCtrd_pos", 3, 4, "eo"); 
  hphiTPCtrd_pos->Rebin(20);
  hphiTPCtrd_neg->Rebin(20);
  TH1D * rphiTPCtrd_negOverPos = (TH1D *) hphiTPCtrd_neg->Clone("rphiTPCtrd_negOverPos");
  rphiTPCtrd_negOverPos->Divide(hphiTPCtrd_pos);
  hphiTPCtrd_neg->SetTitle("negative matched, TRD");
  hphiTPCtrd_pos->SetTitle("positive matched, TRD");
  rphiTPCtrd_negOverPos->SetTitle("negative/positive matched, TRD");
  fout->cd();
  hphiTPCtrd_neg->Write();
  hphiTPCtrd_pos->Write();
  rphiTPCtrd_negOverPos->Write();
  
 // trd / primaries (for pos and neg)
  TH1D * rphiTPCtrd_pos = (TH1D *) hphiTPCtrd_pos->Clone("rphiTPCtrd_pos");
  rphiTPCtrd_pos->Divide(hphiTPC_pos);
  TH1D * rphiTPCtrd_neg = (TH1D *) hphiTPCtrd_neg->Clone("rphiTPCtrd_neg");
  rphiTPCtrd_neg->Divide(hphiTPC_neg);
  rphiTPCtrd_pos->SetTitle("matched in TRD/primary, positive");
  rphiTPCtrd_neg->SetTitle("matched in TRD/primary, negative");

  fout->cd();
  rphiTPCtrd_neg->Write();
  rphiTPCtrd_pos->Write();
  
   // no trd distribs - projections
  TH2D * hphiTPCnotrd = (TH2D *) lnotrd->FindObject("cutPi10_-1.0sigma.PhiOuterTPC_phiOuterTPC_charge");
  TH1D * hphiTPCnotrd_neg = (TH1D *) hphiTPCnotrd->ProjectionX("hphiTPCnotrd_neg", 0, 1, "eo"); 
  TH1D * hphiTPCnotrd_pos = (TH1D *) hphiTPCnotrd->ProjectionX("hphiTPCnotrd_pos", 3, 4, "eo");
  hphiTPCnotrd_pos->Rebin(20);
  hphiTPCnotrd_neg->Rebin(20);
  TH1D * rphiTPCnotrd_negOverPos = (TH1D *) hphiTPCnotrd_neg->Clone("rphiTPCnotrd_negOverPos");
  rphiTPCnotrd_negOverPos->Divide(hphiTPCnotrd_pos);
  hphiTPCnotrd_neg->SetTitle("negative matched, no TRD");
  hphiTPCnotrd_pos->SetTitle("positive matched, no TRD");
  rphiTPCnotrd_negOverPos->SetTitle("negative/positive matched, no TRD");
 
  fout->cd();
  hphiTPCnotrd_neg->Write();
  hphiTPCnotrd_pos->Write();
  rphiTPCnotrd_negOverPos->Write();

  // no trd / primaries (for pos and neg)
  TH1D * rphiTPCnotrd_pos = (TH1D *) hphiTPCnotrd_pos->Clone("rphiTPCnotrd_pos");
  rphiTPCnotrd_pos->Divide(hphiTPC_pos);
  TH1D * rphiTPCnotrd_neg = (TH1D *) hphiTPCnotrd_neg->Clone("rphiTPCnotrd_neg");
  rphiTPCnotrd_neg->Divide(hphiTPC_neg);
  rphiTPCnotrd_pos->SetTitle("matched no TRD/primary, positive");
  rphiTPCnotrd_neg->SetTitle("matched no TRD/primary, negative");
  fout->cd();
  rphiTPCnotrd_neg->Write();
  rphiTPCnotrd_pos->Write();
 
  //matched neg/pos
  TH1D * hphiTPCmatch_neg = (TH1D*) hphiTPCnotrd_neg->Clone("hphiTPCmatch_neg");
  hphiTPCmatch_neg->Add(hphiTPCtrd_neg,1.0);
  TH1D * hphiTPCmatch_pos = (TH1D*) hphiTPCnotrd_pos->Clone("hphiTPCmatch_pos");
  hphiTPCmatch_pos->Add(hphiTPCtrd_pos,1.0);  
  TH1D * rphiTPCmatch_negOverPos = (TH1D *) hphiTPCmatch_neg->Clone("rphiTPCmatch_negOverPos");
  rphiTPCmatch_negOverPos->Divide(hphiTPCmatch_pos);
  hphiTPCmatch_neg->SetTitle("negative matched, TRD+no TRD");
  hphiTPCmatch_pos->SetTitle("positive matched/primary, TRD+no TRD");
  rphiTPCmatch_negOverPos->SetTitle("negative/positive matched, TRD+no TRD");
  fout->cd();
  hphiTPCmatch_neg->Write();
  hphiTPCmatch_pos->Write();
  rphiTPCmatch_negOverPos->Write();

  return;
}
