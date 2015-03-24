//
// Splits a TH2F using the myHistSplig_sparse_IM_PT_CENT class
//
Bool_t kComputeEfficiency=1;
const Int_t kNhistosData=4;
TString macroDir = "$ASD/ResonAnT/"; //"/Users/bellini/alice/macro/kstar";

void projectKStarMCtrue
(
   const char *nameData = "train1220.root",
   TString listName = "RsnOut_tpc2s_tof3sveto",
   const char *outName  = "tpc2s_tof3sveto",
   char * cutID = "2424",
   Color_t customColor = kBlack,
   Short_t ARver = 69,
   Bool_t isPP = 0,
   Bool_t isMC = 1,
   Bool_t isOnlyKstar = 0,
   Bool_t considerTRD = 0
)
{
   // initial setup
  gStyle->SetOptStat("1111");
  gStyle->SetTextFont(42);
  //Color_t color[]={kRed+1, kOrange+1, kGreen+2, kBlue, kViolet-6};
  Color_t color[]={kOrange+7, kPink+6, kGreen+1, kAzure+1, kBlue+4};
  Color_t marker[]={21, 22, 32, 28, 24};
  Color_t histoColor = customColor;
  
   // load utility macro
  gROOT->LoadMacro(Form("%s/projectorInvMass_Centrality_Pt.C+g",macroDir.Data()));
  
  TFile *dummy = new TFile(Form("projMC_%s",nameData),"RECREATE");
  TList * out = new TList();
  
  // output lists - define once and then clear them
  TList lTruesPM;
  TList lTruesMP;
  lTruesPM.SetName(Form("%s_TruesPM", (!isMC ? "Data" : "MC")));
  lTruesMP.SetName(Form("%s_TruesMP", (!isMC ? "Data" : "MC")));
  
  TList lMother,lAntiMother;
  lMother.SetName(Form("%s_Mother", (!isMC ? "Data" : "MC")));
  lAntiMother.SetName(Form("%s_AntiMother", (!isMC ? "Data" : "MC")));
  
  // open input file
  TFile *fileData;// = 0x0;//, *fileMC = 0x0;
  TList *listData;
  TList *listDataAnti;// = 0x0;//, *listMC = 0x0;
  TString antiListName = listName.Data();
  
  fileData = TFile::Open(nameData);
  
  if (!(fileData && fileData->IsOpen())) { 
    Printf("ERROR: cannot open file");
    return; 
  }
    
  listData = (TList*)fileData->Get(listName.Data());
  if (!isOnlyKstar) {
    if (ARver<69) antiListName.ReplaceAll("RsnOut_","RsnOut_anti_");  //needed for AR version < v5-05-69-AN 
    listDataAnti = (TList*)fileData->Get(antiListName.Data());
    Printf("AntiMother list = %s", listDataAnti->GetName());
  }
  
  // names of computed histograms (one per settings)
  THnSparse* hInput[kNhistosData]= {0,0,0,0};
  if (listData) {
    
    hInput[ 0] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_TruesPM", (isPP ? "pp" : "PbPb"), (isMC? "MC" : "Data"), cutID));
    if (ARver<69)
      hInput[ 2] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_kstar_Mother", (isPP ? "pp" : "PbPb"), (isMC? "MC" : "Data"), cutID));
    else
      hInput[ 2] = (THnSparse*)listData->FindObject(Form("TOFKStar%s%s_%s_Ks_Mother", (isPP ? "pp" : "PbPb"), (isMC? "MC" : "Data"), cutID));
    
    hInput[ 0]->SetName("TruesPM");
    hInput[ 2]->SetName("Mother");
    
    if (!isOnlyKstar){
      hInput[ 1] = (THnSparse*)listDataAnti->FindObject(Form("TOFKStar%s%s_%s_kstar_TruesMP", (isPP ? "pp" : "PbPb"), (isMC? "MC" : "Data"), cutID));
      if (ARver<69)
	hInput[ 3] = (THnSparse*)listDataAnti->FindObject(Form("TOFKStar%s%s_%s_kstar_Mother", (isPP ? "pp" : "PbPb"), (isMC? "MC" : "Data"), cutID));
      else
	hInput[ 3] = (THnSparse*)listDataAnti->FindObject(Form("TOFKStar%s%s_%s_antiKs_Mother", (isPP ? "pp" : "PbPb"), (isMC? "MC" : "Data"), cutID));
      
      hInput[ 1]->SetName("TruesMP");
      hInput[ 3]->SetName("AntiMother");
    }
  }
  
  //trd check
  //Double_t pt[] = {1.0, 2.00, 3.00, 4.00, 5.0, 7.0};
   
   // define binning in pT    
   //old ana  Double_t pt[] = { 0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.00 };
   //new ana signed Double_t pt[] = { 0.0, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.00 };
   //integrated Double_t pt[] = { 1.0, 10.0}; 
   //old tpc ana Double_t pt[] = { 0.0, 0.5, 1.0, 1.50, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0,10.0};
   
   //trd ana signed 
   //Double_t pt[] = {1.0, 2.0, 3.0, 4.0, 5.0, 7.0};
  
   //pA analysis
  Double_t cent[]={0., 100.0}; 
  //Double_t cent[]={0., 20., 40., 60., 80., 100.};   
   //TOF and TPC standalone analyses
   //Double_t pt[] = {0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 12.0, 14.0, 16.0};
   //Double_t pt[] = { 0.0, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };
   //define binning in pT  - 300MeV bins - binning A 
   //Double_t pt[] = {0.0, 0.15, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };
   //binning B
   Double_t pt[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00, 12.0, 15.0, 20.0};
   //define binning in pT  - 200MeV bins - binning C
   //Double_t pt[] = {0.0, 0.1, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 10.00 };  
   
    //binning of pp 7 TeV published
    //Double_t pt[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10.00, 12.0, 15.0, 20.0};

   Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;
   TString centLabel;
   
   //define binning in centrality
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
   for (Int_t icentbin=0; icentbin<ncent;icentbin++){
     if (ncent<5) {
       c[icentbin]= new TCanvas(Form("c_%i",icentbin),Form("cent_%i",icentbin),1200,600);
     }
   } 
   //set title for legends
   TString titleEff(listName);
   titleEff.ReplaceAll("RsnOut_","");
   titleEff.ReplaceAll("_"," ");
   

   //save separate files for each centrality bin
   for (Int_t icentbin=0; icentbin<ncent;icentbin++){
     centLabel=Form("(%3.0f-%3.0f%%)",cent[icentbin], cent[icentbin+1]);
     if (customColor == kBlack) histoColor = color[icentbin];
     
     c[icentbin]->Divide(5,5);
     lTruesPM.Clear();
     lMother.Clear();
     FillListForCentrality(out,&lTruesPM,icentbin);
     FillListForCentrality(out,&lMother,icentbin);
     if (!isOnlyKstar){
       lTruesMP.Clear();
       lAntiMother.Clear();
       FillListForCentrality(out,&lTruesMP,icentbin);
       FillListForCentrality(out,&lAntiMother,icentbin);
     }
     
     //output file
     TString efffilename = Form("efficiency_%s_cent%03.0f-%03.0f.root",listName.Data(),cent[icentbin], cent[icentbin+1]);
     TFile *efffile = TFile::Open(efffilename.Data(), "RECREATE");
     //Create histos for efficiency
     TH1F * hTrueCountsPM = new TH1F("hTrueCountsPM","True Counts K+#pi-", npt, pt);
     TH1F * hTrueCountsMP = new TH1F("hTrueCountsMP","True Counts K-#pi+", npt, pt);
     TH1F * hMotherCounts = new TH1F("hMotherCounts","hMotherCounts", npt, pt);
     TH1F * hAntiMotherCounts = new TH1F("hAntiMotherCounts","hAntiMotherCounts", npt, pt);
     
     TH1F * hEffVsPt = new TH1F("hEffVsPt", "; p_{T} (GeV/#it{c}); #epsilon = reco true / generated", npt, pt);
     hEffVsPt->SetTitle(Form("#epsilon(K*+#bar{K*}), %s %s", titleEff.Data(), isPP? "(min. bias)": centLabel.Data()));
     hEffVsPt->SetLineColor(histoColor);
     hEffVsPt->SetMarkerColor(histoColor);
     hEffVsPt->SetLineWidth(2);
     hEffVsPt->SetMarkerStyle(20);
     hEffVsPt->SetMarkerSize(0.8);

     TH1F * hEffVsPtKstar = new TH1F("hEffVsPtKstar", "Efficiency for K*; p_{T} (GeV/#it{c}); #epsilon = reco true / generated", npt, pt);
     hEffVsPtKstar->SetTitle(Form("#epsilon(K*), %s %s", titleEff.Data(), centLabel.Data()));
     hEffVsPtKstar->SetLineColor(histoColor);
     hEffVsPtKstar->SetMarkerColor(histoColor);
     hEffVsPtKstar->SetMarkerStyle(20);

    TH1F * hEffVsPtAntiKstar = new TH1F("hEffVsPtAntiKstar", "Efficiency for #bar{K*}; p_{T} (GeV/#it{c}); #epsilon = reco true / generated", npt, pt);
    hEffVsPtAntiKstar->SetTitle(Form("#epsilon(#bar{K*}), %s %s", titleEff.Data(), centLabel.Data()));
     hEffVsPtAntiKstar->SetLineColor(histoColor);
     hEffVsPtAntiKstar->SetMarkerColor(histoColor);
     hEffVsPtAntiKstar->SetMarkerStyle(24);
     
     TH1F* ratioAKstarOverKstar = new TH1F("hRatioAKstarOverKstar", "ratio; p_{T} (GeV/#it{c}); ratio = #epsilon(#bar{K*}) / #epsilon(K*)", npt, pt); 
      
     hTrueCountsPM->Sumw2();
     hTrueCountsMP->Sumw2();
     hMotherCounts->Sumw2();
     hAntiMotherCounts->Sumw2();
     hEffVsPt->Sumw2();
     hEffVsPtKstar->Sumw2();
     hEffVsPtAntiKstar->Sumw2();

     //cosmetics
     hEffVsPt->GetYaxis()->SetRangeUser(0.0,1.0);
     hEffVsPt->SetLineColor(histoColor);
     hEffVsPt->SetMarkerColor(histoColor);

     //make-up for true M spectra
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

       ((TH1F*)lMother.At(ih))->SetLineColor(kCyan-2);
       ((TH1F*)lMother.At(ih))->SetMarkerColor(kCyan-2);
       ((TH1F*)lMother.At(ih))->SetMarkerStyle(1);
       ((TH1F*)lMother.At(ih))->SetLineWidth(2);
       Int_t mcounts = ((TH1F*)lMother.At(ih))->Integral();

       // scale mothers counts for TRD coverage fraction if enabled
       if (considerTRD && (listName.Contains("Trd"))) {
	 if (listName.Contains("NoTrd")) {
	   mcounts=mcounts*19./36.;
	 } else {
	   mcounts=mcounts*13./36.;
	 }
       }
       //set mothers counts
       hMotherCounts->SetBinContent(ih+1, mcounts);
       hMotherCounts->SetBinError(ih+1, TMath::Sqrt(mcounts));

       if (!isOnlyKstar){
	 ((TH1F*)lTruesMP.At(ih))->SetLineColor(kOrange-2);
	 ((TH1F*)lTruesMP.At(ih))->SetMarkerColor(kOrange-2);
	 ((TH1F*)lTruesMP.At(ih))->SetMarkerStyle(21);
	 ((TH1F*)lTruesMP.At(ih))->SetMarkerSize(0.5); 
	 ((TH1F*)lTruesMP.At(ih))->SetLineWidth(1);
	 ((TH1F*)lTruesMP.At(ih))->SetFillColor(kYellow+10);
	 ((TH1F*)lTruesMP.At(ih))->SetFillStyle(3002);
	 //((TH1F*)lTruesMP.At(ih))->SetDrawOption("bar");
	 Int_t acounts = ((TH1F*)lTruesMP.At(ih))->Integral();
	 hTrueCountsMP->SetBinContent(ih+1, acounts);
	 hTrueCountsMP->SetBinError(ih+1, TMath::Sqrt(acounts));

	 ((TH1F*)lAntiMother.At(ih))->SetLineColor(kBlue+2);
	 ((TH1F*)lAntiMother.At(ih))->SetMarkerColor(kBlue+2);
	 ((TH1F*)lAntiMother.At(ih))->SetMarkerStyle(1);
	 ((TH1F*)lAntiMother.At(ih))->SetLineWidth(2);
	 Int_t amcounts = ((TH1F*)lAntiMother.At(ih))->Integral();
	 
	 // scale anti-mothers counts for TRD coverage fraction if enabled
	 if (considerTRD && (listName.Contains("Trd"))) {
	   if (listName.Contains("NoTrd")) {
	     amcounts=amcounts*19./36.;
	   } else {
	     amcounts=amcounts*13./36.;
	   }
	 }
	 //set antimothers counts
	 hAntiMotherCounts->SetBinContent(ih+1, amcounts);
	 hAntiMotherCounts->SetBinError(ih+1, TMath::Sqrt(amcounts));
       }       

       c[icentbin]->cd(ih+1);
       gPad->SetLogy();
       ((TH1F*)lMother.At(ih))->GetYaxis()->SetRangeUser(0.1,5.e5);
       ((TH1F*)lMother.At(ih))->Draw();
       ((TH1F*)lTruesPM.At(ih))->Draw("same");
       if (!isOnlyKstar){
	 ((TH1F*)lAntiMother.At(ih))->Draw("same");
	 ((TH1F*)lTruesMP.At(ih))->Draw("same");
       }       
     }
     
     //efficiency estimation
     TH1D* halltrue = (TH1D*)hTrueCountsPM->Clone();  
     if (!isOnlyKstar) halltrue->Add(hTrueCountsMP,1);
     
     TH1D* hallmother = (TH1D*)hMotherCounts->Clone();  
     if (!isOnlyKstar) hallmother->Add(hAntiMotherCounts,1);
     
     for (Int_t ii = 1; ii<npt+1; ii++){
       Double_t kstrue = hTrueCountsPM->GetBinContent(ii);
       Double_t kstrueerr = hTrueCountsPM->GetBinError(ii);
       Double_t ksmum = hMotherCounts->GetBinContent(ii);
       Double_t ksmumerr = hMotherCounts->GetBinError(ii);
       
       Double_t akstrue=0.0, akstrueerr=0.0, aksmum=0.0, aksmumerr=0.0;
       if (!isOnlyKstar){
	 akstrue = hTrueCountsMP->GetBinContent(ii);
	 akstrueerr = hTrueCountsMP->GetBinError(ii);
         aksmum = hAntiMotherCounts->GetBinContent(ii);
         aksmumerr = hAntiMotherCounts->GetBinError(ii);
       }
       
       Double_t alltrue = halltrue->GetBinContent(ii);
       Double_t errtrue = halltrue->GetBinError(ii);
       Double_t mothers = hMotherCounts->GetBinContent(ii);
       if (!isOnlyKstar) mothers+=hAntiMotherCounts->GetBinContent(ii);       
       
       Double_t errMothers = hMotherCounts->GetBinError(ii);
       Double_t errAntiMothers = 0.0;
       if (!isOnlyKstar) hAntiMotherCounts->GetBinError(ii);
       Double_t errmum = TMath::Sqrt(errAntiMothers*errAntiMothers+errMothers*errMothers);
       Double_t ratio, ratioks,ratioaks;
       if (mothers>0) 
	 ratio = alltrue/mothers;
       else ratio = 0.0;
       Double_t err = ratio*TMath::Sqrt((errtrue/alltrue)*(errtrue/alltrue)+(errmum/mothers)*(errmum/mothers));
       hEffVsPt->SetBinContent(ii, ratio);
       hEffVsPt->SetBinError(ii, err);

       if (ksmum) ratioks = kstrue/ksmum;
       else ratioks = 0.0;
       Double_t errks = ratioks*TMath::Sqrt((kstrueerr/kstrue)*(kstrueerr/kstrue)+(ksmumerr/ksmum)*(ksmumerr/ksmum));
       if (aksmum) ratioaks = akstrue/aksmum;
       else ratioaks = 0.0;
      Double_t erraks = ratioaks*TMath::Sqrt((akstrueerr/akstrue)*(akstrueerr/akstrue)+(aksmumerr/aksmum)*(aksmumerr/aksmum));
       
       hEffVsPtKstar->SetBinContent(ii, ratioks);
       hEffVsPtKstar->SetBinError(ii, errks);
       if (!isOnlyKstar){
	 hEffVsPtAntiKstar->SetBinContent(ii, ratioaks);
	 hEffVsPtAntiKstar->SetBinError(ii, erraks);
       }
       Printf("Efficiency  %3.1f < pt < %3.1f = %6.5f +/- %6.5f", pt[ii-1], pt[ii], ratio, err);
     }

     //ratio for anti-K* over K* 
     if (!isOnlyKstar){
       ratioAKstarOverKstar = (TH1F*) hEffVsPtAntiKstar->Clone("hRatioAKstarOverKstar");
       ratioAKstarOverKstar->Divide(hEffVsPtKstar);
     }
     
     //Print efficiency plot on canvas
     TCanvas * ceffo = new TCanvas("ceffo","efficiency", 800,600);
     ceffo->cd();
     hEffVsPt->Draw();
     ceffo->Print(Form("%s.png",efffilename.Data()));
     ceffo->Delete();

     // save into a file
     efffile->cd();
     ptbins->Write("ptbins");
     centbins->Write("centbins");
     hMotherCounts->Write();
     hTrueCountsPM->Write();
     hEffVsPtKstar->Write();
     if (!isOnlyKstar) {
       hAntiMotherCounts->Write();
       hTrueCountsMP->Write();
       hEffVsPtAntiKstar->Write();

       //makeup for the ratios
       ratioAKstarOverKstar->SetLineColor(histoColor);
       ratioAKstarOverKstar->SetMarkerColor(histoColor);
       ratioAKstarOverKstar->SetMarkerStyle(0);
       ratioAKstarOverKstar->Write();
     }
     halltrue->SetName("AllTrue");
     halltrue->Write();
     hallmother->SetName("AllMother");
     hallmother->Write();
     hEffVsPt->Write();
     efffile->Close();


     //output file
     TFile *fout = TFile::Open(Form("%s_centBin%02d.root",outName,icentbin), "RECREATE");
     fout->cd();   
     ptbins->Write("ptbins");
     centbins->Write("centbins");
     lTruesPM.Write(hInput[0]->GetName(), TObject::kSingleKey);
     lMother.Write(hInput[2]->GetName(), TObject::kSingleKey);
     if (!isOnlyKstar) {
       lTruesMP.Write(hInput[1]->GetName(), TObject::kSingleKey);
       lAntiMother.Write(hInput[3]->GetName(), TObject::kSingleKey);
     }
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
  if (listName.Contains("AntiMother")) type = "AntiMother";
  
  Int_t nhisto = listIn->GetEntries();
  for (Int_t ih=0;ih<nhisto;ih++){
    TH1F * dummy = listIn->At(ih);
    dummy->GetXaxis()->SetRangeUser(0.6, 1.2);
    TString hName = dummy->GetName();
    if (hName.Contains(type.Data())){
      if (hName.Contains(Form("centBin%02d",icentBin))){
	listOut->AddLast(dummy);
	Printf("Histo found with name: %s",hName.Data());
      }
    }
  } 
  return;
}
