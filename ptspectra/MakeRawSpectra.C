/*
fbellini@cern.ch, 25 october 2012
*/
enum EQuantity{ kYields,
		kMass,
		kWidth};

enum EFitFunction{ kPOLY2,
		   kPOLY3,
		   kLandau,
		   kPOLY1,
		   kEXP,
		   kData};

Color_t color[3][6]={kOrange+7, kPink+6, kGreen+1, kAzure+1, kBlue+4, kBlack, //combined
		     kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack,//tof
		     kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack };  //tpc
		    
Int_t marker[3][6]={21, 22, 32, 28, 24, 20, //combined
		    25, 26, 32, 28, 27, 24,//tof
		    21, 22, 23, 34, 33, 20}; //tpc
//color combined/best spectra
//Color_t color[]={kRed+1, kOrange-2, kSpring+9, kCyan-3, kAzure-6, , kBlack};

Color_t colorFunc[]={kBlue, kRed, kGreen+2, kMagenta+1, kYellow+2, kBlack};
Char_t funcName[5][10]={"BW+POLY2","BW+POLY3", "BW+LAND","BW+POLY1", "BW+EXP"};
//Char_t centLabel[5][6]={"0-20%","20-40%","40-60%","60-80%","80-100%"};

void MakeRawSpectra(Bool_t save =1,
                    TString filein="best_fit_poly2.root",
                    TString fileproj="../_sum_EMnorm1.30-1.50_cut1818_tof2s_train215-216.root",
                    TString suffix="_BWpoly2",
                    TString histoTitle = "fit: BW(#Gamma=#Gamma_{PDG})+poly2",
                    Float_t ptmin=0.0, Float_t ptmax=15.0, Float_t cutChi2=10.0, Bool_t skipLastBin=0)
{
  //graphics
  gROOT->LoadMacro("${ASD}/SetGraphicStyle.C");
  SetGraphicStyle(0);
  gStyle->SetTextFont(42);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  Int_t legendEntryCounter = 1;

  //get bins
  TFile * f= TFile::Open(fileproj.Data());
  if (!f) return;
  //get bins
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  Int_t npt = ptbins->GetNbins();
  const Int_t dimpt = npt+1;
  Double_t pt[dimpt];
  for (Int_t k=0; k<dimpt;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
    //Printf("%5.2f",pt[k]);
  }
  TAxis *centbins = (TAxis*)f->Get("centbins");
  Int_t ncent = centbins->GetNbins();
  const Int_t dimcent = ncent+1;
  Double_t cent[dimcent]; 
  for (Int_t k=0; k<dimcent;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
    Printf("%5.2f",cent[k]);
  }
  
  // if tof analysis show only points in 1.0-8.0 GeV/c range, use 0.5-8.0 for TPC
  Bool_t isTOF = fileproj.Contains("tof");
  //if (isTOF) ptmin=1.0; 

  TFile * fin=TFile::Open(filein.Data());
  Double_t fitParams[11], SoverB=0.0, significance=0.0,normfactorCopy=0.0;
  Int_t centBinID,ptBinID, funcID;
  Double_t rangeInfCopy, rangeSupCopy, ptinfCopy, ptsupCopy, centinfCopy, centsupCopy,  fitrange[2];
  TString*fitfunction;
  TTree *tree= (TTree*) fin->Get("tree");//get fit parameters tree
  if (!tree) {
    Printf("ERROR: invalid tree");
    return;
  }
  tree->SetBranchAddress("signalMass",&fitParams[0]);
  tree->SetBranchAddress("signalMassErr",&fitParams[1]);
  tree->SetBranchAddress("signalWidth",&fitParams[2]);
  tree->SetBranchAddress("signalWidthErr",&fitParams[3]);
  tree->SetBranchAddress("nSignal",&fitParams[4]);
  tree->SetBranchAddress("nSignalErr",&fitParams[5]);
  tree->SetBranchAddress("nBack",&fitParams[6]);
  tree->SetBranchAddress("nBackErr",&fitParams[7]);
  tree->SetBranchAddress("chi2",&fitParams[8]);  
  tree->SetBranchAddress("centBin",&centBinID);
  tree->SetBranchAddress("ptBin",&ptBinID);
  tree->SetBranchAddress("SoverB",&SoverB);
  tree->SetBranchAddress("significance",&significance);
  tree->SetBranchAddress("norm_factor", &normfactorCopy);
  tree->SetBranchAddress("pt_inf", &ptinfCopy);
  tree->SetBranchAddress("pt_sup", &ptsupCopy);
  tree->SetBranchAddress("fitrange_inf",&fitrange[0]);
  tree->SetBranchAddress("fitrange_sup",&fitrange[1]);
 
  TH1F*frameMassVsPt=new TH1F("MassVsPt","Mass vs. p_{T}; p_{T} (GeV/c); Mass (GeV/c^{2})", npt, pt);
  TH1F*frameWidthVsPt=new TH1F("WidthVsPt","Width vs. p_{T}; p_{T} (GeV/c); Width (GeV/c^{2})", npt, pt);
  TH1F*frameYieldVsPt=new TH1F("YieldVsPt","Raw yield vs. p_{T}; p_{T} (GeV/c); dN/dp_{T}", npt, pt);
  TH1F*frameChi2VsPt=new TH1F("Chi2VsPt","#chi^{2} vs. p_{T}; p_{T} (GeV/c); #chi^{2}", npt, pt);
  TH1F*frameSoverBVsPt=new TH1F("SoverBVsPt","S/B; p_{T} (GeV/c); S/B", npt, pt);
  TH1F*frameSignificanceVsPt=new TH1F("SignificanceVsPt","S/#sqrt{S+B}; p_{T} (GeV/c); significance", npt, pt);

  TLine * pdgmass = new TLine(0.,0.896,10.,0.896);
  pdgmass->SetLineStyle(2);
  pdgmass->SetLineWidth(3);
  pdgmass->SetLineColor(kBlack);
  TLine * pdgwidth = new TLine(0.,0.0505,10.,0.0505);
  pdgwidth->SetLineStyle(2);
  pdgwidth->SetLineWidth(3);
  pdgwidth->SetLineColor(kBlack);

  TLegend *mleg=new TLegend(0.8,0.65,0.92,0.88);
  mleg->SetHeader(" Centrality:");
  mleg->SetFillColor(kWhite);
  
  TCanvas *cfit=new TCanvas("MassVsPt","Mass vs p_{T}", 700,600);
  cfit->cd();
  frameMassVsPt->SetLineColor(kWhite);
  frameMassVsPt->SetMarkerColor(kWhite);
  frameMassVsPt->GetYaxis()->SetRangeUser(0.86, 0.92);
  frameMassVsPt->GetYaxis()->SetTitleOffset(1.5);
  //  frameMassVsPt->Draw();

  TCanvas *cwi=new TCanvas("WidthVsPt","Width vs p_{T}", 700,600);
  cwi->cd();
  frameWidthVsPt->SetLineColor(kWhite);
  frameWidthVsPt->SetMarkerColor(kWhite);
  frameWidthVsPt->GetYaxis()->SetRangeUser(0.0, 0.10);
  frameWidthVsPt->GetYaxis()->SetTitleOffset(1.5);
  //  frameWidthVsPt->Draw();

  TCanvas *cry=new TCanvas("cry","Raw yield vs p_{T}", 600,700);
  cry->cd();
  frameYieldVsPt->SetLineColor(kWhite);
  frameYieldVsPt->SetMarkerColor(kWhite);
  frameYieldVsPt->GetYaxis()->SetRangeUser(1e2, 1e7);
  frameYieldVsPt->GetYaxis()->SetTitleOffset(1.1);
  // frameYieldVsPt->Draw();

  TCanvas *cc2=new TCanvas("cc2","#Chi^{2} vs p_{T}", 700,600);
  cc2->cd();
  frameChi2VsPt->SetLineColor(kWhite);
  frameChi2VsPt->SetMarkerColor(kWhite);
  frameChi2VsPt->GetYaxis()->SetRangeUser(0, 5);
  frameChi2VsPt->GetYaxis()->SetTitleOffset(1.1);
  // frameChi2VsPt->Draw();

  TCanvas *cSoB=new TCanvas("cSoB","S/B vs p_{T}", 700,600);
  //cfit->Divide(3,2);
  cSoB->cd();
  frameSoverBVsPt->SetLineColor(kWhite);
  frameSoverBVsPt->SetMarkerColor(kWhite);
  frameSoverBVsPt->GetYaxis()->SetRangeUser(0.001, 1.);
  //->GetYaxis()->SetTitleOffset(1.1);
  //  frameSoverBVsPt->Draw();
  
  TCanvas *csign=new TCanvas("csign","Significance vs p_{T}", 700,600);
  csign->cd();
  frameSignificanceVsPt->SetLineColor(kWhite);
  frameSignificanceVsPt->SetMarkerColor(kWhite);
  frameSignificanceVsPt->GetYaxis()->SetRangeUser(0.1, 150.);
  frameSignificanceVsPt->GetYaxis()->SetTitleOffset(1.1);
  // frameSignificanceVsPt->Draw();

  TH1F *hMassVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hWidthVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hRawYieldVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hChi2VsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hSoverBVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hSignificanceVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TString foutname = Form("RAW_%s", filein.Data());
  TFile * fout = new TFile(foutname.Data(),"recreate");
  Printf("NCENT = %i", ncent);

  /*loop on centrality bins*/
  for (Int_t icentbin=0;icentbin<ncent;icentbin++){
    TString centLabel=Form("%3.0f-%3.0f %",cent[icentbin],cent[icentbin+1]);
    Printf("*******************\n %s \n*******************", centLabel.Data());
    /* mass vs pT plots */
    hMassVsPt[icentbin] = new TH1F(Form("hMassVsPt_%i",icentbin),"Mass vs p_{T}; p_{T} (GeV/c); Mass (GeV/c^{2})", npt, pt);
    hMassVsPt[icentbin]->SetTitle(Form("Fitted K* mass Vs p_{T} (%s)",centLabel.Data()));

    /* mass vs pT plots */
    hWidthVsPt[icentbin] = new TH1F(Form("hWidthVsPt_%i",icentbin),"Width vs p_{T}; p_{T} (GeV/c); Width (GeV/c^{2})", npt, pt);
    hWidthVsPt[icentbin]->SetTitle(Form("Fitted K* width Vs p_{T} (%s)",centLabel.Data()));
    
    /* yield vs pT plots */
    hRawYieldVsPt[icentbin] = new TH1F(Form("hRawYieldVsPt_%i",icentbin),"Raw yield vs p_{T}; p_{T} (GeV/c); dN/dp_{T})", npt, pt);
    hRawYieldVsPt[icentbin]->SetTitle(Form("Raw yields vs p_{T} (%s)",centLabel.Data()));

   /* chi2 vs pT plots */
    hChi2VsPt[icentbin] = new TH1F(Form("hChi2VsPt_%i",icentbin),"#chi^{2} vs p_{T}; p_{T} (GeV/c); dN/dp_{T})", npt, pt);
    hChi2VsPt[icentbin]->SetTitle(Form("#chi^{2} vs p_{T} (%s)",centLabel.Data()));

 /* S/B vs pT plots */
    hSoverBVsPt[icentbin] = new TH1F(Form("hSoverBVsPt_%i",icentbin),"S/B vs p_{T}; p_{T} (GeV/c); S/B)", npt, pt);
    hSoverBVsPt[icentbin]->SetTitle(Form("S/B vs p_{T} (%s)",centLabel.Data()));

   /* signif vs pT plots */
    hSignificanceVsPt[icentbin] = new TH1F(Form("hSignificanceVsPt_%i",icentbin),"Significance vs p_{T}; p_{T} (GeV/c); Significance)", npt, pt);
    hSignificanceVsPt[icentbin]->SetTitle(Form("Significance vs p_{T} (%s)",centLabel.Data()));

    //stop at pt= 6 GeV7c for 80-100%
    if (icentbin==4) skipLastBin = kTRUE;
    if (skipLastBin) ptmax = 6.0;

    for (Int_t ientry=0;ientry<tree->GetEntries();ientry++){
      tree->GetEntry(ientry);
      if ((ptinfCopy<ptmin) || (ptsupCopy>ptmax)) continue;
      if (centBinID==icentbin || centBinID==100){
      	//apply chi2 cut
       	if (fitParams[8]<cutChi2){
	  //check range to be shown
	  /*fill mass vs pt*/
	  //	  Printf("ptBinID = %i  -  pt[x] = %4.2f", ptBinID, pt[ptBinID]);
	  //tree->Show(ientry);
	  hMassVsPt[icentbin]->SetBinContent(ptBinID+1,fitParams[0]);
	  hMassVsPt[icentbin]->SetBinError(ptBinID+1,fitParams[1]);
	  // fin->cd();
	  // hMassVsPt->Write();
	  hMassVsPt[icentbin]->SetMarkerStyle(marker[isTOF][icentbin]);
	  hMassVsPt[icentbin]->SetMarkerColor(color[isTOF][icentbin]);
	  hMassVsPt[icentbin]->SetLineColor(color[isTOF][icentbin]);    
	  hMassVsPt[icentbin]->SetLineWidth(1);
	  hMassVsPt[icentbin]->GetYaxis()->SetRangeUser(0.86, 0.92);
	  hMassVsPt[icentbin]->GetYaxis()->SetTitleOffset(1.5);
	  hMassVsPt[icentbin]->GetYaxis()->SetNdivisions(515);
	  hMassVsPt[icentbin]->GetYaxis()->SetDecimals();
	  hMassVsPt[icentbin]->SetTitle(Form("%s, %s",histoTitle.Data(),centLabel.Data()));
	  /*fill width vs pt*/
	  hWidthVsPt[icentbin]->SetBinContent(ptBinID+1,fitParams[2]);
	  hWidthVsPt[icentbin]->SetBinError(ptBinID+1,fitParams[3]);
	  hWidthVsPt[icentbin]->SetMarkerStyle(marker[isTOF][icentbin]);
	  hWidthVsPt[icentbin]->SetMarkerColor(color[isTOF][icentbin]);
	  hWidthVsPt[icentbin]->SetLineColor(color[isTOF][icentbin]);    
	  hWidthVsPt[icentbin]->SetLineWidth(1);   
	  hWidthVsPt[icentbin]->GetYaxis()->SetRangeUser(0.0, 0.10);
	  hWidthVsPt[icentbin]->SetTitle(Form("%s, %s",histoTitle.Data(),centLabel.Data()));

	  /*fill yield vs pt */
	  //Printf("nsignal = %f    error = %f",fitParams[4],fitParams[5]);
	  Float_t dpt = ptbins->GetBinWidth(ptBinID + 1);
	  Printf("ptBIN = %i --> dpt = %5.2f ", ptBinID, dpt);
	  hRawYieldVsPt[icentbin]->SetBinContent(ptBinID+1,fitParams[4]);
	  hRawYieldVsPt[icentbin]->SetBinError(ptBinID+1,fitParams[5]);//how to estimate error on the integral????
	  hRawYieldVsPt[icentbin]->SetMarkerStyle(marker[isTOF][icentbin]);
	  hRawYieldVsPt[icentbin]->SetMarkerColor(color[isTOF][icentbin]);
	  hRawYieldVsPt[icentbin]->SetLineColor(color[isTOF][icentbin]);  
	  hRawYieldVsPt[icentbin]->SetLineWidth(1);
	  hRawYieldVsPt[icentbin]->GetYaxis()->SetRangeUser(1e1, 1e7);
	  hRawYieldVsPt[icentbin]->GetYaxis()->SetTitleOffset(1.5);
	  hRawYieldVsPt[icentbin]->GetYaxis()->SetTitle("dN/dp_{T}");
	  hRawYieldVsPt[icentbin]->SetTitle(Form("%s, %s",histoTitle.Data(),centLabel.Data()));
	  /*fill chi2 vs pt */
	  //Printf("chi2 = %f",fitParams[8]);
	  hChi2VsPt[icentbin]->SetBinContent(ptBinID+1,fitParams[8]);
	  hChi2VsPt[icentbin]->SetMarkerStyle(marker[isTOF][icentbin]);
	  hChi2VsPt[icentbin]->SetMarkerColor(color[isTOF][icentbin]);
	  hChi2VsPt[icentbin]->SetLineColor(color[isTOF][icentbin]);       
	  hChi2VsPt[icentbin]->SetLineWidth(1);
	  hChi2VsPt[icentbin]->GetYaxis()->SetRangeUser(0., 5.);
	  hChi2VsPt[icentbin]->SetTitle(Form("%s, %s",histoTitle.Data(),centLabel.Data()));
	  /*fill S/B vs pt */
	  hSoverBVsPt[icentbin]->SetBinContent(ptBinID+1,SoverB);
	  hSoverBVsPt[icentbin]->SetMarkerStyle(marker[isTOF][icentbin]);
	  hSoverBVsPt[icentbin]->SetMarkerColor(color[isTOF][icentbin]);
	  hSoverBVsPt[icentbin]->SetLineColor(color[isTOF][icentbin]);     
	  hSoverBVsPt[icentbin]->SetLineWidth(1);
	  hSoverBVsPt[icentbin]->GetYaxis()->SetRangeUser(0., 1.);
	  hSoverBVsPt[icentbin]->SetTitle(Form("%s, %s",histoTitle.Data(),centLabel.Data()));
	  /*fillsignif vs pt */
	  //Printf("significance = %f",significance);
	  hSignificanceVsPt[icentbin]->SetBinContent(ptBinID+1,significance);
	  hSignificanceVsPt[icentbin]->SetMarkerStyle(marker[isTOF][icentbin]);
	  hSignificanceVsPt[icentbin]->SetMarkerColor(color[isTOF][icentbin]);
	  hSignificanceVsPt[icentbin]->SetLineColor(color[isTOF][icentbin]);     
	  hSignificanceVsPt[icentbin]->SetLineWidth(1);
	  hSignificanceVsPt[icentbin]->GetYaxis()->SetRangeUser(0., 150.);
	  hSignificanceVsPt[icentbin]->SetTitle(Form("%s, %s",histoTitle.Data(),centLabel.Data()));
	}//end cut on chi2
      }
    } // loop over tree entries
    mleg->AddEntry(hMassVsPt[icentbin],centLabel.Data(),"lpf");
    if (hMassVsPt[icentbin] && hMassVsPt[icentbin]->GetEntries()>0) {
      cfit->cd();
      if (icentbin>0)
	hMassVsPt[icentbin]->Draw("same");
      else 
	hMassVsPt[icentbin]->Draw();
      fout->cd();
      if (save) hMassVsPt[icentbin]->Write();
    }
    if (hWidthVsPt[icentbin] && hWidthVsPt[icentbin]->GetEntries()>0) {
      cwi->cd();
      if (icentbin>0)
	hWidthVsPt[icentbin]->Draw("same");
      else
	hWidthVsPt[icentbin]->Draw();
      fout->cd();
      if (save) hWidthVsPt[icentbin]->Write();
    }
    cry->cd();
    gPad->SetLogy();
    if (hRawYieldVsPt[icentbin] && hRawYieldVsPt[icentbin]->GetEntries()>0){
      if (icentbin>0) 
	hRawYieldVsPt[icentbin]->Draw("same");
      else 
	hRawYieldVsPt[icentbin]->Draw();
      fout->cd();
      if (save) hRawYieldVsPt[icentbin]->Write();
      legendEntryCounter+=1;
    }
    // cfit->Update();
    // cry->Update();
    if (hSignificanceVsPt[icentbin] && hSignificanceVsPt[icentbin]->GetEntries()>0){
      cc2->cd();
      hChi2VsPt[icentbin]->Draw("same");
      fout->cd();
      if (save) hChi2VsPt[icentbin]->Write();
    }
    if (hSignificanceVsPt[icentbin] && hSignificanceVsPt[icentbin]->GetEntries()>0){
      csign->cd();
      hSignificanceVsPt[icentbin]->Draw("same");
      fout->cd();
      if (save) hSignificanceVsPt[icentbin]->Write();
    }
    if (hSoverBVsPt[icentbin] && hSoverBVsPt[icentbin]->GetEntries()>0){
      cSoB->cd();
      hSoverBVsPt[icentbin]->Draw("same");
      fout->cd();
      if (save) hSoverBVsPt[icentbin]->Write();
    }
  }
  return;
  gROOT->LoadMacro("$ASD/AddPaveText.C");
  cry->cd();
  AddPaveTextStatOnly();
  //  gPad->SetGridy();
  TLegend * autolegry = (TLegend*) cry->BuildLegend(0.4,(0.88-legendEntryCounter*0.03),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegry->SetFillColor(kWhite);
  autolegry->SetLineColor(kWhite);
  //mleg->Draw();
  
  cfit->cd();
  AddPaveTextStatOnly();
  TLegend * autolegm = (TLegend*) cfit->BuildLegend(0.5, 0.15, 0.88, 0.15+(legendEntryCounter*0.04)/*,frameMassVsPt->GetTitle()*/);
  autolegm->SetFillColor(kWhite);
  autolegm->SetLineColor(kWhite);
  //  mleg->Draw();
  //gPad->SetGridy();
  pdgmass->Draw("same");
  autolegm->AddEntry(pdgmass, "PDG value", "l");
  
  cwi->cd();
  AddPaveTextStatOnly();
  TLegend * autolegwi = (TLegend*)cwi->BuildLegend(0.5,(0.88-legendEntryCounter*0.04),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegwi->SetFillColor(kWhite);
  autolegwi->SetLineColor(kWhite);
  //mleg->Draw();
  pdgwidth->Draw("same");
  autolegwi->AddEntry(pdgwidth, "PDG value", "l");

  cc2->cd();
  TLegend * autolegc2 = (TLegend*) cc2->BuildLegend(0.5,(0.88-legendEntryCounter*0.04),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegc2->SetFillColor(kWhite);
  autolegc2->SetLineColor(kWhite);
  gPad->SetGridy();
  //mleg->Draw();
  
  csign->cd();
  TLegend * autolegsign = (TLegend*)csign->BuildLegend(0.5,(0.88-legendEntryCounter*0.04),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegsign->SetFillColor(kWhite);
  autolegsign->SetLineColor(kWhite);
  gPad->SetGridy();
  // mleg->Draw();
  
  cSoB->cd();
  TLegend * autolegSoB = (TLegend*) cSoB->BuildLegend(0.5,(0.88-legendEntryCounter*0.04),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegSoB->SetFillColor(kWhite);
  autolegSoB->SetLineColor(kWhite);
  gPad->SetGridy();
  // gPad->SetLogy();
  // frameSoverBVsPt->GetYaxis()->Unzoom();
  //  mleg->Draw();
  
  TString fimgname(filein.Data());
  fimgname.ReplaceAll(".root",".png");
  cry->SaveAs(Form("RAW_%s",fimgname.Data()));
  cfit->SaveAs(Form("MASS_%s",fimgname.Data()));
  cwi->SaveAs(Form("WIDTH_%s",fimgname.Data()));
  cc2->SaveAs(Form("CHI2_%s",fimgname.Data()));
  cSoB->SaveAs(Form("SoB_%s",fimgname.Data()));
  csign->SaveAs(Form("Signif_%s",fimgname.Data()));
  return;
}
