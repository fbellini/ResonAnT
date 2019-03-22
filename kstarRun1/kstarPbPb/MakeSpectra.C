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

Color_t color[]={kRed+1, kOrange+1, kGreen+2, kBlue+2, kMagenta, kBlack};
Color_t colorFunc[]={kBlue, kRed, kGreen+2, kMagenta+1, kYellow+2, kBlack};
Int_t marker[]={20, 21, 28, 26, 23};//1,1,1,1,1 };
Char_t funcName[5][10]={"BW+POLY2","BW+POLY3", "BW+LAND","BW+POLY1", "BW+EXP"};
//Char_t centLabel[5][6]={"0-20%","20-40%","40-60%","60-80%","80-90%"};

void MakeSpectra10jun13(){
  MakeSpectra(1,"roofit_freeW/","tree_kstar_EXE.root", "_kstar_EMnorm1.30-1.50_train81.root", "_EXE_freeW", "fit: BWxPS+ExE");
  MakeSpectra(1,"roofit_freeW/", "tree_kstar_BWPS3.root", "_kstar_EMnorm1.30-1.50_train81.root", "_BWPS3_freeW", "fit: BWxPS+cheby3");
  MakeSpectra(1, "roofit_freeW/", "tree_antikstar_EXE.root", "_antikstar_EMnorm1.30-1.50_train81.root","_EXE_freeW", "fit: BWxPS+ExE");
  MakeSpectra(1, "roofit_freeW/", "tree_antikstar_BWPS3.root", "_antikstar_EMnorm1.30-1.50_train81.root", "_BWPS3_freeW", "fit: BWxPS+cheby3");
  return;
}

void MakeSpectra31may13_improved(){
  MakeSpectra(1,"roofit_improved/","tree_kstar_EXE.root", "_kstar_EMnorm1.30-1.50_train81.root", "_EXE_improved", "(improved) fit: BWxPS+ExE");
  MakeSpectra(1,"roofit_improved/", "tree_kstar_BWPS3.root", "_kstar_EMnorm1.30-1.50_train81.root", "_BWPS3_improved", "(improved) fit: BWxPS+cheby3");
  MakeSpectra(1, "roofit_improved/", "tree_antikstar_EXE.root", "_kstar_EMnorm1.30-1.50_train81.root","_EXE_improved", "(improved) fit: BWxPS+ExE");
  MakeSpectra(1, "roofit_improved/", "tree_antikstar_BWPS3.root", "_kstar_EMnorm1.30-1.50_train81.root", "_BWPS3_improved", "(improved) fit: BWxPS+cheby3");
  return;
}

void MakeSpectra31may13(){
  MakeSpectra(1,"roofit/","kstar_EXE_31may13.root", "_kstar_EMnorm1.30-1.50_train81.root", "_EXE", "fit: BWxPS+ExE");
  MakeSpectra(1,"roofit/", "kstar_BWPS3_31may13.root", "_kstar_EMnorm1.30-1.50_train81.root", "_BWPS3", "fit: BWxPS+cheby3");
  MakeSpectra(1, "roofit/", "antikstar_EXE_31may13.root", "_kstar_EMnorm1.30-1.50_train81.root","_EXE", "fit: BWxPS+ExE");
  MakeSpectra(1, "roofit/", "antikstar_BWPS3_31may13.root", "_kstar_EMnorm1.30-1.50_train81.root", "_BWPS3", "fit: BWxPS+cheby3");
  return;
}


void MakeSpectra(Bool_t save, Char_t * dirIn="roofit/", TString filein="kstar_31may13.root", Char_t* fileproj="_kstar_EMnorm1.30-1.50_train81.root", TString suffix="_BWPS3", TString histoTitle = "fit: BWxPS+cheby3", Float_t cutChi2=2.0)
{
  
  gROOT->LoadMacro("SetGraphicStyle.C");
  SetGraphicStyle(0);
  gStyle->SetTextFont(42);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  Int_t legendEntryCounter = 1;
  //Char_t projectionFile[100]="sub_kstar.root";
  TFile * f= TFile::Open(fileproj);
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
    //Printf("%5.2f",cent[k]);
  }

  TFile * fin=TFile::Open(Form("%s/%s",dirIn,filein.Data()));

  Double_t fitParams[9], SoverB=0.0, significance=0.0,normfactorCopy=0.0;
  Int_t centBinID,ptBinID, funcID;
  Double_t rangeInfCopy, rangeSupCopy, ptinfCopy, ptsupCopy, centinfCopy, centsupCopy,  fitrange[2];
  TString*fitfunction;
  TTree *tree= (TTree*) fin->Get("tree");//get fit parameters tree
  if (!tree) return;
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
  // tree->SetBranchAddress("norm_inf", &rangeInfCopy);
  // tree->SetBranchAddress("norm_sup", &rangeSupCopy);
  tree->SetBranchAddress("pt_inf", &ptinfCopy);
  tree->SetBranchAddress("pt_sup", &ptsupCopy);
  // tree->SetBranchAddress("cent_inf", &centinfCopy);
  // tree->SetBranchAddress("cent_sup", &centsupCopy);
  //tree->SetBranchAddress("functionID", &funcID);
  tree->SetBranchAddress("fitrange_inf",&fitrange[0]);
  tree->SetBranchAddress("fitrange_sup",&fitrange[1]);
 
  TH1F*frameMassVsPt=new TH1F("MassVsPt","Mass vs. p_{t}; p_{t} (GeV/c); M (GeV/c^{2})", npt, pt);
  TH1F*frameWidthVsPt=new TH1F("WidthVsPt","Width vs. p_{t}; p_{t} (GeV/c); #Gamma (GeV/c^{2})", npt, pt);
  TH1F*frameYieldVsPt=new TH1F("YieldVsPt","Raw yield vs. p_{t}; p_{t} (GeV/c); dN/dp_{t}", npt, pt);
  TH1F*frameChi2VsPt=new TH1F("Chi2VsPt","#chi^{2} vs. p_{t}; p_{t} (GeV/c); #chi^{2}", npt, pt);
  TH1F*frameSoverBVsPt=new TH1F("SoverBVsPt","S/B; p_{t} (GeV/c); S/B", npt, pt);
  TH1F*frameSignificanceVsPt=new TH1F("SignificanceVsPt","S/#sqrt{S+B}; p_{t} (GeV/c); significance", npt, pt);

  TLegend *mleg=new TLegend(0.8,0.65,0.92,0.88);
  mleg->SetHeader(" Centrality:");
  mleg->SetFillColor(kWhite);
  
  TCanvas *cfit=new TCanvas("MassVsPt","Mass vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  cfit->cd();
  frameMassVsPt->SetLineColor(kWhite);
  frameMassVsPt->SetMarkerColor(kWhite);
  frameMassVsPt->GetYaxis()->SetRangeUser(0.86, 0.92);
  frameMassVsPt->GetYaxis()->SetTitleOffset(1.4);
  frameMassVsPt->Draw();

  TCanvas *cwi=new TCanvas("WidthVsPt","Width vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  cwi->cd();
  frameWidthVsPt->SetLineColor(kWhite);
  frameWidthVsPt->SetMarkerColor(kWhite);
  frameWidthVsPt->GetYaxis()->SetRangeUser(0.01, 0.10);
  frameWidthVsPt->GetYaxis()->SetTitleOffset(1.5);
  frameWidthVsPt->Draw();

  TCanvas *cry=new TCanvas("cry","Raw yield vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  cry->cd();
  frameYieldVsPt->SetLineColor(kWhite);
  frameYieldVsPt->SetMarkerColor(kWhite);
  frameYieldVsPt->GetYaxis()->SetRangeUser(1e2, 1e7);
  frameYieldVsPt->GetYaxis()->SetTitleOffset(1.1);
  frameYieldVsPt->Draw();

  TCanvas *cc2=new TCanvas("cc2","#Chi^{2} vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  cc2->cd();
  frameChi2VsPt->SetLineColor(kWhite);
  frameChi2VsPt->SetMarkerColor(kWhite);
  frameChi2VsPt->GetYaxis()->SetRangeUser(0, 10);
  frameChi2VsPt->GetYaxis()->SetTitleOffset(1.1);
  frameChi2VsPt->Draw();

  TCanvas *cSoB=new TCanvas("cSoB","S/B vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  cSoB->cd();
  frameSoverBVsPt->SetLineColor(kWhite);
  frameSoverBVsPt->SetMarkerColor(kWhite);
  frameSoverBVsPt->GetYaxis()->SetRangeUser(0.001, 10.);
  //->GetYaxis()->SetTitleOffset(1.1);
  frameSoverBVsPt->Draw();
  
  TCanvas *csign=new TCanvas("csign","Significance vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  csign->cd();
  frameSignificanceVsPt->SetLineColor(kWhite);
  frameSignificanceVsPt->SetMarkerColor(kWhite);
  frameSignificanceVsPt->GetYaxis()->SetRangeUser(0.1, 400.);
  frameSignificanceVsPt->GetYaxis()->SetTitleOffset(1.1);
  frameSignificanceVsPt->Draw();

  TH1F *hMassVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hWidthVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hRawYieldVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hChi2VsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hSoverBVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hSignificanceVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  if (filein.Contains("fitEM_")) filein.ReplaceAll("fitEM_","");
  if (filein.Contains("aod049_")) filein.ReplaceAll("aod049_","");
  filein.ReplaceAll(".root",Form("%s.root",suffix.Data()));
  TString foutname = Form("rawYields/rawYields_%s", filein.Data());
  TFile * fout = new TFile(foutname.Data(),"recreate");
  
  /*loop on centrality bins*/
  for (Int_t icentbin=0;icentbin<ncent-1;icentbin++){
    TString centLabel=Form("%2.0f-%2.0f%%",cent[icentbin],cent[icentbin+1]);
    
    /* mass vs pT plots */
    hMassVsPt[icentbin] = new TH1F(Form("hMassVsPt_%i",icentbin),"Mass vs p_{t}; p_{t} (GeV/c); M (GeV/c^{2})", npt, pt);
    hMassVsPt[icentbin]->SetTitle(Form("Fitted K* mass Vs p_{t} (%s)",centLabel.Data()));

    /* mass vs pT plots */
    hWidthVsPt[icentbin] = new TH1F(Form("hWidthVsPt_%i",icentbin),"Width vs p_{t}; p_{t} (GeV/c); #Gamma (GeV/c^{2})", npt, pt);
    hWidthVsPt[icentbin]->SetTitle(Form("Fitted K* width Vs p_{t} (%s)",centLabel.Data()));

    /* yield vs pT plots */
    hRawYieldVsPt[icentbin] = new TH1F(Form("hRawYieldVsPt_%i",icentbin),"Raw yield vs p_{t}; p_{t} (GeV/c); dN/dp_{t})", npt, pt);
    hRawYieldVsPt[icentbin]->SetTitle(Form("Raw yields vs p_{t} (%s)",centLabel.Data()));

   /* chi2 vs pT plots */
    hChi2VsPt[icentbin] = new TH1F(Form("hChi2VsPt_%i",icentbin),"#chi^{2} vs p_{t}; p_{t} (GeV/c); dN/dp_{t})", npt, pt);
    hChi2VsPt[icentbin]->SetTitle(Form("#chi^{2} vs p_{t} (%s)",centLabel.Data()));

 /* S/B vs pT plots */
    hSoverBVsPt[icentbin] = new TH1F(Form("hSoverBVsPt_%i",icentbin),"S/B vs p_{t}; p_{t} (GeV/c); S/B)", npt, pt);
    hSoverBVsPt[icentbin]->SetTitle(Form("S/B vs p_{t} (%s)",centLabel.Data()));

   /* signif vs pT plots */
    hSignificanceVsPt[icentbin] = new TH1F(Form("hSignificanceVsPt_%i",icentbin),"Significance vs p_{t}; p_{t} (GeV/c); Significance)", npt, pt);
    hSignificanceVsPt[icentbin]->SetTitle(Form("Significance vs p_{t} (%s)",centLabel.Data()));

    for (Int_t ientry=0;ientry<tree->GetEntries();ientry++){
      tree->GetEntry(ientry);
      if (centBinID==icentbin){
	//apply chi2 cut
	if (fitParams[8]<cutChi2){
	  /*fill mass vs pt*/
	  //	  Printf("ptBinID = %i  -  pt[x] = %4.2f", ptBinID, pt[ptBinID]);
	  hMassVsPt[icentbin]->SetBinContent(ptBinID+1,fitParams[0]);
	  hMassVsPt[icentbin]->SetBinError(ptBinID+1,fitParams[1]);
	  // fin->cd();
	  // hMassVsPt->Write();
	  hMassVsPt[icentbin]->SetMarkerStyle(marker[icentbin]);
	  hMassVsPt[icentbin]->SetMarkerColor(color[icentbin]);
	  hMassVsPt[icentbin]->SetLineColor(color[icentbin]);    
	  hMassVsPt[icentbin]->SetLineWidth(2);
	  hMassVsPt[icentbin]->GetYaxis()->SetRangeUser(0.86, 0.92);
	  hMassVsPt[icentbin]->SetTitle(Form("%s, %s",histoTitle.Data(),centLabel.Data()));

	  /*fill width vs pt*/
	  hWidthVsPt[icentbin]->SetBinContent(ptBinID+1,fitParams[2]);
	  hWidthVsPt[icentbin]->SetBinError(ptBinID+1,fitParams[3]);
	  hWidthVsPt[icentbin]->SetMarkerStyle(marker[icentbin]);
	  hWidthVsPt[icentbin]->SetMarkerColor(color[icentbin]);
	  hWidthVsPt[icentbin]->SetLineColor(color[icentbin]);    
	  hWidthVsPt[icentbin]->SetLineWidth(2);   
	  hWidthVsPt[icentbin]->GetYaxis()->SetRangeUser(0.86, 0.94);
	  hWidthVsPt[icentbin]->SetTitle(Form("%s, %s",histoTitle.Data(),centLabel.Data()));

	  /*fill yield vs pt */
	  //Printf("nsignal = %f    error = %f",fitParams[4],fitParams[5]);
	  Float_t dpt = ptbins->GetBinWidth(ptBinID + 1);
	  //Printf("ptBIN = %i --> dpt = %5.2f ", ptBinID, dpt);
	  hRawYieldVsPt[icentbin]->SetBinContent(ptBinID+1,fitParams[4]/dpt);
	  hRawYieldVsPt[icentbin]->SetBinError(ptBinID+1,fitParams[5]/dpt);//how to estimate error on the integral????
	  hRawYieldVsPt[icentbin]->SetMarkerStyle(marker[icentbin]);
	  hRawYieldVsPt[icentbin]->SetMarkerColor(color[icentbin]);
	  hRawYieldVsPt[icentbin]->SetLineColor(color[icentbin]);  
	  hRawYieldVsPt[icentbin]->SetLineWidth(2);
	  hRawYieldVsPt[icentbin]->GetYaxis()->SetRangeUser(1, 1e7);
	  hRawYieldVsPt[icentbin]->SetTitle(Form("%s, %s",histoTitle.Data(),centLabel.Data()));
	  /*fill chi2 vs pt */
	  //Printf("chi2 = %f",fitParams[8]);
	  hChi2VsPt[icentbin]->SetBinContent(ptBinID+1,fitParams[8]);
	  hChi2VsPt[icentbin]->SetMarkerStyle(marker[icentbin]);
	  hChi2VsPt[icentbin]->SetMarkerColor(color[icentbin]);
	  hChi2VsPt[icentbin]->SetLineColor(color[icentbin]);       
	  hChi2VsPt[icentbin]->SetLineWidth(2);
	  hChi2VsPt[icentbin]->GetYaxis()->SetRangeUser(0., 15.);
	  hChi2VsPt[icentbin]->SetTitle(Form("%s, %s",histoTitle.Data(),centLabel.Data()));
	  /*fill S/B vs pt */
	  hSoverBVsPt[icentbin]->SetBinContent(ptBinID+1,SoverB);
	  hSoverBVsPt[icentbin]->SetMarkerStyle(marker[icentbin]);
	  hSoverBVsPt[icentbin]->SetMarkerColor(color[icentbin]);
	  hSoverBVsPt[icentbin]->SetLineColor(color[icentbin]);     
	  hSoverBVsPt[icentbin]->SetLineWidth(2);
	  hSoverBVsPt[icentbin]->GetYaxis()->SetRangeUser(0., 10.);
	  hSoverBVsPt[icentbin]->SetTitle(Form("%s, %s",histoTitle.Data(),centLabel.Data()));
	  /*fillsignif vs pt */
	  //Printf("significance = %f",significance);
	  hSignificanceVsPt[icentbin]->SetBinContent(ptBinID+1,significance);
	  hSignificanceVsPt[icentbin]->SetMarkerStyle(marker[icentbin]);
	  hSignificanceVsPt[icentbin]->SetMarkerColor(color[icentbin]);
	  hSignificanceVsPt[icentbin]->SetLineColor(color[icentbin]);     
	  hSignificanceVsPt[icentbin]->SetLineWidth(2);
	  hSignificanceVsPt[icentbin]->GetYaxis()->SetRangeUser(0., 300.);
	  hSignificanceVsPt[icentbin]->SetTitle(Form("%s, %s",histoTitle.Data(),centLabel.Data()));
	}//end cut on chi2
      }
    }
    mleg->AddEntry(hMassVsPt[icentbin],centLabel.Data(),"lpf");
    if (hMassVsPt[icentbin] && hMassVsPt[icentbin]->GetEntries()>0) {
      cfit->cd();
      hMassVsPt[icentbin]->Draw("same");
      fout->cd();
      if (save) hMassVsPt[icentbin]->Write();
    }
    if (hWidthVsPt[icentbin] && hWidthVsPt[icentbin]->GetEntries()>0) {
      cwi->cd();
      hWidthVsPt[icentbin]->Draw("same");
      fout->cd();
      if (save) hWidthVsPt[icentbin]->Write();
    }
    cry->cd();
    gPad->SetLogy();
    if (hRawYieldVsPt[icentbin] && hRawYieldVsPt[icentbin]->GetEntries()>0){
      if (icentbin<4) 
	hRawYieldVsPt[icentbin]->Draw("same");
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
  cry->cd();
  gPad->SetGridy();
  TLegend * autolegry = (TLegend*) gPad->BuildLegend(0.5,(0.88-legendEntryCounter*0.04),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegry->SetFillColor(kWhite);
  autolegry->SetLineColor(kWhite);
  //mleg->Draw();
  
  cfit->cd();
  TLegend * autolegm = (TLegend*) gPad->BuildLegend(0.5,(0.88-legendEntryCounter*0.04),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegm->SetFillColor(kWhite);
  autolegm->SetLineColor(kWhite);
  //  mleg->Draw();
  gPad->SetGridy();
  
  cwi->cd();
  TLegend * autolegwi = (TLegend*)gPad->BuildLegend(0.5,(0.88-legendEntryCounter*0.04),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegwi->SetFillColor(kWhite);
  autolegwi->SetLineColor(kWhite);
  //mleg->Draw();

  cc2->cd();
  TLegend * autolegc2 = (TLegend*) gPad->BuildLegend(0.5,(0.88-legendEntryCounter*0.04),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegc2->SetFillColor(kWhite);
  autolegc2->SetLineColor(kWhite);
  gPad->SetGridy();
  //mleg->Draw();
  
  csign->cd();
  TLegend * autolegsign = (TLegend*)gPad->BuildLegend(0.5,(0.88-legendEntryCounter*0.04),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegsign->SetFillColor(kWhite);
  autolegsign->SetLineColor(kWhite);
  gPad->SetGridy();
  // mleg->Draw();
  
  cSoB->cd();
  TLegend * autolegSoB = (TLegend*) gPad->BuildLegend(0.5,(0.88-legendEntryCounter*0.04),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegSoB->SetFillColor(kWhite);
  autolegSoB->SetLineColor(kWhite);
  gPad->SetGridy();
  // gPad->SetLogy();
  // frameSoverBVsPt->GetYaxis()->Unzoom();
  //  mleg->Draw();
  
  TString fimgname(foutname.Data());
  fimgname.ReplaceAll("rawYields/","rawYields/img/");
  fimgname.ReplaceAll(".root",".png");
  cry->SaveAs(fimgname.ReplaceAll("rawYields_","RAWYIELDS_"));
  cfit->SaveAs(fimgname.ReplaceAll("RAWYIELDS_","MASS_"));
  cwi->SaveAs(fimgname.ReplaceAll("MASS_","WIDTH_"));
  cc2->SaveAs(fimgname.ReplaceAll("WIDTH_","CHI2_"));
  cSoB->SaveAs(fimgname.ReplaceAll("CHI2_","SoverB_"));
  csign->SaveAs(fimgname.ReplaceAll("SoverB_","Significance_"));
  return;
}

//-------------------------------------------------------------------------
TH1F* MakeSpectraCent(Int_t quantityID=EQuantity::kYields,Int_t icentbin = 0, TString filein="fitEM_analysisAOD_0-80.root", Char_t* fileproj="sub_analysisAOD_0-80.root", Float_t cutChi2=9999., Bool_t divideByDpt = kFALSE)
{
  gROOT->LoadMacro("SetGraphicStyle.C");
  SetGraphicStyle(0);
  gStyle->SetOptStat(0);
  gStyle->SetTextFont(42);
  
  //Char_t projectionFile[100]="sub_kstar.root";
  TFile * f= TFile::Open(fileproj);
  if (!f) return;
  //get bins
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  Int_t npt = ptbins->GetNbins();
  const Int_t dimpt = npt+1;
  Double_t pt[dimpt];
  for (Int_t k=0; k<dimpt;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
    //    Printf("%5.2f",pt[k]);
  }

  TAxis *centbins = (TAxis*)f->Get("centbins");
  Int_t ncent = centbins->GetNbins();
  const Int_t dimcent = ncent+1;
  Double_t cent[dimcent]; 
  for (Int_t k=0; k<dimcent;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
    //  Printf("%5.2f",cent[k]);
  }

  TFile * fin=TFile::Open(filein.Data());

  Double_t fitParams[9], SoverB=0.0, significance=0.0,normfactorCopy=0.0;
  Int_t centBinID,ptBinID, funcID;
  Double_t rangeInfCopy, rangeSupCopy, ptinfCopy, ptsupCopy, centinfCopy, centsupCopy;
  TString*fitfunction;
  TTree *tree= (TTree*) fin->Get("tree");//get fit parameters tree
  if (!tree)return;
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
  // tree->SetBranchAddress("norm_inf", &rangeInfCopy);
  // tree->SetBranchAddress("norm_sup", &rangeSupCopy);
  tree->SetBranchAddress("pt_inf", &ptinfCopy);
  tree->SetBranchAddress("pt_sup", &ptsupCopy);
  // tree->SetBranchAddress("cent_inf", &centinfCopy);
  // tree->SetBranchAddress("cent_sup", &centsupCopy);
  // tree->SetBranchAddress("functionID", &funcID);

  TH1F*frameMassVsPt=new TH1F("MassVsPt","Mass vs. p_{t}; p_{t} (GeV/c); M (GeV/c^{2})", npt, pt);
  TH1F*frameWidthVsPt=new TH1F("WidthVsPt","Width vs. p_{t}; p_{t} (GeV/c); #Gamma (GeV/c^{2})", npt, pt);
  TH1F*frameYieldVsPt=new TH1F("YieldVsPt","Raw yield vs. p_{t}; p_{t} (GeV/c); dN/dp_{t}", npt, pt);
  TH1F*frameChi2VsPt=new TH1F("Chi2VsPt","#chi^{2} vs. p_{t}; p_{t} (GeV/c); dN/dp_{t}", npt, pt);
  TH1F*frameSoverBVsPt=new TH1F("SoverBVsPt","S/B; p_{t} (GeV/c); S/B", npt, pt);
  TH1F*frameSignificanceVsPt=new TH1F("SignificanceVsPt","S/#sqrt{S+B}; p_{t} (GeV/c); significance", npt, pt);

  TH1F *hMassVsPt;//[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hWidthVsPt;//[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hRawYieldVsPt;//[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hChi2VsPt;//[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hSoverBVsPt;//[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1F *hSignificanceVsPt;//[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};

  TString centLabel=Form("%2.0f-%2.0f %%",cent[icentbin],cent[icentbin+1]);
    
  /* mass vs pT plots */
  hMassVsPt = new TH1F(Form("hMassVsPt_%i",icentbin),"Mass vs p_{t}; p_{t} (GeV/c); M (GeV/c^{2})", npt, pt);
  hMassVsPt->SetTitle(Form("Fitted K* mass Vs p_{t} (%s)",centLabel.Data()));

  /* mass vs pT plots */
  hWidthVsPt = new TH1F(Form("hWidthVsPt_%i",icentbin),"Width vs p_{t}; p_{t} (GeV/c); #Gamma (GeV/c^{2})", npt, pt);
  hWidthVsPt->SetTitle(Form("Fitted K* width Vs p_{t} (%s)",centLabel.Data()));

  /* yield vs pT plots */
  hRawYieldVsPt = new TH1F(Form("hRawYieldVsPt_%i",icentbin),"Raw yield vs p_{t}; p_{t} (GeV/c); dN/dp_{t})", npt, pt);
  hRawYieldVsPt->SetTitle(Form("Raw yields vs p_{t} (%s)",centLabel.Data()));

  /* chi2 vs pT plots */
  hChi2VsPt = new TH1F(Form("hChi2VsPt_%i",icentbin),"#chi^{2} vs p_{t}; p_{t} (GeV/c); dN/dp_{t})", npt, pt);
  hChi2VsPt->SetTitle(Form("#chi^{2} vs p_{t} (%s)",centLabel.Data()));

  /* S/B vs pT plots */
  hSoverBVsPt = new TH1F(Form("hSoverBVsPt_%i",icentbin),"S/B vs p_{t}; p_{t} (GeV/c); S/B)", npt, pt);
  hSoverBVsPt->SetTitle(Form("S/B vs p_{t} (%s)",centLabel.Data()));

  /* signif vs pT plots */
  hSignificanceVsPt = new TH1F(Form("hSignificanceVsPt_%i",icentbin),"Significance vs p_{t}; p_{t} (GeV/c); Significance)", npt, pt);
  hSignificanceVsPt->SetTitle(Form("Significance vs p_{t} (%s)",centLabel.Data()));

  for (Int_t ientry=0;ientry<tree->GetEntries();ientry++){
    tree->GetEntry(ientry);
    if (centBinID==icentbin){
      //apply chi2 cut
      if (fitParams[8]<cutChi2){

	/*fill mass vs pt*/
	//Printf("ptBinID = %i  -  pt[x] = %4.2f", ptBinID, pt[ptBinID]);
	hMassVsPt->SetBinContent(ptBinID+1,fitParams[0]);
	hMassVsPt->SetBinError(ptBinID+1,fitParams[1]);
	hMassVsPt->SetMarkerStyle(marker[icentbin]);
	hMassVsPt->SetMarkerColor(color[icentbin]);
	hMassVsPt->SetLineColor(color[icentbin]);       
	hMassVsPt->GetYaxis()->SetRangeUser(0.86, 0.92);

	/*fill width vs pt*/
	hWidthVsPt->SetBinContent(ptBinID+1,fitParams[2]);
	hWidthVsPt->SetBinError(ptBinID+1,fitParams[3]);
	hWidthVsPt->SetMarkerStyle(marker[icentbin]);
	hWidthVsPt->SetMarkerColor(color[icentbin]);
	hWidthVsPt->SetLineColor(color[icentbin]);       
	hWidthVsPt->GetYaxis()->SetRangeUser(0.01, 0.10);

	/*fill yield vs pt */
	//Printf("nsignal = %f    error = %f",fitParams[4],fitParams[5]);
	Float_t dpt = 1.0;
	if (divideByDpt) dpt = ptbins->GetBinWidth(ptBinID + 1);
	//Printf("ptBIN = %i --> dpt = %5.2f ", ptBinID, dpt);
	hRawYieldVsPt->SetBinContent(ptBinID+1,fitParams[4]/dpt);
	hRawYieldVsPt->SetBinError(ptBinID+1,fitParams[5]/dpt);//how to estimate error on the integral????
	hRawYieldVsPt->SetMarkerStyle(marker[icentbin]);
	hRawYieldVsPt->SetMarkerColor(color[icentbin]);
	hRawYieldVsPt->SetLineColor(color[icentbin]);       
	hRawYieldVsPt->GetYaxis()->SetRangeUser(1, 1e7);

	/*fill chi2 vs pt */
	//Printf("chi2 = %f",fitParams[8]);
	hChi2VsPt->SetBinContent(ptBinID+1,fitParams[8]);
	hChi2VsPt->SetMarkerStyle(marker[icentbin]);
	hChi2VsPt->SetMarkerColor(color[icentbin]);
	hChi2VsPt->SetLineColor(color[icentbin]);       
	hChi2VsPt->GetYaxis()->SetRangeUser(0., 15.);

	/*fill S/B vs pt */
	hSoverBVsPt->SetBinContent(ptBinID+1,SoverB);
	hSoverBVsPt->SetMarkerStyle(marker[icentbin]);
	hSoverBVsPt->SetMarkerColor(color[icentbin]);
	hSoverBVsPt->SetLineColor(color[icentbin]);       
	hSoverBVsPt->GetYaxis()->SetRangeUser(0., 10.);

	/*fillsignif vs pt */
	//Printf("significance = %f",significance);
	hSignificanceVsPt->SetBinContent(ptBinID+1,significance);
	hSignificanceVsPt->SetMarkerStyle(marker[icentbin]);
	hSignificanceVsPt->SetMarkerColor(color[icentbin]);
	hSignificanceVsPt->SetLineColor(color[icentbin]);       
	hSignificanceVsPt->GetYaxis()->SetRangeUser(0., 300.);
      }//end cut on chi2
    }    
  }//loop on tree entries
  if (quantityID==EQuantity::kMass) return hMassVsPt;
  if (quantityID==EQuantity::kYields) return hRawYieldVsPt; 
  if (quantityID=EQuantity::kWidth) return hWidthVsPt;
}


//---------------------------------------------------------
void compareFitFuncPerCent(Int_t centBin=0, TString filefit="roofit/fitEM_BW_POLY2_0.74-1.10.root",Float_t cutChi2=9999.)
{
  //for the selected centrality bins it compares the result of the fit from different functions
 
  gROOT->LoadMacro("SetGraphicStyle.C");
  SetGraphicStyle(0);
  gStyle->SetOptStat(0);
  gStyle->SetTextFont(42);
  
  //Char_t projectionFile[100]="sub_kstar.root";
  TFile * f= TFile::Open(filefit);
  if (!f) return;
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  Int_t nPtBins = ptbins->GetNbins();
  
  Char_t nomefile[5][300];
  sprintf(nomefile[0],"%s", filefit.Data());
  sprintf(nomefile[1],"%s", filefit.ReplaceAll("POLY2","POLY3"));
  sprintf(nomefile[2],"%s", filefit.ReplaceAll("POLY2","LAND"));
  //sprintf(nomefile[3],"roofit/fitEM_BW_EXP_0.74-1.10.root");
  //sprintf(nomefile[4],"roofit/fitEM_BW_POLY1_0.74-1.10.root");
  
  Double_t fitParams[9];
  Int_t ptBinID=-1,centBinID=-1;
  TString centLabel=Form("Centrality %2.0f-%2.0f %%",cent[centBin],cent[centBin+1]);
  TPaveText *centText=new TPaveText(0.25,0.83,0.55,0.88,"NDC");
  centText->AddText(centLabel);

  TFile * fout = new TFile(Form("roofit/systematics_func_cent%i.root",centBin),"recreate");
  TH1F*frameMassVsPt=new TH1F("MassVsPt","Mass vs. p_{t}; p_{t} (GeV/c); M (GeV/c^{2})", npt, pt);
  TH1F*frameWidthVsPt=new TH1F("WidthVsPt","Width vs. p_{t}; p_{t} (GeV/c); #Gamma (GeV/c^{2})", npt, pt);
  TH1F*frameYieldVsPt=new TH1F("YieldVsPt","Raw yield vs. p_{t}; p_{t} (GeV/c); dN/dp_{t}", npt, pt);
  TH1F*frameChi2VsPt=new TH1F("Chi2VsPt","#chi^{2} vs. p_{t}; p_{t} (GeV/c); dN/dp_{t}", npt, pt);

  //sistematics
  TH1F*frameMassVsPtSys=new TH1F("MassVsPtSys","Systematics mass; p_{t} (GeV/c); M (GeV/c^{2})", npt, pt);
  TH1F*frameWidthVsPtSys=new TH1F("WidthVsPtSys","Systematics width; p_{t} (GeV/c); #Gamma (GeV/c^{2})", npt, pt);
  TH1F*frameYieldVsPtSys=new TH1F("YieldVsPtSys","Systematics raw yields; p_{t} (GeV/c); dN/dp_{t}", npt, pt);
 
  TLegend *mleg=new TLegend(0.7,0.65,0.92,0.88);
  mleg->SetHeader(" Fit function:");
  mleg->SetFillColor(kWhite);
  
  TCanvas *cfit=new TCanvas("MassVsPt","Mass vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  cfit->cd();
  frameMassVsPt->GetYaxis()->SetRangeUser(0.86, 0.94);
  frameMassVsPt->Draw();
    
  TCanvas *cwi=new TCanvas("WidthVsPt","Width vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  cwi->cd();
  frameWidthVsPt->GetYaxis()->SetRangeUser(0.01, 0.1);
  frameWidthVsPt->Draw();

  TCanvas *cry=new TCanvas("cry","Raw yield vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  cry->cd();
  frameYieldVsPt->GetYaxis()->SetRangeUser(1, 1e7);
  frameYieldVsPt->Draw();

  TCanvas *cc2=new TCanvas("cc2","#Chi^{2} vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  cc2->cd();
  frameChi2VsPt->GetYaxis()->SetRangeUser(-1, 10);
  frameChi2VsPt->Draw();
  
  TCanvas *cfitSys=new TCanvas("MassVsPtSys","Sys. Mass vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  cfitSys->cd();
  frameMassVsPtSys->GetYaxis()->SetRangeUser(0.,10.);
  frameMassVsPtSys->GetYaxis()->SetTitle("systematic error (%)");
  frameMassVsPtSys->Draw();
    
  TCanvas *cwiSys=new TCanvas("WidthVsPtSys","Sys. Width vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  cwiSys->cd();
  frameWidthVsPtSys->GetYaxis()->SetRangeUser(0., 100.);
  frameWidthVsPtSys->GetYaxis()->SetTitle("systematic error (%)");
  frameWidthVsPtSys->Draw();

  TCanvas *crySys=new TCanvas("crySys","Sys. Raw yield vs p_{t}", 800,600);
  //cfit->Divide(3,2);
  crySys->cd();
  frameYieldVsPtSys->GetYaxis()->SetRangeUser(0.0, 100.0);
  frameYieldVsPtSys->GetYaxis()->SetTitle("systematic error (%)");
  frameYieldVsPtSys->Draw();

  TH1F *hMassVsPt[5] =  {0x0,0x0,0x0,0x0, 0x0};
  TH1F *hWidthVsPt[5] =  {0x0,0x0,0x0,0x0, 0x0};
  TH1F *hRawYieldVsPt[5] =  {0x0,0x0,0x0,0x0, 0x0};
  TH1F *hChi2VsPt[5] =  {0x0,0x0,0x0,0x0, 0x0};
 
  //systematics
  TH1F *hMassVsPtSys[5] =  {0x0,0x0,0x0,0x0, 0x0};
  TH1F *hWidthVsPtSys[5] =  {0x0,0x0,0x0,0x0, 0x0};
  TH1F *hRawYieldVsPtSys[5] =  {0x0,0x0,0x0,0x0, 0x0};
  //central values
  Double_t centralMass[5][15];
  Double_t centralWidth[5][15];
  Double_t centralRawYield[5][15];
  
  for (Int_t ic=0;ic<ncent;ic++){
    Printf("================= CENTRALITY BIN %i =====================", ic);
    Printf("Reference file for central value:  %s", nomefile[2]);
    for (Int_t ipt=0;ipt<npt;ipt++){
      centralMass[ic][ipt]=GetCentralValue(ic, ipt, "m", nomefile[0]);
      centralWidth[ic][ipt]=GetCentralValue(ic, ipt, "w", nomefile[0]);
      centralRawYield[ic][ipt]=GetCentralValue(ic, ipt, "y", nomefile[0]);
      Printf("Pt bin %i: mass = %e, width = %e, rawYield = %e", ipt,centralMass[ic][ipt],centralWidth[ic][ipt],centralRawYield[ic][ipt]);
    }
  }

  for (Int_t ifunc=0;ifunc<3;ifunc++){
    TFile *file = TFile::Open(nomefile[ifunc]);
    
    TTree *tree= (TTree*) file->Get("tree");//get fit parameters tree
    if (!tree)return;
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
    
  // Float_t ptbinwidth[15];
  // for (Int_t j=0;j<15;j++){
    //if (j>=npt) ptbinwidth[j]=0.0;
  //   else ptbinwidth[j]=pt[j+1]-pt[j];
  // }    
    
    /* mass vs pT plots */
    hMassVsPt[ifunc] = new TH1F(Form("hMassVsPt_%i",ifunc),"Mass vs p_{t}; p_{t} (GeV/c); M (GeV/c^{2})", npt, pt);
    hMassVsPt[ifunc]->SetTitle(Form("Fitted K* mass Vs p_{t} (%s)",centLabel.Data()));

    /* mass vs pT plots */
    hWidthVsPt[ifunc] = new TH1F(Form("hWidthVsPt_%i",ifunc),"Width vs p_{t}; p_{t} (GeV/c); #Gamma (GeV/c^{2})", npt, pt);
    hWidthVsPt[ifunc]->SetTitle(Form("Fitted K* width Vs p_{t} (%s)",centLabel.Data()));

    /* yield vs pT plots */
    hRawYieldVsPt[ifunc] = new TH1F(Form("hRawYieldVsPt_%i",ifunc),"Raw yield vs p_{t}; p_{t} (GeV/c); dN/dp_{t})", npt, pt);
    hRawYieldVsPt[ifunc]->SetTitle(Form("Raw yields vs p_{t} (%s)",centLabel.Data()));

   /* chi2 vs pT plots */
    hChi2VsPt[ifunc] = new TH1F(Form("hChi2VsPt_%i",ifunc),"#chi^{2} vs p_{t}; p_{t} (GeV/c); dN/dp_{t})", npt, pt);
    hChi2VsPt[ifunc]->SetTitle(Form("#chi^{2} vs p_{t} (%s)",centLabel.Data()));

    /* syst mass vs pT plots */
    hMassVsPtSys[ifunc] = new TH1F(Form("hMassVsPtSys_%i",ifunc),"Mass vs p_{t}; p_{t} (GeV/c); M-M_{cent} (GeV/c^{2})", npt, pt);
    hMassVsPtSys[ifunc]->SetTitle(Form("Systematics of fitted K* mass Vs p_{t} (%s)",centLabel.Data()));
    
    /* syst mass vs pT plots */
    hWidthVsPtSys[ifunc] = new TH1F(Form("hWidthVsPtSys_%i",ifunc),"Width vs p_{t}; p_{t} (GeV/c); #Gamma-#Gamma_{cent} (GeV/c^{2})", npt, pt);
    hWidthVsPtSys[ifunc]->SetTitle(Form("Systematics of fitted K* width Vs p_{t} (%s)",centLabel.Data()));

    /* syst yield vs pT plots */
    hRawYieldVsPtSys[ifunc] = new TH1F(Form("hRawYieldVsPtSys_%i",ifunc),"Raw yield vs p_{t}; p_{t} (GeV/c); dN/dp_{t}-(dN/dp_{t})_{cent})", npt, pt);
    hRawYieldVsPtSys[ifunc]->SetTitle(Form("Systematics of raw yields vs p_{t} (%s)",centLabel.Data()));

    for (Int_t ientry=0;ientry<tree->GetEntries();ientry++){
      tree->GetEntry(ientry);
      if (centBinID==centBin){	
	//Printf("ptBins[x]=%5.2f    pt[x]=%5.2f",ptbins->GetBinLowEdge(ptBinID+1),pt[ptBinID]);
	if (fitParams[8]<cutChi2){
	/*fill mass vs pt*/
	//  Printf("BIN PT = %i",ptBinID+1);
	hMassVsPt[ifunc]->SetBinContent(ptBinID+1,fitParams[0]);
	hMassVsPt[ifunc]->SetBinError(ptBinID+1,fitParams[1]);
	// fin->cd();
	// hMassVsPt->Write();
	hMassVsPt[ifunc]->SetMarkerStyle(marker[ifunc]);
	hMassVsPt[ifunc]->SetMarkerColor(colorFunc[ifunc]);
	hMassVsPt[ifunc]->SetLineColor(colorFunc[ifunc]);       
	hMassVsPt[ifunc]->GetYaxis()->SetRangeUser(0.86, 0.92);

	/*fill width vs pt*/
	hWidthVsPt[ifunc]->SetBinContent(ptBinID+1,fitParams[2]);
	hWidthVsPt[ifunc]->SetBinError(ptBinID+1,fitParams[3]);
	hWidthVsPt[ifunc]->SetMarkerStyle(marker[ifunc]);
	hWidthVsPt[ifunc]->SetMarkerColor(colorFunc[ifunc]);
	hWidthVsPt[ifunc]->SetLineColor(colorFunc[ifunc]);       
	hWidthVsPt[ifunc]->GetYaxis()->SetRangeUser(0.01, 0.1);

	/*fill yield vs pt */
	Float_t dpt = ptbins->GetBinWidth(ptBinID+1);
	//Printf("#####  nsignal = %f    error = %f",fitParams[4]/dpt,fitParams[5]);
	hRawYieldVsPt[ifunc]->SetBinContent(ptBinID+1,fitParams[4]/dpt);
	hRawYieldVsPt[ifunc]->SetBinError(ptBinID+1,fitParams[5]/dpt);//how to estimate error on the integral????
	hRawYieldVsPt[ifunc]->SetMarkerStyle(marker[ifunc]);
	hRawYieldVsPt[ifunc]->SetMarkerColor(colorFunc[ifunc]);
	hRawYieldVsPt[ifunc]->SetLineColor(colorFunc[ifunc]);       
	hRawYieldVsPt[ifunc]->GetYaxis()->SetRangeUser(1, 1e7);
	
	/*fill chi2 vs pt */
	//Printf("chi2 = %f",fitParams[8]);
	hChi2VsPt[ifunc]->SetBinContent(ptBinID+1,fitParams[8]);
	hChi2VsPt[ifunc]->SetMarkerStyle(marker[ifunc]);
	hChi2VsPt[ifunc]->SetMarkerColor(colorFunc[ifunc]);
	hChi2VsPt[ifunc]->SetLineColor(colorFunc[ifunc]);       
	hChi2VsPt[ifunc]->GetYaxis()->SetRangeUser(0., 5.);

	/*fill syst mass vs pt*/
	// if (pt[ptBinID]==(hMassVsPtSys[ifunc]->GetXaxis()->GetBinLowEdge(ptBinID+1))) Printf("ciccio");
	if (centralMass[centBin][ptBinID]>0){
	  Double_t systM = TMath::Abs(fitParams[0]-centralMass[centBin][ptBinID])/centralMass[centBin][ptBinID];
	  hMassVsPtSys[ifunc]->SetBinContent(ptBinID+1,systM*100.);
	  hMassVsPtSys[ifunc]->SetMarkerStyle(marker[ifunc]);
	  hMassVsPtSys[ifunc]->SetMarkerColor(colorFunc[ifunc]);
	  hMassVsPtSys[ifunc]->SetLineColor(colorFunc[ifunc]); 
	  hMassVsPtSys[ifunc]->GetYaxis()->SetTitle("systematic error (%)");      
	}
	  /*fill syst width vs pt*/
	if (centralWidth[centBin][ptBinID]>0){
	  Double_t systW = TMath::Abs(fitParams[2]-centralWidth[centBin][ptBinID])/centralWidth[centBin][ptBinID];
	  hWidthVsPtSys[ifunc]->SetBinContent(ptBinID+1,systW*100.0);
	  hWidthVsPtSys[ifunc]->SetMarkerStyle(marker[ifunc]);
	  hWidthVsPtSys[ifunc]->SetMarkerColor(colorFunc[ifunc]);
	  hWidthVsPtSys[ifunc]->SetLineColor(colorFunc[ifunc]);       
	  hWidthVsPtSys[ifunc]->GetYaxis()->SetTitle("systematic error (%)");
	}
	/*fill syst yield vs pt */
	if (centralRawYield[centBin][ptBinID]>0){
	  Double_t systRY = TMath::Abs(fitParams[4]-centralRawYield[centBin][ptBinID])/ (centralRawYield[centBin][ptBinID]);
	  Printf("ptbin = %i, cent bin = %i, central raw y = %e  \n fitparams[4]/dpt = %e ----> systRY = %5.2f ",ptBinID, centBin, centralRawYield[centBin][ptBinID]/dpt,fitParams[4]/dpt, systRY*100.);
	  hRawYieldVsPtSys[ifunc]->SetBinContent(ptBinID+1,systRY*100.);
	  hRawYieldVsPtSys[ifunc]->SetMarkerStyle(marker[ifunc]);
	  hRawYieldVsPtSys[ifunc]->SetMarkerColor(colorFunc[ifunc]);
	  hRawYieldVsPtSys[ifunc]->SetLineColor(colorFunc[ifunc]);  
	  hRawYieldVsPtSys[ifunc]->GetYaxis()->SetTitle("systematic error (%)");
	}
      }//end if on chi2
      }//end cent selection
    }//end tree entries
 
    //Draw on canvas
    cfit->cd();
    mleg->AddEntry(hMassVsPt[ifunc], funcName[ifunc],"lpf");
    hMassVsPt[ifunc]->Draw("same");
    cwi->cd();
    hWidthVsPt[ifunc]->Draw("same");
    cry->cd();
    gPad->SetLogy();
    hRawYieldVsPt[ifunc]->Draw("same");
    cc2->cd();
    hChi2VsPt[ifunc]->Draw("same");
    
    cfitSys->cd();
    //hMassVsPtSys[ifunc]->GetYaxis()->SetRangeUser(0.,100.);
    hMassVsPtSys[ifunc]->Draw("same");
    centText->Draw("same");
    cwiSys->cd();
    //hWidthVsPtSys[ifunc]->GetYaxis()->SetRangeUser(0.,100.);
    hWidthVsPtSys[ifunc]->Draw("same");
    centText->Draw("same");
    crySys->cd();
    //hRawYieldVsPtSys[ifunc]->GetYaxis()->SetRangeUser(0.,100.);
    hRawYieldVsPtSys[ifunc]->Draw("same");
    centText->Draw("same");
    //Save on file
    fout->cd();
    if (save){
      hMassVsPt[ifunc]->Write();
      hWidthVsPt[ifunc]->Write();
      hRawYieldVsPt[ifunc]->Write();
      hRawYieldVsPtSys[ifunc]->Write();
      hMassVsPtSys[ifunc]->Write();
      hWidthVsPtSys[ifunc]->Write();
      hChi2VsPt[ifunc]->Write();
      Printf("Opened file %s", file->GetName());
      Printf("Histos saved in file %s", fout->GetName());
    }
  }//end loop on files

  cry->cd();
  mleg->Draw();
  centText->Draw("same");
  gPad->SetGridy();
  cfit->cd();
  mleg->Draw();
  centText->Draw("same");
  gPad->SetGridy();
  cwi->cd();
  mleg->Draw();
  centText->Draw("same");
  gPad->SetGridy();
  cc2->cd();
  mleg->Draw();
  centText->Draw("same");
  gPad->SetGridy();
  cry->SaveAs(Form("rawYields/img/syst_func_RAWYIELDS_cent%i.png",centBin));
  cfit->SaveAs(Form("rawYields/img/syst_func_MASS_cent%i.png",centBin));
  cwi->SaveAs(Form("rawYields/img/syst_func_WIDTH_cent%i.png",centBin));
  cc2->SaveAs(Form("rawYields/img/syst_func_CHI2_cent%i.png",centBin));
  cfitSys->SaveAs(Form("rawYields/img/systPerc_goodChi1.5_func_MASS_cent%i.png",centBin));
  cwiSys->SaveAs(Form("rawYields/img/systPerc_goodChi1.5_func_WIDTH_cent%i.png",centBin));
  crySys->SaveAs(Form("rawYields/img/systPerc_goodChi1.5_func_RAWYIELDS_cent%i.png",centBin));
  return;
}

 
 

//----------------------------------------------------
Double_t GetCentralValue(Int_t icent=0, Int_t ipt=0, TString value = "y", TString nomefile="")
{
  //returns the central value of the specified quantities - from chosen file  
  //need to specify the file with the fit for a given function (poly2?)
  TFile *file = TFile::Open(nomefile.Data());
  if (!file) return 0.0;
  
  TTree *tree= (TTree*) file->Get("tree");//get fit parameters tree
  if (!tree)return;
  
  Double_t fitParams[9];
  Int_t ptBinID,centBinID;

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
  
  Double_t result=0.0;
  for (Int_t ientry=0;ientry<tree->GetEntries();ientry++){
      tree->GetEntry(ientry);
      if ((centBinID==icent)&&(ptBinID==ipt)){
	  if (value.Contains("m")) result = fitParams[0];
	  if (value.Contains("w")) result = fitParams[2];
	  if (value.Contains("y")) result = fitParams[4];
	}
  }
  file->Close();
  return result;  
}
