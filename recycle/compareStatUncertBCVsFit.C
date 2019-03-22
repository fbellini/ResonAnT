void compareStatUncertBCVsFit(TString filein="best_fit_poly2.root", TString fileproj="../_sum_EMnorm1.30-1.50_cut1818_tof2s_train215-216.root")
{
  //graphics
  gROOT->LoadMacro("SetGraphicStyle.C");
  SetGraphicStyle(0);
  gStyle->SetTextFont(42);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  Color_t color[2][6]={kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack, //tpc
		       kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof
  Int_t marker[2][6]={21, 22, 23, 34, 33, 20, //tpc
		      25, 26, 32, 28, 27, 24}; //tof
  
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
  }
  
  // if tof analysis show only points in 1.0-8.0 GeV/c range, use 0.5-8.0 for TPC
  Bool_t isTOF = fileproj.Contains("tof");
  
  TFile * fin=TFile::Open(filein.Data());
  Double_t signal, signalErr, nrawbc, nrawbcErr, nhistoErr;
  Int_t centBinID,ptBinID;
  
  TTree *tree= (TTree*) fin->Get("tree");//get fit parameters tree
  if (!tree) return;
  tree->SetBranchAddress("nSignal",&signal);
  tree->SetBranchAddress("nSignalErr",&signalErr);
  tree->SetBranchAddress("centBin",&centBinID);
  tree->SetBranchAddress("ptBin",&ptBinID);
  tree->SetBranchAddress("nrawbc",&nrawbc);
  tree->SetBranchAddress("nrawbcErr",&nrawbcErr);
  tree->SetBranchAddress("nHistoErr",&nhistoErr);

  TLegend *mleg=new TLegend(0.8,0.65,0.92,0.88);
  mleg->SetHeader(" Centrality:");
  mleg->SetFillColor(kWhite);
  
  TCanvas *cry=new TCanvas("cry","Raw yield vs p_{T}", 600,700);
  TCanvas *cstat=new TCanvas("cstat","Statistical uncertainties", 600,700);
  TCanvas *cstat2=new TCanvas("cstat2","Statistical uncertainties", 700,400);
  
  TH1F *hRawYieldVsPt[5] =  {0x0,0x0,0x0,0x0,0x0};
  TH1F *hRawYieldStatVsPt[5] =  {0x0,0x0,0x0,0x0,0x0};
  TH1F *hBCVsPt[5] =  {0x0,0x0,0x0,0x0,0x0};
  TH1F *hBCStatVsPt[5] =  {0x0,0x0,0x0,0x0,0x0};
  TH1F *hHistoStatVsPt[5] =  {0x0,0x0,0x0,0x0,0x0};
  TH1F *hBCSqrt[5] =  {0x0,0x0,0x0,0x0,0x0};
  TH1F *hBCStat2FitStat[5] =  {0x0,0x0,0x0,0x0,0x0};
  Int_t  legendEntryCounter = 0;  

  /*loop on centrality bins*/
  for (Int_t icentbin=0;icentbin<ncent;icentbin++){
    TString centLabel=Form("%2.0f-%3.0f%%",cent[icentbin],cent[icentbin+1]);
    
    /* yield vs pT plots */
    hRawYieldVsPt[icentbin] = new TH1F(Form("hRawYieldVsPt_%i",icentbin),"; p_{T} (GeV/c); dN/dp_{T})", npt, pt);
    hRawYieldVsPt[icentbin]->SetTitle(Form("Raw yields from fit (%s)",centLabel.Data()));

    hRawYieldStatVsPt[icentbin] = new TH1F(Form("hRawYieldStatVsPt_%i",icentbin),"; p_{T} (GeV/c); stat. uncert.)", npt, pt);
    hRawYieldStatVsPt[icentbin]->SetTitle(Form("Stat. uncert. from fit (%s)",centLabel.Data()));
    
    /* yield from BC vs pT plots */
    hBCVsPt[icentbin] = new TH1F(Form("hBCVsPt_%i",icentbin),"; p_{T} (GeV/c); dN/dp_{T})", npt, pt);
    hBCVsPt[icentbin]->SetTitle(Form("Raw yields from bin counting (%s)",centLabel.Data()));
    
    hBCStatVsPt[icentbin] = new TH1F(Form("hBCStatVsPt_%i",icentbin),"; p_{T} (GeV/c); stat. uncert.)", npt, pt);
    hBCStatVsPt[icentbin]->SetTitle(Form("Stat. uncert. from bin counting (%s)",centLabel.Data()));
  
    /* ratio stat. uncert */
    hBCSqrt[icentbin] = new TH1F(Form("hBCSqrt_%i",icentbin),"; p_{T} (GeV/c); ratio)", npt, pt);
    hBCSqrt[icentbin]->SetTitle(Form(" #sqrt{dN_{bin count}/dp_{T}}/#sqrt{dN_{fit}/dp_{T}} (%s)",centLabel.Data()));
    
    hBCStat2FitStat[icentbin] = new TH1F(Form("hBCStat2FitStat_%i",icentbin),"; p_{T} (GeV/c); ratio)", npt, pt);
    hBCStat2FitStat[icentbin]->SetTitle(Form("Stat. uncert. from bin counting / fit (%s)",centLabel.Data()));
    
    // hHistoStatVsPt[icentbin] = new TH1F(Form("hHistoStatVsPt_%i",icentbin),"; p_{T} (GeV/c); stat. uncert.)", npt, pt);
    // hHistoStatVsPt[icentbin]->SetTitle(Form("#sqrt(S+resB) from histo (%s)",centLabel.Data()));

    for (Int_t ientry=0;ientry<tree->GetEntries();ientry++){

      tree->GetEntry(ientry);
      
      if (centBinID==icentbin){
	/*fill yield vs pt */
	hRawYieldVsPt[icentbin]->SetBinContent(ptBinID+1,signal);
	hRawYieldVsPt[icentbin]->SetBinError(ptBinID+1,signalErr);

 	/*fill yield stat uncert vs pt */
	hRawYieldStatVsPt[icentbin]->SetBinContent(ptBinID+1,signalErr);
	hRawYieldStatVsPt[icentbin]->SetBinError(ptBinID+1,0.0);
	
	/*fill yield vs pt */
	hBCVsPt[icentbin]->SetBinContent(ptBinID+1,nrawbc);
	hBCVsPt[icentbin]->SetBinError(ptBinID+1, nrawbcErr);
	
   	/*fill yield stat uncert vs pt */
	hBCStatVsPt[icentbin]->SetBinContent(ptBinID+1, nrawbcErr);
	hBCStatVsPt[icentbin]->SetBinError(ptBinID+1,0.0);

   	/*fill yield stat uncert ratio vs pt */
	hBCStat2FitStat[icentbin]->SetBinContent(ptBinID+1, nrawbcErr/signalErr);
	hBCStat2FitStat[icentbin]->SetBinError(ptBinID+1,0.0);
 
	/*fill sqrt(yield) vs pt */
	hBCSqrt[icentbin]->SetBinContent(ptBinID+1, signalErr/signal);
	hBCSqrt[icentbin]->SetBinError(ptBinID+1,0.0);

	/*fill sqrt(S+res.bg)_HISTO vs pt */
	// hHistoStatVsPt[icentbin]->SetBinContent(ptBinID+1, nhistoErr);
	// hHistoStatVsPt[icentbin]->SetBinError(ptBinID+1,0.0);
	
     }
    }
    
    hRawYieldVsPt[icentbin]->SetMarkerStyle(marker[isTOF][icentbin]);
    hRawYieldVsPt[icentbin]->SetMarkerColor(color[isTOF][icentbin]);
    hRawYieldVsPt[icentbin]->SetLineColor(color[isTOF][icentbin]);  
    hRawYieldVsPt[icentbin]->SetLineWidth(1);
    hRawYieldVsPt[icentbin]->GetYaxis()->SetRangeUser(1e1, 1e7);
    hRawYieldVsPt[icentbin]->GetYaxis()->SetTitleOffset(1.5);
    hRawYieldVsPt[icentbin]->GetYaxis()->SetTitle("dN/dp_{T}");
	
    hRawYieldStatVsPt[icentbin]->SetMarkerStyle(marker[isTOF][icentbin]);
    hRawYieldStatVsPt[icentbin]->SetMarkerColor(color[isTOF][icentbin]);
    hRawYieldStatVsPt[icentbin]->SetLineColor(color[isTOF][icentbin]);  
    hRawYieldStatVsPt[icentbin]->SetLineWidth(1);
    hRawYieldStatVsPt[icentbin]->GetYaxis()->SetRangeUser(1e-1, 1e7);
    hRawYieldStatVsPt[icentbin]->GetYaxis()->SetTitleOffset(1.5);
    hRawYieldStatVsPt[icentbin]->GetYaxis()->SetTitle("stat. uncert.");
	
    hBCVsPt[icentbin]->SetMarkerStyle(marker[!isTOF][icentbin]);
    hBCVsPt[icentbin]->SetMarkerColor(color[isTOF][icentbin]-5);
    hBCVsPt[icentbin]->SetLineColor(color[isTOF][icentbin]-5);  
    hBCVsPt[icentbin]->SetLineWidth(1);
    hBCVsPt[icentbin]->GetYaxis()->SetRangeUser(1e1, 1e7);
    hBCVsPt[icentbin]->GetYaxis()->SetTitleOffset(1.5);
    hBCVsPt[icentbin]->GetYaxis()->SetTitle("dN/dp_{T}");
	
    // hHistoStatVsPt[icentbin]->SetMarkerStyle(marker[!isTOF][icentbin]);
    // hHistoStatVsPt[icentbin]->SetMarkerColor(color[isTOF][icentbin]-5);
    // hHistoStatVsPt[icentbin]->SetLineColor(color[isTOF][icentbin]-5);  
    // hHistoStatVsPt[icentbin]->SetLineWidth(1);
    // hHistoStatVsPt[icentbin]->GetYaxis()->SetRangeUser(1e1, 1e7);
    // hHistoStatVsPt[icentbin]->GetYaxis()->SetTitleOffset(1.5);
    // hHistoStatVsPt[icentbin]->GetYaxis()->SetTitle("dN/dp_{T}");
	
    hBCStatVsPt[icentbin]->SetMarkerStyle(marker[!isTOF][icentbin]);
    hBCStatVsPt[icentbin]->SetMarkerColor(color[isTOF][icentbin]-5);
    hBCStatVsPt[icentbin]->SetLineColor(color[isTOF][icentbin]-5);  
    hBCStatVsPt[icentbin]->SetLineWidth(1);
    hBCStatVsPt[icentbin]->GetYaxis()->SetRangeUser(1e1, 1e7);
    hBCStatVsPt[icentbin]->GetYaxis()->SetTitleOffset(1.5);
    hBCStatVsPt[icentbin]->GetYaxis()->SetTitle("stat. uncert.");

    hBCStat2FitStat[icentbin]->SetMarkerStyle(marker[isTOF][icentbin]);
    hBCStat2FitStat[icentbin]->SetMarkerColor(color[isTOF][icentbin]);
    hBCStat2FitStat[icentbin]->SetLineColor(color[isTOF][icentbin]);  
    hBCStat2FitStat[icentbin]->SetLineWidth(2);
    hBCStat2FitStat[icentbin]->GetYaxis()->SetRangeUser(0.5, 1.3);
    hBCStat2FitStat[icentbin]->GetYaxis()->SetTitleOffset(1.2);
    hBCStat2FitStat[icentbin]->GetYaxis()->SetTitle("stat. uncert. ratio");

    hBCSqrt[icentbin]->SetMarkerStyle(1);
    hBCSqrt[icentbin]->SetMarkerColor(color[isTOF][icentbin]-5);
    hBCSqrt[icentbin]->SetLineColor(color[isTOF][icentbin]-5);  
    hBCSqrt[icentbin]->SetLineWidth(3);
    hBCSqrt[icentbin]->SetLineStyle(3);
    hBCSqrt[icentbin]->GetYaxis()->SetRangeUser(0.5, 1.3);
    hBCSqrt[icentbin]->GetYaxis()->SetTitleOffset(1.2);
	
    cry->cd();
    gPad->SetLogy();
    if (hRawYieldVsPt[icentbin] && hRawYieldVsPt[icentbin]->GetEntries()>0){
      if (icentbin>0) 
	hRawYieldVsPt[icentbin]->Draw("same");
      else 
	hRawYieldVsPt[icentbin]->Draw();
      hBCVsPt[icentbin]->Draw("same");
      legendEntryCounter+=1;
    }

    cstat->cd();
    gPad->SetLogy();
    if (hRawYieldStatVsPt[icentbin] && hRawYieldStatVsPt[icentbin]->GetEntries()>0){
      if (icentbin>0) 
	hRawYieldStatVsPt[icentbin]->Draw("histosame");
      else 
       	hRawYieldStatVsPt[icentbin]->Draw();
      hBCStatVsPt[icentbin]->Draw("histosame");
    }

    cstat2->cd();
    if (icentbin>0) 
    hBCStat2FitStat[icentbin]->Draw("same");
    else 
      hBCStat2FitStat[icentbin]->Draw();
    //hBCSqrt[icentbin]->Draw("histsame");
  }
  
  gROOT->LoadMacro("$ASD/AddPaveText.C");
  cry->cd();
  TLegend * autolegry = (TLegend*) gPad->BuildLegend(0.4,(0.88-legendEntryCounter*0.03),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegry->SetFillColor(kWhite);
  autolegry->SetLineColor(kWhite);
  
  cstat->cd();
  TLegend * autolegstat = (TLegend*) gPad->BuildLegend(0.4,(0.88-legendEntryCounter*0.03),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegstat->SetFillColor(kWhite);
  autolegstat->SetLineColor(kWhite);
 
  cstat2->cd();
  TLegend * autolegstat2 = (TLegend*) gPad->BuildLegend(0.4,(0.88-legendEntryCounter*0.03),0.88,0.88/*,frameMassVsPt->GetTitle()*/);
  autolegstat2->SetFillColor(kWhite);
  autolegstat2->SetLineColor(kWhite);
  
 TString fimgname(filein.Data());
  fimgname.ReplaceAll(".root",".png");
  //    cry->SaveAs(Form("RAW_%s",fimgname.Data()));
  return;
}
