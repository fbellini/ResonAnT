//gROOT->LoadMacro("$ASD/kstar/MakeSpectra.C");
//gROOT->LoadMacro("$ASD/makeTreeChain.C");
//gROOT->LoadMacro("$ASD/GetPlotRatio.C");
//gROOT->LoadMacro("$ASD/AddPaveText.C");
enum ESystType { kHalfDiff,
		 kMaxDiff,
		 kRatio,
		 kRMS,
		 kUniform,
		 kBinaryTypes, //ways to get syst between 2 measurements
		 // kMaxNDiff,
		 // kUniformN,
		 // kGauss,
		 kNtypes}; //ways to get syst between >2 measurements

  enum EQuantity{ kYields,
		  kMass,
		  kWidth};


void systBinary_pA_PID(Bool_t enaBarlow = 0/*no Barlow check*/,Bool_t display=0, TString deteco = "")
{    
  //21 march 2014
  TString fPath = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/0to100_binB/systPID");  
  TString fCentral = "CORRECTED_br_best_fit_poly2.root"; 
  TString fOther   = "CORRECTED_tpc3s_tof3sveto_br_best_fit_poly2.root"; 
  TString binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/0to100_binB/proj_2323_tpc2s_tof4sveto.root";
  systBinaryFromHisto(fPath, 
		      fCentral, 
		      fOther, 
		      "_tpc3s_tof3sveto.root",
		      0, //all cents
		      EQuantity::kYields, 
		      ESystType::kRMS, 
		      5.0, 
		      binsFile.Data(),
		      enaBarlow, 
		      display,
		      kTRUE, //use corrected spectra
		      deteco.Data());  
  return;
}


void systBinary_pA_PIDstdalone(TString deteco = "TOF",Bool_t enaBarlow = 0/*no Barlow check*/,Bool_t display=0)
{    
  //02 september 2013
  TString fPath = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/systematics",deteco.Data());  
  TString fCentral = "CORRECTED_best_fit_poly2.root"; 
  TString fOther   = "CORRECTED_best_fit_poly2_pid25s.root"; 
  TString binsFile;
  if (deteco.Contains("TOF")) {
    binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1818_tof2s_train215-216.root";
  }  else {
    binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TPC_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1717_tpc2s_train215-216.root";
  }
  systBinaryFromHisto(fPath, 
		      fCentral, 
		      fOther, 
		      "_pid25s.root",
		      -1, //all cents
		      EQuantity::kYields, 
		      ESystType::kRMS, 
		      5.0, 
		      binsFile.Data(),
		      enaBarlow, 
		      display,
		      kTRUE //use corrected spectra
		      );  
  return;
}

void systBinary_pA_BinCount(TString deteco = "TOF",Bool_t enaBarlow = 0/*no Barlow check*/,Bool_t display=0)
{    
  //02 september 2013
  TString fPath = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/systematics",deteco.Data());  
  TString fCentral = "RAW_best_fixedW_poly2.root"; 
  TString fOther; 
  TString binsFile;
  if (deteco.Contains("TOF")) {
    fOther = "RAW_binCount_E4s_P2s.root";
    binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1818_tof2s_train215-216.root";
  }  else {
    binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TPC_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1717_tpc2s_train215-216.root";
    fOther = "RAW_best_binCount_E4s_P2s.root";
  }
  systBinaryFromHisto(fPath, 
		      fCentral, 
		      fOther, 
		      "_binCount_E4s_P2s.root",
		      -1, //all cents
		      EQuantity::kYields, 
		      ESystType::kRMS, 
		      5.0, 
		      binsFile.Data(),
		      enaBarlow, 
		      display);  
  return;
}

void systBinary_pA_freeW(TString deteco = "TOF",Bool_t enaBarlow = 0/*no Barlow check*/,Bool_t display=0)
{    
  //02 september 2013
  TString fPath = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/systematics",deteco.Data());  
  TString fCentral = "RAW_best_fixedW_poly2.root"; 
  TString fOther; 
  TString binsFile;
  if (deteco.Contains("TOF")) {
    fOther = "RAW_best_freeW08_poly2.root";
    binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1818_tof2s_train215-216.root";
  }  else {
    binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TPC_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1717_tpc2s_train215-216.root";
    fOther = "RAW_best_freeW08_poly2.root";
  }
  systBinaryFromHisto(fPath, 
		      fCentral, 
		      fOther, 
		      "_freeW08-15.root",
		      -1, //all cents
		      EQuantity::kYields, 
		      ESystType::kRMS, 
		      5.0, 
		      binsFile.Data(),
		      enaBarlow, 
		      display);  
  return;
}

void systBinary_pA_ResBg(TString deteco = "TOF",Bool_t enaBarlow = 0/*no Barlow check*/,Bool_t display=0)
{    
  //02 september 2013
  TString fPath = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/systematics",deteco.Data());  
  TString fCentral = "RAW_best_fixedW_poly2.root"; 
  TString fOther; 
  TString binsFile;
  if (deteco.Contains("TOF")) {
    fOther = "RAW_best_fixedW_poly1.root";
    binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1818_tof2s_train215-216.root";
  }  else {
    binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TPC_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1717_tpc2s_train215-216.root";
    fOther = "RAW_best_fixedW_poly3.root";
  }
  systBinaryFromHisto(fPath, 
		      fCentral, 
		      fOther, 
		      "_ResBg.root",
		      -1, //all cents
		      EQuantity::kYields, 
		      ESystType::kRMS, 
		      5.0, 
		      binsFile.Data(),
		      enaBarlow, 
		      display);  
  return;
}


void systBinary_pA_LSnorm(TString deteco = "TOF",Bool_t enaBarlow = 0/*no Barlow check*/,Bool_t display=0)
{    
  //02 september 2013
  TString fPath = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/systematics",deteco.Data());  
  TString fCentral = "RAW_best_fixedW_poly2.root"; 
  TString fOther = "RAW_lsNorm13_best_fixedW_poly2.root";
  TString binsFile;
  if (deteco.Contains("TOF"))
    binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1818_tof2s_train215-216.root";
  else 
    binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TPC_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1717_tpc2s_train215-216.root";
  
  systBinaryFromHisto(fPath, 
		      fCentral, 
		      fOther, 
		      "_LSnorm.root",
		      -1, //all cents
		      EQuantity::kYields, 
		      ESystType::kRMS, 
		      5.0, 
		      binsFile.Data(),
		      enaBarlow, 
		      display);  
  return;
}

//-----------------------------------------------------------------
void systBinary_pA_EM(TString deteco = "TOF",Bool_t enaBarlow = 0/*no Barlow check*/,Bool_t display=0)
{    
  //02 september 2013
  TString fPath = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/systematics",deteco.Data());  
  TString fCentral = "RAW_best_fixedW_poly2.root"; 
  TString fOther = "RAW_em_best_fixedW_poly2.root";
  TString binsFile;
  if (deteco.Contains("TOF"))
    binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1818_tof2s_train215-216.root";
  else 
    binsFile = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TPC_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1717_tpc2s_train215-216.root";
  
  systBinaryFromHisto(fPath, 
		      fCentral, 
		      fOther, 
		      "_em.root",
		      -1, //all cents
		      EQuantity::kYields, 
		      ESystType::kRMS, 
		      5.0, 
		      binsFile.Data(),
		      enaBarlow, 
		      display);  
  return;
}

//------------------------------------------------------------------
void systBinaryFromHisto(TString fPath = "normYields", 
			 TString fCentral="cent0_central_spectra.root", 
			 TString fOther="cent0_other_spectra.root",  
			 TString fSyst="kstar.root", 
			 Int_t centrality = -1, 
			 Int_t quantityID = EQuantity::kYields, 
			 Int_t strategy = ESystType::kRMS, 
			 Float_t chi2cut=5.0, 
			 TString filebins="cent0_central_spectra.root",
			 Bool_t enaBarlow = kTRUE,
			 Bool_t display=1,
			 Bool_t useCorrectedSpectra=kFALSE,
			 TString deteco = "")
{
  Color_t color[2][6]={kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack, //tpc
		       kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof
  Int_t marker[2][6]={21, 22, 23, 34, 33, 20, //tpc
		      25, 26, 32, 28, 27, 24}; //tof 
  
  Bool_t isTOF = (fPath.Contains("TOF") || fCentral.Contains("TOF") || fOther.Contains("TOF")) ;
  if (isTOF) deteco = "TOF";
  TString hRawYieldsName;
  if (useCorrectedSpectra) hRawYieldsName = Form("h%sCorrected_",deteco.Data());
  else hRawYieldsName = "hRawYieldVsPt_";
  TString hMassName("hMassVsPt_");
  TString hWidthName("hWidthVsPt_");
  
  //get bins
  Int_t npt_axis = 0, ncent_axis=0; 
  TFile *f=TFile::Open(filebins.Data());
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  TAxis *centbins = (TAxis*)f->Get("centbins");
  if (!ptbins || !centbins) return;
  npt_axis = ptbins->GetNbins();  
  ncent_axis = centbins->GetNbins();
  f->Close();
  const Int_t npt = npt_axis+1;
  const Int_t ncent = ncent_axis+1;
  Double_t pt[npt];
  for (Int_t k=0; k<npt;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
  }
  Double_t cent[ncent]; 
  for (Int_t k=0; k<ncent;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
  }  

  //define output
  gSystem->Exec("mkdir -p systUncert");
  TString outname(fSyst.Data());  
  if (centrality>-1) outname.Prepend(Form("cent%i",centrality));
  else outname.Prepend("allCents");
  if (quantityID == EQuantity::kYields) outname.Prepend("yields_");
  if (quantityID == EQuantity::kMass) outname.Prepend("mass_");
  if (quantityID == EQuantity::kWidth) outname.Prepend("width_");
  if (!enaBarlow) outname.Prepend("noB_");
  outname.Prepend("syst_");
  //  outname.Prepend(Form("syst%i_",strategy));
  outname.Prepend("systUncert/");
  
  TFile * outfile = new TFile(outname.Data(),"recreate"); 

  //Get central spectra from files 
  TFile * finCentral= TFile::Open(Form("%s/%s", fPath.Data(),fCentral.Data()),"read");
  if (!finCentral) {
    Printf("Something wrong in input file 1 (central). Check if it is there!");
    return;
  }
  
  //Get other histo's file
  TFile * finOther= TFile::Open(Form("%s/%s", fPath.Data(),fOther.Data()),"read");
  if (!finOther) {
    Printf("Something wrong in input file 2 (other). Check if it is there!");
      return;
  }
  
  //Get name according to quantity to be displayed
  TString histoName;
  if (quantityID == EQuantity::kYields) histoName = hRawYieldsName;
  if (quantityID == EQuantity::kMass) histoName = hMassName;
  if (quantityID == EQuantity::kWidth) histoName = hWidthName;
  
  //loop on centrality  
  for (Int_t selCent=0; selCent<ncent-1; selCent++) {
    
    //centrality selection if specified
    if ((centrality>-1) && (selCent!=centrality)) continue;
    TString centLabel=Form("%2.0f-%3.0f%%",cent[selCent],cent[selCent+1]);
 
    //Get histo with central values - skip to next cent bin if not found
    TH1F* hCentral = (TH1F*) finCentral->Get(Form("%s%i",histoName.Data(),selCent)); 
    if (!hCentral) {
      Printf("error: cannot get central points histo for cent %i", selCent); 
      continue;
    }
    
    TString centralTitle=fCentral;
    centralTitle.ReplaceAll(".root", ""); //centralTitle.ReplaceAll("_", " ");
    hCentral->SetName(Form("central%i",selCent));
    hCentral->SetTitle(centralTitle.Data());
    hCentral->SetMarkerStyle(20);
    hCentral->SetLineWidth(1);
    
    //Get other histo
    TH1F* hOther = (TH1F*) finOther->Get(Form("%s%i",histoName.Data(),selCent));
    if (!hOther) {
      Printf("error: cannot get other points histo for cent %i",selCent);
      continue;
    }
    
    TString secondTitle=fOther; 
    secondTitle.ReplaceAll(".root", ""); //secondTitle.ReplaceAll("_", " ");
    hOther->SetName(Form("other%i",selCent));
    hOther->SetTitle(secondTitle.Data());
    hOther->SetMarkerStyle(25);
    hOther->SetLineWidth(1);
    
    TH1F* ratio = (TH1F*) hOther->Clone(Form("ratioO2C_%i",selCent));
    ratio->Divide(hCentral);
    ratio->SetTitle("other/central");
    ratio->GetYaxis()->SetRangeUser(0.0,2.0);
    ratio->SetMarkerStyle(28);
    ratio->SetLineWidth(2);
    //save original plots in file
    outfile->cd();
    hCentral->Write();
    hOther->Write();
    ratio->Write();
  
    //define output
    TString titleQvsPt;
    TString titleQvsPt_yaxis;
    if (quantityID==EQuantity::kMass) {
      titleQvsPt = Form("Mass %s", centLabel.Data());
      titleQvsPt_yaxis = Form("M (GeV/c^{2})");
    }
    if (quantityID==EQuantity::kWidth) {
      titleQvsPt = Form("Width %s", centLabel.Data());
      titleQvsPt_yaxis = Form("#Gamma (GeV/c^{2})");
    }
    if (quantityID==EQuantity::kYields) {
      titleQvsPt = Form("Raw yields %s",centLabel.Data());
      titleQvsPt_yaxis = Form("1/N_{ev}*d^{2}N/dp_{t}dy (|y|<0.5)");
    }
        
    TH1F* hSystVsPt = new TH1F(Form("hSystVsPt_%i",selCent),"; p_{t} (GeV/c); syst. uncert.", npt_axis, pt);
    hSystVsPt->SetLineColor(color[isTOF][selCent]);	       
    hSystVsPt->SetMarkerColor(color[isTOF][selCent]);	       
    hSystVsPt->SetMarkerStyle(marker[isTOF][selCent]);	       
    hSystVsPt->SetLineWidth(1);
    
    TH1F* hStatVsPt = new TH1F(Form("hStatVsPt_%i",selCent),"; p_{t} (GeV/c); stat. uncert.", npt_axis, pt);
    hStatVsPt->SetLineColor(color[isTOF][selCent]);	       
    hStatVsPt->SetMarkerColor(color[isTOF][selCent]);	       
    hStatVsPt->SetMarkerStyle(marker[isTOF][selCent]);	       
    hStatVsPt->SetLineWidth(1);
    
    TH1F* hSystVsPtPercentageOfCentral = new TH1F(Form("hSystVsPtPercentageOfCentral_%i",selCent),"; p_{t} (GeV/c); syst. uncert. (% of central value)", npt_axis, pt);
    hSystVsPtPercentageOfCentral->SetTitle(Form("centrality %s", centLabel.Data()));
    hSystVsPtPercentageOfCentral->SetLineColor(color[isTOF][selCent]);	       
    hSystVsPtPercentageOfCentral->SetMarkerColor(color[isTOF][selCent]);	       
    hSystVsPtPercentageOfCentral->SetMarkerStyle(marker[isTOF][selCent]);
    hSystVsPtPercentageOfCentral->SetLineWidth(1);
    
    TH1F* hStatVsPtPercentageOfCentral = new TH1F(Form("hStatVsPtPercentageOfCentral_%i",selCent),"; p_{t} (GeV/c); stat. uncert. (% of central value)", npt_axis, pt);
    hStatVsPtPercentageOfCentral->SetTitle(Form("centrality %s", centLabel.Data()));
    hStatVsPtPercentageOfCentral->SetLineColor(color[isTOF][selCent]);	       
    hStatVsPtPercentageOfCentral->SetMarkerColor(color[isTOF][selCent]);	       
    hStatVsPtPercentageOfCentral->SetMarkerStyle(marker[isTOF][selCent]);
    hStatVsPtPercentageOfCentral->SetLineWidth(1);
    
    TH1F* hYieldVsPt = new TH1F(Form("hYieldVsPt_%i",selCent),"stat. uncert. only; p_{t} (GeV/c);", npt_axis, pt);
    hYieldVsPt->SetTitle(titleQvsPt.Data());
    hYieldVsPt->GetYaxis()->SetTitle(titleQvsPt_yaxis.Data());
    hYieldVsPt->SetMarkerColor(color[isTOF][selCent]+1);
    hYieldVsPt->SetLineColor(color[isTOF][selCent]+1);
    hYieldVsPt->SetMarkerStyle(1);
    hYieldVsPt->SetLineWidth(2);
    
    TH1F* hYieldVsPtWsyst = new TH1F(Form("hYieldVsPtWsyst_%i",selCent),"#sqrt{#sigma_{stat}^{2}+#sigma_{syst}^{2}}; p_{t} (GeV/c);", npt_axis, pt);
    hYieldVsPtWsyst->GetYaxis()->SetTitle(titleQvsPt_yaxis.Data());
    // hYieldVsPtWsyst->SetLineColor(color[isTOF][selCent]+2);	       
    hYieldVsPtWsyst->SetMarkerColor(color[isTOF][selCent]-9);	       
    // hYieldVsPtWsyst->SetMarkerStyle(1);
    hYieldVsPtWsyst->SetFillColor(color[isTOF][selCent]-9);
    hYieldVsPtWsyst->SetLineColor(color[isTOF][selCent]-9);
    hYieldVsPtWsyst->SetMarkerStyle(0);
    hYieldVsPtWsyst->SetOption("E2");
    
    //read input 
    Double_t yield1, yield2, erry1, erry2, syst;//fitParams1[9],fitParams2[9];
    Int_t centBinID1,ptBinID1,centBinID2,ptBinID2;
    
    Int_t nbinsx = hCentral->GetNbinsX()+1;
    for (Int_t ientry = 1;ientry<nbinsx;ientry++){
      
      //index 1 for central value
      yield1=hCentral->GetBinContent(ientry);
      erry1=hCentral->GetBinError(ientry);
      //index 2 for other value
      yield2=hOther->GetBinContent(ientry);
      erry2=hOther->GetBinError(ientry);
      
      //Get systematics according to specified strategy
      if (enaBarlow) {
	if ( !CheckIfStatisticallyCompatible(yield1, erry1, yield2, erry2, kTRUE)) {
	  syst = GetSystUncertainty(strategy,yield1, erry1, yield2, erry2);
	  if (strategy==ESystType::kRatio) 
	    syst = syst*yield1; //ratio = 1-m2/m1 => this is already a percentage!
	} else  {
	  Printf("Error compatible with stat - syst error negligible");
	  syst=1e-5;
	}
      } else {
	syst = GetSystUncertainty(strategy,yield1, erry1, yield2, erry2);
	if (strategy==ESystType::kRatio) 
	  syst = syst*yield1; //ratio = 1-m2/m1 => this is already a percentage!
      }
      
      Printf("Syst ucert for bin %i set to %6.4f\n=======================================", ientry, syst);
      // when fill histo add +2 because definition of bins start
      // norm for bin width already performed at spectra building level
      // Float_t dpt = ptbins->GetBinWidth(ientry);
      if (yield1>0) {
	hYieldVsPt->SetBinContent(ientry, yield1);
	hYieldVsPt->SetBinError(ientry, erry1);
	hYieldVsPtWsyst->SetBinContent(ientry, yield1);
	hYieldVsPtWsyst->SetBinError(ientry, TMath::Sqrt(erry1*erry1+syst*syst));
	hSystVsPt->SetBinContent(ientry, syst);
	hStatVsPt->SetBinContent(ientry, erry1);
	hSystVsPtPercentageOfCentral->SetBinContent(ientry, (syst*100)/yield1);
	hStatVsPtPercentageOfCentral->SetBinContent(ientry, (erry1*100)/yield1);
      }
    }//loop on bins
    
    if (display){
      gStyle->SetOptTitle(0);
      gStyle->SetOptStat(0);
      
      //superposition of the two versions of the spectra
      TCanvas *cspectra=new TCanvas(Form("cspectra_%i",selCent),"Raw yield vs p_{T}", 750,600);
      cspectra->cd();
      hCentral->GetYaxis()->SetRangeUser(1e-5, 15.);
      if (quantityID==EQuantity::kMass)     
	hCentral->GetYaxis()->SetRangeUser(0.87, 0.92);
      if (quantityID==EQuantity::kWidth)     
	hCentral->GetYaxis()->SetRangeUser(0.01, 0.200);
      if (quantityID==EQuantity::kYields) 
	hCentral->GetYaxis()->SetRangeUser(1e-5, 15.); 
      
      hCentral->GetYaxis()->SetTitleOffset(1.2);
      hCentral->Draw();
      hOther->Draw("same");
      if (quantityID==EQuantity::kYields) gPad->SetLogy();
      TLegend * autolegspectra = (TLegend*)gPad->BuildLegend(0.3,0.75,0.88,0.88, Form("Centrality %s",centLabel.Data()));
      autolegspectra->SetFillColor(kWhite);
      autolegspectra->SetLineColor(kWhite);
      autolegspectra->SetTextFont(42);

      TPaveText * labelkind = new TPaveText(0.20,0.91,0.88,0.98,"NDC");
      labelkind->InsertText(secondTitle.Data());

      //central spectrum with errors
      TCanvas *cry=new TCanvas(Form("cry_%i",selCent),"Raw yield vs p_{T} with systematic errors", 750,600);
      cry->cd(1);
      if (quantityID==EQuantity::kYields) 
	hYieldVsPtWsyst->GetYaxis()->SetRangeUser(1e-5, 15.); 
      if (quantityID==EQuantity::kMass)     
	hYieldVsPtWsyst->GetYaxis()->SetRangeUser(0.87, 0.92);
      if (quantityID==EQuantity::kWidth)     
	hYieldVsPtWsyst->GetYaxis()->SetRangeUser(0.01, 0.200);
      hYieldVsPtWsyst->GetYaxis()->SetTitleOffset(1.2);
      hYieldVsPtWsyst->Draw();
      hYieldVsPt->Draw("same");
      labelkind->Draw("same");
      gPad->SetLogy();
      TLegend * autolegry = (TLegend*)gPad->BuildLegend(0.6,0.7,0.88,0.88, Form("Centrality %s",centLabel.Data()));
      autolegry->SetFillColor(kWhite);
      autolegry->SetLineColor(kWhite);
      autolegry->SetTextFont(42);
    
      //systematic uncertainty plot
      TCanvas *cs=new TCanvas(Form("cs_%i",selCent),"Systematic error vs p_{T}", 750,600);
      // cs->Divide(1,2);
      // cs->cd(1);
      // frameSystVsPt->GetYaxis()->SetRangeUser(0.1, 1e6);
      // gPad->SetLogy();
      // frameSystVsPt->Draw();
      // hSystVsPt->Draw("same");
      // cs->cd(2);
      cs->cd();
      hSystVsPtPercentageOfCentral->SetLineWidth(2);
      hSystVsPtPercentageOfCentral->GetYaxis()->SetRangeUser(0.1, 100);
      hSystVsPtPercentageOfCentral->Draw();
      labelkind->Draw("same");
    }
    outfile->cd();
    if (display && (centrality>-1)) {
      gStyle->SetOptTitle(0);
      cry->Write();
      cspectra->Write();
      cs->Write();
      // outname.ReplaceAll(".root",".png");
      // cry->SaveAs(outname.Data());
      // cspectra->SaveAs(outname.Data());
      // cs->SaveAs(outname.Data());
    }
    hYieldVsPt->Write();
    hYieldVsPtWsyst->Write();
    hSystVsPt->Write();
    hStatVsPt->Write();
    hSystVsPtPercentageOfCentral->Write();
    hStatVsPtPercentageOfCentral->Write();
  }
  
  //if all cents together save summary canvas
  if (centrality<0) {
    gROOT->LoadMacro("$ASD/AddPaveText.C");
    TCanvas * centSummary = new TCanvas("centSummary","centSummary",550,600);
    TH1D* systtmp[5];
    for (Int_t ic=0; ic<ncent; ic++) { 
      outfile->cd();
      systtmp[ic] = (TH1D*) outfile->Get(Form("hSystVsPtPercentageOfCentral_%i",ic));
      systtmp[ic]->GetXaxis()->SetRangeUser((isTOF?1.0:0.5),8.0);
      systtmp[ic]->GetYaxis()->SetTitleOffset(1.4);
      centSummary->cd();
      if (ic>0) systtmp[ic]->Draw("same");
      else systtmp[ic]->Draw();
    }
    gStyle->SetOptTitle(0);
    TLegend * legSummary = (TLegend *) centSummary->BuildLegend();
    legSummary->SetFillColor(kWhite);
    legSummary->SetBorderSize(0);
    TLine * line = new TLine(1.0,10,8,10);
    line->SetLineStyle(7);
    line->SetLineWidth(1);
    line->SetLineColor(kBlack);
    centSummary->cd();
    line->Draw("same");
    AddPaveTextDetAnalysis(isTOF?"TOF":"TPC");
    //    if (outname.Contains(".root"))
    outname.ReplaceAll(".root",".png");
    centSummary->SaveAs(outname.Data());
  }
  return;  
  
}

//------------------------------------------------------------------------------------
Bool_t CheckIfStatisticallyCompatible(Double_t m1, Double_t e1, Double_t m2, Double_t e2, Bool_t verbose=0) 
{
  //Barlow method
  //to check if two measurements m1+/-e1 and m2+/-e2 are statistically compatible
  //Returns kFalse if not --> systematic uncertainty (or possible mistake!)
 
  Double_t y1,y2, s1,s2;
  if (e1>e2) {
    y1=m2; s1=e2;
    y2=m1; s2=e1;
  } else {
    y1=m1; s1=e1;
    y2=m2; s2=e2;
  }
  Double_t diffy          = TMath::Abs(y1-y2);
  Double_t deltaSigma2    = TMath::Sqrt(s2*s2-s1*s1);
  if (deltaSigma2==0.0) return kTRUE;
  Double_t numberOfDeltas = diffy / deltaSigma2;

  if (numberOfDeltas<=1.0) return kTRUE;
  if (verbose) Printf("deltaSigma2 = %e , diffy = %e , numberOfDeltas = %5.2f", TMath::Sqrt(deltaSigma2), diffy, numberOfDeltas);
  if (numberOfDeltas>=5.0) Printf("Barlow Check warning: number of deltas = %5.2f > 5.0 - YOU BETTER CHECK THE MEASUREMENT!", numberOfDeltas);
  return kFALSE;
}

//------------------------------------------------------------------------------------
Double_t GetSystUncertainty(Int_t strategy, Double_t m1, Double_t e1, Double_t m2, Double_t e2)
{
  Double_t uncert = 0.0;
  switch (strategy) {
    
  case ESystType::kHalfDiff :
    uncert = TMath::Abs(m2-m1)*0.5;
    break;
    
  case ESystType::kMaxDiff :
    uncert = TMath::Abs(m2-m1);
    break;
    
  case ESystType::kRatio :
    //m1 = the central value
    if (m1>0.0) uncert = TMath::Abs(1-m2/m1);
    else uncert = 0.0;
    break;
    
  case ESystType::kRMS :
    //m1 = the central value
    Double_t valueTmp4rms[2]={m1,m2};
    uncert = TMath::RMS(2, valueTmp4rms);
    break;
    
  case ESystType::kUniform :
    uncert = TMath::Abs(m2-m1)*TMath::Sqrt(12.0);
    break;
    
  default :
    Printf("Invalid systematic uncertainty strategy chosen - returning err_syst = |x-central|");
    uncert = TMath::Abs(m2-m1);
    break;
  }
  
  return uncert;
}

