//---------------------------------------------------------------------------------------
//            fbellini@cern.ch     18 mar 2014
//---------------------------------------------------------------------------------------
  enum EFitFunction{ kPOLY2,
		     kPOLY3,
		     kLandau,
		     kPOLY1,
		     kEXP,
		     kData};

  enum EQuantity{ kYields=0,
		  kMass,
		  kWidth};

void systematics_pA( Char_t * inlist=NULL, 
		     Bool_t is0to100 = kTRUE,
		     Float_t chi2cut= 2.0, 
		     TString filebins="/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/proj_2424_tpc2s_tof3sveto.root", 
		     TString suffix = "poly2", 
		     Int_t quantityID=EQuantity::kYields,
		     TString myCentralFileName = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/tpc2s_tof3sveto/fitEM_norm1_BWpoly2_fixedW/best_fit_poly2.root")
		     //"/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s/fit_syst_poly2_fixedW_l1_r0/best_fit_poly2.root")
{
  Color_t color[2][6]={kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack, //tpc
		       kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof

  Color_t fillcolor[2][6]={kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack, //tpc
		       kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof


  Int_t marker[2][6]={21, 22, 23, 34, 33, 20, //tpc
		      25, 26, 32, 28, 27, 24}; //tof
  Bool_t isTOF=kFALSE;

  //get bins
  Int_t npt_axis = 0, ncent_axis=0; 
  TFile *f=TFile::Open(filebins.Data());
  if (!ptbins || !centbins) return;
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  npt_axis = ptbins->GetNbins();  
  TAxis *centbins = (TAxis*)f->Get("centbins");
  ncent_axis = centbins->GetNbins();
  f->Close();
  const Int_t npt = npt_axis+1;
  const Int_t ncent = ncent_axis+1;
  Double_t pt[npt];
  for (Int_t k=0; k<npt;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
    // Printf("k = %i pt[k]=%3.2f",k,pt[k]);
  }
  Double_t cent[ncent]; 
  for (Int_t k=0; k<ncent;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
    // Printf("k = %i cent[k]=%3.2f",k,cent[k]);
  }

 //define output
 TString outname = Form("%s_chiMax%2.1f.root",suffix.Data(),chi2cut);
 if (quantityID == EQuantity::kYields) outname.Prepend("yields_");
  if (quantityID == EQuantity::kMass) outname.Prepend("mass_");
  if (quantityID == EQuantity::kWidth) outname.Prepend("width_");
  outname.Prepend("syst_");//"%i_",strategy
  //  outname.Prepend("systUncert/");
  TFile * outfile = new TFile(outname.Data(),"recreate");
  
  //chain trees
  TString infile;
  Int_t ifiles=0; 
  FILE * files = fopen(inlist, "r") ; 
  TChain *tree = new TChain("tree");
  while ( infile.Gets(files) ){
    if (infile.Contains("TOF")) isTOF=kTRUE;
    if (infile.Contains("SKIP")==0) {
      tree->Add(infile.Data());
      Printf("Adding tree from %s",infile.Data());
      ifiles++;
#if 0      
      //gROOT->LoadMacro("$ASD/kstar/MakeSpectra.C"); --> copied in this macro instead
      TH1F* spectrum = (TH1F*) MakeSpectraCent(quantityID, selCent, infile.Data(), filebins.Data(), chi2cut);
      outfile->cd();
      //      TString newname = Form("%s",infile.Data());
      TString newname = Form("cent%i_alt%i",selCent,ifiles);
      if (quantityID==EQuantity::kYields) newname.Prepend("yields_");
      if (quantityID==EQuantity::kMass) newname.Prepend("mass_");
      if (quantityID==EQuantity::kWidth) newname.Prepend("width_");
      
      if (newname.Contains("POLY1")) spectrum->SetTitle(Form("BW+Poly1, %s",centLabel.Data()));
      if (newname.Contains("POLY2")) spectrum->SetTitle(Form("BW+Poly2, %s",centLabel.Data()));
      if (newname.Contains("POLY3")) spectrum->SetTitle(Form("BW+Poly3, %s",centLabel.Data()));
      if (newname.Contains("LAND")) spectrum->SetTitle(Form("BW+Landau, %s",centLabel.Data()));
      if (newname.Contains("BOLTZ")) spectrum->SetTitle(Form("RBWxPS+Poly2, %s",centLabel.Data()));
      if (newname.Contains("BWPS3")) spectrum->SetTitle(Form("RBWxPS+Poly3, %s",centLabel.Data()));
      if (newname.Contains("EXE")) spectrum->SetTitle(Form("RBWxPS+Exp, %s",centLabel.Data()));
      spectrum->SetName(newname.Data());
      spectrum->SetMarkerStyle(1);
      spectrum->Write();
#endif
    }
  } 
  Printf("Number of files to be merged = %i\n",ifiles);
  
  //read tree from chain
  Double_t fitParams[9], SoverB=0.0, significance=0.0,normfactorCopy=0.0;
  Int_t centBinID,ptBinID, funcID;
  Double_t rangeInfCopy, rangeSupCopy, ptinfCopy, ptsupCopy, centinfCopy, centsupCopy, massRangeMin,massRangeMax;
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
  tree->SetBranchAddress("fitrange_inf", &massRangeMin);
  tree->SetBranchAddress("fitrange_sup", &massRangeMax);

  for (Int_t selCent=0;selCent<5;selCent++){
    
    TString centLabel=Form("%i-%i%%",(is0to100 ? 0 : cent[selCent]),(is0to100 ? 100 : cent[selCent+1]));
    if (is0to100 && selCent>0) continue;
    
    //define titles
    TString titleHistoY, nameHisto, titleHisto;  
    
    if (quantityID==EQuantity::kYields){
      titleHistoY="dN/dp_{T}";
      nameHisto="hYieldVsPt" ;
      titleHisto=Form("Raw yield %s",centLabel.Data());
    }
    if (quantityID==EQuantity::kMass){
      titleHistoY="M (GeV/c^{2})";
      nameHisto="hMassVsPt";
      titleHisto=Form("BW mass %s",centLabel.Data());
    }
    if (quantityID==EQuantity::kWidth){
      titleHistoY="#Gamma (GeV/c^{2})";
      nameHisto="hWidthVsPt";
      titleHisto=Form("BW width %s",centLabel.Data());
    }
  
    TH1F* hXVsPt = new TH1F(Form("%s_%i",nameHisto.Data(), selCent),";p_{T} (GeV/c);", npt_axis, pt);
    hXVsPt->SetTitle(titleHisto.Data());
    hXVsPt->GetYaxis()->SetTitle(titleHistoY.Data());
    hXVsPt->SetLineColor(color[isTOF][selCent]);	       
    hXVsPt->SetMarkerColor(color[isTOF][selCent]);	       
    hXVsPt->SetMarkerStyle(marker[isTOF][selCent]);	       
  
    TH1F* hXVsPtWsyst = new TH1F(Form("%sWsyst_%i",nameHisto.Data(),selCent),"; p_{T} (GeV/c);", npt_axis, pt);
    hXVsPtWsyst->SetTitle(titleHisto.Data());
    hXVsPtWsyst->GetYaxis()->SetTitle(titleHistoY.Data());
    hXVsPtWsyst->SetLineColor(color[isTOF][selCent]);	       
    hXVsPtWsyst->SetMarkerColor(color[isTOF][selCent]);	       
    hXVsPtWsyst->SetMarkerStyle(marker[isTOF][selCent]);
  
    TH1F* hSystVsPt = new TH1F(Form("hSystVsPt_%i",selCent),"syst. uncert.; p_{T} (GeV/c); syst. uncert.", npt_axis, pt);
    hSystVsPt->SetTitle(Form("syst. uncert. cent. %s", centLabel.Data()));
    hSystVsPt->SetLineColor(color[isTOF][selCent]);	       
    hSystVsPt->SetMarkerColor(color[isTOF][selCent]);	       
    hSystVsPt->SetMarkerStyle(marker[isTOF][selCent]);	       
  
    TH1F* hSystVsPtPercentageOfCentral = new TH1F(Form("hSystVsPtPercentageOfCentral_%i",selCent),"relative syst. uncert. Vs p_{T}; p_{T} (GeV/c); relative syst. uncert.(%)", npt_axis, pt);
    hSystVsPtPercentageOfCentral->SetTitle(Form("cent. %s", centLabel.Data()));
    hSystVsPtPercentageOfCentral->SetLineColor(color[isTOF][selCent]);	       
    hSystVsPtPercentageOfCentral->SetMarkerColor(color[isTOF][selCent]);	       
    hSystVsPtPercentageOfCentral->SetMarkerStyle(marker[isTOF][selCent]);	       
  
    TH1F* hStatVsPt = new TH1F(Form("hStatVsPt_%i",selCent),"stat. uncert.; p_{T} (GeV/c); stat. uncert.", npt_axis, pt);
    hStatVsPt->SetTitle(Form("stat. uncert. cent. %s", centLabel.Data()));
    hStatVsPt->SetLineColor(color[isTOF][selCent]);	       
    hStatVsPt->SetMarkerColor(color[isTOF][selCent]);	       
    hStatVsPt->SetMarkerStyle(marker[isTOF][selCent]);	       
  
    TH1F* hStatVsPtPercentageOfCentral = new TH1F(Form("hStatVsPtPercentageOfCentral_%i",selCent),"relative stat. uncert.; p_{T} (GeV/c); relative stat. uncert. (%)", npt_axis, pt);
    hStatVsPtPercentageOfCentral->SetTitle(Form("cent. %s", centLabel.Data()));
    hStatVsPtPercentageOfCentral->SetLineColor(color[isTOF][selCent]);	       
    hStatVsPtPercentageOfCentral->SetMarkerColor(color[isTOF][selCent]);	       
    hStatVsPtPercentageOfCentral->SetMarkerStyle(marker[isTOF][selCent]);	      
  
    //define central values arrays
    Double_t central[npt],  centralErr[npt];
    const Int_t nfunc=ifiles;
    Double_t yield[npt][nfunc],  yield_statErr[npt][nfunc],  yield_systErr[npt];
    Int_t nused[npt]; 
  
    for (Int_t ipt=0;ipt<npt;ipt++){  
      central[ipt]=-1.;
      centralErr[ipt]=-1.;
      yield_systErr[ipt]=-1.;
      nused[ipt]=0;
      for (Int_t jj=0;jj<nfunc;jj++){    
	yield[ipt][jj]=-1.0;
	yield_statErr[ipt][jj]=-1.0;
      }
    }
  
    //----------------------------------------------
    // assign values from tree
    //----------------------------------------------
    for (Int_t ientry = 0; ientry < tree->GetEntries(); ientry++){
      tree->GetEntry(ientry);
      if ( centBinID!=100 && centBinID!=selCent) continue;
      if ( (fitParams[8]>0.0) && (fitParams[8]<chi2cut) ) {     
	Int_t ipt = ptBinID;
	Int_t ifunc = nused[ipt];
	if (quantityID==EQuantity::kMass){
	  yield[ipt][ifunc]=fitParams[0];
	  yield_statErr[ipt][ifunc]=fitParams[1];  
	} else {
	  if (quantityID==EQuantity::kWidth){
	    yield[ipt][ifunc]=fitParams[2];
	    yield_statErr[ipt][ifunc]=fitParams[3];  
	  } else {
	    if (quantityID==EQuantity::kYields){
	      yield[ipt][ifunc]=fitParams[4];
	      yield_statErr[ipt][ifunc]=fitParams[5];  
	    }
	  }
	}
	nused[ipt]++;
      }//end check on chi2 
    }//end loop on tree entries
  
    for (Int_t ipt=0;ipt<npt-1;ipt++){  
      //divide by bin width
      Float_t dpt = ptbins->GetBinWidth(ipt+1);
    
      //get central value
      Double_t central_tmp[2]={-1.,-1}; 
      Int_t chosenValueID=-1;
      Printf(" \n \n °°°°°°°°°°°°°°°°°°°°°°°° \n°°°°°°°°°°°°°°°°°°°°°°°°°\n pt bin %i, dpt = %4.2f \n °°°°°°°°°°°°°°°°°°°°°°°°°", ipt,dpt);

      if (myCentralFileName.Contains(".root")) {
	//if enabled, use externally set central value for relative errors evaluation 
	GetCentralValueFromFile(myCentralFileName.Data(),quantityID, (is0to100? 100 : selCent), ipt, central_tmp);
	central[ipt]=central_tmp[0]; 
	centralErr[ipt]=central_tmp[1]; 
	Printf("::::: Central value set from external file:\n %s \n::::: Central value = %e +/- %e", myCentralFileName.Data(), central[ipt], centralErr[ipt] );
      } else {
	chosenValueID=GetCentralValue(yield[ipt], yield_statErr[ipt], nused[ipt], central_tmp);
	central[ipt]=central_tmp[0]; 
	centralErr[ipt]=central_tmp[1]; 
	Printf(":::: Cent = %i, pt = %i @ file in list: %i", (is0to100? 100 : selCent), ipt, chosenValueID+1);
	Printf(":::: Central value = %e +/- %e", central[ipt], centralErr[ipt]);
      }
      
      //get syst. uncert.
      yield_systErr[ipt] = GetSystRMS(yield[ipt], yield_statErr[ipt], central[ipt], centralErr[ipt], nused[ipt], chosenValueID);
      
      Printf(":::: Syst. (RMS) = %e ", yield_systErr[ipt]);
      //fill histos for mass or width (no dpt division)
      hXVsPt->SetBinContent(ipt+1, central[ipt]);
      hXVsPt->SetBinError(ipt+1, centralErr[ipt]);
      
      //fill yield histos, syst error only
      hXVsPtWsyst->SetBinContent(ipt+1, central[ipt]);
      hXVsPtWsyst->SetBinError(ipt+1, yield_systErr[ipt]);
      
      //fill systematic error histos
      hSystVsPt->SetBinContent(ipt+1, yield_systErr[ipt]);
      if (central[ipt]>0.0)
	hSystVsPtPercentageOfCentral->SetBinContent(ipt+1, yield_systErr[ipt]*100./(central[ipt]));
      else hSystVsPtPercentageOfCentral->SetBinContent(ipt+1, 0.0);
      Printf(":::: Relative syst. (RMS) = %4.2f %% ", yield_systErr[ipt]*100./(central[ipt]));
      
      //fill statistical error histos
      hStatVsPt->SetBinContent(ipt+1, centralErr[ipt]);
      if (central[ipt]>0.0)
	hStatVsPtPercentageOfCentral->SetBinContent(ipt+1, centralErr[ipt]*100./(central[ipt]));
      else hStatVsPtPercentageOfCentral->SetBinContent(ipt+1, 0.0);
    }//loop on pt
    
    hXVsPt->SetMarkerColor(color[isTOF][selCent]+1);
    hXVsPt->SetMarkerStyle(0);
    hXVsPt->SetFillStyle(0);
    hXVsPtWsyst->SetFillColor(fillcolor[isTOF][selCent]);
    hXVsPtWsyst->SetLineColor(fillcolor[isTOF][selCent]);
    hXVsPtWsyst->SetMarkerStyle(0);
  
    TCanvas *cry=new TCanvas("cry","Raw yield vs p_{T}", 600,700);
    // cry->Divide(2,1);
    cry->cd();
    hXVsPtWsyst->Draw("E2");
    hXVsPt->Draw("same");

    TLegend*errleg = new TLegend(0.65,0.65,0.8,0.88, Form("centrality %s",centLabel.Data()));
    errleg->SetBorderSize(0);
    errleg->SetFillColor(kWhite);
    errleg->AddEntry(hXVsPt, "stat. uncert.","lp");
    errleg->AddEntry(hXVsPtWsyst, "syst. uncert.","f");
    errleg->Draw();
  
    TCanvas *cs=new TCanvas("cs","Systematic error vs p_{T}", 500,600);
    cs->Divide(1,2);
    cs->cd(1);
    hSystVsPt->Draw();
    cs->cd(2);
    hStatVsPtPercentageOfCentral->SetFillColor(kGray);
    hStatVsPtPercentageOfCentral->GetYaxis()->SetRangeUser(0.1, 100.0);
    hStatVsPtPercentageOfCentral->Draw();
    hSystVsPtPercentageOfCentral->SetLineWidth(2);
    hSystVsPtPercentageOfCentral->GetYaxis()->SetRangeUser(0.1, 100.0);
    hSystVsPtPercentageOfCentral->Draw();
  
    outfile->cd();
    cry->Write();
    cs->Write();
    hXVsPt->Write();
    hXVsPtWsyst->Write();
    hSystVsPt->Write();
    hSystVsPtPercentageOfCentral->Write();
    hStatVsPt->Write();
    hStatVsPtPercentageOfCentral->Write();
  }
  return;
}

//------------------------------------------------------------------------------------
Int_t GetCentralValue(Double_t *var_inArray, Double_t *varStatErr_inArray, Int_t nfunc, Double_t *centralValue)
{
  if (!var_inArray || (nfunc<1)) {
    Printf("GetCentralValue:: invalid arguments, doing nothing!");
    return;
  }
  Int_t chosenValueID=-1;
  Int_t usedFunc = nfunc;
  Double_t weightedMean = 0.0, weight = 0.0, deltaValue = 1e20;
  for (Int_t k=0;k<nfunc;k++){
    if ( (var_inArray[k]<=0.0) || (varStatErr_inArray[k]<=0.0) ) {     //skip values that do not pass chi2 cut
      usedFunc--;
    } else {
      //Printf("---- func=%i     value=%e    err=%e    ", k,var_inArray[k],varStatErr_inArray[k]);
      Double_t w2=varStatErr_inArray[k]*varStatErr_inArray[k];
      weightedMean+= var_inArray[k]/w2;
      weight+= 1./w2;
    }
  }//loop func
  
  if (weight>0){
    weightedMean/=weight;
  } else {
    return 0.0;
  } 
  
  centralValue[0] = weightedMean;    
  //Printf("weighted mean = %e used func=%i",weightedMean, usedFunc);
  
  for (Int_t k=0;k<nfunc;k++){
    if ( (var_inArray[k]<=0.0) || (varStatErr_inArray[k]<=0.0) ) continue;     //skip values that do not pass chi2 cut
    Double_t diff = TMath::Abs(centralValue[0]-var_inArray[k]);
    //Printf("var[%i]=%e    central=%e    diff=%e  delta=%e",k, var_inArray[k],centralValue[0], diff, deltaValue );
    if ( diff < deltaValue) {
      deltaValue = diff;
      centralValue[0]=var_inArray[k];
      centralValue[1]=varStatErr_inArray[k];
      chosenValueID=k; 
    }
  }
  Printf("GetCentralValue:: Used functions =%i central value = %e +/- %e", usedFunc, centralValue[0], centralValue[1] );
  return chosenValueID;
}

//-----------------------------------------------------------
Int_t GetCentralValueFromFile(TString myCentralFileName, Int_t quantityID=EQuantity::kYields, Int_t selCent=-1, Int_t selPt=-1, Double_t *myCentralValue)
{
  /* 
     get central value for a single bin from tree in filename, specified by the user
  */
  if (selCent<0 || selPt<0) return -1;

  TFile * myCentralFile = 0x0; 
  if (myCentralFileName) myCentralFile = TFile::Open(myCentralFileName.Data());
  if (!myCentralFile) return -1;
  
  TTree *ctree = (TTree *)myCentralFile->Get("tree");
  Double_t fitParams[9];
  Int_t centBinID,ptBinID;
  ctree->SetBranchAddress("signalMass",&fitParams[0]);
  ctree->SetBranchAddress("signalMassErr",&fitParams[1]);
  ctree->SetBranchAddress("signalWidth",&fitParams[2]);
  ctree->SetBranchAddress("signalWidthErr",&fitParams[3]);
  ctree->SetBranchAddress("nSignal",&fitParams[4]);
  ctree->SetBranchAddress("nSignalErr",&fitParams[5]);
  ctree->SetBranchAddress("nBack",&fitParams[6]);
  ctree->SetBranchAddress("nBackErr",&fitParams[7]);
  ctree->SetBranchAddress("chi2",&fitParams[8]);  
  ctree->SetBranchAddress("centBin",&centBinID);
  ctree->SetBranchAddress("ptBin",&ptBinID);
  
  for (Int_t ientry = 0; ientry < ctree->GetEntries(); ientry++){
    ctree->GetEntry(ientry);
    if ((centBinID!=selCent) || (ptBinID!=selPt)) continue;
    if (quantityID==EQuantity::kMass){
      myCentralValue[0]=fitParams[0];
      myCentralValue[1]=fitParams[1];  
    } else {
      if (quantityID==EQuantity::kWidth){
	myCentralValue[0]=fitParams[2];
	myCentralValue[1]=fitParams[3];  
      } else {
	if (quantityID==EQuantity::kYields){
	  myCentralValue[0]=fitParams[4];
	  myCentralValue[1]=fitParams[5];  
	}
      }
   }
  }//end loop on tree entries
  return 1;
}

//------------------------------------------------------------------------------------
Double_t GetSystRMS(Double_t *var_pt_inArray, Double_t *varStatErr_inArray, Double_t centralValue =0.0, Double_t centralValueErr=0.0, Int_t nfunc=0, Int_t centralFuncID=-1)
{
  /*
    Get syst error as the RMS
  */
  
  if (!var_pt_inArray || !varStatErr_inArray) {
    Printf("GetSystematicError:: invalid arguments, doing nothing!");
    return -1e10;
  }
  if ( (centralValue<=0.0) || (nfunc<1) ){
    Printf("GetSystematicError:: Not enough measurements to estimate systematics, doing nothing");
    return 0.0;
  }
  
  Int_t usedFunc = 0;
  const Int_t nfunctmp =nfunc;
  Bool_t keep[nfunctmp];
  for (Int_t k=0;k<nfunc;k++){
    keep[k]=kFALSE;
    if ((var_pt_inArray[k]<=0.0)||(varStatErr_inArray[k]<=0.0)) continue;
    //if ((k==centralFuncID) || (!IsStatisticallyCompatible(centralValue, centralValueErr, var_pt_inArray[k], varStatErr_inArray[k], 0 /*verbose*/))){
      keep[k]=1;
      usedFunc++;
      //}
  }
  Printf("GetSystematicError:: Systematic from %i points, including central value %e +/- %e \n", usedFunc, centralValue, centralValueErr);
  
  /* RMS from Root */
  const Int_t dim=usedFunc;
  Int_t iused=0;
  //save values to be used for rsm in temp array (includes all alternative + central)
  Double_t valueTmp4rms[dim]; //+1 to include the central value 
  Double_t errorsTmp4rms[dim];
  for (Int_t vv=0;vv<nfunc;vv++){
    if (keep[vv]){ 
      valueTmp4rms[iused]=var_pt_inArray[vv];
      errorsTmp4rms[iused]=varStatErr_inArray[vv];   
      iused++;
    }
  }
  /*Double_t scarto=0.0;
    Double_t mean=TMath::Mean(dim, valueTmp4rms, errorsTmp4rms);
    for (Int_t i=0;i<iused;i++){
    scarto += TMath::Power((valueTmp4rms[i]-mean),2);
    }
    Double_t rms2 = scarto/dim;
    Double_t rms = TMath::Sqrt(rms2);
  */
  //  Printf("DIM = %i   mean = %f rms = %f", dim, mean, rms);
  if (iused<=1) return 0.0;
  Double_t rms = TMath::RMS(iused, valueTmp4rms);
  return rms;
}


//------------------------------------------------------------------------------------
Double_t  GetSystematicError(Double_t *var_pt_inArray, Double_t *varStatErr_inArray, Double_t centralValue =0.0, Double_t centralValueErr=0.0, Int_t nfunc=0, Int_t centralFuncID=-1)
{
  /*
    The sistematic error is taken as the difference between the central and extreme value,
    divided by sqrt(12) (assuming uniform distribution).
    Alternatively it can be defined as the sigma of the gaussian distribution of the difference (value-central_value) if the number of available measurements is >10 (not yet implemented). 
  */
  
  if (!var_pt_inArray || !varStatErr_inArray) {
    Printf("GetSystematicError:: invalid arguments, doing nothing!");
    return;
  }
  if ( (centralValue<=0.0) || (nfunc<1) ){
    Printf("GetSystematicError:: Not enough measurements to estimate systematics, doing nothing");
  }
  
  Int_t usedFunc = 0;
  Double_t minDiff = 1e20, maxDiff = 0.0;
  // TH1D* hDiff = new TH1D("hDiff","hDiff", 2e2, -1e4, 1e4);
  for (Int_t k=0;k<nfunc;k++){
    if ((var_pt_inArray[k]<=0.0)||(varStatErr_inArray[k]<=0.0))   continue;
    usedFunc++;
  }
  
  Printf("GetSystematicError:: Systematic from %i points, for central value %e +/- %e \n",usedFunc-1, centralValue, centralValueErr);
  const Int_t nvalues=usedFunc;//exclude the central from the values
  Double_t usedValues[nvalues], sigmas[nvalues];
  Int_t iused=0;
  
  for (Int_t k=0;k<nfunc;k++){
    //check valid values
    if ((var_pt_inArray[k]<=0.0)||(varStatErr_inArray[k]<=0.0)) continue;
    
    //count the used values besides the central one and exclude those that are statistically compatible
    if ( (k!=centralFuncID) && 
	 (!IsStatisticallyCompatible(centralValue, centralValueErr, var_pt_inArray[k], varStatErr_inArray[k], 0 /*verbose*/)) )
      {
	Printf("=== ifunc = %i --> %e +/- %e is syst @ file %i", usedValues[iused], sigmas[iused], k, k+1);
	usedValues[iused]=var_pt_inArray[k];
	sigmas[iused]=varStatErr_inArray[k];
	iused++;
      } else {
      if (k==centralFuncID) 
	Printf("=== ifunc = %i is CENTRAL = %e +/- %e @@@ file %i",k, centralValue, centralValueErr, k+1); 
      else 
	Printf("=== ifunc = %i --> NOT syst @ file %i",k, k+1);
    }
  }//loop func
  
  if (iused==0) return 0.0;
  else Printf(":::: Available values (+central) for syst. uncert. = %i", iused);
  
  /* RMS from Root */
  const Int_t dim=iused+1;
  //save values to be used for rsm in temp array (includes all alternative + central)
  Double_t valueTmp4rms[dim]; //+1 to include the central value 
  Double_t errorsTmp4rms[dim];
  valueTmp4rms[iused]=centralValue;
  errorsTmp4rms[iused]=centralValueErr;
  for (Int_t vv=0;vv<iused+1;vv++){
    valueTmp4rms[vv]=usedValues[vv];
    errorsTmp4rms[vv]=sigmas[vv];   
  }
  Double_t scarto=0.0;
  Double_t mean=TMath::Mean(dim, valueTmp4rms, errorsTmp4rms);
  for (Int_t i=0;i<iused;i++){
    scarto += TMath::Power((valueTmp4rms[i]-mean),2);
    //rms = TMath::RMS(iused,usedValues[i]);
  }
  Double_t rms2 = scarto/dim;
  Double_t rms = TMath::Sqrt(rms2);
  Printf("DIM = %i   mean = %f rms = %f", dim, mean, rms);
  syst = rms;  

  //Double_t syst=TMath::RMS(dim,valueTmp4rms);
  /*/Uncomment if want ot use uniform distrib
  Double_t max = -1.e10;
  Int_t deltaerrID=-1;
  for (Int_t i=0;i<iused;i++){
    if (TMath::Abs(usedValues[i]-centralValue)>max){
      max=TMath::Abs(usedValues[i]-centralValue);
      deltaerrID=i;
    } 
    Printf("usedValue = %e max = %e",usedValues[i], max);
  }
  //syst=(max-1)*deltaerr[deltaerrID];
  //Printf("GetSystematicErrors:: nfunc <= %i points, syst = (N-1)*Delta = %e",iused-1,syst);  
  //below: syst = Max(y_i-central)/sqrt(12)   
  syst = max/TMath::Sqrt(12);
  Printf("GetSystematicErrors:: nfunc = %i points, max = %e,  syst = |max-central|/sqrt(12) = %e", iused, max, syst);  
  */
  return syst;
}


//------------------------------------------------------------------------------------
Bool_t IsStatisticallyCompatible(Double_t m1, Double_t e1, Double_t m2, Double_t e2, Bool_t verbose=0) 
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
  if (deltaSigma2==0.0) {
    Printf("Barlow check: deltaSigma = 0.0");
    return kTRUE;
  }
  Double_t numberOfDeltas = diffy / deltaSigma2;
  if (numberOfDeltas<=1.0) {
    Printf("Barlow check: measurements are statistically compatible");
    return kTRUE;
  }
  if (verbose) Printf("deltaSigma2 = %e , diffy = %e , numberOfDeltas = %5.2f", TMath::Sqrt(deltaSigma2), diffy, numberOfDeltas);
  if (numberOfDeltas>=5.0) Printf("Barlow Check warning: number of deltas = %5.2f > 5.0 - YOU BETTER CHECK THE MEASUREMENT!", numberOfDeltas);
  return kFALSE;
}

//------------------------------------------------------------------------------------
Double_t syst(Int_t strategy, Double_t m1, Double_t e1, Double_t m2, Double_t e2)
{
  //several ways for computing syst uncert - not used
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


//-------------------------------------------------------------------------
TH1F* MakeSpectraCent(Int_t quantityID=EQuantity::kYields,
		      Int_t icentbin = 0, 
		      TString filein="best_fit_poly2.root", 
		      Char_t* fileproj="/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/proj_2424_tpc2s_tof3sveto.root", 
		      Float_t cutChi2=10., 
		      Bool_t divideByDpt = kFALSE)
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
