//---------------------------------------------------------------------------------------
//            fbellini@cern.ch     12 nov 2012
//---------------------------------------------------------------------------------------
  enum EFitFunction{ kPOLY2,
		     kPOLY3,
		     kLandau,
		     kPOLY1,
		     kEXP,
		     kData};

  enum EQuantity{ kYields,
		  kMass,
		  kWidth};

void runSystAllCent(Char_t * inlist=NULL, 
		    Float_t chi2cut=2.0,
		    TString filebins="$HOME/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TPC_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1717_tpc2s_train215-216.root", 
		    TString suffix = "rangeWidth", 
		    Int_t quantityID=EQuantity::kYields)
{
  for (Int_t j=0;j<5;j++){
    systematics(inlist,
		j, 
		chi2cut,
		filebins.Data(),
		suffix.Data(), 
		quantityID);
  }
  return;
}
//---------------------------------------------------------------------------------------
void systematics( Char_t * inlist=NULL, 
		  Int_t selCent = 0, 
		  Float_t chi2cut=2.0, 
		  TString filebins="", 
		  TString suffix = "poly2", 
		  Int_t quantityID=EQuantity::kYields)
{
  
  if (quantityID==0) quantityID = EQuantity::kYields;

  Color_t color[]={kRed, kOrange, kGreen+2, kBlue, kMagenta, kBlack};
  Int_t marker[]={20, 21, 28, 22, 23};
  
  gROOT->LoadMacro("kstar/MakeSpectra.C");

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
  
  TString centLabel=Form("%3.0f-%3.0f%%",cent[selCent],cent[selCent+1]);
  
  //define output
  TString outname = Form("%s_cent%i_chiMax%2.1f.root",suffix.Data(),selCent,chi2cut);
  if (quantityID == EQuantity::kYields) outname.Prepend("yields_");
  if (quantityID == EQuantity::kMass) outname.Prepend("mass_");
  if (quantityID == EQuantity::kWidth) outname.Prepend("width_");
  outname.Prepend("syst_");//"%i_",strategy
  outname.Prepend("systUncert/");
  TFile * outfile = new TFile(outname.Data(),"recreate");
  
  //chain trees
  TString infile;
  Int_t ifiles=0; 
  FILE * files = fopen(inlist, "r") ; 
  TChain *tree = new TChain("tree");
  while ( infile.Gets(files) ){
    if (infile.Contains("SKIP")==0) {
      tree->Add(infile.Data());
      //Printf("Adding tree from %s",infile.Data());
      ifiles++;
      
      TH1F* spectrum = (TH1F*) MakeSpectraCent(quantityID, selCent, infile.Data(), filebins.Data(), chi2cut);
      outfile->cd();
      TString newname = Form("%s",infile.Data());
      newname.ReplaceAll("_aod049_kstar.root","");
      if (quantityID==EQuantity::kYields) newname.ReplaceAll("fit","rawYields_");
      if (quantityID==EQuantity::kMass) newname.ReplaceAll("fit","mass_");
      if (quantityID==EQuantity::kWidth) newname.ReplaceAll("fit","width_");
  
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
    }
  } 
  Printf("Number of files to be merged = %i\n",ifiles);
  
  //read tree from chain
  Double_t fitParams[9], SoverB=0.0, significance=0.0,normfactorCopy=0.0;
  Int_t centBinID,ptBinID, funcID;
  Float_t rangeInfCopy, rangeSupCopy, ptinfCopy, ptsupCopy, centinfCopy, centsupCopy, massRangeMin,massRangeMax;
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
  tree->SetBranchAddress("norm_inf", &rangeInfCopy);
  tree->SetBranchAddress("norm_sup", &rangeSupCopy);
  tree->SetBranchAddress("pt_inf", &ptinfCopy);
  tree->SetBranchAddress("pt_sup", &ptsupCopy);
  tree->SetBranchAddress("cent_inf", &centinfCopy);
  tree->SetBranchAddress("cent_sup", &centsupCopy);
  tree->SetBranchAddress("functionID", &funcID);
  tree->SetBranchAddress("fitrange_inf", &massRangeMin);
  tree->SetBranchAddress("fitrange_sup", &massRangeMax);

  //define titles
  TString titleHistoY, nameHisto, titleHisto;  

  if (quantityID==EQuantity::kYields){
    titleHistoY="dN/dp_{t}";
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
  
  //define histos for output
  TH1F*frameXVsPt=new TH1F("frameXVsPt","; p_{t} (GeV/c);", npt_axis, pt);
  frameXVsPt->SetLineColor(kWhite);
  frameXVsPt->SetMarkerColor(kWhite);
  frameXVsPt->SetTitle(titleHisto.Data());//"K^{*0}+#bar{K^{*0}} raw yield");
  frameXVsPt->GetYaxis()->SetTitle(titleHistoY.Data());
  
  TH1F*frameSystVsPt=new TH1F("frameSystVsPt","Systematic error vs. p_{t};p_{t} (GeV/c);syst error", npt_axis, pt);

  TH1F* hXVsPt = new TH1F(Form("%s_%i",nameHisto.Data(), selCent),";p_{t} (GeV/c);", npt_axis, pt);
  hXVsPt->SetTitle(titleHisto.Data());
  hXVsPt->GetYaxis()->SetTitle(titleHistoY.Data());
  hXVsPt->SetLineColor(color[selCent]+1);	       
  hXVsPt->SetMarkerColor(color[selCent]+1);	       
  hXVsPt->SetMarkerStyle(marker[selCent]);	       
  
  TH1F* hXVsPtWsyst = new TH1F(Form("%sWsyst_%i",nameHisto.Data(),selCent),"; p_{t} (GeV/c);", npt_axis, pt);
  hXVsPtWsyst->SetTitle(titleHisto.Data());
  hXVsPtWsyst->GetYaxis()->SetTitle(titleHistoY.Data());
  hXVsPtWsyst->SetLineColor(color[selCent]+2);	       
  hXVsPtWsyst->SetMarkerColor(color[selCent]+2);	       
  hXVsPtWsyst->SetMarkerStyle(marker[selCent]);
  
  TH1F* hSystVsPt = new TH1F(Form("hSystVsPt_%i",selCent),"systematic error; p_{t} (GeV/c); systematic error", npt_axis, pt);
  hSystVsPt->SetTitle(Form("syst. uncert. cent. %s", centLabel.Data()));
  hSystVsPt->SetLineColor(color[selCent]);	       
  hSystVsPt->SetMarkerColor(color[selCent]);	       
  hSystVsPt->SetMarkerStyle(marker[selCent]);	       
  
  TH1F* hSystVsPtPercentageOfCentral = new TH1F(Form("hSystVsPtPercentageOfCentral_%i",selCent),"systematic error Vs p_{t}; p_{t} (GeV/c); systematic error (% of central value)", npt_axis, pt);
  hSystVsPtPercentageOfCentral->SetTitle(Form("cent. %s", centLabel.Data()));
  hSystVsPtPercentageOfCentral->SetLineColor(color[selCent]);	       
  hSystVsPtPercentageOfCentral->SetMarkerColor(color[selCent]);	       
  hSystVsPtPercentageOfCentral->SetMarkerStyle(marker[selCent]);	       
  
  TH1F* hStatVsPt = new TH1F(Form("hStatVsPt_%i",selCent),"statistical error; p_{t} (GeV/c); statistical error", npt_axis, pt);
  hStatVsPt->SetTitle(Form("stat. uncert. cent. %s", centLabel.Data()));
  hStatVsPt->SetLineColor(color[selCent]);	       
  hStatVsPt->SetMarkerColor(color[selCent]);	       
  hStatVsPt->SetMarkerStyle(marker[selCent]);	       
  
  TH1F* hStatVsPtPercentageOfCentral = new TH1F(Form("hStatVsPtPercentageOfCentral_%i",selCent),"statistical error Vs p_{t}; p_{t} (GeV/c); statistical error (% of central value)", npt_axis, pt);
  hStatVsPtPercentageOfCentral->SetTitle(Form("cent. %s", centLabel.Data()));
  hStatVsPtPercentageOfCentral->SetLineColor(color[selCent]);	       
  hStatVsPtPercentageOfCentral->SetMarkerColor(color[selCent]);	       
  hStatVsPtPercentageOfCentral->SetMarkerStyle(marker[selCent]);	      
  
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
    if (!(centBinID==selCent)) continue;
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
      //Printf("#%i: cent=%i pt=%i index=%i --> value = %e +/- %e", ientry, centBinID, ptBinID, ifunc, yield[ipt][ifunc], yield_statErr[ipt][ifunc]);
    }//end check on chi2 
  }//end loop on tree entries
  
  for (Int_t ipt=1;ipt<npt-1;ipt++){  
    //divide by bin width
    Float_t dpt = ptbins->GetBinWidth(ipt+1);
    //get central value
    Double_t central_tmp[2]={-1.,-1};
    Int_t chosenValueID=-1;
    Printf(" \n \n °°°°°°°°°°°°°°°°°°°°°°°° pt bin %i °°°°°°°°°°°°°°°°°°°°°°°°° \n °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°dpt = %4.2f", ipt,dpt);
    chosenValueID=GetCentralValue(yield[ipt], yield_statErr[ipt], nused[ipt], central_tmp);
    Printf("------------------------------------------------ cent = %i, pt = %i @ file in list: %i", selCent, ipt, chosenValueID+1);
    central[ipt]=central_tmp[0]; 
    centralErr[ipt]=central_tmp[1]; 
    yield_systErr[ipt]=GetSystematicError(yield[ipt], yield_statErr[ipt], central[ipt], centralErr[ipt], nused[ipt],chosenValueID);

    //to be divided by dpt
    //fill yield histos, stat error only
    if (quantityID==EQuantity::kYields){
      hXVsPt->SetBinContent(ipt+1, central[ipt]/dpt);
      hXVsPt->SetBinError(ipt+1, centralErr[ipt]/dpt);
  
      //fill yield histos, syst error only
      hXVsPtWsyst->SetBinContent(ipt+1, central[ipt]/dpt);
      hXVsPtWsyst->SetBinError(ipt+1, yield_systErr[ipt]/dpt);
  
      //fill systematic error histos
      hSystVsPt->SetBinContent(ipt+1, yield_systErr[ipt]/dpt);
      if (central[ipt]>0)
	hSystVsPtPercentageOfCentral->SetBinContent(ipt+1, yield_systErr[ipt]*100/(dpt*central[ipt]));

    } else {
      //fill histos for mass or width (no dpt division)
      hXVsPt->SetBinContent(ipt+1, central[ipt]);
      hXVsPt->SetBinError(ipt+1, centralErr[ipt]);
      
      //fill yield histos, syst error only
      hXVsPtWsyst->SetBinContent(ipt+1, central[ipt]);
      hXVsPtWsyst->SetBinError(ipt+1, yield_systErr[ipt]);
  
      //fill systematic error histos
      hSystVsPt->SetBinContent(ipt+1, yield_systErr[ipt]);
      if (central[ipt]>0.0)
	hSystVsPtPercentageOfCentral->SetBinContent(ipt+1, yield_systErr[ipt]*100/(central[ipt]));

      //fill statistical error histos
      hStatVsPt->SetBinContent(ipt+1, centralErr[ipt]);
      if (central[ipt]>0.0)
	hStatVsPtPercentageOfCentral->SetBinContent(ipt+1, centralErr[ipt]*100/(central[ipt]));
    }
  }//loop on pt
  
  hXVsPt->SetMarkerColor(color[selCent]+1);
  hXVsPt->SetMarkerStyle(0);
  hXVsPt->SetFillStyle(0);
  hXVsPtWsyst->SetFillColor(color[selCent]-9);
  hXVsPtWsyst->SetLineColor(color[selCent]-9);
  hXVsPtWsyst->SetMarkerStyle(0);

  TCanvas *cry=new TCanvas("cry","Raw yield vs p_{T}", 1000,600);
  // cry->Divide(2,1);
  cry->cd();
  if (quantityID==EQuantity::kYields){
    frameXVsPt->GetYaxis()->SetRangeUser(9e2, 5e6);
    gPad->SetLogy();
  }
  if (quantityID==EQuantity::kMass){
    frameXVsPt->GetYaxis()->SetRangeUser(0.86,0.92);
  }
  if (quantityID==EQuantity::kWidth){
    frameXVsPt->GetYaxis()->SetRangeUser(0.01,0.1);
  }
  frameXVsPt->Draw();
  hXVsPtWsyst->Draw("E2 same");
  hXVsPt->Draw("same");
  TLegend*errleg = new TLegend(0.65,0.65,0.8,0.88, Form("centrality %s",centLabel.Data()));
  errleg->SetBorderSize(0);
  errleg->SetFillColor(kWhite);
  errleg->AddEntry(hXVsPt, "stat. uncert.","lp");
  errleg->AddEntry(hXVsPtWsyst, "syst. uncert.","f");
  errleg->Draw();
  // cry->cd(2);
  // frameXVsPt->GetYaxis()->SetRangeUser(1, 1e7);
  // gPad->SetLogy();
  //frameXVsPt->Draw();
  
  TCanvas *cs=new TCanvas("cs","Systematic error vs p_{T}", 800,600);
  cs->Divide(1,2);
  cs->cd(1);
  if (quantityID==EQuantity::kYields){
    frameSystVsPt->GetYaxis()->SetRangeUser(0.1, 1e6);
    gPad->SetLogy();
  } 
  frameSystVsPt->Draw();
  hSystVsPt->Draw("same");
  cs->cd(2);
  hSystVsPtPercentageOfCentral->GetYaxis()->SetRangeUser(0.1, 50.0);
  hSystVsPtPercentageOfCentral->Draw();
  
  outfile->cd();
  cry->Write();
  cs->Write();
  frameXVsPt->Write();
  hXVsPt->Write();
  frameSystVsPt->Write();
  hXVsPtWsyst->Write();
  hSystVsPt->Write();
  hSystVsPtPercentageOfCentral->Write();
  hStatVsPt->Write();
  hStatVsPtPercentageOfCentral->Write();
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
  Double_t syst = 0.0; 
  
  for (Int_t k=0;k<nfunc;k++){
    if ((var_pt_inArray[k]<=0.0)||(varStatErr_inArray[k]<=0.0))   continue;
    //copy into other array only good values and increase index
    //usedValues[iused]=var_pt_inArray[k];
    // if (k==centralFuncID){
    //   Printf("=== ifunc = %i is CENTRAL = %e +/- %e @@@ file %i",k, centralValue, centralValueErr, k+1);
    //   usedValues[k]=0.0; //set to zero even if value used for syst
    //   sigmas[k]=0.0; //set to zero if value not used for syst
    // } else {
	// Double_t delta = var_pt_inArray[k]-centralValue;
	// Double_t deltaerr = TMath::Sqrt(TMath::Abs(centralValueErr*centralValueErr-varStatErr_inArray[k]*varStatErr_inArray[k]));
	// //Printf("delta = %e \n deltaerr = %e", delta, deltaerr);
	// Double_t systcheck;
	// if (deltaerr>0.0) systcheck =TMath::Abs(delta/deltaerr);
	// else systcheck = 0.0;
	// if (systcheck>1.0) {
      //      usedValues[k]=systcheck;
      // 	sigmas[k]=deltaerr;
      // 	Printf("=== ifunc= %i --> systcheck = Delta/sigma_delta = %5.3f --> used for syst @ file %i",k,systcheck, k+1);
      // 	iused++;
      // } else {
      // 	Printf("=== ifunc= %i --> systcheck = Delta/sigma_delta = %5.3f --> NOT syst @ file %i",k,systcheck, k+1);
      // 	usedValues[k]=0.0; //set to zero even if value used for syst
      // }
    if ( (k!=centralFuncID) && 
	 (!IsStatisticallyCompatible(centralValue, centralValueErr, var_pt_inArray[k], varStatErr_inArray[k], 0 /*verbose*/)) ) {
      usedValues[iused]=var_pt_inArray[k];
      sigmas[iused]=varStatErr_inArray[k];
      Printf("=== ifunc = %i --> %e +/- %e is syst @ file %i", usedValues[iused], sigmas[iused], k, k+1);
      iused++;
    } else {
      if (k==centralFuncID) 
	Printf("=== ifunc = %i is CENTRAL = %e +/- %e @@@ file %i",k, centralValue, centralValueErr, k+1); 
      else 
	Printf("=== ifunc = %i --> NOT syst @ file %i",k, k+1);
    }
  }//loop func
  
  if (iused==0) return 0.0;
  else Printf("====================================\n Available values (!=central) for syst. uncert. = %i \n====================================", iused);
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
  if (deltaSigma2==0.0) return kTRUE;
  Double_t numberOfDeltas = diffy / deltaSigma2;
  if (numberOfDeltas<=1.0) return kTRUE;
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

