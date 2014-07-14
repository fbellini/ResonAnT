void kstarSystematics(Char_t * outfile=NULL, Char_t * list=NULL, TString syst_type="FRN", Float_t chi2cut=50.)
{
  
  Int_t npt_axis = 0, ncent_axis=0; 
  TFile *f=TFile::Open("/Users/bellini/alice/resonances/myKstar/pwglf_train_out/data/good49_syst/sub_EMnorm1.30-1.50_aod049_kstar.root");
  if (!ptbins || !centbins) return;
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  npt_axis = ptbins->GetNbins();  
  TAxis *centbins = (TAxis*)f->Get("centbins");
  ncent_axis = centbins->GetNbins();
  f->Close();
  const Int_t npt = npt_axis+1, ncent = ncent_axis+1;
  Double_t pt[npt];
  for (Int_t k=0; k<npt;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
  }
  Double_t cent[ncent]; 
  for (Int_t k=0; k<ncent;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
  }
  Printf("npt = %i",npt);
 const Int_t nfunc=3;

  //chain trees
  TString infile;
  Int_t ifiles=0;
  //get bins
  
  FILE * files = fopen(list, "r") ; 
  TChain *tree = new TChain("tree");
  while ( infile.Gets(files) ){
    tree->Add(infile.Data());
    //Printf("Adding tree from %s",infile.Data());
    ifiles++;
  } 
  Printf("Number of files to be merged = %i\n",ifiles);
  
  //read tree from chain
  Double_t fitParams[9], SoverB=0.0, significance=0.0,normfactorCopy=0.0;
  Int_t centBinID,ptBinID, funcID;
  Float_t rangeInfCopy, rangeSupCopy, ptinfCopy, ptsupCopy, centinfCopy, centsupCopy, massRangeMin,massRangeMax;
  
TH1F*frameYieldVsPt=new TH1F("YieldVsPt","Raw yield vs. p_{t}; p_{t} (GeV/c); dN/dp_{t}", npt_axis, pt);
  TH1F*frameSystVsPt=new TH1F("SystVsPt","Systematic error vs. p_{t}; p_{t} (GeV/c); syst error", npt_axis, pt);
  TH1F* hYieldVsPt[ncent];
  TH1F* hSystVsPt[ncent];
 
  for (Int_t icentbin=0;icentbin<ncent;icentbin++){
    hYieldVsPt[icentbin] = new TH1F(Form("hYieldVsPt_%i",icentbin),"Yield Vs p_{t}; p_{t} (GeV/c); M (GeV/c^{2})", npt_axis, pt);
    hSystVsPt[icentbin] = new TH1F(Form("hSystVsPt_%i",icentbin),"systematic error Vs p_{t}; p_{t} (GeV/c); systematic error", npt_axis, pt);
  }

  //TString*fitfunction;
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
  //  tree->SetBranchAddress("function", &fitfunction);
  tree->SetBranchAddress("functionID", &funcID);
  tree->SetBranchAddress("fitrange_inf", &massRangeMin);
  tree->SetBranchAddress("fitrange_sup", &massRangeMax);
 
  Double_t yield_cent0[ncent][npt][nfunc], yieldStatErr_cent0[ncent][npt][nfunc], yield_central_cent0[ncent][npt], yield_centralErr_cent0[ncent][npt], yield_systErr_cent0[ncent][npt];
  Double_t yield_cent1[npt][nfunc], yieldStatErr_cent1[npt][nfunc], yield_central_cent1[npt], yield_centralErr_cent1[npt], yield_systErr_cent1[npt];
  Double_t yield_cent2[npt][nfunc], yieldStatErr_cent2[npt][nfunc], yield_central_cent2[npt], yield_centralErr_cent2[npt], yield_systErr_cent2[npt];
  Double_t yield_cent3[npt][nfunc], yieldStatErr_cent3[npt][nfunc], yield_central_cent3[npt], yield_centralErr_cent3[npt], yield_systErr_cent3[npt];
  Double_t yield_cent4[npt][nfunc], yieldStatErr_cent4[npt][nfunc], yield_central_cent4[npt], yield_centralErr_cent4[npt], yield_systErr_cent4[npt];

  for (Int_t ipt=0;ipt<npt;ipt++){  
    for (Int_t k=0;k<nfunc;k++){
      for (Int_t icent=0;icent<ncent;icent++){
	yield_cent0[icent][ipt][k] = 0.0;
	yield_central_cent0[icent][ipt]=0.0;
	yield_systErr_cent0[icent][ipt]=0.0;
      }
      yield_cent1[ipt][k]= 0.0;
      yield_central_cent1[ipt]=0.0;
      yield_systErr_cent1[ipt]=0.0;
      
      yield_cent2[ipt][k] = 0.0;
      yield_central_cent2[ipt]=0.0;
      yield_systErr_cent2[ipt]=0.0;
      
      yield_cent3[ipt][k] = 0.0;
      yield_central_cent3[ipt]=0.0;
      yield_systErr_cent3[ipt]=0.0;
      
      yield_cent4[ipt][k] = 0.0;
      yield_central_cent4[ipt]=0.0;
      yield_systErr_cent4[ipt]=0.0;
      //      yieldStatErr_cent0[ipt][k] = 0.0;
      yieldStatErr_cent1[ipt][k] = 0.0;
      yieldStatErr_cent2[ipt][k] = 0.0;
      yieldStatErr_cent3[ipt][k] = 0.0;
      yieldStatErr_cent4[ipt][k] = 0.0;  
      //yield_centralErr_cent0[ipt][k] = 0.0;
      yield_centralErr_cent1[ipt]= 0.0;
      yield_centralErr_cent2[ipt] = 0.0;
      yield_centralErr_cent3[ipt] = 0.0;
      yield_centralErr_cent4[ipt] = 0.0;  
    }//loop func
  }
  
  //sistematics for fit function
  TString funcLabel;
  Printf("Entries: %i",tree->GetEntries());
  for (Int_t ientry = 0; ientry < tree->GetEntries(); ientry++){
    tree->GetEntry(ientry);

    //set desired index
    Int_t index = -1;
    Int_t index_fitFunc = funcID;
    Int_t index_fitMinArr=-1; Float_t fitMinArr[4] = {0.70, 0.74, 0.78, 0.82};  
    Int_t index_fitMaxArr=-1; Float_t fitMaxArr[3] = {1.14, 1.18, 1.20};  
    Int_t index_normMinArr=-1; Float_t normMinArr[3] = {1.30, 1.35, 1.40};  
    
    for (Int_t j=0;j<4;j++){
      if (fitMinArr[j]==massRangeMin) index_fitMinArr=j;
      if (j<3){
	if (fitMaxArr[j]==massRangeMax) index_fitMaxArr=j;
	if (normMinArr[j]==rangeInfCopy) index_normMinArr=j;
      }
    }
    
    if (syst_type.Data()=="R") index = index_fitMaxArr+index_fitMinArr*3;
    if (syst_type.Data()=="N") index = index_normMinArr;
    if (syst_type.Data()=="F") index = index_fitFunc;
    if ( syst_type.Contains("F") && syst_type.Contains("R") && syst_type.Contains("N")) index = tree->GetEntries()/ifiles;

    //assign values of the variable only if the cut on chi2 is passed, for different centralities
    //it remains 0.0 if the chi2 does not pass the cut 
    if ((fitParams[8]<chi2cut) ) {
      yield_cent0[centBinID][ptBinID][index]=fitParams[4];   
      yieldStatErr_cent0[centBinID][ptBinID][k]=fitParams[5];
    }    
    if ( (centBinID==1) && (fitParams[8]<chi2cut) ) {
      yield_cent1[ptBinID][index]=fitParams[4];
      yieldStatErr_cent1[ptBinID][k]=fitParams[5];   
    }   
    if ( (centBinID==2) && (fitParams[8]<chi2cut) ) {
      yield_cent2[ptBinID][index]=fitParams[4];   
      yieldStatErr_cent2[ptBinID][k]=fitParams[5];
    }
    if ( (centBinID==3) && (fitParams[8]<chi2cut) ) {
      yield_cent3[ptBinID][index]=fitParams[4];   
      yieldStatErr_cent3[ptBinID][k]=fitParams[5];
    }
    if ( (centBinID==4) && (fitParams[8]<chi2cut) ) {
      yield_cent4[ptBinID][index]=fitParams[4];   
      yieldStatErr_cent4[ptBinID][k]=fitParams[5];
    }    
    //Printf("%i %i %i %i %e",ientry,centBinID, ptBinID, funcID, fitParams[4]);
  }
  
  //loop on pt bins
  for (Int_t ipt=0;ipt<npt;ipt++){  
    for (Int_t icent=0;icent<ncent;icent++){
      Double_t centralValue[2]={0.0, 0.0};
      GetCentralValue(yield_cent0[icent][ipt], yieldStatErr_cent0[icent][ipt], nfunc, centralValue);
      yield_central_cent0[icent][ipt]=centralValue[0];
      yield_centralErr_cent0[icent][ipt]=centralValue[1];
      Printf("icent=%i ipt=%i central=%e error=%e", icent, ipt, yield_central_cent0[icent][ipt], yield_centralErr_cent0[icent][ipt]);
      hYieldVsPt[icent]->SetBinContent(ipt+1, yield_central_cent0[icent][ipt]);
      hYieldVsPt[icent]->SetBinError(ipt+1, yield_centralErr_cent0[icent][ipt]);
      yield_systErr_cent0[icent][ipt]=(Double_t)GetSystematicError(yield_cent0[icent][ipt], yield_central_cent0[icent][ipt], nfunc);
      hSystVsPt[icent]->SetBinContent(ipt+1, yield_systErr_cent0[icent][ipt]);
    }
    /*
    Double_t centralValue[2]={0.0, 0.0};
    GetCentralValue(yield_cent1[ipt], yieldStatErr_cent1[ipt], nfunc, centralValue);
    yield_central_cent1[ipt]=centralValue[0];
    yield_centralErr_cent1[ipt]=centralValue[1];
    hYieldVsPt[1]->SetBinContent(ipt, yield_central_cent1[ipt]);
    hYieldVsPt[1]->SetBinError(ipt, yield_centralStatErr_cent1[ipt]);

    //yield_central_cent2[ipt]=(Double_t)GetCentralValue(yield_cent2[ipt], yieldStatErr_cent2[ipt], nfunc);
    //  yield_central_cent3[ipt]=(Double_t)GetCentralValue(yield_cent3[ipt], yieldStatErr_cent3[ipt], nfunc);
    
    //yield_systErr_cent0[ipt]=(Double_t)GetSystematicError(yield_cent0[ipt], yield_central_cent0[ipt],nfunc);
    yield_systErr_cent1[ipt]=(Double_t)GetSystematicError(yield_cent1[ipt],  yield_central_cent1[ipt], nfunc);
    hSystVsPt[1]->SetBinContent(ipt, yield_systErr_cent1[ipt]);
    //yield_systErr_cent2[ipt]=(Double_t)GetSystematicError(yield_cent2[ipt],  yield_central_cent2[ipt], nfunc);
    // yield_systErr_cent3[ipt]=(Double_t)GetSystematicError(yield_cent3[ipt],  yield_central_cent3[ipt], nfunc);
  
    // Printf("pt=%i cent=0 central=%9.3f",ipt, ptBinID, yield_central_cent0[ipt]);
    Printf("pt=%i cent=1 central=%9.3f",ipt, ptBinID, yield_central_cent1[ipt]);
    //Printf("#####  pt=%i cent=2 central=%9.3f",ipt, ptBinID, yield_central_cent2[ipt]);
    //  Printf("pt=%i cent=3 central=%9.3f",ipt, ptBinID, yield_central_cent3[ipt]); 
    */
  }
  
    
  TCanvas *cry=new TCanvas("cry","Raw yield vs p_{T}", 800,600);
  //cfit->Divide(3,2);
  cry->cd();
  frameYieldVsPt->GetYaxis()->SetRangeUser(1, 1e7);
  frameYieldVsPt->Draw();
  for (Int_t icentbin=0;icentbin<ncent;icentbin++){
    hYieldVsPt[icentbin]->Draw("same");
  }

  TCanvas *cs=new TCanvas("cs","Raw yield vs p_{T}", 800,600);
  //cfit->Divide(3,2);
  cs->cd();
  frameSystVsPt->GetYaxis()->SetRangeUser(0., 1e4);
  frameSystVsPt->Draw();
  for (Int_t icentbin=0;icentbin<ncent;icentbin++){
    hSystVsPt[icentbin]->Draw("same");
  }

  return;
  
}

//------------------------------------------------------------------------------------
Double_t GetCentralValue(Double_t *var_inArray, Double_t *varStatErr_inArray, Int_t nfunc, Double_t *centralValue)
{
  if (!var_inArray || (nfunc<1)) {
    Printf("GetCentralValuesArray: invalid arguments, doing nothing!");
    return;
  }
  //return;
  Int_t usedFunc = nfunc;
  Double_t weightedMean = 0.0, weight = 1.0, deltaValue = 1e20;
  for (Int_t k=0;k<usedFunc;k++){
    if ( (var_inArray[k]<=0.0) || (varStatErr_inArray[k]<=0.0) ) {     //skip values that do not pass chi2 cut
      usedFunc--;
    } else {
      weightedMean+= var_inArray[k]/(varStatErr_inArray[k]*varStatErr_inArray[k]);
      weight+= 1/(varStatErr_inArray[k]*varStatErr_inArray[k]);
    }
  }//loop func
  
  if (weight>0){
    weightedMean/=weight;
  } else {
    return 0.0;
  }
  
  centralValue[0] = weightedMean;    
  for (Int_t k=0;k<nfunc;k++){
    if ( (var_inArray[k]==0.0)|| (varStatErr_inArray[k]<=0.0) ) continue;
    Double_t diff = TMath::Abs(centralValue[0]-var_inArray[k]);
    if ( diff < deltaValue) {
      deltaValue = diff;
      centralValue[0]=var_inArray[k];
      centralValue[1]=varStatErr_inArray[k];
    }
  }
  Printf("----- Used functions =%i : central value = %e +/- %e", usedFunc, centralValue[0], centralValue[1] );
  return centralValue[0];
}

//------------------------------------------------------------------------------------
Double_t  GetSystematicError(Double_t *var_pt_inArray, Double_t centralValue =0.0, Int_t nfunc=0)
{
  /*
    The sistematic error is estimated as the sigma of the gaussian distribution of the difference (value-central_value) 
    if the number of available measurements is >10. 
    Otherwise it is taken as the difference between the two extreme values,
    divided by sqrt(12) (assuming uniform distribution).
  */

  if (!var_pt_inArray || (centralValue<=0.0) || (nfunc<1)) {
    Printf("GetSystematicErrorsArray: invalid arguments, doing nothing!");
    return;
  }
  //skip values that do not pass chi2 cut
  if ( centralValue == 0.0) return 0.0; 


  Double_t minDiff = 1e20, maxDiff = 0.0, minValue = 0.0, maxValue=1e20;
  TH1D* hDiff = new TH1D("hDiff","hDiff", 2e5, -1e5, 1e5);
  Printf("Gaussian fit over nfunc points = %i",nfunc);

  for (Int_t k=0;k<nfunc;k++){
    if (var_pt_inArray[k]<=0.0) continue;
    Double_t dummydiff = var_pt_inArray[k]-centralValue;
    Printf("%i -> diff = %9.4f", k, dummydiff);
    hDiff->Fill(dummydiff);
    if (nfunc<10) {
      if (dummydiff<minDiff) {
	minDiff = dummydiff;
	minValue = var_pt_inArray[k];
      }
      if (dummydiff>maxDiff) {
	maxDiff = dummydiff;
	maxValue= var_pt_inArray[k];
      }
    }
  }//loop func
  if (nfunc<10) {
    return (maxValue-minValue)/TMath::Sqrt(12);
  }else{
    hDiff->Fit("gaus");
    Double_t peakDiff=(hDiff->GetFunction("gaus"))->GetParameter(1);
    Double_t spreadDiff=(hDiff->GetFunction("gaus"))->GetParameter(2);
    Double_t peakDiffErr=(hDiff->GetFunction("gaus"))->GetParError(1);
    Double_t spreadDiffErr=(hDiff->GetFunction("gaus"))->GetParError(2);
    Printf("##################################################################################");
    Printf("Main peak diff (gaus): mean = %f +- %f\n",peakDiff,peakDiffErr );
    Printf("                       sigma = %f +- %f\n",spreadDiff,spreadDiffErr );
    return spreadDiff;    
  }
}



 // switch funcID {
    //   case EFitFunction::kPOLY2 :
    // 	funcLabel = "POLY2";
    // 	break;
    //   case EFitFunction::kPOLY3 :
    // 	funcLabel = "POLY3";
    // 	break;
    //   case EFitFunction::kPOLY1 :
    // 	funcLabel = "POLY1";
    // 	break;
    //   case EFitFunction::kLandau :
    // 	funcLabel = "LAND";
    // 	break;
    //   case EFitFunction::EXP :
    // 	funcLabel = "EXP";
    // 	break;
    //   default:
    // 	funcLabel = "none";
    // 	break;
    //   } 
