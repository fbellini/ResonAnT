//---------------------------------------------------------------------------------------
//            fbellini@cern.ch     09 gen 2013
//---------------------------------------------------------------------------------------

void runSystAllCent(Char_t * inlist=NULL, Float_t chi2cut=2.)
{
  for (Int_t j=0;j<4,j++){
    systematics(inlist, chi2cut);
  }
  return;
}

//---------------------------------------------------------------------------------------
TString makeTreeChain(Char_t * inlist=NULL){
  //chain trees
  TString infile;
  Int_t ifiles=0; 
  FILE * files = fopen(inlist, "r") ; 
  TChain *tree = new TChain("tree");
  while ( infile.Gets(files) ){
    tree->Add(infile.Data());
    //Printf("Adding tree from %s",infile.Data());
    ifiles++;
  } 
  Printf("Number of merged trees = %i\n",ifiles);
  TString name = Form("tree_%s.root",inlist);
  tree->SaveAs(name.Data());
  return name.Data();
}

//---------------------------------------------------------------------------------------
void systesiMass(/*TString path = "",*/Char_t * inlist=NULL,  Char_t * inlist2=NULL, Int_t selCent = 0, Float_t chi2cut=2.0, TString filebins="sub_analysisAOD_0-80.root", Bool_t display=0, TString suffix)
{
  enum EFitFunction{ kPOLY2,
		     kPOLY3,
		     kLandau,
		     kPOLY1,
		     kEXP,
		     kData};
  
  enum EQuantity{ kYields,
		  kMass,
		  kWidth};
  
  Color_t color[]={kRed, kOrange, kGreen+2, kBlue, kMagenta, kBlack};
  Int_t marker[]={20, 21, 28, 22, 23};
  

  gROOT->LoadMacro("kstar/MakeSpectra.C");

  //get bins
  Int_t npt_axis = 0, ncent_axis=0; 
  TFile *f=TFile::Open(filebins.Data());//"/Users/bellini/alice/resonances/myKstar/pwglf_train_out/data/aod049_28nov_rebinHighPt/sub_EMnorm1.30-1.50_aod049_kstar.root");
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

  TString centLabel=Form("%2.0f-%2.0f%%",cent[selCent],cent[selCent+1]);
  
  //define output
  TString outname(inlist);
  outname.ReplaceAll(".lst",".root");
  outname.Prepend(Form("syst_cent%i_chi%2.1f_",selCent,chi2cut));
  //  TFile * outfile = new TFile(outname.Data(),"recreate"); 
  TFile * outfile = new TFile(Form("tesi%s_cent0%i_chiMax%2.1f.root",suffix.Data(), selCent,chi2cut),"recreate");
    
  //chain trees
  TString infname1 = Form("%s",inlist);
  TString infname2 = Form("%s",inlist2);
  TString name1 = makeTreeChain(infname1.Data());
  TString name2 = makeTreeChain(infname2.Data());
  
  //display pt spectrum from list 1
  TH1F* spectrum1 = (TH1F*) MakeSpectraCent(EQuantity::kMass,selCent, name1.Data(), filebins.Data(), chi2cut);
  TString newname1 = Form("%s",inlist);
  newname1.ReplaceAll(".root","");
  newname1.ReplaceAll("fitEM","raw_EM");
  // if (newname1.Contains("POLY2")) spectrum1->SetTitle(Form("Poly 2, %s",centLabel.Data()));
  // if (newname1.Contains("POLY3")) spectrum1->SetTitle(Form("Poly 3, %s",centLabel.Data()));
  // if (newname1.Contains("LAND")) spectrum1->SetTitle(Form("Landau, %s",centLabel.Data()));
  spectrum1->SetName("h1");
  spectrum1->SetTitle(newname1.Data());
  spectrum1->SetMarkerStyle(20);
  outfile->cd();
  spectrum1->Write();
  
  //display pt spectrum from list 2
  TH1F* spectrum2 = (TH1F*) MakeSpectraCent(EQuantity::kMass,selCent, name2.Data(), filebins.Data(), chi2cut);
  TString newname2 = Form("%s",inlist2);
  newname2.ReplaceAll(".root","");
  newname2.ReplaceAll("fitEM","raw_EM");
  // if (newname2.Contains("POLY2")) spectrum2->SetTitle(Form("Poly 2, %s",centLabel.Data()));
  // if (newname2.Contains("POLY3")) spectrum2->SetTitle(Form("Poly 3, %s",centLabel.Data()));
  // if (newname2.Contains("LAND")) spectrum2->SetTitle(Form("Landau, %s",centLabel.Data()));
  spectrum2->SetName("h2");
  spectrum2->SetTitle(newname2.Data());
  spectrum2->SetMarkerStyle(25);
  outfile->cd();
  spectrum2->Write();
  
  TH1F* ratio = (TH1F*) spectrum1->Clone("h1overh2");
  ratio->Divide(h2);
  ratio->SetTitle("h1/h2");
  ratio->SetMarkerStyle(28);
  outfile->cd();
  ratio->Write();

  //define output
  TH1F*frameMassVsPt=new TH1F("frameMassVsPt","Raw yield vs. p_{t}; p_{t} (GeV/c); dN/dp_{t}", npt_axis, pt);
  TH1F*frameSystVsPt=new TH1F("frameSystVsPt","; p_{t} (GeV/c); syst. error", npt_axis, pt);

  TH1F* hMassVsPt = new TH1F(Form("hMassVsPt_%i",selCent),"stat. err. only; p_{t} (GeV/c); M (GeV/c^{2})", npt_axis, pt);
  hMassVsPt->SetTitle(Form("centrality %s", centLabel.Data()));
  hMassVsPt->SetLineColor(kBlack/*color[selCent]*/);	       
  hMassVsPt->SetMarkerColor(kBlack/*color[selCent]*/);	       
  hMassVsPt->SetMarkerStyle(marker[selCent]);	       
  
  TH1F* hMassVsPtWsyst = new TH1F(Form("hMassVsPtWsyst_%i",selCent),"#sqrt{#sigma_{stat}^{2}+#sigma_{syst}^{2}}; p_{t} (GeV/c); M (GeV/c^{2})", npt_axis, pt);
  hMassVsPtWsyst->SetTitle(Form("centrality %s", centLabel.Data()));
  hMassVsPtWsyst->SetLineColor(color[selCent]+2);	       
  hMassVsPtWsyst->SetMarkerColor(color[selCent]+2);	       
  hMassVsPtWsyst->SetMarkerStyle(marker[selCent]);
  
  TH1F* hSystVsPt = new TH1F(Form("hSystVsPt_%i",selCent),"; p_{t} (GeV/c); systematic error", npt_axis, pt);
  hSystVsPt->SetTitle(Form("centrality %s", centLabel.Data()));
  hSystVsPt->SetLineColor(color[selCent]);	       
  hSystVsPt->SetMarkerColor(color[selCent]);	       
  hSystVsPt->SetMarkerStyle(marker[selCent]);	       
  
  TH1F* hSystVsPtPercentageOfCentral = new TH1F(Form("hSystVsPtPercentageOfCentral_%i",selCent),"; p_{t} (GeV/c); systematic error (% of central value)", npt_axis, pt);
  hSystVsPtPercentageOfCentral->SetTitle(Form("centrality %s", centLabel.Data()));
  hSystVsPtPercentageOfCentral->SetLineColor(color[selCent]);	       
  hSystVsPtPercentageOfCentral->SetMarkerColor(color[selCent]);	       
  hSystVsPtPercentageOfCentral->SetMarkerStyle(marker[selCent]);


  //read input trees
  Double_t yield1, yield2, erry1, erry2;//fitParams1[9],fitParams2[9];
  Int_t centBinID1,ptBinID1,
    centBinID2,ptBinID2;
  TFile * f1 = TFile::Open(name1.Data());
  TChain *tree1 = (TTree*) f1->Get("tree");
  tree1->SetBranchAddress("nSignal",&yield1);
  tree1->SetBranchAddress("nSignalErr",&erry1);
  tree1->SetBranchAddress("centBin",&centBinID1);
  tree1->SetBranchAddress("ptBin",&ptBinID1);

  TFile * f2 = TFile::Open(name2.Data());
  TChain *tree2 = (TTree*) f2->Get("tree");
  tree2->SetBranchAddress("nSignal",&yield2);
  tree2->SetBranchAddress("nSignalErr",&erry2);
  tree2->SetBranchAddress("centBin",&centBinID2);
  tree2->SetBranchAddress("ptBin",&ptBinID2);

  Int_t nentries = TMath::Max(tree1->GetEntries(),tree2->GetEntries());
  Printf("Max n. entries = %i", nentries);
  for (Int_t ientry = 0;ientry<nentries;ientry ++){
    if ( (tree1->GetEntry(ientry)>0) && (tree2->GetEntry(ientry)>0) ) {
      printf("-- Pt1=%i, Pt2=%i ===> ", ptBinID1,ptBinID2);
      Double_t diffy=0.0, delta2=0.0, syst=0.0, nsigma=0.0; 
      diffy = TMath::Abs(yield1-yield2);
      delta2 = TMath::Abs(erry2*erry2-erry1*erry1);
      nsigma = diffy/TMath::Sqrt(delta2);
      Printf("Delta=%e , diffy=%e , nsigma=%5.2f", TMath::Sqrt(delta2), diffy, nsigma);
      syst = diffy/2.0;// (nsigma-1)*TMath::Sqrt(delta2);
      if (syst<=0) {
	Printf("Error compatible with stat - syst error negligible");
	syst=1e-5;
      }
      //when fill histo add +2 because definiion of bins start
      Float_t dpt = ptbins->GetBinWidth(ptBinID1 + 1);
      hMassVsPt->SetBinContent(ptBinID1+1, yield1/dpt);
      hMassVsPt->SetBinError(ptBinID1+1, erry1/dpt);
      hMassVsPtWsyst->SetBinContent(ptBinID1+1, yield1/dpt);
      hMassVsPtWsyst->SetBinError(ptBinID1+1, TMath::Sqrt(erry1*erry1/(dpt*dpt)+syst*syst));
      hSystVsPt->SetBinContent(ptBinID1+1, syst/dpt);
      hSystVsPtPercentageOfCentral->SetBinContent(ptBinID1+1, (syst*100)/yield1);
    }  
    else {
      if  (tree1->GetEntry(ientry)>0) {Printf("Nothing done because tree 2 does not have entry %i",ientry);}
      else { Printf("Nothing done because tree 1 does not have entry %i",ientry);}
    }
  }//loop on entries

  if (display){
  TCanvas *cry=new TCanvas("cry","Raw yield vs p_{T}", 1000,600);
  cry->Divide(2,1);
  cry->cd(1);
  frameMassVsPt->GetYaxis()->SetRangeUser(9e2, 5e6);
  gPad->SetLogy();
  frameMassVsPt->Draw();
  hMassVsPt->Draw("same");
  cry->cd(2);
  frameMassVsPt->GetYaxis()->SetRangeUser(1, 1e7);
  gPad->SetLogy();
  frameMassVsPt->Draw();
  hMassVsPtWsyst->Draw("same");

  TCanvas *cs=new TCanvas("cs","Systematic error vs p_{T}", 800,600);
  cs->Divide(1,2);
  cs->cd(1);
  frameSystVsPt->GetYaxis()->SetRangeUser(0.1, 1e6);
  gPad->SetLogy();
  frameSystVsPt->Draw();
  hSystVsPt->Draw("same");
  cs->cd(2);
  hSystVsPtPercentageOfCentral->GetYaxis()->SetRangeUser(0.1, 1e2);
  hSystVsPtPercentageOfCentral->Draw();
  }
  outfile->cd();
  if (display) {
    cry->Write();
    cs->Write();
  }
  frameMassVsPt->Write();
  hMassVsPt->Write();
  frameSystVsPt->Write();
  hMassVsPtWsyst->Write();
  hSystVsPt->Write();
  hSystVsPtPercentageOfCentral->Write();
  return;  
  
}

//------------------------------------------------------------------------------------
Double_t GetCentralValueWM(const Double_t *varIn, const Double_t *varIn_StatErr, Int_t nfunc, Double_t *centralValue)
{
  //set centralValue to the min difference wrt to the weighted mean
  //the weighted mean is returned by the function
  
  if (!varIn || (nfunc<1)) {
    //    Printf("GetCentralValue:: invalid arguments, doing nothing!");
    return;
  }
  
  Double_t weightedMean = TMath::Mean(nfunc, varIn, varIn_StatErr);
  
  centralValue[0] = weightedMean;    
  Printf("weighted mean = %e",weightedMean);
  for (Int_t k=0;k<nfunc;k++){
    Double_t diff = TMath::Abs(centralValue[0]-varIn[k]);
    //   Printf("var[%i]=%e    central=%e    diff=%e  delta=%e",k, varIn[k],centralValue[0], diff, deltaValue );
    if ( diff < deltaValue) {
      deltaValue = diff;
      centralValue[0]=varIn[k];
      centralValue[1]=varIn_StatErr[k];
    }
  }
  Printf("GetCentralValue:: Used functions =%i central value = %e +/- %e", nfunc, centralValue[0], centralValue[1] );
  return weightedMean;
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
      Printf("---- func=%i     value=%e    err=%e    ", k,var_inArray[k],varStatErr_inArray[k]);
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
  Printf("weighted mean = %e used func=%i",weightedMean, usedFunc);
  
  for (Int_t k=0;k<nfunc;k++){
    if ( (var_inArray[k]<=0.0) || (varStatErr_inArray[k]<=0.0) ) continue;     //skip values that do not pass chi2 cut
    Double_t diff = TMath::Abs(centralValue[0]-var_inArray[k]);
    Printf("var[%i]=%e    central=%e    diff=%e  delta=%e",k, var_inArray[k],centralValue[0], diff, deltaValue );
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
Double_t  GetSystematicError(Double_t *var_pt_inArray, Double_t *varStatErr_inArray, Double_t centralValue =0.0, Double_t centralValueErr=0.0, Int_t nfunc=0)
{
  /*
    The sistematic error is estimated as the sigma of the gaussian distribution of the difference (value-central_value) 
    if the number of available measurements is >10. 
    Otherwise it is taken as the difference between the two extreme values,
    divided by sqrt(12) (assuming uniform distribution).
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
  
  Printf("=== GetSystematicError:: Systematic from %i points, for central value %e +/- %e",usedFunc, centralValue, centralValueErr);
  const Int_t nvalues=usedFunc;
  Double_t usedValues[nvalues], sigmas[nvalues];
  Int_t iused=0;
  Double_t syst = 0.0; 
  for (Int_t k=0;k<nfunc;k++){
    if ((var_pt_inArray[k]<=0.0)||(varStatErr_inArray[k]<=0.0))   continue;
    //copy into other array only good values and increase index
    //usedValues[iused]=var_pt_inArray[k];
    iused++;
    //save diff into th1 for gaussian fit
    Double_t delta = var_pt_inArray[k]-centralValue;
    Double_t deltaerr = TMath::Sqrt(TMath::Abs(varStatErr_inArray[k]*varStatErr_inArray[k]-centralValueErr*centralValueErr));
    Printf("delta = %e \n deltaerr = %e", delta, deltaerr);
    Double_t systcheck=TMath::Abs(delta/deltaerr);
    if (systcheck>1.0) {
      usedValues[k]=systcheck;
      sigmas[k]=deltaerr;
      Printf("ifunc= %i --> systcheck = Delta/sigma_delta = %5.3f --> used for syst",k,systcheck);
    } else {
      Printf("ifunc= %i --> systcheck = Delta/sigma_delta = %5.3f --> NOT syst",k,systcheck);
      usedValues[k]=0.0;
    }
  }//loop func
  
  if (iused==0) return 0.0;

  Double_t max = -1e10;
  Int_t deltaerrID=-1;
  //Double_t min =  1e10;
  for (Int_t i=0;i<iused;i++){
    // if (usedValues[i]<min)
    //   min=usedValues[i];
    if (usedValues[i]>max){
      max=usedValues[i];
      deltaerrID=i;
    }
  }
  syst=(max-1)*deltaerr[deltaerrID];
  //Printf("min= %e, max=%e",min,max);
  //old   syst = TMath::Abs(max-min)/TMath::Sqrt(12);
  
  Printf("GetSystematicErrors:: Systematic error with nfunc= %i points (unif) = %e",iused,syst);  
  return syst;
}



