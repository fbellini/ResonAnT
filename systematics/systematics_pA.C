  enum EQuantity{ kYields=0,
		  kCorrYields,
		  kMass,
		  kWidth};

void systematics(Int_t selCent = 0,
		 Float_t checkBarlow = 0.0,
		 Bool_t usetree = kFALSE,
		 TString inlist = "list.lst",
		 Int_t quantityID = EQuantity::kYields,
		 TString centralFileName = 
		 "/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/fitEM_norm1_BWpoly2_fixedW_default/best_fit_poly2.root",
		 //"/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/fitEM_norm1_BWpoly2_fixedW/RAW_best_fit_poly2.root",
		 //"/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/0to100_binB/systPID/tpc2s_tof3sveto/fitEM_norm1_BWpoly2_fixedW/best_fit_poly2.root",
		 TString filebins = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/proj_2424_tpc2s_tof3sveto.root",
		 //"/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/0to100_binB/proj_2424_tpc2s_tof3sveto.root",
		 Bool_t verbose = 1,
		 Int_t isTOF = 0)
{

  Color_t color[3][6]={kOrange+7, kPink+6, kGreen+1, kAzure+1, kBlue+3, kBlack, //combined
		       kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-4, kBlack, //tpc
		       kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof
 
  Color_t fillcolor[3][6]={kOrange+7, kPink+6, kGreen+1, kAzure+1, kBlue+3, kBlack, //combined
			   kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack, //tpc
			   kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof
  
  Int_t marker[3][6]={21, 22, 32, 28, 24, 20, //combined
		      21, 22, 23, 34, 33, 20, //tpc
		      25, 26, 32, 28, 27, 24}; //tof
  
  //get bins
  Int_t npt_axis = 0, ncent_axis=0; 
  TFile *f=TFile::Open(filebins.Data());
  if (!ptbins || !centbins) return;
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  npt_axis = ptbins->GetNbins();  
  TAxis *centbins = (TAxis*)f->Get("centbins");
  ncent_axis = centbins->GetNbins();
  f->Close();
  const Int_t npt = npt_axis;
  const Int_t ncent = ncent_axis;
  Double_t pt[npt+1];
  for (Int_t k=0; k<npt+1;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
  }
  Double_t cent[ncent+1]; 
  for (Int_t k=0; k<ncent+1;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
  }

  //define arrays with values
  Double_t val[20]; //values per pT bin 
  Double_t valErr[20]; //error on value per pT bin
  Double_t central = 0.0; //central value for a single pT bin
  Double_t centralErr = 0.0; //error on central value for a single pT bin
  Int_t nused = 0; //number of used measurements
  Int_t navail= 0; //number of available measurments
  TTree * tree = 0x0; TTree * centralTree = 0x0; //input tree
  TList * list = new TList();
  TList * centralList = new TList();//input list

  //define titles
  TString titleHistoY, nameHisto, titleHisto;  
  
  if (quantityID==EQuantity::kCorrYields){
    titleHistoY="1/N_{evts}*1/B.R.*dN/dydp_{T}*1/#epsilon";
    nameHisto=Form("hCorrected_%i",selCent);//"hYieldVsPt" ;
    // titleHisto=Form("Corrected yield %s",centLabel.Data());
  }
  if (quantityID==EQuantity::kYields){
    titleHistoY="dN/dp_{T}";
    nameHisto=Form("hRawYieldVsPt_%i",selCent);// "raw";//
    // titleHisto=Form("Raw yield %s",centLabel.Data());
  }
  if (quantityID==EQuantity::kMass){
    titleHistoY="M (GeV/c^{2})";
    nameHisto="mass";//"hMassVsPt";
    //  titleHisto=Form("BW mass %s",centLabel.Data());
  }
  if (quantityID==EQuantity::kWidth){
    titleHistoY="#Gamma (GeV/c^{2})";
    nameHisto="gamma";//"hWidthVsPt";
    //  titleHisto=Form("BW width %s",centLabel.Data());
  }
  
  
  TFile * centralFile = TFile::Open(centralFileName.Data());
  if (!centralFile) { Printf("----- Central file not valid"); return; }
  //get input tree or histogram list
  if (usetree) { 
    if (verbose) Printf("::::: Getting input values from tree");
    tree = (TTree*) GetMergedTree(inlist);  
    if (!tree) { printf("----- Main tree not found!"); return; }
    centralTree = (TTree*) centralFile->Get("tree");
    if (!centralTree) { printf("----- Central Tree not found!"); return; }  
  } else {
    if (verbose) Printf("::::: Getting input values from histos %s",  nameHisto.Data());
    GetListOfHistos(inlist, nameHisto.Data(), list);
    list->SetName("input");
    centralHisto = (TH1*) centralFile->Get(nameHisto.Data());
    if (!centralHisto) { printf("----- Central Histo not found!"); return; }
    centralList = new TList(); centralList->Add(centralHisto);
    if (!centralList) { printf("----- Central List not found!"); return; }
    centralList->SetName("central");
  }
  
  TCanvas * c = new TCanvas(Form("cNiDelta_c%i",selCent),Form("cNiDelta_c%",selCent), 1200,900);
  c->Divide(6,4);

  //TODO add here definition of histograms
  Int_t indexLS = 0;
  TString newfilename = Form("systematics");
  if (centralFileName.Contains("fitLS")) {
    newfilename.Append("LS");
    indexLS = 1;
  }  else newfilename.Append("EM");
  if  (checkBarlow>0) newfilename.Append(Form("_Barlow%2.1f", checkBarlow));
  newfilename.Append(Form("_%i.root", selCent));
  TFile * fout = new TFile(newfilename.Data(),"recreate");
  
  Double_t systRel = 0.0; Double_t syst = 0.0;

  //loop on centrality and pT -- each pt bin processed separately
  //for (Int_t selCent = 100; selCent<101; selCent++){    
    TString centLabel=Form("%3.0f-%3.0f%%",(selCent==100 ? 0 : cent[selCent]),(selCent==100 ? 100 : cent[selCent+1]));
    
    Int_t index = ((selCent==100)? 5 : selCent);
    TH1F* hXVsPt = new TH1F(Form("%s_%i",nameHisto.Data(), selCent),";p_{T} (GeV/c);", npt_axis, pt);
    hXVsPt->SetTitle(titleHisto.Data());
    hXVsPt->GetYaxis()->SetTitle(titleHistoY.Data());
    hXVsPt->SetLineColor(color[isTOF][index]);	       
    hXVsPt->SetMarkerColor(color[isTOF][index]);	       
    hXVsPt->SetMarkerStyle(marker[isTOF][index]);	       
  
    TH1F* hXVsPtWsyst = new TH1F(Form("%sWsyst_%i",nameHisto.Data(),selCent),"; p_{T} (GeV/c);", npt_axis, pt);
    hXVsPtWsyst->SetTitle(titleHisto.Data());
    hXVsPtWsyst->GetYaxis()->SetTitle(titleHistoY.Data());
    hXVsPtWsyst->SetLineColor(color[isTOF][index]);	       
    hXVsPtWsyst->SetMarkerColor(color[isTOF][index]);	       
    hXVsPtWsyst->SetMarkerStyle(marker[isTOF][index]);
  
    TH1F* hSystVsPt = new TH1F(Form("hSystVsPt_%i",selCent),"syst. uncert.; p_{T} (GeV/c); syst. uncert.", npt_axis, pt);
    hSystVsPt->SetTitle(Form("syst. uncert. cent. %s", centLabel.Data()));
    hSystVsPt->SetLineColor(color[isTOF][index]);	       
    hSystVsPt->SetMarkerColor(color[isTOF][index]);	       
    hSystVsPt->SetMarkerStyle(marker[isTOF][index]);	       
  
    TH1F* hSystVsPtPercentageOfCentral = new TH1F(Form("hSystVsPtPercentageOfCentral_%i",selCent),"relative syst. uncert. Vs p_{T}; p_{T} (GeV/c); relative syst. uncert.(%)", npt_axis, pt);
    hSystVsPtPercentageOfCentral->SetTitle(Form("cent. %s", centLabel.Data()));
    hSystVsPtPercentageOfCentral->SetLineColor(color[isTOF+indexLS][index]);	       
    hSystVsPtPercentageOfCentral->SetMarkerColor(color[isTOF+indexLS][index]+indexLS);	       
    hSystVsPtPercentageOfCentral->SetMarkerStyle(marker[isTOF+indexLS][index]);	       
  
    TH1F* hStatVsPt = new TH1F(Form("hStatVsPt_%i",selCent),"stat. uncert.; p_{T} (GeV/c); stat. uncert.", npt_axis, pt);
    hStatVsPt->SetTitle(Form("stat. uncert. cent. %s", centLabel.Data()));
    hStatVsPt->SetLineColor(color[isTOF][index]);	       
    hStatVsPt->SetMarkerColor(color[isTOF][index]);	       
    hStatVsPt->SetMarkerStyle(marker[isTOF][index]);	       
  
    TH1F* hStatVsPtPercentageOfCentral = new TH1F(Form("hStatVsPtPercentageOfCentral_%i",selCent),"relative stat. uncert.; p_{T} (GeV/c); relative stat. uncert. (%)", npt_axis, pt);
    hStatVsPtPercentageOfCentral->SetTitle(Form("%s", centLabel.Data()));
    hStatVsPtPercentageOfCentral->SetLineColor(color[isTOF][index]);	       
    hStatVsPtPercentageOfCentral->SetMarkerColor(color[isTOF][index]);	       
    hStatVsPtPercentageOfCentral->SetMarkerStyle(marker[isTOF][index]);	      

    for (Int_t selPt=0; selPt<npt; selPt++){
      //reset values
      nused = 0;
      for (Int_t j=0;j<20;j++) { 
	val[j]=-1.0; 
	valErr[j]=-1.0;
      }
      central = 0.0;  
      centralErr = 0.0;
      
      //get central value and its error +
      //get all available values for each pT bin 
      //and their number
      if (usetree) {
	nused = GetValuesFromTree(tree, selCent, selPt, quantityID, val, valErr);
	GetValuesFromTree(centralTree, selCent, selPt, quantityID, &central, &centralErr);
      } else {
	nused = GetValuesFromHisto(list, selCent, selPt, val, valErr);
	GetValuesFromHisto(centralList, selCent, selPt, &central, &centralErr);
      }
      
      //if centrality is 80-100% stop at 6 GeV/c
      if (selCent==4 && selPt>=17) {
	central = 0.0;
	centralErr = 0.0;
      } 
      //some printouts
      if (verbose) {
	Printf("::::: Pt bin %i - Central value = %e +/- %e", selPt, central, centralErr); 
	for (Int_t j=0;j<nused;j++) { 
	  Printf("x[%i] = %e +/- %e", j, val[j], valErr[j]);
	}
	Printf("::::: Pt bin %i - Values available for systematics = %i", selPt, nused);
      }
      
      //calculate systematic uncertainty as RMS
      TH1D * hVal = new TH1D(Form("hc%ipt%i",selCent,selPt),Form("hc%ipt%i",selCent,selPt), 200, central-10.0*centralErr, central+10.0*centralErr);
      hVal->Fill(central);
      for (Int_t i=0; i<nused; i++){
	if (val[i]<=0) { Printf("Skipping invalid value %i", i); continue; }
	if ( (val[i]<central-10.0*centralErr) || (val[i]>central+10.0*centralErr) ) {
	  Printf("CAREFUL: value %i is outside range of +/-10sigma from central!!!",i);
	  continue;
	}
	if ((checkBarlow>0.0) && IsStatisticallyCompatible(central, centralErr, val[i], valErr[i], 0, checkBarlow)) {
	  if (verbose) Printf("Barlow: value %i is statistically compatible (Nb=%f) with central",i, checkBarlow);
	  continue;
	}	
	hVal->Fill(val[i]);
	Printf("-- Using value %i", i);
      }
      Printf("::::: Pt bin %i: histogram entries = %3.0f", selPt, hVal->GetEntries());
      syst = hVal->GetRMS(1);/*x-axis*/
      if (central>0.0) systRel = syst*100./central;
      else  systRel = 0.0;
      Printf("::::: Systematic uncert = %e (%6.4f%%)", syst, systRel);	
      fout->cd();
      hVal->Write();

      //fill histos for mass or width (no dpt division)
      hXVsPt->SetBinContent(selPt+1, central);
      hXVsPt->SetBinError(selPt+1, centralErr);
      
      //fill yield histos, syst error only
      hXVsPtWsyst->SetBinContent(selPt+1, central);
      hXVsPtWsyst->SetBinError(selPt+1, syst);
      
      //fill systematic error histos
      hSystVsPt->SetBinContent(selPt+1, syst);
      hSystVsPtPercentageOfCentral->SetBinContent(selPt+1, systRel);
      
      //fill statistical error histos
      hStatVsPt->SetBinContent(selPt+1, centralErr);
      if (central>0) hStatVsPtPercentageOfCentral->SetBinContent(selPt+1, centralErr*100./central);
      else hStatVsPtPercentageOfCentral->SetBinContent(selPt+1, 0.0);

      //plot Ni (Barlow) and get median and sigma --> save in output file
      TH1D * hNiDelta = new TH1D(Form("hNiDelta_c%ipt%i",selCent,selPt),Form("hNiDelta_c%ipt%i",selCent,selPt), 160, -20., 20.);
      for (Int_t i=0; i<nused; i++){
	Float_t nDeltaBarlow = GetNiDeltaBarlow(central, centralErr, val[i], valErr[i]);
	hNiDelta->Fill(nDeltaBarlow);
      }
      Float_t median = Median(hNiDelta);
      Float_t mean = hNiDelta->GetMean();
      Float_t sigma = hNiDelta->GetRMS();
      TH1D * hNiDeltaCheck = new TH1D(Form("hNiDeltaStat_c%ipt%i",selCent,selPt),Form("hNiDeltaStat_c%_ipt%i",selCent,selPt), 3, 0., 3.);
      hNiDeltaCheck->GetXaxis()->SetBinLabel(1,"median"); hNiDeltaCheck->SetBinContent(1, median); 
      hNiDeltaCheck->GetXaxis()->SetBinLabel(2,"mean"); hNiDeltaCheck->SetBinContent(2, mean); 
      hNiDeltaCheck->GetXaxis()->SetBinLabel(3,"sigma"); hNiDeltaCheck->SetBinContent(3, sigma); 
      // TPaveText * pave = new TPaveText(0.5, 0.4, 0.9, 0.8,"NDC");
      // pave->SetFillColor(kWhite);
      // pave->AddText(Form("median = %4.2f", median));
      // pave->AddText(Form("mean   = %4.2f", mean));
      // pave->AddText(Form("RMS    = %4.2f", sigma));
      //TCanvas * c = new TCanvas(Form("cNiDelta_c%ipt%i",selCent,selPt),Form("cNiDelta_c%ipt%i",selCent,selPt), 600,600);
      c->cd(selPt+1);
      gStyle->SetOptStat("rems");
      hNiDelta->Draw();
      gPad->Update();
      TPaveStats * stats = (TPaveStats * ) hNiDelta->FindObject("stats");
      stats->SetX1NDC(0.5);
      stats->SetX2NDC(0.95);
      stats->SetY1NDC(0.5);
      stats->SetY2NDC(0.90);
      
   // pave->Draw("same");
      fout->cd();
      hNiDelta->Write();
      hNiDeltaCheck->Write();
      //      c->Write();
    }
    
    //prepare canvases
    TCanvas *cry=new TCanvas("cry","Raw yield vs p_{T}", 600,700);
    cry->cd();
    hXVsPtWsyst->Draw("E2");
    hXVsPt->Draw("same");
    
    TLegend*errleg = new TLegend(0.65,0.65,0.8,0.88, Form("%s",centLabel.Data()));
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
    hStatVsPtPercentageOfCentral->GetYaxis()->SetRangeUser(0.1, 50.0);
    hStatVsPtPercentageOfCentral->Draw();
    hSystVsPtPercentageOfCentral->SetLineWidth(2);
    hSystVsPtPercentageOfCentral->GetYaxis()->SetRangeUser(0.1, 50.0);
    hSystVsPtPercentageOfCentral->Draw();

    //save plot per centrality in file
    fout->cd();
    cry->Write();
    cs->Write();
    hXVsPt->Write();
    hXVsPtWsyst->Write();
    hSystVsPt->Write();
    hSystVsPtPercentageOfCentral->Write();
    hStatVsPt->Write();
    hStatVsPtPercentageOfCentral->Write();
  return;
}

//------------------------------------------------------------------------------------
Double_t GetDeltaBarlow(Double_t m1, Double_t e1, Double_t m2, Double_t e2)
{
  //returns Delta_i = sqrt(s_i^2 - s_c^2)
  Double_t y1,y2, s1,s2;
  if (e1>e2) {
    y1=m2; s1=e2;
    y2=m1; s2=e1;
  } else {
    y1=m1; s1=e1;
    y2=m2; s2=e2;
  }
  return TMath::Sqrt(s2*s2-s1*s1);
}

//------------------------------------------------------------------------------------
Double_t GetNiDeltaBarlow(Double_t m1, Double_t e1, Double_t m2, Double_t e2)
{
  //returns N_i = |y_i - y_c| / Delta_i 
  Double_t y1,y2, s1,s2;
  Int_t sign = 1;
  if (e1>e2) {
    y1=m2; s1=e2;
    y2=m1; s2=e1;
  } else {
    sign = -1;
    y1=m1; s1=e1;
    y2=m2; s2=e2;
  }
  Double_t diffy          = TMath::Abs(y1-y2);
  Double_t deltaSigma2    = TMath::Sqrt(s2*s2-s1*s1);
  if (deltaSigma2 == 0) return 0;
  return sign*diffy/deltaSigma2;
}
//------------------------------------------------------------------------------------
Bool_t IsStatisticallyCompatible(Double_t m1, Double_t e1, Double_t m2, Double_t e2, Bool_t verbose=0, Float_t Nb = 1.0) 
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
  if (numberOfDeltas<=Nb) {
    Printf("Barlow check (Nb=%f): measurements are statistically compatible", Nb);
    return kTRUE;
  }
  if (verbose) Printf("deltaSigma2 = %e , diffy = %e , numberOfDeltas = %5.2f", TMath::Sqrt(deltaSigma2), diffy, numberOfDeltas);
  if (numberOfDeltas>=5.0) Printf("Barlow Check warning: number of deltas = %5.2f > 5.0 - YOU BETTER CHECK THE MEASUREMENT!", numberOfDeltas);
  return kFALSE;
}



Int_t GetValuesFromTree(TTree * tree, Int_t selCent, Int_t selPt, Int_t quantityID, Double_t * val, Double_t * valErr)
{
  if (!tree) return -1;
  if (!val || !valErr) return -1;
  
  Float_t chi2cut = 5.0;
  //read tree from chain
  Double_t fitParams[9], SoverB=0.0;
  Int_t centBinID, ptBinID;
  tree->SetBranchAddress("signalMass",&fitParams[0]);
  tree->SetBranchAddress("signalMassErr",&fitParams[1]);
  tree->SetBranchAddress("signalWidth",&fitParams[2]);
  tree->SetBranchAddress("signalWidthErr",&fitParams[3]);
  tree->SetBranchAddress("nSignal",&fitParams[4]);
  tree->SetBranchAddress("nSignalErr",&fitParams[5]);
  tree->SetBranchAddress("chi2",&fitParams[8]);  
  tree->SetBranchAddress("centBin",&centBinID);
  tree->SetBranchAddress("ptBin",&ptBinID);
  
  Int_t used=0;
  for (Int_t ientry = 0; ientry < tree->GetEntries(); ientry++){
    tree->GetEntry(ientry);
    if (centBinID!=selCent) continue;
    if (ptBinID!=selPt) continue;
    if ( (fitParams[8]>0.0) && (fitParams[8]<chi2cut) ) {
      if (quantityID==EQuantity::kMass){
	val[used]=fitParams[0];
	valErr[used]=fitParams[1];  
      } else {
	if (quantityID==EQuantity::kWidth){
	  val[used]=fitParams[2];
	  valErr[used]=fitParams[3];  
	} else {
	  if (quantityID==EQuantity::kYields){
	    val[used]=fitParams[4];
	    valErr[used]=fitParams[5];  
	  }
	}
      }
      used++;
    }//end check on chi2 
  }//end loop on tree entries
  return used;
}

TTree * GetMergedTree(Char_t * inlist = "tree_raw.lst")
{
  //chain trees
  TString infile;
  Int_t ifiles=0; 
  FILE * files = fopen(inlist, "r") ; 
  TChain *tree = new TChain("tree");
  while ( infile.Gets(files) ){
    if (infile.Contains("SKIP")==0) {
      tree->Add(infile.Data());
      Printf("Adding tree from %s",infile.Data());
      ifiles++;
    }
  }
  Printf("Number of files merged = %i\n",ifiles);
  return tree;
} 

void GetListOfHistos(Char_t * inlist = "histo_raw.lst", TString nomeHisto = "raw", TList * outList)
{
  //add all histos in a unique list
  if (!outList) {
    Printf("----- Invalid list passed as input");
    return;
  }
  TH1D * histo;
  //chain trees
  TString infile;
  Int_t ifiles=0; 
  FILE * files = fopen(inlist, "r") ; 
  while ( infile.Gets(files) ){
    if (infile.Contains("SKIP")==0) {
      TFile * fin = TFile::Open(infile.Data());
      if (!fin) {Printf("Cannot open input file"); return;}
      histo = (TH1D*) fin->Get(nomeHisto.Data());
      if (!histo) { Printf("something wrong here"); return;} 
      outList->AddLast(histo);
      Printf("Adding histo from %s",infile.Data());
      ifiles++;
    }
  }
  Printf("Number of histogram added to the list = %i\n",ifiles);
  outList->ls();
  return;
}

Int_t * GetValuesFromHisto(TList * list, Int_t selCent, Int_t selPt, Double_t * val,Double_t *valErr)
{
  //add all histos in a unique list
  Int_t nhistos = list->GetEntries();
  for (Int_t j=0; j<nhistos; j++){
    TH1D* histo = (TH1D*) list->At(j);
    val[j] = histo->GetBinContent(selPt+1);
    valErr[j] = histo->GetBinError(selPt+1);
  }
  return nhistos;
}

double Median(const TH1D * h1) 
{ 
  //Get Median from histo
   int n = h1->GetXaxis()->GetNbins();  
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray(); 
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]); 
}

