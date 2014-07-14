enum EFitFunction{ kPOLY2,
		   kPOLY3,
		   kLandau,
		   kPOLY1,
		   kEXP,
		   kData};

Color_t color[]={kRed, kOrange, kGreen+2, kBlue, kMagenta, kBlack};
Color_t colorFunc[]={kBlue, kRed, kGreen+2, kMagenta+1, kYellow+2, kBlack};
Int_t marker[]={20, 21, 28, 22, 23};
Char_t funcName[5][10]={"BW+POLY2","BW+POLY3", "BW+LAND","BW+POLY1", "BW+EXP"};


void MakeTableFromFitTree(TString filein="fitEM_analysisAOD_0-80.root", Char_t* fileproj="sub_analysisAOD_0-80.root")
{
  //style
  gROOT->LoadMacro("SetGraphicStyle.C");
  SetGraphicStyle(0);
  gStyle->SetOptStat(0);
  gStyle->SetTextFont(42);

  //get bins
  TFile * f= TFile::Open(fileproj);
  if (!f) return;
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  Int_t npt = ptbins->GetNbins();
  const Int_t dimpt = npt+1;
  Double_t pt[dimpt];
  for (Int_t k=0; k<dimpt;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
    Printf("%5.2f",pt[k]);
  }

  TAxis *centbins = (TAxis*)f->Get("centbins");
  Int_t ncent = centbins->GetNbins();
  const Int_t dimcent = ncent+1;
  Double_t cent[dimcent]; 
  for (Int_t k=0; k<dimcent;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
    Printf("%5.2f",cent[k]);
  }

  //get tree
  TFile * fin=TFile::Open(filein.Data());
  Double_t fitParams[9], SoverB=0.0, significance=0.0,normfactorCopy=0.0;
  Int_t centBinID,ptBinID, funcID;
  Double_t rangeInfCopy, rangeSupCopy, ptinfCopy, ptsupCopy,massRangeInfCopy,massRangeSupCopy;
  char fitfunction[5];
  Int_t centinfCopy, centsupCopy;
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
  tree->SetBranchAddress("fitrange_inf", &massRangeInfCopy);
  tree->SetBranchAddress("fitrange_sup", &massRangeSupCopy);
  
  ofstream ftxt;
  TString nametxt = Form("%s",filein.Data());
  nametxt.ReplaceAll(".root",".txt");
  nametxt.ReplaceAll("roofit/","");
  ftxt.open(nametxt.Data(),ios::out);
  for (Int_t j=0;j<tree->GetEntries();j++){
    tree->GetEntry(j);
    if (funcID==EFitFunction::kPOLY2) sprintf(fitfunction,"poly2");
    if (funcID==EFitFunction::kPOLY3) sprintf(fitfunction, "poly3");
    if (funcID==EFitFunction::kLandau) sprintf(fitfunction,"Landau");

    centinfCopy = centBinID*20;
    centsupCopy = (centBinID+1)*20;
  
    //showpoint
    //width(int)
    //precision(int)
    //ftxt.width(3);
    //ftxt <<  " & ";
    ftxt << " " << centinfCopy <<   "-" << centsupCopy << "  &  ";   // cent bin
    ftxt.precision(1);
    ftxt << showpoint << ptinfCopy      <<   "-" <<  showpoint << ptsupCopy;  // pt bin
    //ftxt <<  " & " << showpoint << rangeInfCopy   <<   "-" << showpoint << rangeSupCopy;    // bg normalization range
    
    // ftxt.precision(3);
    // ftxt <<   " & " << normfactorCopy; // normalization factor
   
    // ftxt.precision(2);
    // ftxt << " & " << showpoint << massRangeInfCopy   << "-" << showpoint << massRangeSupCopy; //fit range
    //ftxt << " & " << fitfunction; //fit function    
    ftxt.precision(5);
    // ftxt <<  " & " << fixed << fitParams[0];//   << " $\\pm$ " 
    // ftxt <<  " & " << fixed << fitParams[1]; // fitted mass 
    ftxt <<  "  &  " << scientific << fitParams[4];// " $\\pm$ " 
    ftxt <<  "  &  " << scientific << fitParams[5];  //yield
    ftxt <<  "  &  " << scientific << fitParams[6];//  << " $\\pm$ " 
    ftxt <<  "  &  " << scientific << fitParams[7];  //background      
    ftxt.precision(2);
    ftxt <<  "  &  " << fixed << fitParams[8];  //chi_2
    ftxt.precision(2);
    // ftxt <<  " & " << fixed << SoverB  <<  //S/B
    //   " & " << fixed << significance  <<   //S/sqrt(S+B)
    ftxt << "  \\\\" << endl;
  }
  ftxt.close();
  return;
}
