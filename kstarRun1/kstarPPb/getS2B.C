void getS2B(TString outname = "table_s2b.tmp",
	    Float_t nsigma=5.0,
	    TString file_fit_name = "/Users/bellini/alice/resonances/myKstar/pwglf_train_out/train81/roofit/antikstar_BWPS3_31may13.root",
	    TString file_distr_name = "/Users/bellini/alice/resonances/myKstar/pwglf_train_out/train81/_antikstar_EMnorm1.30-1.50_train81.root",
	    Bool_t isVerbose = 1)
{
  ofstream ftxt;
  outname.ReplaceAll(".tmp",Form("_%2.1fs.tmp",nsigma));
  ftxt.open(outname.Data(),ios::out);
  ftxt <<"Fit file: "<<file_fit_name.Data() << endl;
  ftxt <<"Histo file: "<< file_distr_name.Data() << endl;
  ftxt << "cent    pt    signal    signalErr      Bg     ResBg    S/B   S/resB   S/sqrt(S+B)    S/sqrt(S+resB)"<<endl;
  for (Int_t ic = 0;ic<4;ic++){
    for (Int_t ip=1; ip<11;ip++){      
      Double_t result[8]={0.,0.,0.,0.,0.,0.,0.,0.};
      getS2B(result,
	     ic, 
	     ip, 
	     nsigma,
	     file_fit_name.Data(),
	     file_distr_name.Data(),
	     isVerbose);
      ftxt.precision(2);
      ftxt  << ic << "     " << ip << "   ";  // pt bin
      ftxt.precision(6);
      ftxt << result[0]<<"    " << result[1] << "    " << result[2] << "    " << result[3] << "    "<< result[4]<<"    "<< result[5]<<"    "<< result[6]<<"    "<< result[7] << endl; 
    }
  } 
  ftxt.close();
  return;
}

//---------------------------------------------
Double_t getS2B(Double_t * result,
		Int_t centbin = 0, 
	  	Int_t ptbin = 1, 
		Float_t nsigma=5.0,
		TString file_fit_name = "/Users/bellini/alice/resonances/myKstar/pwglf_train_out/train81/roofit/kstar_BWPS3_31may13.root",
		TString file_distr_name = "/Users/bellini/alice/resonances/myKstar/pwglf_train_out/train81/_kstar_EMnorm1.30-1.50_train81.root",
		Bool_t isVerbose = kTRUE)
{ 

  if (!result) return 0.0;

  //from fit to S+res.bg
  Double_t treesig, sig = 0.0;      //integral of signal from fit
  Double_t treesigerr, sigErr = 0.0;   //error on integral of signal from fit
  // Double_t resbg = 0.0;    //integral of res. bg from fit
  // Double_t resbgErr = 0.0; //error on integral of res. bg from fit
  Int_t centBinID,ptBinID;

  //from histos
  Double_t embg = 0.0;     //integral of em bg histo
  Double_t usint = 0.0;    //integral of US distrib
  Double_t subint = 0.0;   // integral of em-bg subtracted US distrib.
  Double_t bg = 0.0;       //integral of US distrib - integral of signal from fit
  Double_t resbg = 0.0;   // integral of em-bg subtracted US distrib. - integral of signal from fit

  //Get fit result from fit output file
  TFile * file_fit = TFile::Open(file_fit_name.Data());
  TTree * tree = (TTree*) file_fit->Get("tree");
  if (!tree) return;
  tree->SetBranchAddress("nSignal",&treesig);
  tree->SetBranchAddress("nSignalErr",&treesigerr);
  tree->SetBranchAddress("centBin",&centBinID);
  tree->SetBranchAddress("ptBin",&ptBinID);
  for (Int_t ientry=0;ientry<tree->GetEntries();ientry++){
    tree->GetEntry(ientry);
    if ((centBinID==centbin) && (ptBinID==ptbin)) {
      sig = treesig;
      sigErr = treesigerr;
    }
  }
  
  //get histo from sub* file
  TFile * file_distr = TFile::Open(file_distr_name.Data());
  TH1D * hus = (TH1D *) file_distr->Get(Form("Signal_ptBin%02i_centBin%02i", ptbin, centbin));
  TH1D * hembg = (TH1D *) file_distr->Get(Form("norm_Mixing_ptBin%02i_centBin%02i", ptbin, centbin));
  TH1D * hsubint = (TH1D *) file_distr->Get(Form("sub_norm_Mixing_ptBin%02i_centBin%02i", ptbin, centbin));
  
  //Get x_low and x_up, aka range for integrals
  Double_t xlow = 0.89595-nsigma*(0.0505)/2.35;
  Double_t xup =  0.89595+nsigma*(0.0505)/2.35;
  Int_t ilow = ((TAxis*) hus->GetXaxis())->FindBin(xlow);
  Int_t iup = ((TAxis*) hus->GetXaxis())->FindBin(xup);
  
  //Get integrals
  embg = hembg->Integral(ilow,iup);
  usint = hus->Integral(ilow,iup);
  subint = hsubint->Integral(ilow,iup);
  
  //Get Bg integrals
  bg = usint-sig;
  resbg = subint-sig;
  
  Double_t s2b=0.0;
  Double_t s2resb=0.0;
  Double_t s2sqrt=0.0;
  Double_t s2sqrtres=0.0;
  
  if (nsigma==1.0) sig = sig*0.683;
  if (nsigma==2.0) sig = sig*0.9545;
  if (nsigma==3.0) sig = sig*0.997;

  if (bg>0) s2b = sig/bg;
  else Printf("Warning: US-sig <=0!");

  if (resbg>0) s2resb = sig/resbg;
  else Printf("Warning: sub-sig <=0!");
  
  if ((usint)>0) s2sqrt = sig/TMath::Sqrt(usint);
  else Printf("Warning: bg+sig <=0!");

  if ((usint)>0) s2sqrtres = sig/TMath::Sqrt(subint);
  else Printf("Warning: resbg+sig <=0!");

  if (isVerbose){
    Printf("==============================================\n========================================== S/B \n Open file with distributions: %s \n Open file with fit result: %s",file_distr_name.Data(),file_fit_name.Data() );
    Printf("Integral range: %5.4f (%02i) < M < %5.4f (%i)", xlow, ilow, xup, iup);
    Printf("-------------------------------------------------------------------------------------");
    Printf (" US integral = %e \n EM integral = %e \n Sub integral = %e", usint, embg, subint);
    Printf(" Bg integral = %e \n Res. Bg integral = %e \n Signal = %e", bg, resbg, sig);
    Printf("-------------------------------------------------------------------------------------");
    Printf(" S/B = %e \n S/resB = %e", s2b, s2resb);
    Printf(" S/Sqrt(S+B) = %e \n S/sqrt(S+resB) = %e", s2sqrt, s2sqrtres);
  }
  result[0] = sig;
  result[1] = sigErr;
  result[2] = bg;
  result[3] = resbg;
  result[4] = s2b;
  result[5] = s2resb;
  result[6] = s2sqrt;
  result[7] = s2sqrtres;
  
  return s2b;
}
