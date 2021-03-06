#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit;

enum EFitFunction{ kPOLY2,
		   kPOLY3,
		   kLandau,
		   kPOLY1,
		   kEXP,
		   kData};

Color_t color[]={kRed, kOrange, kGreen+2, kBlue, kMagenta, kBlack};
Color_t colorfit[]={kBlue, kRed, kGreen+2, kMagenta+1, kYellow+2, kBlack};
Int_t marker[]={20, 21, 22, 23, 28};

enum EbackgndType{ kUnlikeSE, kUnlikeME,  kLikeSE, knBackgndTypes};

// retrieve reference values in the database PDG
Int_t PDG = 313;
TDatabasePDG *pdg     = TDatabasePDG::Instance();
TParticlePDG *part    = pdg->GetParticle(PDG);
Double_t      pdgMass    = part->Mass(); // const Double_t pdgMass = 0.89594;
Double_t      pdgWidth   = part->Width(); // 0.0487;// const Double_t pdgWidth = 0.0487;

Double_t xrangem, xrangeM;
Float_t  nWidthRangeFit;
Double_t fitCutMassMin,fitCutMassMax,fitCutWidthMin, fitCutWidthMax;

//gROOT->SetStyle("Plain");
gSystem->Load("libRooFit");


//--------------------------------------------------
void fitKStar(TString infilename="sub_aod49_kstar.root", 
	      Bool_t isRebin=kFALSE,
	      Int_t backgndType=EbackgndType::kUnlikeME, 
	      TString function="POLY2",
	      Int_t selCentBin=-1,
	      Int_t istartptbin=2, 
	      Int_t istopptbin= 14, 
	      Int_t canvasSplit = 4, 
	      Double_t massRangeMin = 0.74, 
	      Double_t massRangeMax = 1.1, 
	      Bool_t enaFitResPave=1,
	      Bool_t printAll = 0,
	      Float_t nsigmaPeak= 7.0, 
	      Float_t cutMassMin=-1.0,
	      Float_t cutMassMax=-1.0,
	      Float_t cutWidthMin=-1.0,
	      Float_t cutWidthMax=-1.0
	      )
{
  
  gROOT->LoadMacro("/Users/bellini/alice/macro/SetGraphicStyle.C");
  SetGraphicStyle();
  TGaxis::SetMaxDigits(3);
  // //open input file
  if (!infilename){
    printf("ERROR: invalid file name provided. \n");
    return;
  }
  TFile* fin=TFile::Open(infilename.Data());
  if (!fin) {
    printf("ERROR: invalid or corrupted input file. Exiting.\n");
    return;
  }

  //get bins
  TAxis *ptbins = (TAxis*)fin->Get("ptbins");
  Int_t npt = ptbins->GetNbins();
  TAxis *centbins = (TAxis*)fin->Get("centbins");
  Int_t ncent = centbins->GetNbins();

  const Int_t dimpt = npt+1;
  Double_t pt[dimpt];
  for (Int_t k=0; k<dimpt;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
    Printf("%5.2f",pt[k]);
  }
  const Int_t dimcent = ncent+1;
  Double_t cent[dimcent]; 
  for (Int_t k=0; k<dimcent;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
    Printf("%5.2f",cent[k]);
  }
  fin->Close();

  Int_t nbinstodisplay = istopptbin-istartptbin;

  //define output
  TString fileout="";
  if (backgndType==EbackgndType::kUnlikeME){
    //fileout.Form("roofit/%sfitEM_peak%2.1f_%s_%3.2f-%3.2f_%s",(isRebin? "rebin_" : ""), nsigmaPeak,function.Data(),massRangeMin,massRangeMax,infilename.Data());
    fileout.Form("roofit/%sfitEM_%s_%3.2f-%3.2f_%s",(isRebin? "rebin_" : ""),function.Data(),massRangeMin,massRangeMax,infilename.Data());
    fileout.ReplaceAll("sub_","");
  } 
  if (backgndType==EbackgndType::kLikeSE){
    //fileout.Form("roofit/%sfitLS_peak%2.1f_%s_%3.2f-%3.2f_%s",(isRebin? "rebin_" : ""), nsigmaPeak,function.Data(),massRangeMin,massRangeMax,infilename.Data());
    fileout.Form("roofit/%sfitLS_%s_%3.2f-%3.2f_%s",(isRebin? "rebin_" : ""),function.Data(),massRangeMin,massRangeMax,infilename.Data());
    fileout.ReplaceAll("sub_","");
 } 

  if (selCentBin>-1) fileout.ReplaceAll(".root",Form("_cent%i.root",selCentBin));
  if (nbinstodisplay==1) fileout.ReplaceAll(".root",Form("_Pt%i.root",istartptbin));
  if (fileout.Contains("../")) fileout.ReplaceAll("../","");
  TFile * fout=new TFile(fileout.Data(),"recreate");
  Double_t fitParams[9], SoverB=0.0, significance=0.0, normfactorCopy;
  Int_t centBinID,ptBinID, funcID;
  Float_t rangeInfCopy, rangeSupCopy, ptinfCopy, ptsupCopy, centinfCopy, centsupCopy, massRangeInfCopy,massRangeSupCopy;
  
  TString fitfunction;
  TTree *tree=new TTree("tree","fit parameters tree");
  tree->Branch("signalMass",&fitParams[0],"signalMass/D");
  tree->Branch("signalMassErr",&fitParams[1],"signalMassErr/D");
  tree->Branch("signalWidth",&fitParams[2],"signalWidth/D");
  tree->Branch("signalWidthErr",&fitParams[3],"signalWidthErr/D");
  tree->Branch("nSignal",&fitParams[4],"nSignal/D");
  tree->Branch("nSignalErr",&fitParams[5],"nSignalErr/D");
  tree->Branch("nBack",&fitParams[6],"nBack/D");
  tree->Branch("nBackErr",&fitParams[7],"nBackErr/D");
  tree->Branch("chi2",&fitParams[8],"chi2/D");  
  tree->Branch("centBin",&centBinID,"centBin/I");
  tree->Branch("ptBin",&ptBinID,"ptBin/I");
  tree->Branch("SoverB",&SoverB,"SoverB/D");
  tree->Branch("significance",&significance,"significance/D");
  tree->Branch("norm_factor", &normfactorCopy, "norm_factor/D");
  tree->Branch("norm_inf", &rangeInfCopy, "norm_inf/F");
  tree->Branch("norm_sup", &rangeSupCopy, "norm_sup/F");
  tree->Branch("pt_inf", &ptinfCopy, "pt_inf/F");
  tree->Branch("pt_sup", &ptsupCopy, "pt_sup/F");
  tree->Branch("cent_inf", &centinfCopy, "cent_inf/F");
  tree->Branch("cent_sup", &centsupCopy, "cent_sup/F");
  tree->Branch("function", &fitfunction, "function", 1000);
  tree->Branch("functionID", &funcID, "functionID/I");
  tree->Branch("fitrange_inf", &massRangeInfCopy, "fitrange_inf/F");
  tree->Branch("fitrange_sup", &massRangeSupCopy, "fitrange_sup/F");
 

  TString centBinLabel;
  if (istopptbin==-1) istopptbin=npt;
  RooPlot *plotframe;
  TLegend * leg = new TLegend(0.6,0.55,0.88,0.88);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetBorderSize(0.0);

  //set fit range in tree
  massRangeInfCopy=massRangeMin; 
  massRangeSupCopy=massRangeMax;
  //set function in tree
  fitfunction=function;
  if (fitfunction.Contains("POLY2")) funcID=EFitFunction::kPOLY2;
  if (fitfunction.Contains("POLY3")) funcID=EFitFunction::kPOLY3;
  if (fitfunction.Contains("LAND")) funcID=EFitFunction::kLandau;
  if (fitfunction.Contains("EXP")) funcID=EFitFunction::kEXP;
  if (fitfunction.Contains("POLY1")) funcID=EFitFunction::kPOLY1;
  //loop on bins
  for (Int_t icentbin=0; icentbin<ncent;icentbin++){   
    if ( (selCentBin>-1) && (!(icentbin==selCentBin)) ) continue;
    
    char cutstring[100];
    sprintf(cutstring,"");
    centBinID=icentbin;
    Int_t canvasDim = 400;
    if (canvasSplit<4) canvasDim = 600;
    TCanvas *cfit=new TCanvas(Form("cfit_%i",icentbin),Form("cfit_%i",icentbin), canvasDim*canvasSplit,(canvasDim*0.7*nbinstodisplay)/canvasSplit);
    cfit->Divide(canvasSplit,(nbinstodisplay)/canvasSplit);
 
    //define histos
    //TH2F*hMassVsPt2=new TH2F(Form("hMassVsPt_%i",icentbin),Form("K* mass (%2.0f-%2.0f %); p_{t} (GeV/c); M (GeV/c^{2})", cent[icentbin], cent[icentbin+1]), npt, pt, 350, 0.75, 1.1);
    TH1F*hMassVsPt=new TH1F(Form("hMassVsPt_%i",icentbin),Form("K* mass (%2.0f-%2.0f %); p_{t} (GeV/c); M (GeV/c^{2})",cent[icentbin],cent[icentbin+1]), npt, pt);

    for (Int_t iptbin=istartptbin; iptbin<istopptbin;iptbin++){   
      //open input file
      if (!infilename){
	printf("ERROR: invalid file name provided. \n");
	return;
      }
      TFile* fin=TFile::Open(infilename.Data());
      if (!fin) {
	printf("ERROR: invalid or corrupted input file. Exiting.\n");
	return;
      }
      
      //read normalization tree
      Float_t rangeInf, rangeSup, ptinf, ptsup, centinf, centsup;
      Double_t normfactor;
      Int_t treept, treecent;
      TTree *ntree = (TTree*)fin->Get("ntree");
      ntree->SetBranchAddress("factor", &normfactor[0]);
      ntree->SetBranchAddress("inf", &rangeInf);
      ntree->SetBranchAddress("sup", &rangeSup);
      ntree->SetBranchAddress("pt_inf", &ptinf);
      ntree->SetBranchAddress("pt_sup", &ptsup);
      ntree->SetBranchAddress("pt_bin", &treept);
      ntree->SetBranchAddress("cent_inf", &centinf);
      ntree->SetBranchAddress("cent_sup", &centsup);
      ntree->SetBranchAddress("cent_bin", &treecent);
      
   
      cfit->cd(iptbin-istartptbin+1);
      ptBinID=iptbin;
      
      //loop on normalization tree to get normalization info and copy in fit tree
      for (Int_t ientry=0;ientry<ntree->GetEntries();ientry++){
	ntree->GetEntry(ientry);
	if(centBinID==treecent && ptBinID==treept){
	  rangeInfCopy=rangeInf; 
	  rangeSupCopy=rangeSup;
	  ptinfCopy=ptinf;
	  ptsupCopy=ptsup;
	  centinfCopy=centinf; 
	  centsupCopy=centsup;
	  normfactorCopy=normfactor;
	}
      }      
      
      //draw fit result
      TString testoback="";
      TString hSignalName="";
      TString testopt;
      testopt.Form(" %3.1f#leqp_{t}#leq%3.1f GeV/c ",pt[iptbin],pt[iptbin+1]);
      TString testocent;
      testocent.Form("(%2.0f-%2.0f%%)",cent[icentbin],cent[icentbin+1]);
      
      if (backgndType==EbackgndType::kUnlikeME){
	hSignalName.Form("sub_norm_Mixing_ptBin%02i_centBin%02i",iptbin,icentbin);
	testoback.Form("(EM backgnd)");
      }
      if (backgndType==EbackgndType::kLikeSE){
	hSignalName.Form("sub_norm_Like_ptBin%02i_centBin%02i",iptbin,icentbin);
	testoback.Form("(LS backgnd)");
      }
      TH1D *hnsig = (TH1D*)fin->Get(hSignalName.Data());   
      if (!hnsig){
	printf("fitKstar - ERROR: input histogram %s not found in file %s . returning...\n", hSignalName.Data(),infilename.Data() );
	return;
      }
      plotframe=(RooPlot*)fitKStar(hnsig,
				   fitParams,
				   function.Data(),
				   massRangeMin, 
				   massRangeMax, 
				   nsigmaPeak, 
				   cutMassMin, 
				   cutMassMax, 
				   cutWidthMin, 
				   cutWidthMax,
				   isRebin);
      //draw legend
      leg->Clear();
      Char_t names[10][100];
      for (int i=0; i<plotframe->numItems(); i++) {
	TString obj_name=plotframe->nameOf(i); 
	if (obj_name=="") continue;
	cout << Form("%d. '%s'\n",i,obj_name.Data());
	sprintf(names[i], "%s",obj_name.Data());
      }
      Char_t caption[10][30] = {
	"Data",
	"BW+poly2",
	"poly2",
	"BW+poly3",
	"poly3",
	"BW+Landau",
	"Landau",
	"fitted signal only",
	"fitted signal only",
       	"fitted signal only"
      };
      leg->AddEntry(hnsig,"data","lp");
      for (Int_t i=1;i<3;i++) {
	TObject *obj = plotframe->findObject(names[i]);
	if (!obj) {
	  Warning("fitKstar",Form("Can't find item='%s' in the frame2!\n",names[i]));
	  continue;
	}
	Printf("%s",caption);
	leg->AddEntry(obj,caption[i],"lp");
      }
      //============== end Draw legend
      TString frameTitle="";//K*+#bar{K*}";
      //frameTitle.Append(testoback);
      frameTitle.Append(testopt);
      frameTitle.Append(testocent);
      plotframe->SetTitle(frameTitle.Data());
      plotframe->SetTitleOffset(0.7);
      if (isRebin) plotframe->SetYTitle("counts / (20 MeV/c^{2})");
	else {
	  plotframe->SetYTitle("counts / (10 MeV/c^{2})");
	    }
      TPaveText *textFit = 0x0;
      if (printAll){
	textFit = new TPaveText(0.55,0.5,0.95,0.87,"NDC");
	textFit->SetBorderSize(1);
	textFit->SetFillColor(kWhite);
      } else {
	if (selCentBin<0) textFit = new TPaveText(0.55,0.65,0.90,0.89,"NDC");
        if (selCentBin==0) textFit = new TPaveText(0.12,0.13,0.75,0.43,"NDC");
        if (selCentBin==1) textFit = new TPaveText(0.12,0.13,0.65,0.38,"NDC");
        if (selCentBin>=2) textFit = new TPaveText(0.55,0.6,0.89,0.89,"NDC");
	//	if (selCentBin==3) textFit = new TPaveText(0.12,0.13,0.75,0.43,"NDC");
        textFit->SetBorderSize(0);
	textFit->SetFillStyle(0);
	textFit->SetTextAlign(12); //32 = vert middle, horiz right
      } 
      if (printAll)textFit->AddText(Form("%s, %4.2f#leqM_{inv}#leq%4.2f",function.Data(),massRangeMin,massRangeMax));
      textFit->AddText(Form("M(GeV/c^{2}) =  %6.4f #pm %6.4f",fitParams[0],  fitParams[1]));
      if (printAll) textFit->AddText(Form("#Gamma(K*) =  %6.4f #pm %6.4f GeV/c^{2} \n",fitParams[2],  fitParams[3]));
      textFit->AddText(Form("N_{raw} =  %8.0f #pm %8.0f \n",fitParams[4],  fitParams[5]));
      textFit->AddText(Form("N_{Bg} = %8.0f #pm %8.0f \n",fitParams[6],  fitParams[7]));
      textFit->AddText(Form("#chi^{2} = %6.4f \n",fitParams[8]));
      textFit->SetTextFont(42);
      textFit->SetTextColor(kBlack);
      // if (function.Contains("POLY1")) textFit->SetTextColor(colorfit[EFitFunction::kPOLY1]);
      // if (function.Contains("POLY3")) textFit->SetTextColor(colorfit[EFitFunction::kPOLY3]);
      // if (function.Contains("EXP")) textFit->SetTextColor(colorfit[EFitFunction::kEXP]);
      // if (function.Contains("POLY2")) textFit->SetTextColor(colorfit[EFitFunction::kPOLY2]);
      // if (function.Contains("LAND")) textFit->SetTextColor(colorfit[EFitFunction::kLandau]);
      SoverB = fitParams[4]/fitParams[6];
      significance = fitParams[4]/TMath::Sqrt(fitParams[4]+fitParams[6]);
      //---------- end new
     
      //plotframe=(RooPlot*)fitKStar(iptbin,icentbin,infilename,backgndType,fitParams);
      //hMassVsPt2->Fill(icentbin,fitParams[0]);
      hMassVsPt->SetBinContent(iptbin+1,fitParams[0]);
      hMassVsPt->SetBinError(iptbin+1,fitParams[1]);
      hnsig->SetMarkerSize(0.7);
      hnsig->SetMarkerColor(colorfit[EFitFunction::kData]);
      hnsig->SetLineColor(colorfit[EFitFunction::kData]);
      hnsig->SetLineWidth(1);
      hnsig->GetXaxis()->SetRangeUser(0.6,1.3);
      hnsig->GetXaxis()->SetTitleFont(42);
      hnsig->GetXaxis()->SetTitleOffset(1.2);
      
      //Draw in canvas
      if (isRebin) plotframe->GetYaxis()->SetTitle("counts / (20 MeV/c^{2})");//RangeUser(0.01,7e5);
      else  plotframe->GetYaxis()->SetTitle("counts / (10 MeV/c^{2})");
      plotframe->GetYaxis()->SetTitleOffset(1.0);//RangeUser(0.01,7e5);
      plotframe->GetYaxis()->SetTitleSize(0.05);//RangeUser(0.01,7e5);
      plotframe->GetYaxis()->SetLabelSize(0.045);//RangeUser(0.01,7e5);
      plotframe->GetXaxis()->SetTitleOffset(0.9);//RangeUser(0.01,7e5);
      plotframe->GetXaxis()->SetLabelOffset(0.006);//RangeUser(0.01,7e5);
      plotframe->GetXaxis()->SetTitleSize(0.05);//RangeUser(0.01,7e5);
      plotframe->GetXaxis()->SetLabelSize(0.045);//RangeUser(0.01,7e5);
      plotframe->Draw(); 
      
      //leg->Draw("");
      hnsig->Draw("same");
      if (enaFitResPave) textFit->Draw("same");
      tree->Fill();
    }//end loop on pt
    hMassVsPt->SetMarkerStyle(marker[icentbin]);
    hMassVsPt->SetMarkerColor(color[icentbin]);
    hMassVsPt->SetLineColor(color[icentbin]);                 
    fout->cd();
    //hMassVsPt->Write();
    hMassVsPt->Write();
    if (nbinstodisplay==1) cfit->SaveAs(Form("roofit/fit%s_%s_range%3.2f-%3.2f_cent%i_Pt%i.png", (backgndType==EbackgndType::kUnlikeME ? "EM" : "LS"), function.Data(), massRangeMin,massRangeMax,icentbin,istartptbin));
    else cfit->SaveAs(Form("roofit/fit%s_%s_range%3.2f-%3.2f_canvas%ix%i_cent%i_startPt%i.png", (backgndType==EbackgndType::kUnlikeME ? "EM" : "LS"), function.Data(), massRangeMin,massRangeMax,canvasSplit,(nbinstodisplay)/canvasSplit,icentbin,istartptbin));
    
  }//END LOOP on cent
  fout->cd();
  tree->Write();
  centbins->Write();
  ptbins->Write();
  fout->Close();
    
  // printf("**************************************************************\n************************ SUMMARY *****************************\n
// Resonance PDG params:
//                    pdgMass = %6.5f \n
//                    pdgWidth = %6.5f \n
// range of S+B fit:
//                    xrangem = %6.2f \n
//                    xrangeM = %6.2f \n
// range of S mass fit:
//                    nWidthRangeFit = %6.3f \n
//                    fitCutMassMin  = %6.3f \n
//                    fitCutMassMax  = %6.3f \n
//                    fitCutWidthMin = %6.3f \n
//                    fitCutWidthMax = %6.3f \n ", 
// 	 pdgMass, pdgWidth, xrangem, xrangeM, nWidthRangeFit,
// 	 fitCutMassMin, fitCutMassMax, fitCutWidthMin,fitCutWidthMax);
  return;
}


//------------------------------------------------------------------------------------------------
RooPlot* fitKStar(TH1D*hnsig, 
		  Double_t* fitParams, 
		  TString function, 
		  Double_t massRangeMin = 0.0, 
		  Double_t massRangeMax=0.0, 
		  Float_t nsigmaPeak= 7.0, 
		  Float_t cutMassMin=-1.0,
		  Float_t cutMassMax=-1.0,
		  Float_t cutWidthMin=-1.0,
		  Float_t cutWidthMax=-1.0,
		  Bool_t isRebin=kFALSE
		  )
{
  /*
    performs fit on K* plot histo
    where it is assumed that already (norm+)subtracted distribution is given as input
  */

  if (!hnsig){
    printf("fitKstar - ERROR: input histogram %s not found in file . returning...\n", hSignalName.Data() );
    return;
  }
  Double_t histo_integral = hnsig->Integral();
  if (isRebin) hnsig->Rebin(2);
  SetFitParams(massRangeMin, massRangeMax, nsigmaPeak, cutMassMin, cutMassMax, cutWidthMin, cutWidthMax);
  
  Double_t signalMass,signalMassErr,
    signalWidth,signalWidthErr,
    chi2;
  Int_t nrange = 0;
  Double_t nBack=0;
  Double_t nBackErr=0;
  Double_t nSignal=0;
  Double_t nSignalErr=0;
  Double_t showxrangem = 0.7,showxrangeM = 1.3;
  
  RooRealVar x("m","M (GeV/c^{2})", xrangem, xrangeM);//showxrangem, showxrangeM
  //RooRealVar width("width","width",pdgWidth,fitCutWidthMin,fitCutWidthMax);
  //RooRealVar width("width","width",pdgWidth,pdgWidth*0.99,1.5*pdgWidth);
  RooRealVar width("width","width",pdgWidth,pdgWidth,pdgWidth);
  RooRealVar mean("mean","mean",pdgMass,fitCutMassMin,fitCutMassMax);
  RooBreitWigner breit("breit","signal Breit-Wigner",x,mean,width);
  RooDataHist data("data","data",RooArgList(x),hnsig);

  /*-----------------------------------------------------
    residual background 1 - Chebychev polynomials 
  -----------------------------------------------------*/
  
  //Poly 3
   RooRealVar c0("c0","coefficient #0", -1.0, 1.0) ;
   RooRealVar c1("c1","coefficient #1", -0.1, 0.0) ; //-0.5, 0.0 poly2
   RooRealVar c2("c2","coefficient #2",  0.0, 1.0) ;
  //Poly 2
  // RooRealVar a0("a0","coefficient #0", -1.0, 0.0) ; //-1,1
  // RooRealVar a1("a1","coefficient #1", -0.1, 0.1) ; //-0.5, 0.0 poly2
  RooRealVar a0("a0","coefficient #0", -1.0, 0.0) ; //-1,1
  RooRealVar a1("a1","coefficient #1", -0.5, 0.5) ; //-0.5, 0.0 poly2
      
  // RooChebychev cheby1("cheby","residual background Chebychev poli", x, RooArgList(c0));
   RooChebychev cheby2("cheby","residual background Chebychev poli", x, RooArgList(a0,a1));
   RooChebychev cheby3("cheby","residual background Chebychev poli", x, RooArgList(c0,c1,c2)); 
  
  
  /* new 15/12/12 */

  // //Poly 2 f(x)= 1 + p1*x + p2*x*x
  // RooRealVar a1("a1","coefficient #1", -2., 2.) ; //-0.5, 0.0 poly2
  // RooRealVar a2("a2","coefficient #2", -0.5.,0.) ;
  // //Poly 3 = f(x)= 1 + p1*x + p2*x*x + p3*x*x*x 
  // RooRealVar p1("p1","coefficient #1", -1.0, 1.0) ; //-0.5, 0.0 poly2
  // RooRealVar p2("p0","coefficient #0", -1., 1.) ;
  // RooRealVar p3("p1","coefficient #1", -1., 1.0) ; //-0.5, 0.0 poly2
  
  //solo per 60-80
  // RooRealVar a0("a0","coefficient #0", -2, 2.0) ; //-1,1
  // RooRealVar a1("a1","coefficient #1", -10., 10.) ; //-0.5, 0.0 poly2
  
  RooPolynomial cheby1("cheby","residual background Chebychev poli", x, RooArgList(a0));
  //RooPolynomial cheby2("cheby","residual background Chebychev poli", x, RooArgList(a0,a1));
  // RooPolynomial cheby3("cheby","residual background Chebychev poli", x, RooArgList(p1,p2,p3)); 
  
  
  /*-----------------------------------------------------
    residual background 2 - Landau
    -----------------------------------------------------*/
  RooRealVar sigmaL("sigmaL","Landau sigma",50., 0., 100.);
  RooRealVar meanL("meanL","Landau mean", 0.72, 0.1,0.83);  
  RooLandau landau("landau","residual background Landau",x,meanL,sigmaL);
  
  /*-----------------------------------------------------
    residual background 3 - exponential 
    -----------------------------------------------------*/
  RooRealVar pare("pare","par expo ",-35., -50., -0.1);
  RooExponential expo("expo","residual background exponential",x,pare);
  
  /* method 1: use fraction of signal w.r.t. S+B */
  // model(x) = fsig*sig(x) + (1-fsig)*bkg(x)
  //RooRealVar fsig(“fsig","signal fraction",0.5,0.,1.) ;
  //RooAddPdf model(“model","model",RooArgList(breit,cheby),fsig) ;

  /* method 2: use counts of signal and backgnd */
  RooRealVar nsig("nsig","signal fraction",histo_integral*0.5, 0., histo_integral) ;
  RooRealVar nbkg("nbkg","background fraction",histo_integral*0.5, 0., histo_integral) ;

  RooAddPdf model("model","model",RooArgList(breit,cheby2),RooArgList(nsig,nbkg)) ;
  RooAddPdf model1("model1","model1",RooArgList(breit,cheby1),RooArgList(nsig,nbkg)) ;
  RooAddPdf model3("model3","model3",RooArgList(breit,cheby3),RooArgList(nsig,nbkg)) ;
  RooAddPdf modelL("modelL","modelL",RooArgList(breit,landau),RooArgList(nsig,nbkg)) ;
  RooAddPdf modelE("modelE","modelE",RooArgList(breit,expo),RooArgList(nsig,nbkg)) ;
  
  //prin meaning of parameters
  RooArgSet* params = model->getVariables() ;
  params->Print("v") ;
  
  //Display data
  RooPlot* xframe = x.frame(0.6,1.3);//massRangeMin, massRangeMax);
  data.plotOn(xframe, Name("data"), LineColor(colorfit[EFitFunction::kData]), MarkerColor(colorfit[EFitFunction::kData]), MarkerSize(0.2), LineWidth(2), DrawOption("E1X0")) ;//nota: add Name("data") to get chi2
  
  //Perform fit 
  if (function.Contains("POLY1")){
    model1.fitTo(data, Extended(kTRUE), SumW2Error(kFALSE), Save());/* Range(xrangem, xrangeM),*/
    model1.plotOn(xframe, Name("model1"), LineColor(colorfit[EFitFunction::kPOLY1]),LineWidth(2),Range(xrangem, xrangeM));//nota: add Name("model") to get chi2
    //model1.plotOn(xframe, Components(breit),LineColor(colorfit[EFitFunction::kPOLY1]),LineWidth(1),Range(xrangem, xrangeM)) ;
    model1.plotOn(xframe, Components(cheby1),LineStyle(kDashed),LineColor(colorfit[EFitFunction::kPOLY1]),LineWidth(2),Range(xrangem, xrangeM)) ;
    chi2=xframe->chiSquare("model1", "data", 3);
  }
  
  if (function.Contains("POLY3")){
    model3.fitTo(data, Extended(kTRUE), SumW2Error(kFALSE), Save());/* Range(xrangem, xrangeM),*/
    model3.plotOn(xframe, Name("model3"), LineColor(colorfit[EFitFunction::kPOLY3]),LineWidth(2),Range(xrangem, xrangeM));//nota: add Name("model") to get chi2
    model3.plotOn(xframe, Components(cheby3), LineStyle(kDashed),LineColor(colorfit[EFitFunction::kPOLY3]),LineWidth(2),Range(xrangem, xrangeM)) ;
    //model3.plotOn(xframe, Components(breit), LineColor(colorfit[EFitFunction::kPOLY3]),LineWidth(1),Range(xrangem, xrangeM)) ;
    chi2=xframe->chiSquare("model3", "data", 3);
  }
   
  if (function.Contains("EXP")){
    modelE.fitTo(data, Extended(kTRUE), SumW2Error(kFALSE), Save());/* Range(xrangem, xrangeM),*/
    modelE.plotOn(xframe, Name("modelE"), LineColor(colorfit[EFitFunction::kEXP]),LineWidth(2),Range(xrangem, xrangeM));//nota: add Name("model") to get chi2
    modelE.plotOn(xframe, Components(expo), LineStyle(kDashed),LineColor(colorfit[EFitFunction::kEXP]),LineWidth(2),Range(xrangem, xrangeM)) ;
    modelE.plotOn(xframe, Components(breit), LineColor(colorfit[EFitFunction::kEXP]),LineWidth(1),Range(xrangem, xrangeM)) ;
    chi2=xframe->chiSquare("modelE", "data", 3);
  }

  if (function.Contains("POLY2")){
    model.fitTo(data, Extended(kTRUE), SumW2Error(kFALSE));/* Range(xrangem, xrangeM),*/
    model.plotOn(xframe, Name("model"), LineColor(colorfit[EFitFunction::kPOLY2]),LineWidth(2),Range(xrangem, xrangeM));//nota: add Name("model") to get chi2
     model.plotOn(xframe, Components(cheby2),LineStyle(kDashed),LineColor(colorfit[EFitFunction::kPOLY2]),LineWidth(2),Range(xrangem, xrangeM)) ;
     //model.plotOn(xframe, Components(breit),LineColor(colorfit[EFitFunction::kPOLY2]),LineWidth(1),Range(xrangem, xrangeM)) ;
    chi2=xframe->chiSquare("model", "data", 3);
  }

  if (function.Contains("LAND")){
    modelL.fitTo(data, Extended(kTRUE), SumW2Error(kFALSE));/* Range(xrangem, xrangeM),*/  
    modelL.plotOn(xframe, Name("modelL"), LineColor(colorfit[EFitFunction::kLandau]),LineWidth(2),Range(xrangem, xrangeM));//nota: add Name("model") to get chi2
		  modelL.plotOn(xframe, Components(landau),LineStyle(kDashed),LineColor(colorfit[EFitFunction::kLandau]),LineWidth(2),Range(xrangem, xrangeM)) ;
		  //modelL.plotOn(xframe, Components(breit),LineColor(colorfit[EFitFunction::kLandau]),LineWidth(1),Range(xrangem, xrangeM)) ;
    chi2=xframe->chiSquare("modelL", "data", 3);
  }


  signalMass=mean.getVal();
  signalMassErr=mean.getError();
  signalWidth=width.getVal() ;
  signalWidthErr=width.getError() ;
  nSignal=nsig.getVal();
  nSignalErr=nsig.getError();
  nBack=nbkg.getVal();
  nBackErr=nbkg.getError();
    
  if (fitParams){
    fitParams[0]=signalMass;
    fitParams[1]=signalMassErr;
    fitParams[2]=signalWidth;
    fitParams[3]=signalWidthErr;
    fitParams[4]=nSignal;
    fitParams[5]=nSignalErr;
    fitParams[6]=nBack;
    fitParams[7]=nBackErr;
    fitParams[8]=chi2;
  }
  
//   printf("**************************************************************\n
// ************************ FIT RESULT *****************************\n
// Input Resonance PDG params:
// - pdgMass = %6.5f 
// - pdgWidth = %6.5f \n
// Range of S+B fit: %6.4f <= Minv <= %6.4f GeV/c^2\n
// Fit settings:
// - nWidthRangeFit = %4.2f \n
// - %6.4f < M < %6.4f GeV/c^2 \n
// - %6.4f < W < %6.4f GeV/c^2\n 
// \n ************ Fit result
// function : %s
// M(K*) =  %6.4f +/- %6.4f GeV/c^2 \n
// W(K*) =  %6.4f +/- %6.4f GeV/c^2 \n
// raw counts =  %e +/- %e \n
// N_Bg = %e +/- %e \n
// Chi2 = %6.4f \n", 
// 	 pdgMass, pdgWidth, xrangem, xrangeM, nWidthRangeFit,
// 	 fitCutMassMin, fitCutMassMax, fitCutWidthMin,fitCutWidthMax,
// 	 function.Data(),
// 	 fitParams[0],  fitParams[1],  fitParams[2],  fitParams[3],
// 	 fitParams[4],  fitParams[5],  fitParams[6],  fitParams[7], fitParams[8] );

   xframe->Draw("E");
  return xframe;
}

//--------------------------------------------------
void CheckFitRanges(){
  
  if (fitCutMassMin < xrangem) {
    fitCutMassMin=xrangem;
  }
  if (fitCutMassMax<xrangeM) {
    fitCutMassMax=xrangeM;
  }
  return;
}

//--------------------------------------------------
void SetFitParams(
		  Double_t massRangeMin = 0.78,
		  Double_t massRangeMax = 1.26,
		  Float_t nsigmaPeak= 7.0, 
		  Float_t cutMassMin=-1.0,
		  Float_t cutMassMax=-1.0,
		  Float_t cutWidthMin=-1.0,
		  Float_t cutWidthMax=-1.0
		  ){
  xrangem = massRangeMin;
  xrangeM = massRangeMax;
  nWidthRangeFit = nsigmaPeak;
  
  if (cutMassMin<0.) fitCutMassMin = pdgMass - nWidthRangeFit * pdgWidth / 2.35;
  else fitCutMassMin=cutMassMin;
  
  if (cutMassMax<0.)  fitCutMassMax = pdgMass + nWidthRangeFit * pdgWidth / 2.35;
  else fitCutMassMax=cutMassMax;
  
  if (fitCutWidthMin<0.) fitCutWidthMin=pdgWidth;
  else fitCutWidthMin=cutWidthMin;
  
  if (fitCutWidthMax<0.) fitCutWidthMax=1*pdgWidth;//+pdgWidth*0.5;
  else fitCutWidthMax=cutWidthMax;
 
  CheckFitRanges();
  
  return;
}
