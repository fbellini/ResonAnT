#include "TCanvas.h"
#include "RooPlot.h"
#include "RooAbsReal.h"

using namespace RooFit;

enum EFitFunction{ kBW,
		   kPOLY1,
		   kPOLY2,
		   kPOLY3,
		   kLAND,
		   kEXP,
		   kVOIGT,
		   kRELBW,
		   kREL3,
		   kRELE,
		   kRELBWPS,
		   kBWPS3,
		   kEXE,
		   kData};

Color_t color[]={kRed+1, kOrange, kGreen+2, kBlue+1, kMagenta, kBlack};
Color_t colorfit[]={ kTeal-5, kMagenta+2, kBlue, kRed, kMagenta+1, kYellow+2, kSpring+2, kBlue+2, kOrange, kAzure-7, kPink-7, kAzure+8, kGreen-5, kBlack};
Int_t marker[]={20, 21, 22, 23, 28};

enum EbackgndType{ kUnlikeSE, kUnlikeME,  kLikeSE, kTrue, knBackgndTypes};

// retrieve reference values in the database PDG
Int_t PDG = 313;
TDatabasePDG *pdg     = TDatabasePDG::Instance();
TParticlePDG *part    = pdg->GetParticle(PDG);
Double_t      pdgMass    = part->Mass(); // const Double_t pdgMass = 0.89594;
Double_t      pdgWidth   = part->Width(); // 0.0487;// const Double_t pdgWidth = 0.0487;

// Double_t xrangem, xrangeM;
// Float_t  nWidthRangeFit;
// Double_t fitCutMassMin,fitCutMassMax,fitCutWidthMin, fitCutWidthMax;

//gROOT->SetStyle("Plain");
gSystem->Load("libRooFit");


//--------------------------------------------------
void fitKStarPro(
	      TString infilename="sub_aod49_kstar.root", 
	      TString outdirname = "roofit/", /* / mandatory*/
	      Bool_t isRebin=kFALSE,
	      Bool_t useChi2 = kTRUE,
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
	      Bool_t isShowFullRrange = 0,
	      Float_t cutMassMin= 0.0,
	      Float_t cutMassMax= 2.0,
	      Bool_t isWconst = kTRUE,
	      Float_t cutWidthMin= 0.0,
	      Float_t cutWidthMax= 1.0,
	      Bool_t isResVconst = kTRUE,
	      Double_t res = 0.003,
	      Double_t cutResMin=0.001,
	      Double_t cutResMax=0.050
	       )
{
  
  gROOT->LoadMacro("/Users/bellini/alice/macro/SetGraphicStyle.C");
  //  SetGraphicStyle(0,0,0);
#ifdef __CINT__
  gROOT->ProcessLine(".x $ASD/kstar/roofit/RooRelBW.cxx+") ;
  //  gROOT->ProcessLine(".x RooBoltzPS.cxx+") ;
  gROOT->ProcessLine(".x $ASD/kstar/roofit/RooRelBWPS2.cxx+") ;
#endif
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
    //Printf("%5.2f",pt[k]);
  }
  const Int_t dimcent = ncent+1;
  Double_t cent[dimcent]; 
  for (Int_t k=0; k<dimcent;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
    //Printf("%5.2f",cent[k]);
  }
  fin->Close();

  Int_t nbinstodisplay = istopptbin-istartptbin;
  
  //define output
  TString fileout="";
  if (backgndType==EbackgndType::kUnlikeME){
    fileout.Form("%sfitEM_%s_%3.2f-%3.2f_%s",(isRebin? "rebin_" : ""),function.Data(),massRangeMin,massRangeMax,infilename.Data());
    fileout.ReplaceAll("sub_","");
  } 
  if (backgndType==EbackgndType::kLikeSE){
    fileout.Form("%sfitLS_%s_%3.2f-%3.2f_%s", (isRebin? "rebin_" : ""),function.Data(),massRangeMin,massRangeMax,infilename.Data());
    fileout.ReplaceAll("sub_","");
  } 
  
  if (backgndType==EbackgndType::kTrue){
    fileout.Form("%sfitTrue_%s_%3.2f-%3.2f_%s",(isRebin? "rebin_" : ""),function.Data(),massRangeMin,massRangeMax,infilename.Data());
    fileout.ReplaceAll("projMC","");
    fileout.ReplaceAll("proj/","");
  } 
  
  if (selCentBin>-1) fileout.ReplaceAll(".root",Form("_cent%i.root",selCentBin));
  if (nbinstodisplay==1) fileout.ReplaceAll(".root",Form("_Pt%i.root",istartptbin));
  if (fileout.Contains("../")) fileout.ReplaceAll("../","");
  fileout.Prepend(Form("%s/",outdirname.Data()));
  TFile * fout=new TFile(fileout.Data(),"recreate");
  
  Double_t fitParams[11], SoverB=0.0, significance=0.0, normfactorCopy;
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
  tree->Branch("signalRes",&fitParams[9],"signalRes/D");
  tree->Branch("signalResErr",&fitParams[10],"signalResErr/D");
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
  tree->Branch("functionID", &funcID, "functionID/I");
  tree->Branch("fitrange_inf", &massRangeInfCopy, "fitrange_inf/F");
  tree->Branch("fitrange_sup", &massRangeSupCopy, "fitrange_sup/F");
  //new from 24/04 for bin counting
  // tree->Branch("nSignal3s",&nSignal3s,"nSignal3s/D");
  // tree->Branch("nSignal3sErr",&nSignal3sErr,"nSignal3sErr/D");
  // tree->Branch("nBack3s",&nBack3s,"nBack3s/D");
  // tree->Branch("nBack3sErr",&nBack3sErr,"nBack3sErr/D");
  
  TString centBinLabel;
  if (istopptbin==-1) istopptbin=npt;
  RooPlot *plotframe;
  TLegend * leg = new TLegend(0.6,0.7,0.88,0.88);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetBorderSize(0.0);

  //set fit range in tree
  massRangeInfCopy=massRangeMin; 
  massRangeSupCopy=massRangeMax;
  //set function in tree
  fitfunction=function;
  if (fitfunction.Contains("BREIT")) funcID=EFitFunction::kBW;
  if (fitfunction.Contains("POLY1")) funcID=EFitFunction::kPOLY1;
  if (fitfunction.Contains("POLY2")) funcID=EFitFunction::kPOLY2;
  if (fitfunction.Contains("POLY3")) funcID=EFitFunction::kPOLY3;
  if (fitfunction.Contains("LAND")) funcID=EFitFunction::kLAND;
  if (fitfunction.Contains("EXP")) funcID=EFitFunction::kEXP;
  if (fitfunction.Contains("VOIGT")) funcID=EFitFunction::kVOIGT;
  if (fitfunction.Contains("RELBW")) funcID=EFitFunction::kRELBW;
  if (fitfunction.Contains("REL3")) funcID=EFitFunction::kREL3;
  if (fitfunction.Contains("RELE")) funcID=EFitFunction::kRELE;
  if (fitfunction.Contains("BOLTZ")) funcID=EFitFunction::kRELBWPS;
  if (fitfunction.Contains("BWPS3")) funcID=EFitFunction::kBWPS3;
  if (fitfunction.Contains("EXE")) funcID=EFitFunction::kEXE;

  //define array for minuit exit status
  Int_t statusArray[4]={-1,-1,-1,-1};

  //loop on bins
  for (Int_t icentbin=0; icentbin<ncent;icentbin++){   
    if ( (selCentBin>-1) && (!(icentbin==selCentBin)) ) continue;
    
    char cutstring[100];
    sprintf(cutstring,"");
    centBinID=icentbin;
    Int_t canvasDim = 700;
    if (canvasSplit<4) canvasDim = 550;

    TCanvas *cfit=new TCanvas(Form("cfit_%i",icentbin),Form("cfit_%i",icentbin), canvasDim*canvasSplit,(canvasDim*0.7*nbinstodisplay)/canvasSplit);
    cfit->Divide(canvasSplit,(nbinstodisplay)/canvasSplit);
 
    //define histos
    TH1F*hMassVsPt=new TH1F(Form("hMassVsPt_%i",icentbin),Form("K* mass (%2.0f-%2.0f %); p_{t} (GeV/c); M (GeV/c^{2})",cent[icentbin],cent[icentbin+1]), npt, pt);
    TH1F*hWidthVsPt=new TH1F(Form("hWidthVsPt_%i",icentbin),Form("K* width (%2.0f-%2.0f %); p_{t} (GeV/c); #Gamma (GeV/c^{2})",cent[icentbin],cent[icentbin+1]), npt, pt);
    TH1F*hRawVsPt=new TH1F(Form("hRawVsPt_%i",icentbin),Form("K* raw yield (%2.0f-%2.0f %); p_{t} (GeV/c); raw yield, dN/dp_{T}",cent[icentbin],cent[icentbin+1]), npt, pt);

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
      Double_t normfactor_em, normfactor_ls ;
      Int_t treept, treecent;
      TTree *ntree = (TTree*)fin->Get("ntree");
      if (!ntree) {
	Printf("WARNING: no info from norm tree saved.");
      } else {
	ntree->SetBranchAddress("factor_em", &normfactor_em);
	ntree->SetBranchAddress("factor_ls", &normfactor_ls);
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
	    normfactorCopy=normfactor_em;
	  }
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
      if (backgndType==EbackgndType::kTrue){
	hSignalName.Form("TruesPM_ptBin%02i_centBin%02i",iptbin,icentbin);
	testoback.Form("True MC");
      }
      
      TH1D * hnsig = (TH1D*)fin->Get(hSignalName.Data()); 
      
      if (!hnsig){
	printf("fitKstar - ERROR: input histogram %s not found in file %s . returning...\n", hSignalName.Data(),infilename.Data() );
	return;
      }
      //reset minuit status array
      for (Int_t j=0;j<4;j++) {statusArray[j]=-1;}
      //fit
      plotframe=(RooPlot*) fitKStar(hnsig, 	
				    fitParams, 
				    statusArray,
				    funcID, 
				    ptbins->GetBinCenter(iptbin+1),
				    useChi2,
				    isRebin, 	  
				    isShowFullRrange,
				    massRangeMin, 
				    massRangeMax,  
				    cutMassMin, 
				    cutMassMax, 
				    isWconst,
				    cutWidthMin, 
				    cutWidthMax,
				    isResVconst,
				    0.003,
				    cutResMin,
				    cutResMax);

      for (Int_t j=0;j<4;j++) Printf("Status %i = %i", j, statusArray[j]); 
      
      //check status array: if at least one step failed do not save output!     
      if (isMinuitFailedFit(statusArray)) {
	for (Int_t zz=0;zz<5;zz++){
	  Printf("                ||||                  ||||");
	  Printf("                vvvv                  vvvv");
	  Printf("                FIT FAILED cent %i pt %i", icentbin, iptbin);
	}
	continue; //if fit fails, do not fill histos nor tree nor save img if one bin at the time
      }

      // else {
      // 	if (fitParams[8]>5) {
      // 	  Printf("                ||||                  ||||");
      // 	  Printf("                vvvv                  vvvv");
      // 	  Printf("                BAD CHI^2(>5) cent %i pt %i", icentbin, iptbin);
      // 	}
      // 	continue;
      // }

      //draw legend
      leg->Clear();
      Char_t names[10][100];
      for (int i=0; i<plotframe->numItems(); i++) {
	TString obj_name=plotframe->nameOf(i); 
	if (obj_name=="") continue;
	cout << Form("%d. '%s'\n",i,obj_name.Data());
	sprintf(names[i], "%s",obj_name.Data());
      }
      Char_t caption[3][50] = {
	"Data",
	"Breit-Wigner x PS + background",//"BW+poly2",
	"Residual background"
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
	textFit = new TPaveText(0.55,0.65,0.95,0.89,"NDC");
	//textFit = new TPaveText(0.12,0.12,0.55,0.39,"NDC");
	textFit->SetBorderSize(1);
	textFit->SetFillColor(kWhite);
      } else {
	if (selCentBin<0) textFit = new TPaveText(0.55,0.65,0.90,0.89,"NDC");
        //if (selCentBin==0) textFit = new TPaveText(0.12,0.13,0.75,0.43,"NDC");
        if (selCentBin==0) textFit = new TPaveText(0.52,0.7,0.89,0.89,"NDC"); //this for pA
        if (selCentBin==1) textFit = new TPaveText(0.12,0.13,0.65,0.38,"NDC");
        if (selCentBin>=2) textFit = new TPaveText(0.55,0.6,0.89,0.89,"NDC");
	//	if (selCentBin==3) textFit = new TPaveText(0.12,0.13,0.75,0.43,"NDC");
        textFit->SetBorderSize(0);
	textFit->SetFillStyle(0);
	textFit->SetTextAlign(12); //32 = vert middle, horiz right
      } 
      if (printAll) textFit->AddText(Form("%s, %4.2f#leqM_{inv}#leq%4.2f",function.Data(),massRangeMin,massRangeMax));
      textFit->AddText(Form("M(GeV/c^{2}) =  %6.4f #pm %6.4f",fitParams[0],  fitParams[1]));
      if (printAll) textFit->AddText(Form("#Gamma(K*) =  %6.4f #pm %6.4f GeV/c^{2} \n",fitParams[2],  fitParams[3]));
      if (printAll && function.Contains("VOIGT")) textFit->AddText(Form("#sigma(K*) =  %6.4f #pm %6.4f GeV/c^{2} \n",fitParams[9],  fitParams[10]));
      textFit->AddText(Form("N_{raw} =  %8.0f #pm %8.0f \n",fitParams[4],  fitParams[5]));
      textFit->AddText(Form("N_{Bg} = %8.0f #pm %8.0f \n",fitParams[6],  fitParams[7]));
      textFit->AddText(Form("#chi^{2} = %6.4f \n",fitParams[8]));
      
      textFit->SetTextFont(42);
      textFit->SetTextColor(kBlack);
      SoverB = fitParams[4]/fitParams[6];
      significance = fitParams[4]/TMath::Sqrt(fitParams[4]+fitParams[6]);
     
      //makeup for data plot
      hnsig->SetMarkerSize(0.7);
      hnsig->SetMarkerColor(colorfit[EFitFunction::kData]);
      hnsig->SetLineColor(colorfit[EFitFunction::kData]);
      hnsig->SetLineWidth(1);
      hnsig->GetXaxis()->SetRangeUser(0.6,1.2);
      hnsig->GetXaxis()->SetTitleFont(42);
      hnsig->GetXaxis()->SetTitleOffset(1.2);
      hnsig->GetXaxis()->SetTitle("M(K#pi) (GeV/c^{2})");

      //Draw in canvas
      if (isRebin) plotframe->GetYaxis()->SetTitle("counts / (20 MeV/c^{2})");//RangeUser(0.01,7e5);
      else  plotframe->GetYaxis()->SetTitle("counts / (10 MeV/c^{2})");
      plotframe->GetXaxis()->SetTitle("M(K#pi) (GeV/c^{2})");
      plotframe->GetYaxis()->SetTitleOffset(1.1);//RangeUser(0.01,7e5);
      plotframe->GetYaxis()->SetTitleSize(0.045);//RangeUser(0.01,7e5);
      plotframe->GetYaxis()->SetLabelSize(0.045);//RangeUser(0.01,7e5);
      plotframe->GetXaxis()->SetTitleOffset(0.9);//RangeUser(0.01,7e5);
      plotframe->GetXaxis()->SetLabelOffset(0.008);//RangeUser(0.01,7e5);
      plotframe->GetXaxis()->SetTitleSize(0.045);//RangeUser(0.01,7e5);
      plotframe->GetXaxis()->SetLabelSize(0.04);//RangeUser(0.01,7e5);
      plotframe->Draw(); 
      
      leg->Draw("");
      hnsig->Draw("same");
      if (enaFitResPave) textFit->Draw("same");
      //if (IsSaveResult()) {	 
      //fill mass vs pt plot
      hMassVsPt->SetBinContent(iptbin+1,fitParams[0]);
      hMassVsPt->SetBinError(iptbin+1,fitParams[1]);
      //fill width vs pt plot
      hWidthVsPt->SetBinContent(iptbin+1,fitParams[2]);
      hWidthVsPt->SetBinError(iptbin+1,fitParams[3]);
      //fill raw yield vs pt plot
      Double_t binwidth = ptbins->GetBinWidth(iptbin+1);
      hRawVsPt->SetBinContent(iptbin+1,fitParams[4]/binwidth);
      hRawVsPt->SetBinError(iptbin+1,fitParams[5]/binwidth);
      //add result to tree
      tree->Fill();	 
       // } else {
       // 	 Printf("Fit result NOT saved in tree nor histos"); 
       // 	 continue;
       // }
    }//end loop on pt
    
    hMassVsPt->SetMarkerStyle(marker[icentbin]);
    hMassVsPt->SetMarkerColor(color[icentbin]);
    hMassVsPt->SetLineColor(color[icentbin]);                 
   
    hWidthVsPt->SetMarkerStyle(marker[icentbin]);
    hWidthVsPt->SetMarkerColor(color[icentbin]);
    hWidthVsPt->SetLineColor(color[icentbin]);                 
   
    hRawVsPt->SetMarkerStyle(marker[icentbin]);
    hRawVsPt->SetMarkerColor(color[icentbin]);
    hRawVsPt->SetLineColor(color[icentbin]);                 

    fout->cd();
    hMassVsPt->Write();    
    hWidthVsPt->Write();
    hRawVsPt->Write();

    //TString nameimg1 = Form("roofit/fit%s_%s_range%3.2f-%3.2f_cent%i_Pt%i.png", (backgndType==EbackgndType::kUnlikeME ? "EM" : "LS"), function.Data(), massRangeMin,massRangeMax,icentbin,istartptbin);
    //TString nameimg2= Form("roofit/fit%s_%s_range%3.2f-%3.2f_canvas%ix%i_cent%i_startPt%i.png", (backgndType==EbackgndType::kUnlikeME ? "EM" : "LS"), function.Data(), massRangeMin,massRangeMax,canvasSplit,(nbinstodisplay)/canvasSplit,icentbin,istartptbin);
    TString nameimg1 = Form("%s_c%i_pt%i_range%3.2f-%3.2f.png",function.Data(),icentbin,istartptbin,massRangeMin,massRangeMax );
    TString nameimg2 = Form("%s_c%i_startpt%i_range%3.2f-%3.2f.png",function.Data(),icentbin,istartptbin, massRangeMin,massRangeMax);
    if (infilename.Contains("_antikstar_")) {
      nameimg1.Prepend("antikstar_");
      nameimg2.Prepend("antikstar_");
    }
    if (infilename.Contains("_kstar_")) {
      nameimg1.Prepend("kstar_");
      nameimg2.Prepend("kstar_");
    }
    if (infilename.Contains("_sum_")) {
      nameimg1.Prepend("sum_");
      nameimg2.Prepend("sum_");
    }
    nameimg1.Prepend(Form("%s/",outdirname.Data()));
    nameimg2.Prepend(Form("%s/",outdirname.Data()));
    
    if (nbinstodisplay==1) cfit->SaveAs(nameimg1.Data());
    else cfit->SaveAs(nameimg2.Data());
    
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
RooPlot* fitKStar(TH1D * hnsig, 
		  Double_t * fitParams,
		  Int_t *statusArray,
		  Int_t funcID = EFitFunction::kBWPS3, 
		  Double_t momt,
		  Bool_t useChi2=kTRUE,
		  Bool_t isRebin = kTRUE,
		  Bool_t isShowFullRrange = kTRUE,
		  Double_t xrangem = 0.6, 
		  Double_t xrangeM = 1.3, 
		  Float_t cutMassMin = 0.0,
		  Float_t cutMassMax = 2.0,
		  Bool_t isWconst = kTRUE,
		  Float_t cutWidthMin = -1.0,
		  Float_t cutWidthMax = -1.0,
		  Bool_t isResVconst = kTRUE,
		  Double_t res = 0.003,
		  Double_t cutResMin=0.001,
		  Double_t cutResMax=0.050
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
  
  if (isRebin) hnsig->Rebin(2);

  if ( (xrangem<0.64)||(xrangeM>1.3) ) {
    Printf("Requested fit range exceeds histo range - setting to M_low = 0.64 GeV/c^2, M_up = 1.3 GeV/c^2");
    xrangem=0.64;
    xrangeM=1.3;
  }

  Double_t histo_integral = hnsig->Integral();  
  Double_t signalMass,signalMassErr, signalRes,signalResErr,
    signalWidth,signalWidthErr,
    chi2;
  Int_t nrange = 0;
  Double_t nBack=0;
  Double_t nBackErr=0;
  Double_t nSignal=0;
  Double_t nSignalErr=0;
  Double_t showxrangem = 0.64,showxrangeM = 1.2;
  Double_t res=0.003,fitResMin=0.001,fitResMax=0.050; //Voigtian resolution param
  
  RooRealVar x("m","M (GeV/c^{2})", xrangem, xrangeM);//showxrangem, showxrangeM
  RooDataHist data("data","data",RooArgList(x), hnsig);
  
  /*-----------------------------------------------------
    peak functions
    -----------------------------------------------------*/
  if (!isWconst && (cutWidthMin<=0.0)) cutWidthMin=0.5*pdgWidth;
  if (!isWconst && (cutWidthMax<=0.0)) cutWidthMax=1.5*pdgWidth;
  RooRealVar width("width","width",pdgWidth,cutWidthMin,cutWidthMax);
  if (isWconst) width.setConstant(kTRUE);
  RooRealVar mean("mean","mean",pdgMass,cutMassMin,cutMassMax);
  RooRealVar vres("res","res", res ,cutResMin,cutResMax);
  if (isResVconst) vres.setConstant(kTRUE);
  RooRealVar pt("pt","pt", momt, momt, momt);
  pt.setConstant(kTRUE); 
  // RooRealVar T("T","temperature", 0.2, 0.140, 0.3); //in GeV/c2
  //  T.setConstant(kTRUE); 
  
  RooBreitWigner breit("breit","non-rel. Breit-Wigner",x, mean, width);
  RooRelBW relbw("relbw","relativistic BW", x, mean, width);
  RooRelBWPS2 relbwps("relbwps","BW (X) PSfactor", x, pt, mean, width) ;
  RooVoigtian voigt("voigt","signal Voigt",x, mean, width, vres, kTRUE);//last arg enables fast algorithm if true
 
  /*-----------------------------------------------------
    residual background 1 - Chebychev polynomials 
    -----------------------------------------------------*/
  //Poly 3
  RooRealVar c0("c0","coefficient #0", -0.5, -1.0, 0.0) ;
  RooRealVar c1("c1","coefficient #1", -0.05, -0.5, 0.0) ; //-0.5, 0.0 poly2
  RooRealVar c2("c2","coefficient #2", -0.01, -0.5, 0.5) ;
  //Poly 2
  // RooRealVar a0("a0","coefficient #0", -1.0, 0.0) ; //-1,1
  // RooRealVar a1("a1","coefficient #1", -0.1, 0.1) ; //-0.5, 0.0 poly2
  RooRealVar a0("a0","coefficient #0", -0.5, 0.0) ; //-1,1
  RooRealVar a1("a1","coefficient #1", -0.1, 0.10) ; //-0.5, 0.0 poly2
  
  //poly1
  // f(x) = sum_i a_i * x^i
  // By default coefficient a_0 is chosen to be 1, as polynomial probability density functions have one degree of freedome less than polynomial functions due to the normalization condition
  RooRealVar p0("p0","coefficient #0", -4.0, 1.0) ; //-1,1
  RooRealVar p1("p1","coefficient #1", -10.0, 10.00) ; //-0.5, 0.0 
  
  RooPolynomial poly1("cheby","residual background Chebychev poli", x, RooArgList(p0, p1));
  // RooChebychev cheby1("cheby","residual background Chebychev poli", x, RooArgList(c0));
  RooChebychev cheby2("cheby","residual background Chebychev poli", x, RooArgList(a0,a1));
  RooChebychev cheby3("cheby","residual background Chebychev poli", x, RooArgList(c0,c1,c2)); 
  
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

  /*-----------------------------------------------------
    residual background 4 - exp (x) erf
    ----------------------------------------------------*/
   RooRealVar p1("p1","p1", 0.6, 0., 1.5);
   RooRealVar p2("p2","p2", 0.6, 0., 1.5);
   RooFormulaVar acc("acc","acc","0.5*(1-TMath::Erf((@1-@0)/@0/@2))", RooArgSet(x,p1,p2));
  
   RooRealVar sigmaG("sigmaG","Gaussian sigma", 0.05, 0.01, 0.10);
   RooRealVar meanG("meanG","Gaussian mean", 0.7, 0.64, 0.8);  
   RooGaussian gauss("gaus","res. bg. gaussian smearing", x, meanG,sigmaG );

   RooRealVar t("t","t", -10.0, -20., 0.);
   RooExponential expo2("expo2","residual background exponential",x,t);
    
   RooEffProd lxg("lxg","exp x acc", expo2, acc);

   //RooEffProd lxg("lxg","exp x gauss", expo2, gauss); 
  // RooFFTConvPdf lxg("lxg","landau (X) gauss",x,landau,gauss);
  // RooNumConvPdf lxg("lxg","landau x gauss",x,landau,gauss) ; 

   /*----------------------------------------------
     Model for S+B - extended likelihood
   ----------------------------------------------*/
   RooRealVar nsig("nsig","signal fraction",histo_integral*0.5, 0., histo_integral) ;
   RooRealVar nbkg("nbkg","background fraction",histo_integral*0.5, 0., histo_integral) ;

      //use non-relativistic bw
   RooAddPdf modelBW("modelBW","modelBW",RooArgList(breit),RooArgList(nsig)) ;
   RooAddPdf model1("model1","model1",RooArgList(breit,poly1),RooArgList(nsig,nbkg)) ;
   RooAddPdf model2("model2","model2",RooArgList(breit,cheby2),RooArgList(nsig,nbkg)) ;
   RooAddPdf model3("model3","model3",RooArgList(breit,cheby3),RooArgList(nsig,nbkg)) ;
   RooAddPdf modelL("modelL","modelL",RooArgList(breit,landau),RooArgList(nsig,nbkg)) ;
   RooAddPdf modelE("modelE","modelE",RooArgList(breit,expo),RooArgList(nsig,nbkg)) ;

   //Use voigtian
   RooAddPdf modelV("modelV","modelV",RooArgList(voigt,cheby3),RooArgList(nsig,nbkg)) ;

   //Use relativistic BW 
   RooAddPdf modelR2("modelR2","modelR2",RooArgList(relbw,cheby2),RooArgList(nsig,nbkg)) ;
   RooAddPdf modelR3("modelR3","modelR3",RooArgList(relbw,cheby3),RooArgList(nsig,nbkg)) ;
   RooAddPdf modelRE("modelRE","modelRE",RooArgList(relbw,lxg),RooArgList(nsig,nbkg)) ;

   //use relativistic BW x Boltzmann factor
   RooAddPdf modelRPS2("modelRPS2","modelRPS2",RooArgList(relbwps,cheby2),RooArgList(nsig,nbkg)) ;
   RooAddPdf modelRPS3("modelRPS3","modelRPS3",RooArgList(relbwps,cheby3),RooArgList(nsig,nbkg)) ;
   RooAddPdf modelEXE("modelEXE","modelEXE",RooArgList(relbwps,lxg),RooArgList(nsig,nbkg)) ;

   // //print meaning of parameters
   // RooArgSet * params = model.getVariables() ;
   // RooArgSet * paramsPS = modelRPS.getVariables() ;
   // RooArgSet * paramsEXE = model.getVariables() ;
   // if (function.Contains("BOLTZ")) paramsPS->Print("v") ;
   // else {
   //   if (function.Contains("EXE")) paramsEXE->Print("v") ;
   //   else params->Print("v") ;
   // } 
   
  /*-------------------------------------------------
     display data
     -------------------------------------------------*/
   if (!isShowFullRrange) {
     showxrangem=xrangem;
     showxrangeM=xrangeM;    
   } else {
     showxrangem=0.64;
     showxrangeM=1.2;    
   }
   RooPlot * xframe = x.frame(showxrangem, showxrangeM);
   data.plotOn(xframe, Name("data"), LineColor(colorfit[EFitFunction::kData]), MarkerColor(colorfit[EFitFunction::kData]), MarkerSize(0.2), LineWidth(2), DrawOption("E1X0")) ;//nota: add Name("data") to get chi2
   
   /*-------------------------------------------------
     perform chi2 or likelihood fit & display result
     -------------------------------------------------*/
   RooFitResult * result;   
   Printf("°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°\n°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° \nFunction ID = %i \n°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°\n°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°", funcID);
   Bool_t useImprove = kFALSE;
   Bool_t useMinos = kTRUE;

   if (funcID == EFitFunction::kBW) {
     RooChi2Var *fitchi2= new RooChi2Var("fitchi2","fitchi2", modelBW, data, kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", modelBW, data);
     if (useChi2)
       result = (RooFitResult*) performChi2Fit(fitchi2, statusArray, useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll, statusArray, useMinos, useImprove);
     modelBW.plotOn(xframe, Name("modelBW"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
     chi2=xframe->chiSquare("modelBW", "data", 3);
   }
   
   if (funcID == EFitFunction::kPOLY1) { 
     RooChi2Var * fitchi2= new RooChi2Var("fitchi2","fitchi2", model1, data, kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", model1, data, Extended(kTRUE),Verbose(kTRUE));
     if (useChi2)
       result = (RooFitResult*) performChi2Fit(fitchi2, statusArray, useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll, statusArray, useMinos, useImprove);
     model1.plotOn(xframe, Name("model1"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
     model1.plotOn(xframe, Components(poly1),LineStyle(kDashed),LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM)) ;
     //model1.plotOn(xframe, Components(breit),LineColor(colorfit[funcID]),LineWidth(1),Range(xrangem, xrangeM)) ;
     chi2=xframe->chiSquare("model1", "data", 3);
   }
   
   if (funcID == EFitFunction::kPOLY2) {
     RooChi2Var *fitchi2= new RooChi2Var("fitchi2","fitchi2", model2, data, kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", model2, data, Extended(kTRUE),Verbose(kTRUE));
     if (useChi2)
       result = (RooFitResult*) performChi2Fit(fitchi2, statusArray, useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll, statusArray, useMinos, useImprove);
     model2.plotOn(xframe, Name("model2"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
     model2.plotOn(xframe, Components(cheby2),LineStyle(kDashed),LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM)) ;
     //model2.plotOn(xframe, Components(breit),LineColor(colorfit[funcID]),LineWidth(1),Range(xrangem, xrangeM)) ;
     chi2=xframe->chiSquare("model2", "data", 3);
   }
     
   if (funcID == EFitFunction::kPOLY3 ) {
     RooChi2Var *fitchi2= new RooChi2Var("fitchi2","fitchi2", model3, data, kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", model3, data, Extended(kTRUE),Verbose(kTRUE));
      if (useChi2)
	result = (RooFitResult*) performChi2Fit(fitchi2, statusArray, useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll, statusArray, useMinos, useImprove);
      model3.plotOn(xframe, Name("model3"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
      model3.plotOn(xframe, Components(cheby3),LineStyle(kDashed),LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM)) ;
      //model3.plotOn(xframe, Components(breit),LineColor(colorfit[funcID]),LineWidth(1),Range(xrangem, xrangeM)) ;
      chi2=xframe->chiSquare("model3", "data", 3);
   };
       
   if (funcID ==  EFitFunction::kLAND ) {
     RooChi2Var *fitchi2= new RooChi2Var("fitchi2","fitchi2", modelL, data, kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", modelL, data, Extended(kTRUE),Verbose(kTRUE));
     if (useChi2)
       result = (RooFitResult*) performChi2Fit(fitchi2, statusArray, useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll,statusArray,  useMinos, useImprove);
     modelL.plotOn(xframe, Name("modelL"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
     modelL.plotOn(xframe, Components(landau),LineStyle(kDashed),LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM)) ;
     //model1.plotOn(xframe, Components(breit),LineColor(colorfit[funcID]),LineWidth(1),Range(xrangem, xrangeM)) ;
     chi2=xframe->chiSquare("modelL", "data", 3);
   };

   if (funcID ==  EFitFunction::kEXP ) {
     RooChi2Var *fitchi2= new RooChi2Var("fitchi2","fitchi2", modelE, data, kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", modelE, data, Extended(kTRUE),Verbose(kTRUE));
      if (useChi2)
       result = (RooFitResult*) performChi2Fit(fitchi2,statusArray,  useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll,statusArray,  useMinos, useImprove);
    modelE.plotOn(xframe, Name("modelE"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
     modelE.plotOn(xframe, Components(expo),LineStyle(kDashed),LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM)) ;
     //modelE.plotOn(xframe, Components(breit),LineColor(colorfit[funcID]),LineWidth(1),Range(xrangem, xrangeM)) ;
     chi2=xframe->chiSquare("modelE", "data", 3);
   };

   if (funcID ==  EFitFunction::kVOIGT ) {
     RooChi2Var *fitchi2= new RooChi2Var("fitchi2","fitchi2", modelV, data, kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", modelV, data, Extended(kTRUE),Verbose(kTRUE));
     if (useChi2)
       result = (RooFitResult*) performChi2Fit(fitchi2,statusArray,  useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll,statusArray, useMinos, useImprove);
     modelV.plotOn(xframe, Name("modelV"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
     modelV.plotOn(xframe, Components(cheby2),LineStyle(kDashed),LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM)) ;
     //model1.plotOn(xframe, Components(breit),LineColor(colorfit[funcID]),LineWidth(1),Range(xrangem, xrangeM)) ;
     chi2=xframe->chiSquare("modelV", "data", 3);
   };

   if (funcID ==  EFitFunction::kRELBW ) {
     RooChi2Var * fitchi2= new RooChi2Var("fitchi2","fitchi2", modelR2, data, kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", modelR2, data, Extended(kTRUE),Verbose(kTRUE));
     if (useChi2)
       result = (RooFitResult*) performChi2Fit(fitchi2, statusArray, useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll, statusArray, useMinos, useImprove);
     modelR2.plotOn(xframe, Name("modelR2"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
     modelR2.plotOn(xframe, Components(cheby2),LineStyle(kDashed),LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM)) ;
     //modelR2.plotOn(xframe, Components(relbw),LineColor(colorfit[funcID]),LineWidth(1),Range(xrangem, xrangeM)) ;
     chi2=xframe->chiSquare("modelR2", "data", 3);
   };

   if (funcID ==  EFitFunction::kREL3 ) {
     RooChi2Var *fitchi2= new RooChi2Var("fitchi2","fitchi2", modelR3, data, kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", modelR3, data, Extended(kTRUE),Verbose(kTRUE));
     if (useChi2)
       result = (RooFitResult*) performChi2Fit(fitchi2, statusArray, useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll,statusArray,  useMinos, useImprove);
     modelR3.plotOn(xframe, Name("modelR3"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
     modelR3.plotOn(xframe, Components(cheby3),LineStyle(kDashed),LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM)) ;
     //modelR3.plotOn(xframe, Components(relbw),LineColor(colorfit[funcID]),LineWidth(1),Range(xrangem, xrangeM)) ;
     chi2=xframe->chiSquare("modelR3", "data", 3);
   };

   if (funcID ==  EFitFunction::kRELE ) {
     RooChi2Var *fitchi2= new RooChi2Var("fitchi2","fitchi2", modelRE, data,kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", modelRE, data, Extended(kTRUE),Verbose(kTRUE));
     if (useChi2)
       result = (RooFitResult*) performChi2Fit(fitchi2, statusArray, useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll,statusArray,  useMinos, useImprove);
     modelRE.plotOn(xframe, Name("modelRE"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
     modelRE.plotOn(xframe, Components(lxg),LineStyle(kDashed),LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM)) ;
     //modelRE.plotOn(xframe, Components(relbw),LineColor(colorfit[funcID]),LineWidth(1),Range(xrangem, xrangeM)) ;
     chi2=xframe->chiSquare("modelRE", "data", 3);
   };

   if (funcID ==  EFitFunction::kRELBWPS ) {
     RooChi2Var *fitchi2= new RooChi2Var("fitchi2","fitchi2", modelRPS2, data,kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", modelRPS2, data, Extended(kTRUE),Verbose(kTRUE));
     if (useChi2)
       result = (RooFitResult*) performChi2Fit(fitchi2, statusArray, useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll,statusArray,  useMinos, useImprove);
     modelRPS2.plotOn(xframe, Name("modelRPS2"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
     modelRPS2.plotOn(xframe, Components(cheby2),LineStyle(kDashed),LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM)) ;
     //modelRPS2.plotOn(xframe, Components(relbwps),LineColor(colorfit[funcID]),LineWidth(1),Range(xrangem, xrangeM)) ;
     chi2=xframe->chiSquare("modelRPS2", "data", 3);
   };

   if (funcID ==  EFitFunction::kBWPS3 ) {
     Printf("FITTING BWPS3");
     RooChi2Var *fitchi2= new RooChi2Var("fitchi2","fitchi2", modelRPS3, data, kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", modelRPS3, data, Extended(kTRUE),Verbose(kTRUE));
     if (useChi2)
       result = (RooFitResult*) performChi2Fit(fitchi2, statusArray, useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll, statusArray, useMinos, useImprove);
     modelRPS3.plotOn(xframe, Name("modelRPS3"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
     modelRPS3.plotOn(xframe, Components(cheby3),LineStyle(kDashed),LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM)) ;
     //modelRPS3.plotOn(xframe, Components(relbwps),LineColor(colorfit[funcID]),LineWidth(1),Range(xrangem, xrangeM)) ;
     chi2=xframe->chiSquare("modelRPS3", "data", 3);
   };

   if (funcID ==  EFitFunction::kEXE ) {
     RooChi2Var *fitchi2= new RooChi2Var("fitchi2","fitchi2", modelEXE, data,kTRUE) ;
     RooNLLVar *nll = new RooNLLVar("nll","nll", modelEXE, data, Extended(kTRUE),Verbose(kTRUE));
     if (useChi2)
       result = (RooFitResult*) performChi2Fit(fitchi2, statusArray, useMinos, useImprove);
     else 
       result = (RooFitResult*) performLikelihoodFit(nll,statusArray,  useMinos, useImprove);
     modelEXE.plotOn(xframe, Name("modelEXE"), LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM));
     modelEXE.plotOn(xframe, Components(lxg),LineStyle(kDashed),LineColor(colorfit[funcID]),LineWidth(2),Range(xrangem, xrangeM)) ;
     //modelEXE.plotOn(xframe, Components(relbwps),LineColor(colorfit[funcID]),LineWidth(1),Range(xrangem, xrangeM)) ;
     chi2=xframe->chiSquare("modelEXE", "data", 3);
   };

   result->Print("v");
 
   signalMass=mean.getVal();
   signalMassErr=mean.getError();
   if (funcID == EFitFunction::kVOIGT){
     signalRes=vres.getVal();
     signalResErr=vres.getError();
   } else { 
     signalRes=-1.0;
     signalResErr=-1.0;
   }
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
     fitParams[9]=signalRes;
     fitParams[10]=signalResErr;
   }
   xframe->Draw("E");
   return xframe;
}

//--------------------------------------------------------------------
RooFitResult * performChi2Fit(RooChi2Var * fitchi2, Int_t *statusArray, Bool_t useMinos = 1, Bool_t useImprove = 0)
{
  //performs chi2 fit and saves status in array passed as second argument at every minuit step
  if (!fitchi2) return 0x0;
  
  RooMinuit * m2 = new RooMinuit(*fitchi2) ;
  m2->setStrategy(2);
  m2->migrad();
  tmpresult = (RooFitResult*) m2->save();
  statusArray[0]=tmpresult->status();
  m2->hesse();
  tmpresult = (RooFitResult*) m2->save();
  statusArray[1]=tmpresult->status();
  
  if (useImprove) {
    m2->improve();
  tmpresult = (RooFitResult*) m2->save();
  statusArray[2]=tmpresult->status();
  } else statusArray[2]=0;
  
  if (useMinos) {
    m2->minos();
    tmpresult = (RooFitResult*) m2->save();
    statusArray[3]=tmpresult->status();
  }  else statusArray[3]=0;
  
  result = (RooFitResult*) m2->save();
  //result->Print("v");
  return result;
}

//--------------------------------------------------------------------
RooFitResult * performLikelihoodFit(RooNLLVar * nll, Int_t *statusArray, Bool_t useMinos = 1, Bool_t useImprove = 0)
{
  //performs maximum likelihood fit and savees status at each step in passed array
  if (!nll) return 0x0;
  RooRealVar * offset = new RooRealVar("offset","offset",-95633288.84);; 
  RooAbsReal * L = new RooFormulaVar("L","L","@0-@1", RooArgSet(*nll,*offset));
  
  //create minuit session
  RooMinuit * m1 = new RooMinuit(*L);
  //m1->setPrintLevel(4);
  m1->setStrategy(1);
  m1->setEps(1e-16);
  m1->migrad();
  m1->setStrategy(2);
  m1->migrad();
  tmpresult = (RooFitResult*) m1->save();
  statusArray[0]=tmpresult->status();
  m1->hesse();
  tmpresult = (RooFitResult*) m1->save();
  statusArray[1]=tmpresult->status();

  if (useImprove) { 
    m1->improve();
    tmpresult = (RooFitResult*) m1->save();
    statusArray[2]=tmpresult->status();
  } else statusArray[2]=0;
  
  if (useMinos) {
    m1->minos();
    tmpresult = (RooFitResult*) m1->save();
    statusArray[3]=tmpresult->status();
  } else statusArray[3]=0;
  
  result = (RooFitResult*) m1->save();
  return result;
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

//------------------------------------------------------------------------------
Double_t NormalizedIntegral(RooAbsPdf& function, RooRealVar& integrationVar, Double_t lowerLimit, Double_t upperLimit)
{
  integrationVar.setRange("integralRange", lowerLimit, upperLimit) ;
  RooAbsReal* realintegral = function.createIntegral(integrationVar, RooArgSet(integrationVar), Range("integralRange")) ;
  Double_t normalizedIntegralValue = realintegral­>getVal();
  printf("\nIntegrating %s from %f to %f\n", function.GetName(), lowerLimit, upperLimit);
  printf("Integral value: %f\n", normalizedIntegralValue);
  return normalizedIntegralValue;
}

//------------------------------------------------------------------------------
Bool_t isMinuitFailedFit(Int_t *statusArray)
{
  if (!statusArray){
    Printf("Invalid status array passed to checker. Return 1.");
    return kTRUE;
  }

  Printf("=============== CHECKING FIT OUTPUT VIA MINUIT STATUS FLAGS ===============");
  for (Int_t j=0;j<4;j++) {
    if (statusArray[j]==4) {
      Printf("Minuit status is 4: migrad did not converge - FIT FAILED");
      return kTRUE;
    } else 
      if (statusArray[j]==-1) Printf("Minuit status unset - check if fit is performed");
      else 
	if (statusArray[j]!=0) Printf("Minuit status different from 0 or 4. More infos needed.");
  }
  return 0;  
}


Bool_t Pause()
{
  char dummy[100];
  cout << "Pausing... Type 'c' to continue " << endl;
  cin >> dummy;
  switch (dummy){
  case 'c':
    return kTRUE;
  default:
    return kFALSE;
  }
  
  // if (!strcmp(dummy, "c")) return kFALSE;
  // if (!strcmp(dummy, "q")) return kTRUE;  
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
   
