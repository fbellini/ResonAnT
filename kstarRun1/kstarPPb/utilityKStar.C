#define DEBUG 1

const Double_t kBigNumber=1E10;
const Double_t kSmallNumber=1E-10;
//define binning in invariant mass
const Double_t invMassInf=0.6;
const Double_t invMassSup=1.5;
const Double_t invMassBinWidth=0.01; //GeV/c^2
const Double_t reductionFactor=1.;//background reduction factor 
const char kXaxisTitle[20]="M_{inv} (GeV/c^{2})";
const char kYaxisTitle[10]="pairs";

enum EpairType { kUnlikePM, kUnlikeMP, kMixingPM, kMixingMP, kLikePP, kLikeMM, knPairTypes};
enum EbackgndType{ kUnlikeSE, kUnlikeME,  kLikeSE, knBackgndTypes};
enum Epart{ kstar, antikstar, knPart};
TString pairsName[6]={"UnlikePM", "UnlikeMP", "MixingPM", "MixingMP","LikePP", "LikeMM" };
TString backgndNames[3]={"UnlikeSE","UnlikeME","LikeSE"};

//define binning in pt
// Double_t pt[] = {    0.00, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 
// 		     2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.20,
// 		     4.40, 4.60, 4.80, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00, 9.00, 10.00 };
Double_t pt[] = {    0.00, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 4.00, 5.00, 10.00 };

//define binning in centrality
Double_t cent[]={ 0.0, 10.0, 20.0, 40.0, 60.0, 80.0, 90.0};

Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;
Int_t ncent  = sizeof(cent) / sizeof(cent[0]) - 1;

TString macroDir = "/Users/bellini/alice/macro/kstar/";

//------------------------------------------------------------
TList * analyseCentralityBin(Int_t icentbin){

  TString listName;
  listName.Form("Centrality_%2.0f-%20f",cent[icentbin],cent[icentbin+1]);
  TList *list=new TList();
  list->SetName(listName.Data());


  return list;
}


//------------------------------------------------------------
void SumKStarAntiKStarAllPt(TString infilename="subAllPt_proj_analysisAOD.root", Bool_t save2file=kFALSE){
  
  /*
    It sums K* and anti-K* spectra after background subraction 
    for mixing event and like-sign backgnd
    Input file must be the file with distributions after backgnd sub, integrated aver all pT, 
    different centrality bins
    
    -> output is saved on file only if save2file=1
       in that case output file name is infilename preceeded by fixed prefix
   */

  gROOT->LoadMacro("/Users/bellini/alice/macro/SetGraphicStyle.C");
  SetGraphicStyle(0);

  if (!infilename){
    printf("ERROR: invalid file name provided. \n");
    return;
  }
  TFile* fin=TFile::Open(infilename.Data());
  if (!fin) {
    printf("ERROR: invalid or corrupted input file. Exiting.\n");
    return;
  }
  TList outlist;
  
  TH1D * hkstarmix[2]={0,0};
  TH1D * hkstarlike[2]={0,0};
  TString hnameKStar[2][2];
  
  TCanvas *dummycanvas2=new TCanvas("dummycanvas2","dummycanvas2",1000,700);
  dummycanvas2->Divide(3,2);
  TCanvas *dummycanvas3=new TCanvas("dummycanvas3","dummycanvas3",1000,700);
  dummycanvas3->Divide(3,2);


  for (Int_t icentbin=0; icentbin<ncent;icentbin++){    
    //event mixing
    // hnameKStar[0][0].Clear();
    // hnameKStar[0][1].Clear();
    hnameKStar[0][0].Form("subEM_Data_KStar_ptBinAll_centBin%02i",icentbin);
    hnameKStar[0][1].Form("subEM_Data_antiKStar_ptBinAll_centBin%02i",icentbin);
    //likesign
    // hnameKStar[1][0].Clear();
    // hnameKStar[1][1].Clear();
    hnameKStar[1][0].Form("subLS_Data_KStar_ptBinAll_centBin%02i",icentbin);
    hnameKStar[1][1].Form("subLS_Data_antiKStar_ptBinAll_centBin%02i",icentbin);
  
    for (Int_t i=0;i<Epart::knPart;i++){      
      hkstarmix[i]=(TH1D*)fin->Get(hnameKStar[0][i]);
      hkstarlike[i]=(TH1D*)fin->Get(hnameKStar[1][i]);      
    }//loop over part.le   
    
    TString dummyname;
    TString dummytitle;
     
    TH1D*hmix=new TH1D("hmix"," ", TMath::Abs(invMassSup-invMassInf)/invMassBinWidth , invMassInf, invMassSup);
    hmix=(TH1D*)hkstarmix[0]->Clone();
    hmix->Add(hkstarmix[1],1.);
    dummyname.Form("subEM_ptAll_centBin%02i",icentbin);
    dummytitle.Form("K*+#bar{K*} - EM backgnd - all p_{T}, (%3.0f - %3.0f) central",cent[icentbin],cent[icentbin+1] ); 
    hmix->SetNameTitle(dummyname.Data(),dummytitle.Data());
    hmix->GetXaxis()->SetTitle(kXaxisTitle);
    hmix->GetYaxis()->SetTitle(kYaxisTitle);

    outlist.Add(hmix);
    
    //draw on canvas
    dummycanvas2->cd(icentbin+1);
    hmix->Draw();
    dummycanvas2->Update();

    TH1D*hlike=new TH1D("hlike"," ", TMath::Abs(invMassSup-invMassInf)/invMassBinWidth , invMassInf, invMassSup);
    hlike=(TH1D*)hkstarlike[0]->Clone();
    hlike->Add(hkstarlike[1],1.);
    dummyname.Form("subLS_ptAll_centBin%02i",icentbin);
    dummytitle.Form("K*+#bar{K*} - LS backgnd - all p_{T}, (%3.0f - %3.0f) central",cent[icentbin],cent[icentbin+1] ); 
    hlike->SetNameTitle(dummyname.Data(),dummytitle.Data()); 
    hlike->GetXaxis()->SetTitle(kXaxisTitle);
    hlike->GetYaxis()->SetTitle(kYaxisTitle);


    outlist.Add(hlike);
    //draw on canvas
    dummycanvas3->cd(icentbin+1);
    hlike->Draw();
    dummycanvas3->Update();

  }//loop ove cent bins
  
  if (save2file){
    TString outfilename;
    outfilename.Form("sumKAntiK_%s",infilename.Data());
    TFile * fout =new TFile(outfilename.Data(),"recreate");
    fout->cd();
    outlist.Write();
    fout->Close();
    printf("INFO: output saved as requested into file %s \n",outfilename.Data());
  }
  
  return;
}

//-------------------------------------------------------
void sumAllPt(TString infilename="proj_analysisAOD.root", Bool_t save2file=kFALSE){
  
  /*
    It sums K* and anti-K* spectra over all pT bins and performs background subraction 
    for mixing event and like-sign backgnd
    Input file must be the file with signal and backgnd distrib, different pT bins, 
    different centrality bins
    
    -> output is saved on file only if save2file=1
       in that case output file name is infilename preceeded by fixed prefix
   */

  gROOT->LoadMacro("/Users/bellini/alice/macro/SetGraphicStyle.C");
  SetGraphicStyle();

  if (!infilename){
    printf("ERROR: invalid file name provided. \n");
    return;
  }
  TFile* fin=TFile::Open(infilename.Data());
  if (!fin) {
    printf("ERROR: invalid or corrupted input file. Exiting.\n");
    return;
  }

  TString outfilename, sumAllPtDirName,sumPartAllPtDirName ;
  outfilename.Form("sumDistribAllPt_%s",infilename.Data());
  sumAllPtDirName="sumAllPt";
  sumAllPtDirName="sumPartAllPt";
  TFile * fout =new TFile(outfilename.Data(),"recreate");
  fout->mkdir(sumAllPtDirName.Data());
  fout->mkdir(sumPartAllPtDirName.Data());
  
  //define output histos
  TH1D * histosKStarSum[6];
  for (Int_t ih=0;ih<6;ih++){
    histosKStarSum[ih]=new TH1D("hdummy"," ", TMath::Abs(invMassSup-invMassInf)/invMassBinWidth , invMassInf, invMassSup);
  }

  TString hnameKStarSum[6];//signal, EM, LS
  TString htitleKStarSum[6];//signal, EM, LS
  
  TH1D * histosPartSum[3];
  TString hnamePartSum[3];//signal, EM, LS
  TString htitlePartSum[3];//signal, EM, LS
  
  //read input histos
  TString hnameKStar[6];//signal, EM, LS
  TH1D *histosKStar[6]; //2 particles * 3 types of pairs
  Bool_t isInitialized=0;
  
  TCanvas *dummycanvas=new TCanvas("dummycanvas","dummycanvas",900,600);
  dummycanvas->Divide(3,2);
  dummycanvas->Update();
  TCanvas *dummycanvas2=new TCanvas("dummycanvas2","dummycanvas2",900,600);
  dummycanvas2->Divide(3,2);
  dummycanvas2->Update();
  
  for (Int_t icentbin=0; icentbin<ncent;icentbin++){   
    isInitialized=0;
    for (Int_t ipair=0; ipair < EpairType::knPairTypes; ipair++){
      histosKStarSum[ipair]->Reset("ICES");
      for (Int_t iptbin=0; iptbin<npt;iptbin++){ 
	hnameKStar[ipair].Form("Data_%s_ptBin%02i_centBin%02i", pairsName[ipair].Data(),iptbin, icentbin);
 	histosKStar[ipair]= (TH1D*)fin->Get(hnameKStar[ipair]);
	
#if DEBUG
	printf("-------------------------------------------------------- ipair = %i:  %i, %i \n",ipair, ipair%2,ipair/2);
	printf("Debug: histo = histosKstar[%i][%i] = %s \n",  ipair%2, ipair/2, hnameKStar[ipair].Data());
#endif	
 	if (!histosKStar[ipair]) {
	  printf("Histo retrieval failed for histo %s! \n", hnameKStar[ipair].Data());	 
	} else {
	  histosKStarSum[ipair]->Add(histosKStar[ipair],1.);
	  printf("adding -> %s ",histosKStar[ipair]->GetName());
	}
      }//loop on pt bins
      hnameKStarSum[ipair].Form("Data_%s_ptBinAll_centBin%02i", pairsName[ipair].Data(), icentbin);
      htitleKStarSum[ipair].Form("K*- all p_{T}, (%3.0f - %3.0f) central",cent[icentbin],cent[icentbin+1] ); 
      histosKStarSum[ipair]->SetNameTitle(hnameKStarSum[ipair].Data(), htitleKStarSum[ipair].Data());
      histosKStarSum[ipair]->GetXaxis()->SetTitle(kXaxisTitle);
      histosKStarSum[ipair]->GetYaxis()->SetTitle(kYaxisTitle);
      
      if (!fout){
	printf("error in output file\n");
	return;
      }
      fout->cd(sumAllPtDirName.Data());	     
      histosKStarSum[ipair]->Write();      
    }//loop on pair types        
    
    for (Int_t j=0;j<3;j++){      
      //sum k* and anti-K*
      histosPartSum[j]=(TH1D*)histosKStarSum[2*j]->Clone();
      histosPartSum[j]->Add(histosKStarSum[2*j+1],1.);
      hnamePartSum[j].Form("Data_%s_ptBinAll_centBin%02i", backgndNames[j].Data(), icentbin);
      htitlePartSum[j].Form("K*+#bar{K*} %s distribution , all p_{T}, (%3.0f - %3.0f) central",backgndNames[j].Data(),cent[icentbin],cent[icentbin+1] );       
      histosPartSum[j]->SetNameTitle(hnamePartSum[j].Data(), htitlePartSum[j].Data());
      printf("sumAllPt - INFO: adding to list histo %s\n",hnamePartSum[j].Data());
      
      switch (j){
      case 0:
	histosPartSum[j]->SetLineColor(kBlack);
	histosPartSum[j]->SetMarkerColor(kBlack);
	break;
      case 1:
	histosPartSum[j]->SetLineColor(kRed);
	histosPartSum[j]->SetMarkerColor(kRed);
	break;
      case 2:
	histosPartSum[j]->SetLineColor(kBlue);
	histosPartSum[j]->SetMarkerColor(kBlue);
	break;
      default:
	break;
      }
      if (save2file){
    	fout->cd(sumPartAllPtDirName.Data());	
    	histosPartSum[j]->Write();
      }          
    }

    dummycanvas->cd(icentbin+1);
    histosPartSum[0]->Draw();
    // TH1D* normLS=normalize2Integral(histosPartSum[2],histosPartSum[0]);
    // TH1D* normEM=normalize2Integral(histosPartSum[1],histosPartSum[0]);
    TH1D* normLS=normValuesInterval(histosPartSum[2],histosPartSum[0]);
    TH1D* normEM=normValuesInterval(histosPartSum[1],histosPartSum[0]);
    normLS->Draw("same");
    normEM->Draw("same");
    dummycanvas->Update();
    TLegend *l=new TLegend(0.45,0.6,0.8,0.8);
    l->AddEntry( histosPartSum[0],"unlike SE","lpf");
    l->AddEntry( normEM,"unlike ME (norm)","lpf");
    l->AddEntry( normLS,"like SE (norm)","lpf");
    l->Draw("same");

    TH1D* hRatioSBoverNormLS=(TH1D*)histosPartSum[0]->Clone();
    TH1D* hRatioSBoverNormEM=(TH1D*)histosPartSum[0]->Clone();
    
    TString titleLS, nameLS;
    titleLS.Form("(S+B) /normalized LS backgnd - cent %2.0f-%2.0f",cent[icentbin],cent[icentbin+1]);
    nameLS.Form("hRatioSBoverNormLS_cent%02i",icentbin);
    hRatioSBoverNormLS->Divide(normLS);    
    hRatioSBoverNormLS->SetNameTitle(nameLS.Data(),titleLS.Data());
    hRatioSBoverNormLS->SetLineColor(kBlue);
    hRatioSBoverNormLS->SetMarkerColor(kBlue);
    hRatioSBoverNormLS->GetYaxis()->SetRangeUser(0.,2.);
    hRatioSBoverNormLS->GetYaxis()->SetTitle("ratio (S+B)/B");
   
    TString titleEM,nameEM;
    titleEM.Form("(S+B) /normalized LS backgnd - cent %2.0f-%2.0f",cent[icentbin],cent[icentbin+1]);
    nameEM.Form("hRatioSBoverNormEM_cent%02i",icentbin);
    hRatioSBoverNormEM->Divide(normEM);
    hRatioSBoverNormEM->SetNameTitle(nameEM.Data(),titleEM.Data());
    hRatioSBoverNormEM->SetLineColor(kRed);
    hRatioSBoverNormEM->SetMarkerColor(kRed);
    hRatioSBoverNormEM->GetYaxis()->SetRangeUser(0.,2.);
    hRatioSBoverNormEM->GetYaxis()->SetTitle("ratio (S+B)/B");

    gStyle->SetOptStat(0);
    dummycanvas2->cd(icentbin+1);
    hRatioSBoverNormLS->Draw();
    hRatioSBoverNormEM->Draw("same");
    dummycanvas2->Update();
 
    TLegend *l2=new TLegend(0.45,0.2,0.8,0.4);
    l2->AddEntry( hRatioSBoverNormEM,"unlike ME (norm)","lpf");
    l2->AddEntry( hRatioSBoverNormLS,"like SE (norm)","lpf");
    l2->Draw("same");x    
    if (save2file){
      hRatioSBoverNormLS->Write();
      hRatioSBoverNormEM->Write();
    }  
  }//loop over cent bins
  fout->cd();
  dummycanvas->Write();
  dummycanvas2->Write();  
  fin->Close();
  fout->Close();
  printf("INFO: output saved as requested into file %s \n",outfilename.Data());
  return;
}

//---------------------------------------------------------------------------------
TH1D* normalize2Integral(TH1D* hist=NULL, TH1D*histRef=NULL, Double_t factor=1.){
  
  if ((!hist)||(!histRef)){
    printf("ERROR: invalid histogram address passed to normalization function\n");
    return 0;
  }
  
  Int_t nbins=hist->GetXaxis()->GetNbins();
  Int_t nbinsRef=histRef->GetXaxis()->GetNbins();
  
  if (nbins!=nbinsRef){
    printf("ERROR: histgrams have different binning. Doing nothing.\n");
    return 0;    
  }
  
  TH1D * cloneH = (TH1D*) hist->Clone();
  Double_t integralS= histRef->Integral();
  Double_t integralB= hist->Integral();
  
  if (integralS==0){
    printf("INFO: signal + backgnd histogram integral = 0. Skipping normalization.\n");
    return 0;
  }
  
  if (integralB==0){
    printf("INFO: backgnd histogram integral = 0. Skipping normalization.\n");
    return 0;
  }
  
  if (factor<=0) factor=1.;
  Double_t normFactor= factor*integralS/integralB;
  cloneH->Scale(normFactor);
  return cloneH;
}

//---------------------------------------------------------------------------------
Float_t GetRangeValuesNormalizationFactor(TH1D* hist=NULL, TH1D* histRef=NULL, Double_t valueMin=1.1, Double_t valueMax=1.5){
  
  /*
    estimate normalization factor to be used on hist
    by using histRef as reference 
    in *values* range given by valueMin - valueMax
    
  */

  if ((!hist)||(!histRef)){
    printf("GetRangeValuesNormalizationFactor - ERROR: invalid histgram address passed to normalization function\n");
    return 0;
  }
  
  Int_t nbins=hist->GetXaxis()->GetNbins();
  Int_t nbinsRef=histRef->GetXaxis()->GetNbins();
  
  if (nbins!=nbinsRef){
    printf("GetRangeValuesNormalizationFactor - ERROR: histgrams have different binning. Doing nothing.\n");
    return 0;    
  }
  
  Int_t ibinmin,ibinmax;
  
  if ((valueMin<=kSmallNumber) || (valueMin>=valueMax)) {
    valueMin = hist->GetXaxis()->GetBinLowEdge(1);
    ibinmin=1;
    printf("GetRangeValuesNormalizationFactor - INFO:  min value used for normalization range as default.\n");
  } else {
    ibinmin= 1 + ( valueMin-(hist->GetXaxis()->GetBinLowEdge(1)) )/invMassBinWidth;
  } 
  
  if ((valueMax>=kBigNumber) || (valueMax<valueMin)) {
    valueMax = hist->GetXaxis()->GetBinUpEdge(nbins);
    ibinmax=nbins;
    printf("GetRangeValuesNormalizationFactor - INFO:  max value used for normalization range as default.\n");  
  } else {    
    ibinmax = nbins - ( (hist->GetXaxis()->GetBinLowEdge(nbins)) - valueMax )/invMassBinWidth;  
  }
  
  printf("--------------- Normalizing histogram in interval [ %5.2f - %5.2f ]\n",valueMin,valueMax);
  printf("                corresponding to (x-axis) bin interval [ %i - %i ]\n",ibinmin,ibinmax);
  
  Double_t integralS= histRef->Integral(ibinmin,ibinmax);
  Double_t integralB= hist->Integral(ibinmin,ibinmax);
  
  if (integralS==0){
    printf("GetRangeValuesNormalizationFactor - INFO: reference histogram integral = 0. Skipping normalization.\n");
    return 0;
  }
  if (integralB==0){
    printf("GetRangeValuesNormalizationFactor - INFO: histogram to be normalized has integral = 0. Skipping normalization.\n");
    return 0;
  }
  
  Double_t normFactor= integralS/integralB;
  printf("GetRangeValuesNormalizationFactor - INFO: Normalization factor estimate = %6.3f\n",normFactor);
  return normFactor;
  
}

//---------------------------------------------------------------------------------
TH1D* normValuesInterval(TH1D* hist=NULL, TH1D* histRef=NULL, Double_t valueMin=1.1, Double_t valueMax=1.5, Double_t factor=1.){
  
  if ((!hist)||(!histRef)){
    printf("ERROR: invalid histgram address passed to normalization function\n");
    return 0;
  }
    printf("normValuesInterval - INFO: Normalizing background in IM interval [ %5.2f - %5.2f ]\n",valueMin,valueMax);
  
  TH1D * cloneH = (TH1D*) hist->Clone();  
  Double_t normFactor= GetRangeValuesNormalizationFactor(hist,histRef,valueMin,valueMax);
  printf("normValuesInterval - INFO: normalized histo %s with norm factor = (norm)%6.3f * (custom)%6.3f = (tot) %6.3f\n",cloneH->GetName(),normFactor,factor,normFactor*factor);
  
  if (factor<=0) factor=1.;
  normFactor=normFactor*factor;
  cloneH->Scale(normFactor);
  return cloneH;
  
}
