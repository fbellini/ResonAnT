/*
fbellini@cern.ch - 25/10/2012

Starting from projected histogram per centrality bins and pt bins:
- normalize background
- compute signal (K*+anti-k*), LS and EM background
- subtract bg
- display

*/
const Double_t kBigNumber=1E10;
const Double_t kSmallNumber=1E-10;
const Double_t invMassBinWidth=0.01; //GeV/c^2

TString macroDir = "/Users/fbellini/alice/macro/ResonAnT";
enum ECompType_t {kSum,
		  kHalfSum,
		  kKStar,
		  kAntiKStar};

TH1D*    BgLike(TH1D* hpp = NULL, TH1D* hmm = NULL);
TH1D*    SumSignal(TH1D* hpm = NULL, TH1D* hmp = NULL);
TH1D*    BgMixing(TH1D *hpm = NULL, TH1D* hmp = NULL);
void     HistoMakeUp(TH1D* histo = NULL, Int_t markerStyle=-1);
Double_t BestNormalization(TH1D *hSig, TH1D *hBg, Double_t normMin, Double_t normMax, Double_t normStep, Int_t firstBin, Int_t lastBin);
TH1D*    normalize2Integral(TH1D* hist=NULL, TH1D* histRef=NULL, Double_t factor=1.);
Float_t  GetRangeValuesNormalizationFactor(TH1D* hist=NULL, TH1D* histRef=NULL, Double_t valueMin=1.3, Double_t valueMax=1.5);
TH1D*    normValuesInterval(TH1D* hist=NULL, TH1D* histRef=NULL, Double_t valueMin=1.3, Double_t valueMax=1.5, Double_t factor=1., Double_t *outputFactor = 0);
TH1D*    subtractBackgnd(TH1D* hist=NULL, TH1D* hist2=NULL);


Int_t phiXeXeAnalysis(	 TString projectionFile="proj_20180123_RsnOut.root", 
			 Float_t emNormInf = 1.050, Float_t emNormSup = 1.150, 
			 Int_t ipt=-1, Int_t icent= -1,
			 Bool_t isRebin = 1, 
			 Float_t scaleEMNormFactor = 1.0, Bool_t doNormLS = 0, 
			 Bool_t enableCanvas = 1, Bool_t saveImg = 0, Bool_t useCorrLS = 0)
{
  TGaxis::SetMaxDigits(3);
 
  gStyle->SetTextFont(42);
  TFile * f= TFile::Open(projectionFile.Data());
  if (!f) return 1;
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  Int_t nPtBins = ptbins->GetNbins();
  TAxis *centbins = (TAxis*)f->Get("centbins");
  Int_t nCentBins = centbins->GetNbins();
      
  Bool_t nextBin;
  Short_t display = 1;
  if (ipt>=0 && icent>=0) display = 1;  
  if (ipt<0 && icent>=0) display = 2;
  if (ipt>=0 && icent<0) display = 3;
  if (!enableCanvas) display = 0;

  TCanvas * cdisplay[2];
  Char_t prefix[2][3]={"EM","LS"};
  
  switch (display)
    {
    case 1:
      for (Int_t j=0;j<2;j++){ 
	cdisplay[j] = new TCanvas(Form("display_%s",prefix[j]), Form("%s bg - pt bin %i - centrality bin %i", prefix[j], ipt,icent),600,600);
	cdisplay[j]->Divide(2,1);
      }
      break;
    case 2:
      for (Int_t j=0;j<2;j++){ 
	cdisplay[j] = new TCanvas(Form("cent%i_%s",icent,prefix[j]), Form("%s bg - centrality bin %i", prefix[j],icent),1000,650);
	cdisplay[j]->Divide(nPtBins/3,3);
      }
      break;
    case 3:
      for (Int_t j=0;j<2;j++){ 
	cdisplay[j] = new TCanvas(Form("pt%i_%s",ipt,prefix[j]), Form("%s bg - pt bin %i", prefix[j], ipt),1000,650);
	cdisplay[j]->Divide(3,2);
      }
      break;
    default:
      cdisplay[0] = 0x0;
      cdisplay[1] = 0x0;
      break;
    }
  
  TString fout_name = Form("sub_norm%3.2f-%3.2f_%s", emNormInf,emNormSup, projectionFile.Data());
  
  TFile * fout = new TFile(fout_name.Data(),"recreate");
  ptbins->Write("ptbins");
  centbins->Write("centbins");
  
  Float_t rangeInf = emNormInf, rangeSup = emNormSup;
  Double_t normfactor_em[1],normfactor_ls[1];
  Float_t ptinf=0.0, ptsup=0.0, centinf=0.0, centsup=0.0;
  Int_t treept=-1, treecent=-1;
  TTree *ntree = new TTree("ntree","normFactors");
  ntree->Branch("factor_em", &normfactor_em[0], "factor/D");
  ntree->Branch("factor_ls", &normfactor_ls[0], "factor/D");
  ntree->Branch("inf", &rangeInf, "inf/F");
  ntree->Branch("sup", &rangeSup, "sup/F");
  ntree->Branch("pt_inf", &ptinf, "pt_inf/F");
  ntree->Branch("pt_sup", &ptsup, "pt_sup/F");
  ntree->Branch("pt_bin", &treept, "pt_bin/I");
  ntree->Branch("cent_inf", &centinf, "cent_inf/F");
  ntree->Branch("cent_sup", &centsup, "cent_sup/F");
  ntree->Branch("cent_bin", &treecent, "cent_bin/I");

  for (Int_t icentbin=0;icentbin<nCentBins;icentbin++){
    //if only one bin selected skip the others
    if ((icent>=0) && (icentbin!=icent)) continue;
    
    for (Int_t iptbin=0;iptbin<nPtBins;iptbin++){
      //if only one bin selected skip the others      
      if ((ipt>0) && (iptbin!=ipt)) continue;
      
      Printf("*************************************************");
      Printf("*********** cent %i - pt %i **************", icentbin, iptbin);
      Printf("*************************************************");

      Double_t lowPt=ptbins->GetBinLowEdge(iptbin+1);
      Double_t upPt=ptbins->GetBinLowEdge(iptbin+2);
      Double_t lowC=centbins->GetBinLowEdge(icentbin+1);
      Double_t upC=centbins->GetBinLowEdge(icentbin+2);
  
      //assign variables for norm factors tree
      ptinf=lowPt;
      ptsup=upPt;
      centinf=lowC;
      centsup=upC;
      treept=iptbin; 
      treecent=icentbin;

      //get input histos
      TString hs1_name=Form("Data_UnlikePM_ptBin%02i_centBin%02i",iptbin,icentbin);
      TString hb1em_name=Form("Data_MixingPM_ptBin%02i_centBin%02i",iptbin,icentbin);
      TString hb1ls_name=Form("Data_LikePP_ptBin%02i_centBin%02i",iptbin,icentbin);
      TString hb2ls_name=Form("Data_LikeMM_ptBin%02i_centBin%02i",iptbin,icentbin);

      // signal
      TH1D * hs1 = (TH1D*) f->Get(hs1_name.Data())->Clone();
      TH1D * hs = (TH1D*)hs1->Clone("");
      Printf("Signal+background SE phi: %s",hs->GetName());
      hs->SetNameTitle(Form("Signal_ptBin%02i_centBin%02i",iptbin,icentbin),
		       Form("S+Bg: %4.2f<p_{T}<%5.2f GeV/c (%2.0f-%2.0f%%)",lowPt,upPt,lowC,upC));

      // mixed-event background
      TH1D * hb1em = (TH1D*) f->Get(hb1em_name.Data())->Clone();
      TH1D * hbem = (TH1D*)hb1em->Clone();
      Printf("EM background for phi: %s",hbem->GetName());
      hbem->SetNameTitle(Form("Mixing_ptBin%02i_centBin%02i",iptbin,icentbin),
			 Form("EM bg: %4.2f<p_{T}<%5.2f GeV/c (%2.0f-%2.0f%%)",lowPt,upPt,lowC,upC));

      // like-sign background
      TH1D * hb1ls = (TH1D*) f->Get(hb1ls_name.Data())->Clone();
      TH1D * hb2ls = (TH1D*) f->Get(hb2ls_name.Data())->Clone();
      TH1D * hbls = (TH1D*) BgLike(hb1ls, hb2ls);
      Printf("LS background for phi: %s",hbls->GetName());
      hbls->SetNameTitle(Form("Like_ptBin%02i_centBin%02i",iptbin,icentbin), 
			 Form("LS bg: %4.2f<p_{T}<%5.2f GeV/c (%2.0f-%2.0f%%)",lowPt,upPt,lowC,upC));
      
      //rebin 
      if (isRebin) {
	hs->Rebin(3);
	hbem->Rebin(3);
	hbls->Rebin(3);
      }
      
      //signal make up
      TH1D * hs_bis = (TH1D*) hs->Clone();
      TH1D * hs_tris = (TH1D*) hs->Clone();
      HistoMakeUp(hs,20);  
      HistoMakeUp(hs_bis, 20);  
      HistoMakeUp(hs_tris, 20);  
      
      Int_t inorm1 =  hbem->GetXaxis()->FindBin(1.1);
      Int_t inorm2 =  hbem->GetXaxis()->FindBin(1.2);
      
      //MEB make up
      TH1D * hbem_norm;
      if ((emNormInf<=0) || (emNormSup<=0)) {
	hbem_norm = (TH1D*) hbem->Clone(Form("norm_%s", hbem->GetName()));
	Double_t norm_best_factor = BestNormalization(hs, hbem_norm, 0.005, 0.500, 0.001, inorm1, inorm2);
	hbem_norm->Scale(norm_best_factor);
	Printf(":::: EM Best normalization by factor = %5.3f", norm_best_factor);
	normfactor_em[0]=norm_best_factor;
      } else {
	hbem_norm = (TH1D*) normValuesInterval(hbem, hs, emNormInf, emNormSup, scaleEMNormFactor, normfactor_em); 
	Printf(":::: EM Normalization in [%3.2f,%3.2f] by factor = %5.3f", emNormInf, emNormSup, normfactor_em);
      }
      hbem_norm->SetName(Form("norm_%s", hbem->GetName()));
      HistoMakeUp(hbem_norm, 1);  
      
      //LSB make up
      TH1D * hbls_norm = (TH1D*) hbls->Clone(Form("norm_%s", hbls->GetName()));
      if (doNormLS) {
	if ((emNormInf<=0) || (emNormSup<=0)) {
	  Double_t norm_best_factor = BestNormalization(hs, hbls_norm, 0.005, 0.500, 0.001, inorm1, inorm2);
	  hbls_norm->Scale(norm_best_factor);
	  Printf(":::: LS Best normalization by factor = %5.3f", norm_best_factor);
	  normfactor_ls[0] = norm_best_factor;      
	} else {
	  hbls_norm = (TH1D*) normValuesInterval(hbls, hs, emNormInf, emNormSup, scaleEMNormFactor, normfactor_ls); 
	  Printf(":::: LS Normalization in [%3.2f,%3.2f] by factor = %5.3f", emNormInf, emNormSup, normfactor_ls);
	}
      }
      hbls_norm->SetName(Form("norm_%s", hbls->GetName()));
      HistoMakeUp(hbls,1);
      HistoMakeUp(hbls_norm,1);
      
      fout->cd();
      hs->Write();
      hbem->Write();
      hbem_norm->Write();
      hbls->Write();
      hbls_norm->Write();
      
      /***************************************************************************
	 EVENT MIXING BACKGROUND
	 background subtraction formula: sub = s-b
	 s = ls_pm + ls_mp       (only if K*+anti-K* is enabled, otherwise only one of the two ls)
	 b = (em_pm + em_mp)*0.5      
      ***************************************************************************/
      TH1D * sub_em =(TH1D*) subtractBackgnd(hs, hbem_norm);
      sub_em->SetName(Form("sub_%s",hbem_norm->GetName()));
      sub_em->SetTitle(Form("S+res.Bg (EM): %4.2f<p_{T}<%5.2fGeV/c (%2.0f-%2.0f%%)",lowPt,upPt,lowC,upC));
      HistoMakeUp(sub_em,28);
      //fill normalization tree
      ntree->Fill();
      
      //save to file
      fout->cd();
      sub_em->Write();

      //draw
      switch (display){
      case 1:
	cdisplay[0]->cd();
	hs_tris->Draw();
	hbem_norm->Draw("same");
	sub_em->Draw("same");
	break;
      case 2:
	cdisplay[0]->cd(iptbin+1);
	hs_tris->Draw();
	hbem_norm->Draw("same");
	sub_em->Draw("same");
	break;
      case 3:
      cdisplay[0]->cd(icentbin+1);
	hs_tris->Draw();
	hbem_norm->Draw("same");
	sub_em->Draw("same");
	break;
      default:
	break;
      }
      /***************************************************************************
	LIKE SIGN BACKGROUND
	background subtraction formula: sub = s-b
	s = ls_pm + ls_mp       (only if K*+anti-K* is enabled, otherwise only one of the two ls)
	b = Sqrt (ls_pp * ls_mm)
      ***************************************************************************/      
       //correct by EM +-/++ and -+/-- factors 
      /*
	TString corrFactorLSfile = "/Users/bellini/alice/resonances/myPhi/check_LSwithMixingBg/correctionLSfromMix.root";
	if (useCorrLS){
	TFile * fcorr = TFile::Open(corrFactorLSfile.Data(),"read");
	Printf("========== CORRECTION FOR LS BG FROM EM USED \n reading correction factors from file %s", corrFactorLSfile.Data());
	TH1D * corrp = (TH1D*) fcorr->Get("corrp");
	TH1D * corrm = (TH1D*) fcorr->Get("corrm");
	hb1ls->Multiply(corrp);
	hb2ls->Multiply(corrm);
	TH1D* hb1ls_copy = (TH1D*)hb1ls->Clone(Form("corrected_%s",hb1ls->GetName()));
	TH1D* hb2ls_copy = (TH1D*)hb2ls->Clone(Form("corrected_%s",hb2ls->GetName()));
	TH1D* hs1_copy = (TH1D*)hs1->Clone(Form("subCorrLS_%s",hs1->GetName()));
	TH1D* hs2_copy = (TH1D*)hs2->Clone(Form("subCorrLS_%s",hs2->GetName()));
	if (isRebin) {
	  hb1ls_copy->Rebin(2);
	  hb2ls_copy->Rebin(2);
	  hs1_copy->Rebin(2);
	  hs2_copy->Rebin(2);
	}
	hs1_copy->Add(hb1ls_copy,-1);
	hs2_copy->Add(hb2ls_copy,-1);
	fout->cd();
	hb1ls_copy->Write();
	hb2ls_copy->Write();
	hs1_copy->Write();
	hs2_copy->Write();
      }
      */
      //subtract
      TH1D * sub_ls =(TH1D*) subtractBackgnd(hs_bis, hbls_norm);
      sub_ls->SetName(Form("sub_%s",hbls_norm->GetName()));
      HistoMakeUp(sub_ls, 20);
      sub_ls->SetTitle(Form("S+Res.bg(LS): %4.2f<p_{T}<%5.2fGeV/c (%2.0f-%2.0f%%)",lowPt,upPt,lowC,upC));
      //save to file
      fout->cd();
      sub_ls->Write();
          
      //draw
      switch (display){
      case 1:
	cdisplay[1]->cd();
	hs_tris->Draw();
	hbls->Draw("same");
	hbls_norm->Draw("same");
	sub_ls->Draw("same");
	break;
      case 2:
	cdisplay[1]->cd(iptbin+1);
	hs_tris->Draw();
	hbls->Draw("same");
	hbls_norm->Draw("same");
	sub_ls->Draw("same");
	break;
      case 3:
      cdisplay[1]->cd(icentbin+1);
	hs_tris->Draw();
	hbls->Draw("same");
	hbls_norm->Draw("same");
	sub_ls->Draw("same");
	break;
      default:
	break;
      }
    }
    
    if (cdisplay[0]){
    fout->cd();
    cdisplay[0]->Write();
    if (saveImg) cdisplay[0]->SaveAs(Form("%s_canvas.png",cdisplay[0]->GetName()));
    }
    if (cdisplay[1]){
      fout->cd();
      cdisplay[1]->Write();
      if (saveImg) cdisplay[1]->SaveAs(Form("%s_canvas.png",cdisplay[1]->GetName()));
    }
    cdisplay->Update();
    cout << "Continue? (1/0)" << endl;
    cin >> nextBin;
     switch (nextBin) 
       {
       case 1:
     	 continue;
       default:
     	 return;
       }
  }
  fout->cd();
  ntree->Write();
  Printf("========================================= Succesfully saved output file %s",fout->GetName());
  return 0;
}

/*****************************************/
TH1D* BgLike(TH1D *hpp, TH1D *hmm)
{
   Int_t i, nbins = hpp->GetNbinsX();
   Double_t y1, y2, e1, e2, bg, err;
   
   TH1D *out = (TH1D*)hpp->Clone();
   
   out->SetName("tmp");
   out->Clear();
   out->SetEntries(hpp->GetEntries() + hmm->GetEntries());
   
   for (i = 1; i <= nbins; i++) {
      y1 = hpp->GetBinContent(i);
      y2 = hmm->GetBinContent(i);
      e1 = hpp->GetBinError(i);
      e2 = hmm->GetBinError(i);
      
      if (y1 < 1.0 || y2 < 1.0) {
         out->SetBinContent(i, 0.0);
         out->SetBinError(i, 0.0);
      } else {
	bg =  2.0*TMath::Sqrt(y1*y2);//multiply by two if signal is not divided by two bg = 2.0 * TMath::Sqrt(y1*y2);
	err = y1*y1*e2*e2 + y2*y2*e1*e1;
	err /= y1*y2;
	out->SetBinContent(i, bg);
	out->SetBinError(i, TMath::Sqrt(err));
      }
   }
   
   return out;
}

/*****************************************/
TH1D* SumSignal(TH1D *hpm, TH1D *hmp)
{
   Int_t i, nbins = hpm->GetNbinsX();
   Double_t y1, y2, e1, e2, sig, err;   
   TH1D *out = (TH1D*)hpm->Clone();  
   out->SetName("tmp");
   out->Clear();
   out->SetEntries(hpm->GetEntries() + hmp->GetEntries());   
   for (i = 1; i <= nbins; i++) {
      y1 = hpm->GetBinContent(i);
      y2 = hmp->GetBinContent(i);
      e1 = hpm->GetBinError(i);
      e2 = hmp->GetBinError(i);
      if (y1 < 1.0 || y2 < 1.0) {
         out->SetBinContent(i, 0.0);
         out->SetBinError(i, 0.0);
      } else {
	sig = (y1+y2); //divide by 2 if kHalfSum: done in the main function
	err = e2*e2 + e1*e1;
	out->SetBinContent(i, sig);
	out->SetBinError(i, TMath::Sqrt(err));
      }
   }   
   return out;
}
/*****************************************/
TH1D* BgMixing(TH1D *hpm, TH1D *hmp)
{
   Int_t i, nbins = hpm->GetNbinsX();
   Double_t y1, y2, e1, e2, bg, err;   
   TH1D *out = (TH1D*)hpm->Clone();  
   out->SetName("tmp");
   out->Clear();
   out->SetEntries(hpm->GetEntries() + hmp->GetEntries());   
   for (i = 1; i <= nbins; i++) {
      y1 = hpm->GetBinContent(i);
      y2 = hmp->GetBinContent(i);
      e1 = hpm->GetBinError(i);
      e2 = hmp->GetBinError(i);
      if (y1 < 1.0 || y2 < 1.0) {
         out->SetBinContent(i, 0.0);
         out->SetBinError(i, 0.0);
      } else {
         bg = (y1+y2);
         err = e2*e2 + e1*e1;
         out->SetBinContent(i, bg);
         out->SetBinError(i, TMath::Sqrt(err));
      }
   }   
   return out;
}

/*****************************************/
void HistoMakeUp(TH1D*histo, Int_t markerStyle)
{
  if (!histo) return;
  TString hname = Form("%s",histo->GetName());

  // histo->Rebin(2);

  TAxis * xaxis = (TAxis*)histo->GetXaxis();
  xaxis->SetTitle("#it{M}_{KK} (GeV/#it{c}^{2})");

  TAxis * yaxis = (TAxis*)histo->GetYaxis();
  yaxis->SetTitle("Counts / (0.001 GeV/#it{c}^{2})");

  Color_t color;
  if (hname.Contains("Signal") || hname.Contains("UnlikePM")){
    color = kBlack;
    markerStyle = 20;
  }
  if (hname.Contains("Mixing")){
    color = kOrange-3;
    markerStyle = 24;
  }
  if (hname.Contains("Like")) {
    color = kAzure+1;
    markerStyle = 25;
  }
  if (hname.Contains("Mixing") && hname.Contains("norm")) {
    color = kRed;
    markerStyle = 24;
  }
  if (hname.Contains("Like") && hname.Contains("norm")) {
    color = kBlue;
    markerStyle = 25;
  }
  if (hname.Contains("Like") && hname.Contains("sub")) {
    color = kBlack;
    markerStyle = 20;
  }
  if (hname.Contains("Mixing") && hname.Contains("sub")) {
    color = kGray;
    markerStyle = 2;
  }
  histo->SetLineColor(color);
  histo->SetLineWidth(1);
  histo->SetMarkerColor(color);
  histo->SetFillColor(kWhite);
  histo->SetFillStyle(0);
  if (markerStyle>=0) histo->SetMarkerStyle(markerStyle);
 
  return;
}

//-----------------------------------------------------------------------
Double_t BestNormalization(TH1D *hSig, TH1D *hBg, Double_t normMin, Double_t normMax, Double_t normStep, Int_t firstBin, Int_t lastBin)
{
//
// search best factor which leaves all subtractions positive within error
//
   Double_t norm, bestNorm = normMin;
   Double_t min, bestMin  = 1E20;
   TH1D *htmp = NULL;
   for (norm = normMax; norm >= normMin; norm -= normStep) {
      htmp = (TH1D*)hSig->Clone("tmp");
      htmp->Add(hBg, -norm);
      min = 1E20;
      for (Int_t k = firstBin; k <= lastBin; k++) {
         Double_t y  = htmp->GetBinContent(k);
         if (htmp->GetBinError(k) < 1E-6) continue;
         y -= htmp->GetBinError(k);
         if (y < min) min = y;
      }
      if (min < 0.0)  
         continue;
      else if (min < bestMin) {
         bestMin  = min;
         bestNorm = norm;
      }
   }
   Double_t IntS = htmp->Integral(firstBin, lastBin);
   Double_t IntB = hBg->Integral(firstBin, lastBin);
   Printf(":::: Normalization: best -> intS = %f, intB = %f, norm = %8.4f", IntS, IntB, bestNorm);
   return bestNorm;
}


//---------------------------------------------------------------------------------
TH1D* normalize2Integral(TH1D* hist, TH1D*histRef, Double_t factor)
{
  
  if ((!hist)||(!histRef)){
    printf("normalize2Integral - ERROR: invalid histogram address passed to normalization function\n");
    return 0;
  }
  
  Int_t nbins=hist->GetXaxis()->GetNbins();
  Int_t nbinsRef=histRef->GetXaxis()->GetNbins();
  
  if (nbins!=nbinsRef){
    printf("normalize2Integral - ERROR: histgrams have different binning. Doing nothing.\n");
    return 0;    
  }
  
  TH1D * cloneH = (TH1D*) hist->Clone();
  Double_t integralS= histRef->Integral();
  Double_t integralB= hist->Integral();
  
  if (integralS==0){
    printf("normalize2Integral - INFO: signal + backgnd histogram integral = 0. Skipping normalization.\n");
    return 0;
  }
  
  if (integralB==0){
    printf("normalize2Integral - INFO: backgnd histogram integral = 0. Skipping normalization.\n");
    return 0;
  }
  
  if (factor<=0) factor=1.;
  Double_t normFactor= factor*integralS/integralB;
  cloneH->Scale(normFactor);
  return cloneH;
}

//---------------------------------------------------------------------------------
Float_t GetRangeValuesNormalizationFactor(TH1D* hist, TH1D* histRef, Double_t valueMin, Double_t valueMax){
  
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
  
  // printf("--------------- Normalizing histogram in interval [ %5.2f - %5.2f ]\n",valueMin,valueMax);
  //  printf("                corresponding to (x-axis) bin interval [ %i - %i ]\n",ibinmin,ibinmax);
  
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
  //  printf("GetRangeValuesNormalizationFactor - INFO: Normalization factor estimate = %6.3f\n",normFactor);
  return normFactor;
  
}

//---------------------------------------------------------------------------------
TH1D* normValuesInterval(TH1D* hist, TH1D* histRef, Double_t valueMin, Double_t valueMax, Double_t factor, Double_t * outputFactor)
{
  
  if ((!hist)||(!histRef)){
    printf("normValuesInterval - ERROR: invalid histgram address passed to normalization function\n");
    return 0;
  }
    printf("normValuesInterval - INFO: Normalizing background in IM interval [ %5.2f - %5.2f ]\n",valueMin,valueMax);
  
  TH1D * cloneH = (TH1D*) hist->Clone();  
  Double_t normFactor= GetRangeValuesNormalizationFactor(hist,histRef,valueMin,valueMax);
  printf("normValuesInterval: normalized histo %s with norm factor = (norm)%6.3f * (custom)%6.3f = (tot) %6.3f\n",cloneH->GetName(),normFactor,factor,normFactor*factor);  
  if (factor<=0) factor=1.;
  normFactor=normFactor*factor;
  cloneH->Scale(normFactor);
  outputFactor[0] = normFactor; 
  return cloneH;
  
}

//---------------------------------------------------------------------------------
TH1D* subtractBackgnd(TH1D* hist, TH1D* hist2)
{
  //hist - hist2
  if ((!hist)||(!hist2)){
    printf("subtractBackgnd - ERROR: invalid histgram address passed to base subtraction function\n");
    return 0x0;
  }
  
  Int_t nbins=hist->GetXaxis()->GetNbins();
  Int_t nbinsRef=hist2->GetXaxis()->GetNbins();
  
  if (nbins!=nbinsRef){
    printf("subtractBackgnd - ERROR: histgrams have different binning. Doing nothing.\n");
    return 0x0;    
  }
  
  hist->Add(hist2,-1);
  return hist;
}
