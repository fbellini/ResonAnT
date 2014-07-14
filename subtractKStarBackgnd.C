#define DEBUG 0

const Double_t kBigNumber=1E10;
const Double_t kSmallNumber=1E-10;
TString macroDir = "/Users/bellini/alice/macro/kstar";
const Double_t invMassBinWidth=0.01; //GeV/c^2

//---------------------------------------------------------------------------------
TString DefineInputHistoName(Int_t ipair=-1, Int_t iptbin=-1, Int_t icentbin=-1)
{

  //builds histogram name according to fixed pattern
  if ((ipair<0) ||(iptbin<0) ||(icentbin<0) ){
    printf("InputHistoName::ERROR - Invalid bin chosen. \n");
    return;
  }
  TString name="";
  //name.Clear();
  name.Form("Data_%s_ptBin%02i_centBin%02i", pairsName[ipair].Data(),iptbin, icentbin);
#if DEBUG 
 printf("Histo name defined: %s\n",name.Data());
#endif
  return name;
}



//---------------------------------------------------------------------------------
TH1D* subtractMixingBackgnd(TH1D* hist=NULL, TH1D* hist2=NULL, Float_t factor=1.0){
  
  if ((!hist)||(!hist2)){
    printf("subtractMixingBackgnd - ERROR: invalid histogram address passed to subtraction function\n");
    return 0x0;
  }
  
  TH1D *backgnd=new TH1D(); 
  if (factor>0.0){
    //  backgnd=(TH1D*) normalize2Integral(hist2, hist, normFactor);    
    backgnd=(TH1D*)normValuesInterval(hist2, hist, 1.1, 1.5, factor);
  } else {
    backgnd=(TH1D*)hist2->Clone(); 
  }
  // backgnd->SaveAs(Form("%s.root",backgnd->GetName()));
  // hist->SaveAs(Form("%s.root",hist->GetName()));
  return subtractBackgnd(hist, backgnd);  
}

//---------------------------------------------------------------------------------
TH1D* subtractLikeSignBackgnd(TH1D* hist=NULL, TH1D* hist2=NULL, Float_t factor=1.0){
  
  if ((!hist)||(!hist2)){
    printf("subtractLikeSignBackgnd - ERROR: invalid histgram address passed to subtraction function\n");
    return 0x0;
  }  
  /*
    use if not normalize
  if (factor!=0.0) hist2->Scale(factor);
  return subtractBackgnd(hist, hist2);  
  */
  TH1D *backgnd=new TH1D(); 
  if (factor>0.0){
    backgnd=(TH1D*)normValuesInterval(hist2, hist, 1.1, 1.5, factor);
  } else {
    backgnd=(TH1D*)hist2->Clone(); 
  }
  return subtractBackgnd(hist, backgnd);  
  
}

//---------------------------------------------------------------------------------
TH1D* subtractBackgnd(TH1D* hist=NULL, TH1D* hist2=NULL){
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

//---------------------------------------------------------------------------------
TH1D* normalize2Integral(TH1D* hist=NULL, TH1D*histRef=NULL, Double_t factor=1.){
  
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
Float_t GetRangeValuesNormalizationFactor(TH1D* hist=NULL, TH1D* histRef=NULL, Double_t valueMin=1.3, Double_t valueMax=1.5){
  
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
    printf("normValuesInterval - ERROR: invalid histgram address passed to normalization function\n");
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

//---------------------------------------------------------------------------------
TH1D* normValuesInterval_old(TH1D* hist=NULL, TH1D* histRef=NULL, Double_t valueMin=1.1, Double_t valueMax=1.5, Double_t factor=1.){
  
  if ((!hist)||(!histRef)){
    printf("ERROR: invalid histgram address passed to normalization function\n");
    return 0;
  }
  
  Int_t nbins=hist->GetXaxis()->GetNbins();
  Int_t nbinsRef=histRef->GetXaxis()->GetNbins();
  
  if (nbins!=nbinsRef){
    printf("ERROR: histgrams have different binning. Doing nothing.\n");
    return 0;    
  }
  
  Int_t ibinmin=-1,ibinmax=-1;

  if ((valueMin<=kSmallNumber) || (valueMin>=valueMax)) {
    valueMin = hist->GetXaxis()->GetBinLowEdge(1);
    ibinmin=1;
    printf("INFO:  min value used for normalization range as default.\n");
  } else {
    ibinmin= 1 + ( ( valueMin - (hist->GetXaxis()->GetBinLowEdge(1)) ) / invMassBinWidth );
  } 
  
  if ((valueMax>=kBigNumber) || (valueMax<valueMin)) {
    valueMax = hist->GetXaxis()->GetBinUpEdge(nbins);
    ibinmax=nbins;
    printf("INFO:  max value used for normalization range as default.\n");  
  } else {    
    ibinmax = nbins - ( (hist->GetXaxis()->GetBinLowEdge(nbins)) - valueMax )/invMassBinWidth;  
  }
  
  printf("--------------- Normalizing background in IM interval [ %5.2f - %5.2f ]\n",valueMin,valueMax);
  printf("                corresponding to IM (x-axis) bins [ %i - %i ]\n",ibinmin,ibinmax);
  
  TH1D * cloneH = (TH1D*) hist->Clone();
  Double_t integralS= histRef->Integral(ibinmin,ibinmax);
  Double_t integralB= hist->Integral(ibinmin,ibinmax);
  
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
  printf("normValuesInterval - INFO: normalized histo %s with norm factor = %6.3f\n",cloneH->GetName(),normFactor);
  cloneH->Scale(normFactor);
  return cloneH;
  
}
//-----------------------------------------------------------------------
Double_t BestNormalization(TH1D *hSig, TH1D *hBg, Double_t normMin, Double_t normMax, Double_t normStep, Int_t firstBin, Int_t lastBin)
{
//
// search best factor which leaves all subtractions positive within error
//
   Double_t norm, bestNorm = normMin;
   Double_t min, bestMin  = 1E20;
   for (norm = normMax; norm >= normMin; norm -= normStep) {
      TH1D *htmp = (TH1D*)hSig->Clone("tmp");
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
   Double_t IntS=htmp->Integral(firstBin, lastBin);
   Double_t IntB = hBg->Integral(firstBin, lastBin);
   Printf("Normalization: best -> intS = %f, intB = %f, norm = %8.4f", IntS, IntB, bestNorm);
   return bestNorm;
}
