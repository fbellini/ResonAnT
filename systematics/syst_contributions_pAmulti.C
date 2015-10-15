void syst_contributions_pAmulti(Int_t icent=0, 
				TString bg = "LS",
				TString date="12oct15", 
				Int_t barlow = 0, 
				Int_t smooth = 2,
				Bool_t doFinal = 1 )
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TString deteco = "Combo";
  Int_t ipid=0;
  //set input name
  TString fPath = Form("$HOME/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/systUncert/");
  //For final analysis in september 2014
  //TString fPathCorr = Form("$HOME/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/fit%s_norm%i_BWpoly2_fixedW/", bg.Data(),bg.Contains("EM"));
  //For check of june 2015 and normalization to visible Xsec
  TString fPathCorr = Form("$HOME/alice/resonances/kstar_pA5.02TeV/output_LF5455/multi/central/fit%s_norm%i_BWpoly2_fixedW/", bg.Data(),bg.Contains("EM"));

  TString corrFile = (doFinal? "CORRCHECK_br_normVXS_best_fit_poly2.root" : "CORRECTED_br_best_fit_poly2.root");
  TString hCorrYieldName = Form("hCorrected_%i",icent);//,deteco.Data());
  
  //pA analysis
  Double_t cent[]={ 0.0, 20.0, 40.0, 60.0, 80.0, 100.0};   
  Double_t pt[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00, 12.0, 15.0};  
  Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;   
  Int_t   ncent  = sizeof(cent) / sizeof(cent[0]) - 1;
  TString centLabel=Form("%3.0f-%3.0f%%",cent[icent], cent[icent+1]);
  
  //cosmetics  
  Color_t color[3][6]={kRed, kPink+6, kGreen+1, kViolet-3, kAzure+7, kBlack, //combined final (as the phi)
		       //kOrange+7, kPink+6, kGreen+1, kAzure+1, kBlue+3, kBlack, //combined prelim
		       kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack, //tpc
		       kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof
   
  Int_t marker[3][6]={20, 21, 32, 28, 24, 20, //combined
		      21, 22, 23, 34, 33, 20, //tpc
		      25, 26, 32, 28, 27, 24}; //tof
  
  //create axis to reproduce the binning
  TAxis *ptbins = new TAxis(npt, pt);
  TAxis *centbins = new TAxis(ncent, cent);
  
  //statistical 
  TH1F * statunc = new TH1F(Form("statunc_%i",icent),Form("Statistical uncertainty",icent), npt, pt);
  statunc->SetLineWidth(3);
  statunc->SetLineColor(kGreen+2);
  statunc->SetMarkerColor(kGreen+2);
  statunc->SetLineStyle(1);
  statunc->SetMarkerStyle(0);

  //total systematics
  TH1F * sum2 = new TH1F(Form("sum2_%i",icent),Form("Sum pt-independent contr.",icent), npt, pt);
  sum2->SetLineWidth(3);
  sum2->SetLineColor(kRed);
  sum2->SetMarkerColor(kRed);
  sum2->SetLineStyle(1);
  sum2->SetMarkerStyle(0); 

  //total uncorrelated with pi
  TH1F * sum2_uncorr = new TH1F(Form("sum2_uncorr_%i",icent),Form("p_{T}-uncorrelated sys. uncert."), npt, pt);
  sum2_uncorr->SetLineWidth(3);
  sum2_uncorr->SetLineColor(kGray+2);
  sum2_uncorr->SetMarkerColor(kGray+2);
  sum2_uncorr->SetLineStyle(1);
  sum2_uncorr->SetMarkerStyle(0); 
  //total uncorrelated with pi
  TH1F * sum2_uncorrPi = new TH1F(Form("sum2_uncorrPi_%i",icent),Form("Sum pt-independent contr., #pi uncorrelated",icent), npt, pt);
  sum2_uncorrPi->SetLineWidth(2);
  sum2_uncorrPi->SetLineColor(kRed+1);
  sum2_uncorrPi->SetMarkerColor(kRed+1);
  sum2_uncorrPi->SetLineStyle(3);
  sum2_uncorrPi->SetMarkerStyle(0); 

  //total uncorrelated with K
  TH1F * sum2_uncorrKa = new TH1F(Form("sum2_uncorrKa_%i",icent),Form("Sum pt-independent contr., K uncorrelated",icent), npt, pt);
  sum2_uncorrKa->SetLineWidth(2);
  sum2_uncorrKa->SetLineColor(kBlue+1);
  sum2_uncorrKa->SetMarkerColor(kBlue+1);
  sum2_uncorrKa->SetLineStyle(4);
  sum2_uncorrKa->SetMarkerStyle(0); 

 //total uncorrelated with p
  TH1F * sum2_uncorrPro = new TH1F(Form("sum2_uncorrPro_%i",icent),Form("Sum pt-independent contr., p uncorrelated",icent), npt, pt);
  sum2_uncorrPro->SetLineWidth(2);
  sum2_uncorrPro->SetLineColor(kGreen+1);
  sum2_uncorrPro->SetMarkerColor(kGreen+1);
  sum2_uncorrPro->SetLineStyle(5);
  sum2_uncorrPro->SetMarkerStyle(0); 

  //pt dependent
  TH1F * material = new TH1F(Form("material_%i",icent), Form("Material Budget"), npt, pt);
  TH1F * hadrint = new TH1F(Form("hadrint_%i",icent),Form("Hadronic inter."), npt, pt);
  TH1F * bincount = new TH1F(Form("bincount_%i",icent),Form("Bin counting"), npt, pt);
  TH1F * evtmix = new TH1F(Form("evtmix_%i",icent),Form("Event mixing"), npt, pt);
  TH1F * bgnorm = new TH1F(Form("bgnorm_%i",icent),Form("Bg. normalisation"), npt, pt);
  TH1F * range = new TH1F(Form("range_%i",icent),Form("Fit range",icent), npt, pt);
  TH1F * function = new TH1F(Form("function_%i",icent),Form("Res. bg. fit function",icent), npt, pt);
  TH1F * width = new TH1F(Form("width_%i",icent),Form("Breit Wigner width & range",icent), npt, pt);
  TH1F * pid = new TH1F(Form("pid_%i",icent),Form("PID"), npt, pt);

  //pt independent
  TH1F * vertexing = new TH1F("vertexing","Vertexing efficiency",npt, pt);
  vertexing->SetLineWidth(2);
  vertexing->SetLineColor(kAzure+3);
  vertexing->SetMarkerColor(kAzure+3);
  vertexing->SetLineStyle(3);
  vertexing->SetMarkerStyle(0); 
 
  TH1F * tracking = new TH1F("tracking","Global tracking",npt, pt);
  tracking->SetLineWidth(2);
  tracking->SetLineColor(kOrange);
  tracking->SetMarkerColor(kOrange);
  tracking->SetLineStyle(9);
  tracking->SetMarkerStyle(0); 
 
  TH1F * trackcuts = new TH1F("trackcuts","Track cuts",npt, pt);
  trackcuts->SetLineWidth(2);
  trackcuts->SetLineColor(kAzure+10);
  trackcuts->SetMarkerColor(kAzure+10);
  trackcuts->SetLineStyle(2);
  trackcuts->SetMarkerStyle(0); 

#if 0 //not a systematic effect
  TH1F * tofmatching = new TH1F("tofmatching","TOF matching",npt, pt);  
  tofmatching->SetLineWidth(3);
  tofmatching->SetLineColor(kGreen+2);
  tofmatching->SetMarkerColor(kGreen+2);
  tofmatching->SetLineStyle(5);
  tofmatching->SetMarkerStyle(0); 

  TH1F * tofmatchingMC = new TH1F("tofmatchingMC","TOF matching data/MC",npt, pt);  
  tofmatchingMC->SetLineWidth(2);
  tofmatchingMC->SetLineColor(kSpring+5);
  tofmatchingMC->SetMarkerColor(kSpring+5);
  tofmatchingMC->SetLineStyle(0);
  tofmatchingMC->SetMarkerStyle(0); 

  TH1F * pidResponse = new TH1F(Form("pidResponse_%i",icent),Form("TOF PID response tune"), npt, pt);
  pidResponse->SetLineWidth(3);
  pidResponse->SetLineColor(kOrange+1);
  pidResponse->SetMarkerColor(kOrange+1);
  pidResponse->SetLineStyle(2);
  pidResponse->SetMarkerStyle(0); 

  //pt dependent    
  TFile * fBinCount = TFile::Open(Form("%s/Syst_fromRatio_SystBinCounting_BCtoFIT.root",fPath.Data()));
  TH1F * dummyBC = (TH1F*) fBinCount->Get(Form("hSystVsPtPercentage_100"));
  bincount = (TH1F*) dummyBC->Clone("binCount"); 
  bincount->SetTitle("Bin Counting");
  bincount->SetLineWidth(3);
  bincount->SetLineColor(kViolet-5);
  bincount->SetMarkerColor(kViolet-5);
  bincount->SetLineStyle(3);
  bincount->SetMarkerStyle(0);
  
  TFile * fEvMix = TFile::Open(Form("%s/Syst_fromRatio_SystBG_EMtoLS.root",fPath.Data()));
  TH1F * dummyEM = (TH1F*) fEvMix->Get(Form("hSystVsPtPercentage_100"));
  evtmix = (TH1F*) dummyEM->Clone("evMix"); 
  evtmix->SetTitle("Event mixing");
  evtmix->SetLineWidth(3);
  evtmix->SetLineColor(kBlue);
  evtmix->SetMarkerColor(kBlue);
  evtmix->SetLineStyle(2);
  evtmix->SetMarkerStyle(0);

  TFile * fWidthSyst = TFile::Open(Form("%s/syst_yields_SystWidth.root",fPath.Data()));
  TH1F * dummyW = (TH1F*) fWidthSyst->Get(Form("hSystVsPtPercentageOfCentral_100"));
  width = (TH1F*) dummyW->Clone("width"); 
  width->SetTitle("Width constraints");
  width->SetLineWidth(3);
  width->SetLineColor(kCyan-2);
  width->SetMarkerColor(kCyan-2);
  width->SetLineStyle(2);
  width->SetMarkerStyle(0); 
#endif

  TFile * fMaterial = TFile::Open(Form("%s/systMaterial.root",fPath.Data()));
  TH1F * dummyMT = (TH1F*) fMaterial->Get(Form("hSystVsPtPercentageOfCentral"));
  material = (TH1F*) dummyMT->Clone("material"); 
  material->SetTitle("Material budget");
  material->SetLineWidth(4);
  material->SetLineColor(kPink+2);
  material->SetMarkerColor(kPink+2);
  material->SetLineStyle(3);
  material->SetMarkerStyle(0); 
  TH1F * hMaterial4ratio2Pi = (TH1F*) fMaterial->Get("hSyst4Ratio2Pi");  
  TH1F * hMaterial4ratio2Ka = (TH1F*) fMaterial->Get("hSyst4Ratio2Ka");  
  
  TFile * fHadrSyst = TFile::Open(Form("%s/systHadrInt.root",fPath.Data()));
  TH1F * dummyHI = (TH1F*) fHadrSyst->Get(Form("hSystVsPtPercentageOfCentral"));
  hadrint = (TH1F*) dummyHI->Clone("hadrint"); 
  hadrint->SetTitle("Hadronic int.");
  hadrint->SetLineWidth(2);
  hadrint->SetLineColor(kOrange+1);
  hadrint->SetMarkerColor(kOrange+1);
  hadrint->SetLineStyle(1);
  hadrint->SetMarkerStyle(0); 
  TH1F * hHadrInt4ratio2Pi = (TH1F*) fHadrSyst->Get("hSyst4Ratio2Pi");  
  TH1F * hHadrInt4ratio2Ka = (TH1F*) fHadrSyst->Get("hSyst4Ratio2Ka");  

  TString normfile;
  TString sysHistoName = Form("hSystVsPtPercentageOfCentral_%i", icent);
  if (smooth>0) {
    normfile = Form("smooth_sysNorm_%s_Norm.root",bg.Data());
    sysHistoName = Form("hSystVsPtPercentageOfCentral_%i_smooth%i", icent,smooth);      
      // normfile = Form("smooth_systematics%s_Norm.root",bg.Data()); 
    // if (barlow>0) normfile.ReplaceAll("_Norm.root",Form("_Barlow%i.0_Norm.root",barlow));
    // sysHistoName = Form("hSystVsPtPercentageOfCentral_%i_smooth%i", icent,smooth);
  } else {
    if (barlow>0) normfile = Form("sysNorm_%s_B%i_%i.root", bg.Data(), barlow, icent);
    else normfile = Form("sysNorm_%s_%i.root", bg.Data(), icent);
  }
  TFile * fBgNorm = TFile::Open(Form("%s/%s",fPath.Data(),normfile.Data()));
  TH1F * dummyBN = (TH1F*) fBgNorm->Get(sysHistoName.Data());
  bgnorm = (TH1F*) dummyBN->Clone("bgNorm"); 
  bgnorm->SetTitle("Bg. norm. range");
  bgnorm->SetLineWidth(2);
  bgnorm->SetLineColor(kTeal+2);
  bgnorm->SetMarkerColor(kTeal+2);
  bgnorm->SetLineStyle(1);
  bgnorm->SetMarkerStyle(0);
  
  TString rangefile;// = Form("sysRange_%s_0.root", bg.Data());
  TString sysHistoName = Form("hSystVsPtPercentageOfCentral_%i", icent);
  if (smooth>0) {
    rangefile = Form("smooth_sysRangeW_%s_RangeW.root",bg.Data());
    sysHistoName = Form("hSystVsPtPercentageOfCentral_%i_smooth%i", icent,smooth);      
    // rangefile = Form("smooth_systematics%s_RangeW.root",bg.Data()); 
    // if (barlow>0) rangefile.ReplaceAll("_RangeW.root",Form("_Barlow%i.0_RangeW.root",barlow));
    // sysHistoName.Append(Form("_smooth%i",smooth));
  } else {
    if (barlow>0) rangefile = Form("sysRangeW_%s_B%i_%i.root", bg.Data(), barlow,icent);
    else rangefile = Form("sysRangeW_%s_%i.root", bg.Data(),icent);
  }

  TFile * fRangeSyst = TFile::Open(Form("%s/%s",fPath.Data(),rangefile.Data()));
  TH1F * dummyR = (TH1F*) fRangeSyst->Get(sysHistoName.Data());
  range = (TH1F*) dummyR->Clone("range"); 
  range->SetTitle("Fit range and free width");
  range->SetLineWidth(3);
  range->SetLineColor(kPink+8);
  range->SetMarkerColor(kPink+8);
  range->SetLineStyle(8);
  range->SetMarkerStyle(0);

  TString fcnfile;// = Form("sysFcn_%s_0.root", bg.Data());
  TString sysHistoName = Form("hSystVsPtPercentageOfCentral_%i", icent);
  if (smooth>0) {
    fcnfile = Form("smooth_sysFcn_%s_Fcn.root",bg.Data());
    sysHistoName = Form("hSystVsPtPercentageOfCentral_%i_smooth%i", icent,smooth);      
    // fcnfile = Form("smooth_systematics%s_Fcn.root",bg.Data()); 
    // if (barlow>0) fcnfile.ReplaceAll("_Fcn.root",Form("_Barlow%i.0_Fcn.root",barlow));
    // sysHistoName.Append(Form("_smooth%i",smooth));
  } else {
    if (barlow>0) fcnfile = Form("sysFcn_%s_B%i_%i.root", bg.Data(), barlow,icent);
    else fcnfile = Form("sysFcn_%s_%i.root", bg.Data(),icent);
  }

  TFile * fFuncSyst = TFile::Open(Form("%s/%s",fPath.Data(),fcnfile.Data()));
  TH1F * dummyF = (TH1F*) fFuncSyst->Get(sysHistoName.Data());
  function = (TH1F*) dummyF->Clone("function"); 
  function->SetTitle("Res. bg. fit function");
  function->SetLineWidth(2);
  function->SetLineColor(kBlue+2);
  function->SetMarkerColor(kBlue+2);
  function->SetLineStyle(1);
  function->SetMarkerStyle(0); 

  TString PIDfile;// = Form("sysPID_%s_0.root", bg.Data());
  TString sysHistoName = Form("hSystVsPtPercentageOfCentral_%i", icent);
  if (smooth>0) {
    PIDfile = Form("smooth_sysPID_%s_PID.root",bg.Data());
    sysHistoName = Form("hSystVsPtPercentageOfCentral_%i_smooth%i", icent,smooth);      
    // PIDfile = Form("smooth_systematics%s_PID.root",bg.Data()); 
    // if (barlow>0) PIDfile.ReplaceAll("_PID.root",Form("_Barlow%i.0_PID.root",barlow));
    // sysHistoName.Append(Form("_smooth%i",smooth));
  } else {
    if (barlow>0) PIDfile = Form("sysPID_%s_B%i_%i.root", bg.Data(), barlow, icent);
    else PIDfile = Form("sysPID_%s_%i.root", bg.Data(),icent);
  }

  TFile * fPidSyst = TFile::Open(Form("%s/%s",fPath.Data(),PIDfile.Data()));
  TH1F * dummyP = (TH1F*) fPidSyst->Get(sysHistoName.Data());
  pid = (TH1F*) dummyP->Clone("PID"); 
  pid->SetTitle("PID");
  pid->SetLineWidth(3);
  pid->SetLineColor(kBlack);
  pid->SetMarkerColor(kBlack);
  pid->SetLineStyle(2);
  pid->SetMarkerStyle(0); 

  Double_t tracking_1trk = 0.03;
  Double_t trackCuts_rsn = 0.025;
  Double_t vertexEff = 0.01;
  for (Int_t ii = 0;ii<npt;ii++){
    Int_t ibin = ii+1;
    //material->SetBinContent(ibin, 4);
    tracking->SetBinContent(ibin, tracking_1trk*2.0*100.0);
    trackcuts->SetBinContent(ibin, trackCuts_rsn*100.0);
    if (icent==4) vertexing->SetBinContent(ibin, vertexEff*100.0);
    else  vertexing->SetBinContent(ibin, 0.0);
    // tofmatching->SetBinContent(ibin, 6);
    // tofmatchingMC->SetBinContent(ibin, 4);
    // pidResponse->SetBinContent(ibin, 4);
    // bgnorm->SetBinContent(ibin, 0.02); //use pt-dep uncert.
    // evtplane->SetBinContent(ibin, 0.3); //use pt-dep uncert.
    // evtmix->SetBinContent(ibin, 0.5); //use pt-dep uncert.  
  }

  //sum all contributions in quadrature
  Double_t syst_ptFullyCorr_sum2 = (tracking_1trk*2.0)*(tracking_1trk*2.0); //tracking^2 

  //estimate uncertainty per each pt bin
  for (Int_t ii = 0;ii<npt;ii++){
    Int_t ibin = ii+1;
    
#if 0 //not a systematic effect    
    Double_t bincount_syst = dummyBC->GetBinContent(ibin)/100.;
    Double_t evmix_syst = dummyEM->GetBinContent(ibin)/100.;
    Double_t width_syst = dummyW->GetBinContent(ibin)/100.;
#endif
    Double_t material_syst = dummyMT->GetBinContent(ibin)/100.;
    Double_t hadrint_syst = dummyHI->GetBinContent(ibin)/100.;
    Double_t range_syst = dummyR->GetBinContent(ibin)/100.;
    Double_t func_syst = dummyF->GetBinContent(ibin)/100.;
    Double_t bgnorm_syst = dummyBN->GetBinContent(ibin)/100.;
    Double_t pid_syst = dummyP->GetBinContent(ibin)/100.;
    Double_t hadrint_syst_r2pi = hHadrInt4ratio2Pi->GetBinContent(ibin)/100.;
    Double_t hadrint_syst_r2ka = hHadrInt4ratio2Ka->GetBinContent(ibin)/100.;
    Double_t material_syst_r2pi = hMaterial4ratio2Pi->GetBinContent(ibin)/100.;
    Double_t material_syst_r2ka = hMaterial4ratio2Ka->GetBinContent(ibin)/100.;

    Double_t ptuncorr2 = 0.0; //default 10% se non c'è
    //correlated and uncorrelated uncertainties for particle ratios
    Double_t uncorrelPi2 = 0.0; //uncert. uncorrelated with pi uncert. ^2
    Double_t uncorrelKa2 = 0.0;  //uncert. uncorrelated with K uncert. ^2
  
    //vertexing
    if (icent==4) {
      ptuncorr2+=(vertexEff*vertexEff);
      uncorrelPi2+=(vertexEff*vertexEff);
      uncorrelKa2+=(vertexEff*vertexEff);
    }
    
    //track cuts
    ptuncorr2+=(trackCuts_rsn*trackCuts_rsn);
    uncorrelPi2+=(trackCuts_rsn*trackCuts_rsn);
    uncorrelKa2+=(trackCuts_rsn*trackCuts_rsn);

    //fit range
    if (range_syst>0.0) {
      ptuncorr2+=range_syst*range_syst;
      
      //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=range_syst*range_syst;
      uncorrelKa2+=range_syst*range_syst;
    }
 
    //norm bg
    if (bgnorm_syst>0.0) {
      ptuncorr2+=bgnorm_syst*bgnorm_syst;
 
     //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=bgnorm_syst*bgnorm_syst;
      uncorrelKa2+=bgnorm_syst*bgnorm_syst;
    } 

    //function res bg
    if (func_syst>0.0){
      ptuncorr2+=func_syst*func_syst;      
      //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=func_syst*func_syst;
      uncorrelKa2+=func_syst*func_syst;
    }

    //PID
    if (pid_syst>0.0) {
      ptuncorr2+=pid_syst*pid_syst;
      
      //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=pid_syst*pid_syst;
      uncorrelKa2+=pid_syst*pid_syst;
    }

   //material budget
    if (material_syst>0.0) {
      ptuncorr2+=material_syst*material_syst;

     //removed contribution from correlated uncertainty for particle ratios to pi and K, not for p
      uncorrelPi2+=material_syst_r2pi*material_syst_r2pi;
      uncorrelKa2+=material_syst_r2ka*material_syst_r2ka;
    }

    //hadronic interaction cross section
    if (hadrint_syst>0.0) {
      ptuncorr2+=hadrint_syst*hadrint_syst;

       //removed contribution to correlated uncertainty for particle ratios to pi and K, not for p
      uncorrelPi2+=hadrint_syst_r2pi*hadrint_syst_r2pi;
      uncorrelKa2+=hadrint_syst_r2ka*hadrint_syst_r2ka;
    }
 

    //skip those for which no syst uncert
#if 0 //not a systematic effect        
    bincount_syst=0.0;
    evmix_syst=0.0;
    bgnorm_syst=0.0;
    bgwidth_syst=0.0;
    
    //bin counting
    if (bincount_syst>0.0) {
      ptuncorr2+=bincount_syst*bincount_syst;

      //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=bincount_syst*bincount_syst;
      uncorrelKa2+=bincount_syst*bincount_syst;
    }
    //event mixing
    if (evmix_syst>0.0) {
      ptuncorr2+=evmix_syst*evmix_syst;
      //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=evmix_syst*evmix_syst;
      uncorrelKa2+=evmix_syst*evmix_syst;
    } 
    //width
    if (width_syst>0.0) {
      ptuncorr2+=width_syst*width_syst;
      //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=width_syst*width_syst;
      uncorrelKa2+=width_syst*width_syst;
    } 
#endif

    Double_t totsyst = TMath::Sqrt(syst_ptFullyCorr_sum2+ptuncorr2);
    sum2->SetBinContent(ibin, totsyst*100.);
    Printf("bin %i tot. perc. %6.4f", ibin, totsyst*100.);

    Double_t totsystUncorr = TMath::Sqrt(ptuncorr2);
    sum2_uncorr->SetBinContent(ibin, totsystUncorr*100.);

    Double_t totsystUncorrPi = TMath::Sqrt(uncorrelPi2);
    sum2_uncorrPi->SetBinContent(ibin, totsystUncorrPi*100.);

    Double_t totsystUncorrKa = TMath::Sqrt(uncorrelKa2);
    sum2_uncorrKa->SetBinContent(ibin, totsystUncorrKa*100.);

    Double_t totsystUncorrPro = TMath::Sqrt(ptuncorr2);
    sum2_uncorrPro->SetBinContent(ibin, totsystUncorrPro*100.);
  }
  
  //assign syst err to data
  TFile * fdata = TFile::Open(Form("%s/%s",fPathCorr.Data(),corrFile.Data()));
  TH1F * data = (TH1F*) fdata->Get(Form("%s",hCorrYieldName.Data()));  
  TH1F * data_Wsyst = (TH1F*) data->Clone(Form("%s%i_syst",hCorrYieldName.Data(),icent));
  TH1F * data_Wsyst_uncorr = (TH1F*) data->Clone(Form("%s%i_syst_uncorr",hCorrYieldName.Data(),icent));
  TH1F * data_Wsyst_Wstat = (TH1F*) data->Clone(Form("%s%i_syst_stat",hCorrYieldName.Data(),icent));
  
  TH1F * data_Wsyst_uncorrPi = (TH1F*) data->Clone(Form("%s%i_syst_uncorrPi",hCorrYieldName.Data(),icent));
  TH1F * data_Wsyst_uncorrKa = (TH1F*) data->Clone(Form("%s%i_syst_uncorrKa",hCorrYieldName.Data(),icent));
  TH1F * data_Wsyst_uncorrPro = (TH1F*) data->Clone(Form("%s%i_syst_uncorrPro",hCorrYieldName.Data(),icent));
  TH1F * data_Wsyst_Wstat_uncorrPi = (TH1F*) data->Clone(Form("%s%i_syst_Wstat_uncorrPi",hCorrYieldName.Data(),icent));
  TH1F * data_Wsyst_Wstat_uncorrKa = (TH1F*) data->Clone(Form("%s%i_syst_Wstat_uncorrKa",hCorrYieldName.Data(),icent));
  TH1F * data_Wsyst_Wstat_uncorrPro = (TH1F*) data->Clone(Form("%s%i_syst_Wstat_uncorrPro",hCorrYieldName.Data(),icent));

  for (Int_t ii = 1;ii<npt;ii++){   
    Int_t ibin = ii;
    Double_t yd = data->GetBinContent(ibin);
    Double_t yd_stat = data->GetBinError(ibin);
    Double_t perc = sum2->GetBinContent(ibin);
    Double_t yd_syst = yd*perc/100.;
    statunc->SetBinContent(ibin, (yd>0? (yd_stat*100.0/yd) : 0.0));
    data_Wsyst->SetBinContent(ibin,yd);
    data_Wsyst->SetBinError(ibin, yd_syst);
    data_Wsyst_Wstat->SetBinContent(ibin,yd);
    data_Wsyst_Wstat->SetBinError(ibin,TMath::Sqrt(yd_syst*yd_syst+yd_stat*yd_stat));
    Printf("bin %i     yield = %e     syst err = %e    stat err = %e", ibin, yd, yd_syst, yd_stat);

    Double_t perc_uncorr = sum2_uncorr->GetBinContent(ibin);
    Double_t yd_syst_uncorr = yd*perc_uncorr/100.;
    data_Wsyst_uncorr->SetBinContent(ibin,yd);
    data_Wsyst_uncorr->SetBinError(ibin, yd_syst_uncorr);


    Double_t perc_uncorrPi = sum2_uncorrPi->GetBinContent(ibin);
    Double_t yd_syst_uncorrPi = yd*perc_uncorrPi/100.;
    data_Wsyst_uncorrPi->SetBinContent(ibin,yd);
    data_Wsyst_uncorrPi->SetBinError(ibin, yd_syst_uncorrPi);
    data_Wsyst_Wstat_uncorrPi->SetBinContent(ibin,yd);
    data_Wsyst_Wstat_uncorrPi->SetBinError(ibin,TMath::Sqrt(yd_syst_uncorrPi*yd_syst_uncorrPi+yd_stat*yd_stat));

    Double_t perc_uncorrKa = sum2_uncorrKa->GetBinContent(ibin);
    Double_t yd_syst_uncorrKa = yd*perc_uncorrKa/100.;
    data_Wsyst_uncorrKa->SetBinContent(ibin,yd);
    data_Wsyst_uncorrKa->SetBinError(ibin, yd_syst_uncorrKa);
    data_Wsyst_Wstat_uncorrKa->SetBinContent(ibin,yd);
    data_Wsyst_Wstat_uncorrKa->SetBinError(ibin,TMath::Sqrt(yd_syst_uncorrKa*yd_syst_uncorrKa+yd_stat*yd_stat));

    Double_t perc_uncorrPro = sum2_uncorrPro->GetBinContent(ibin);
    Double_t yd_syst_uncorrPro = yd*perc_uncorrPro/100.;
    data_Wsyst_uncorrPro->SetBinContent(ibin,yd);
    data_Wsyst_uncorrPro->SetBinError(ibin, yd_syst_uncorrPro);
    data_Wsyst_Wstat_uncorrPro->SetBinContent(ibin,yd);
    data_Wsyst_Wstat_uncorrPro->SetBinError(ibin,TMath::Sqrt(yd_syst_uncorrPro*yd_syst_uncorrPro+yd_stat*yd_stat));
  }
  
  //systematic uncertainty plot
  TCanvas *cs=new TCanvas("cs","Systematic uncertainty vs p_{t}", 750,600);
  cs->cd();
  sum2->SetTitle("Total syst. uncert.");
  sum2->GetYaxis()->SetTitle("relative uncert. (%)");
  sum2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  sum2->GetYaxis()->SetRangeUser(0.1, 35);
  if (icent<4)  sum2->GetXaxis()->SetRangeUser(0.0, 14.9);
  else   sum2->GetXaxis()->SetRangeUser(0.0, 5.9);

  sum2->Draw();
  range->Draw("same");
  function->Draw("same");
  pid->Draw("same");
  //  bincount->Draw("same");
  //  evtmix->Draw("same");
  //width->Draw("same");
  bgnorm->Draw("same");
  if (icent==4) vertexing->Draw("same");
  tracking->Draw("same");
  trackcuts->Draw("same");
  material->Draw("same");
  hadrint->Draw("same");
  
  sum2_uncorr->Draw("same");
  // if (ipid) {
  //   tofmatching->Draw("same");
  //   tofmatchingMC->Draw("same");
  //   pidResponse->Draw("same");
  // }
  
  TLegend * autolegry = (TLegend*)gPad->BuildLegend(0.25,0.60,0.88,0.88, Form("V0A multiplicity class %s",centLabel.Data()));
  autolegry->SetFillColor(kWhite);
  autolegry->SetLineColor(kWhite);
  autolegry->SetTextFont(42);
  autolegry->SetNColumns(2); 
  //color data
  data->SetMarkerColor(color[ipid][icent]+1);
  data->SetLineColor(color[ipid][icent]+1);
  data->SetMarkerStyle(marker[ipid][icent]);
  data->SetMarkerSize(0.7);
  data->SetLineWidth(2);
  data->SetTitle(Form("%s (stat. uncert.)",centLabel.Data()));
  data_Wsyst->SetMarkerColor(color[ipid][icent]);	       
  data_Wsyst->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst->SetLineColor(color[ipid][icent]);
  data_Wsyst->SetMarkerStyle(0);
  data_Wsyst->SetOption("E2");
  data_Wsyst->SetTitle(Form("%s (syst. uncert.)",centLabel.Data()));
  TPaveText * barlowLabel = new TPaveText(0.25,0.57,0.88,0.65,"NDC");
  barlowLabel->SetFillColor(kWhite);
  barlowLabel->SetBorderSize(0);
  barlowLabel->SetTextFont(42);
  if (barlow>0) barlowLabel->InsertText(Form("Stat. compatibility test for N_{b} = #Delta_{i}/#sigma_{i} = %i",barlow));
  //else barlowLabel->InsertText("No stat. compatibility test");
  cs->cd();
  barlowLabel->Draw();

  TString imagefilename;
  if (barlow>0) imagefilename = Form("summaryAllSystUncert_%s_B%i_cent%i_%s", bg.Data(), barlow, icent, date.Data()) ;
  else imagefilename = Form("summaryAllSystUncert_%s_cent%i_%s", bg.Data(), icent, date.Data());
  if (smooth>0) imagefilename.ReplaceAll("summaryAllSystUncert",Form("summaryAllSystUncert_smooth%i",smooth));

  cs->Print(Form("%s.png", imagefilename.Data()));
  cs->Print(Form("%s.C", imagefilename.Data()));

  //make-up
  data_Wsyst_uncorr->SetMarkerColor(color[ipid][icent]);	       
  data_Wsyst_uncorr->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst_uncorr->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst_uncorr->SetLineColor(color[ipid][icent]);
  data_Wsyst_uncorr->SetMarkerStyle(0);
  data_Wsyst_uncorr->SetOption("E2");
  data_Wsyst_uncorr->SetTitle(Form("%s (syst. uncert., p_{T}-uncorr.)",centLabel.Data()));
  
  data_Wsyst_uncorrPi->SetMarkerColor(color[ipid][icent]);	       
  data_Wsyst_uncorrPi->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst_uncorrPi->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst_uncorrPi->SetLineColor(color[ipid][icent]);
  data_Wsyst_uncorrPi->SetMarkerStyle(0);
  data_Wsyst_uncorrPi->SetOption("E2");
  data_Wsyst_uncorrPi->SetTitle(Form("%s (syst. uncert., #pi uncorr.)",centLabel.Data()));

  data_Wsyst_uncorrKa->SetMarkerColor(color[ipid][icent]);	       
  data_Wsyst_uncorrKa->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst_uncorrKa->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst_uncorrKa->SetLineColor(color[ipid][icent]);
  data_Wsyst_uncorrKa->SetMarkerStyle(0);
  data_Wsyst_uncorrKa->SetOption("E2");
  data_Wsyst_uncorrKa->SetTitle(Form("%s (syst. uncert., K uncorr.)",centLabel.Data()));

  data_Wsyst_uncorrPro->SetMarkerColor(color[ipid][icent]);	       
  data_Wsyst_uncorrPro->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst_uncorrPro->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst_uncorrPro->SetLineColor(color[ipid][icent]);
  data_Wsyst_uncorrPro->SetMarkerStyle(0);
  data_Wsyst_uncorrPro->SetOption("E2");
  data_Wsyst_uncorrPro->SetTitle(Form("%s (syst. uncert., p uncorr.)",centLabel.Data()));

  data_Wsyst_Wstat->SetMarkerColor(color[ipid][icent]);	       
  data_Wsyst_Wstat->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst_Wstat->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst_Wstat->SetLineColor(color[ipid][icent]);
  data_Wsyst_Wstat->SetMarkerStyle(0);
  data_Wsyst_Wstat->SetTitle(Form("%s (#sqrt{syst^{2}+stat^{2}})",centLabel.Data()));

  data_Wsyst_Wstat_uncorrPi->SetMarkerColor(color[ipid][icent]);	       
  data_Wsyst_Wstat_uncorrPi->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst_Wstat_uncorrPi->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst_Wstat_uncorrPi->SetLineColor(color[ipid][icent]);
  data_Wsyst_Wstat_uncorrPi->SetMarkerStyle(0);
  data_Wsyst_Wstat_uncorrPi->SetTitle(Form("%s (#sqrt{syst^{2}+stat^{2}}, #pi uncorr.)",centLabel.Data()));

  data_Wsyst_Wstat_uncorrKa->SetMarkerColor(color[ipid][icent]);	       
  data_Wsyst_Wstat_uncorrKa->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst_Wstat_uncorrKa->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst_Wstat_uncorrKa->SetLineColor(color[ipid][icent]);
  data_Wsyst_Wstat_uncorrKa->SetMarkerStyle(0);
  data_Wsyst_Wstat_uncorrKa->SetTitle(Form("%s (#sqrt{syst^{2}+stat^{2}}, K uncorr.)",centLabel.Data()));
  
  data_Wsyst_Wstat_uncorrPro->SetMarkerColor(color[ipid][icent]);	       
  data_Wsyst_Wstat_uncorrPro->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst_Wstat_uncorrPro->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst_Wstat_uncorrPro->SetLineColor(color[ipid][icent]);
  data_Wsyst_Wstat_uncorrPro->SetMarkerStyle(0);
  data_Wsyst_Wstat_uncorrPro->SetTitle(Form("%s (#sqrt{syst^{2}+stat^{2}}, p uncorr.)",centLabel.Data()));

  TCanvas *cunc=new TCanvas("cunc","Summary of uncertainty vs p_{t}", 750,600);
  cunc->cd();
  sum2->Draw();
  sum2_uncorr->Draw("HIST same");
  sum2_uncorrPi->Draw("HIST same");
  sum2_uncorrKa->Draw("HIST same");
  statunc->Draw("HIST same");
  TLegend * autolegry2 = (TLegend*)gPad->BuildLegend(0.25,0.65,0.88,0.88);
  autolegry2->SetFillColor(kWhite);
  autolegry2->SetLineColor(kWhite);
  autolegry2->SetTextFont(42);
  autolegry2->SetNColumns(1); 
  autolegry2->Draw();

  //save to out file
  TString outfilename;
  if (barlow>0) outfilename = Form("finalWsyst_%s_B%i_%s_%i.root", bg.Data(),barlow, date.Data(),icent);
  else outfilename = Form("finalWsyst_%s_%s_%i.root", bg.Data(), date.Data(),icent);
  if (smooth>0) outfilename.ReplaceAll("finalWsyst",Form("finalWsyst_smooth%i",smooth));

  TString imagefilename2;
  if (barlow>0) imagefilename2 = Form("summaryUncert_%s_B%i_cent%i_%s", bg.Data(), barlow, icent, date.Data()) ;
  else imagefilename2 = Form("summaryUncert_%s_cent%i_%s", bg.Data(), icent, date.Data());
  if (smooth>0) imagefilename2.ReplaceAll("summaryUncert", Form("summaryUncert_smooth%i",smooth));
  cunc->SaveAs(Form("%s.png", imagefilename2.Data()));
  cunc->SaveAs(Form("%s.C", imagefilename2.Data()));

  
  TFile * fout = new TFile(outfilename.Data(),"recreate");
  fout->cd();
  material->Write();
  vertexing->Write();
  tracking->Write();
  trackcuts->Write();
  // if (ipid) {
  //   tofmatching->Write();
  //   tofmatchingMC->Write();
  //   pidResponse->Write();
  // }
  //  bincount->Write();
  // evtmix->Write();
  pid->Write();
  range->Write();
  bgnorm->Write();
  function->Write();
  hadrint->Write();
  statunc->Write();
  data->Write();
  sum2->Write();
  sum2_uncorr->Write();
  sum2_uncorrPi->Write();
  sum2_uncorrKa->Write();
  sum2_uncorrPro->Write();
  data_Wsyst->Write();
  data_Wsyst_Wstat->Write();
  data_Wsyst_uncorr->Write();
  data_Wsyst_uncorrPi->Write();
  data_Wsyst_Wstat_uncorrPi->Write();
  data_Wsyst_uncorrKa->Write();
  data_Wsyst_Wstat_uncorrKa->Write();
  data_Wsyst_uncorrPro->Write();
  data_Wsyst_Wstat_uncorrPro->Write();
  
  cunc->Write();
  cs->Write();
  fout->Close(); 
  
  return;
  
}

