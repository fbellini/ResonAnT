void syst_contributions_pAmulti(Int_t icent=0,TString date="05apr14")
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TString deteco = "Combo";
  Int_t ipid=0;
  //set input name
  TString fPath = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/systUncert/");
  TString fPathCorr = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/fitEM_norm1_BWpoly2_fixedW/");
  TString corrFile = "CORRECTED_br_best_fit_poly2.root";
  TString hCorrYieldName = Form("hCorrected_%i",icent);//,deteco.Data());
  
  //pA analysis
  Double_t cent[]={ 0.0, 100.0};   
  Double_t pt[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00, 12.0, 15.0};  
  Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;   
  Int_t   ncent  = sizeof(cent) / sizeof(cent[0]) - 1;
  TString centLabel=Form("%i-%i%%",icent*20,(icent+1)*20);
  
  //cosmetics  
  Color_t color[3][6]={kOrange+7, kPink+6, kGreen+1, kAzure+1, kBlue+4, kBlack, //combined
		       kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack, //tpc
		       kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof
   
  Int_t marker[3][6]={21, 22, 32, 28, 24, 20, //combined
		      21, 22, 23, 34, 33, 20, //tpc
		      25, 26, 32, 28, 27, 24}; //tof
  
  //create axis to reproduce the binning
  TAxis *ptbins = new TAxis(npt, pt);
  TAxis *centbins = new TAxis(ncent, cent);
  
  //total 
  TH1F * sum2 = new TH1F(Form("sum2_%i",icent),Form("Sum pt-independent contr.",icent), npt, pt);
  sum2->SetLineWidth(3);
  sum2->SetLineColor(kRed);
  sum2->SetMarkerColor(kRed);
  sum2->SetLineStyle(1);
  sum2->SetMarkerStyle(0); 

  //total uncorrelated with pi
  TH1F * sum2_uncorrPi = new TH1F(Form("sum2_uncorrPi_%i",icent),Form("Sum pt-independent contr., #pi uncorrelated",icent), npt, pt);
  sum2_uncorrPi->SetLineWidth(3);
  sum2_uncorrPi->SetLineColor(kRed);
  sum2_uncorrPi->SetMarkerColor(kRed);
  sum2_uncorrPi->SetLineStyle(1);
  sum2_uncorrPi->SetMarkerStyle(0); 

  //total uncorrelated with K
  TH1F * sum2_uncorrKa = new TH1F(Form("sum2_uncorrKa_%i",icent),Form("Sum pt-independent contr., K uncorrelated",icent), npt, pt);
  sum2_uncorrKa->SetLineWidth(3);
  sum2_uncorrKa->SetLineColor(kRed);
  sum2_uncorrKa->SetMarkerColor(kRed);
  sum2_uncorrKa->SetLineStyle(1);
  sum2_uncorrKa->SetMarkerStyle(0); 

 //total uncorrelated with p
  TH1F * sum2_uncorrPro = new TH1F(Form("sum2_uncorrPro_%i",icent),Form("Sum pt-independent contr., p uncorrelated",icent), npt, pt);
  sum2_uncorrPro->SetLineWidth(3);
  sum2_uncorrPro->SetLineColor(kRed);
  sum2_uncorrPro->SetMarkerColor(kRed);
  sum2_uncorrPro->SetLineStyle(1);
  sum2_uncorrPro->SetMarkerStyle(0); 

  //pt dependent
  TH1F * material = new TH1F(Form("material_%i",icent), Form("Material Budget"), npt, pt);
  TH1F * hadrint = new TH1F(Form("hadrint_%i",icent),Form("Hadronic inter."), npt, pt);
  TH1F * bincount = new TH1F(Form("bincount_%i",icent),Form("Bin counting"), npt, pt);
  TH1F * evtmix = new TH1F(Form("evtmix_%i",icent),Form("Event mixing"), npt, pt);
  TH1F * bgnorm = new TH1F(Form("bgnorm_%i",icent),Form("LS bg. normalisation"), npt, pt);
  TH1F * range = new TH1F(Form("range_%i",icent),Form("Fit range",icent), npt, pt);
  TH1F * function = new TH1F(Form("function_%i",icent),Form("Res. bg. fit function",icent), npt, pt);
  TH1F * width = new TH1F(Form("width_%i",icent),Form("Breit Wigner width and range",icent), npt, pt);
  TH1F * pid = new TH1F(Form("pid_%i",icent),Form("PID"), npt, pt);

  //pt independent
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

#if 0 //not systematic 
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
  TH1F * hHadrInt4ratio2Pi = (TH1F*) fHadrInt->Get("hSyst4Ratio2Pi");  
  TH1F * hHadrInt4ratio2Ka = (TH1F*) fHadrInt->Get("hSyst4Ratio2Ka");  

  TFile * fBgNorm = TFile::Open(Form("%s/systematicsEM_Norm_%i.root",fPath.Data(),icent));
  TH1F * dummyBN = (TH1F*) fBgNorm->Get(Form("hSystVsPtPercentageOfCentral_%i",icent));
  bgnorm = (TH1F*) dummyBN->Clone("bgNorm"); 
  bgnorm->SetTitle("Bg. norm. range");
  bgnorm->SetLineWidth(2);
  bgnorm->SetLineColor(kTeal+2);
  bgnorm->SetMarkerColor(kTeal+2);
  bgnorm->SetLineStyle(1);
  bgnorm->SetMarkerStyle(0);
  
  TFile * fRangeSyst = TFile::Open(Form("%s/systematicsEM_RangeW_%i.root",fPath.Data(),icent));
  TH1F * dummyR = (TH1F*) fRangeSyst->Get(Form("hSystVsPtPercentageOfCentral_%i",icent)));
  range = (TH1F*) dummyR->Clone("range"); 
  range->SetTitle("Fit range and free width");
  range->SetLineWidth(2);
  range->SetLineColor(kPink+8);
  range->SetMarkerColor(kPink+8);
  range->SetLineStyle(1);
  range->SetMarkerStyle(0);

  TFile * fFuncSyst = TFile::Open(Form("%s/systematicsEM_Fcn_%i.root",fPath.Data(),icent));
  TH1F * dummyF = (TH1F*) fFuncSyst->Get(Form("hSystVsPtPercentageOfCentral_%i",icent));
  function = (TH1F*) dummyF->Clone("function"); 
  function->SetTitle("Res. bg. fit function");
  function->SetLineWidth(2);
  function->SetLineColor(kBlue+2);
  function->SetMarkerColor(kBlue+2);
  function->SetLineStyle(1);
  function->SetMarkerStyle(0); 

  TFile * fPidSyst = TFile::Open(Form("%s/systematicsEM_PID_%i.root",fPath.Data(),icent));
  TH1F * dummyP = (TH1F*) fPidSyst->Get(Form("hSystVsPtPercentageOfCentral_%i",icent));
  pid = (TH1F*) dummyP->Clone("PID"); 
  pid->SetTitle("PID");
  pid->SetLineWidth(3);
  pid->SetLineColor(kBlack);
  pid->SetMarkerColor(kBlack);
  pid->SetLineStyle(2);
  pid->SetMarkerStyle(0); 

  Double_t tracking_1trk = 0.03;
  Double_t trackCuts_rsn = 0.025;
  for (Int_t ii = 0;ii<npt;ii++){
    Int_t ibin = ii+1;
    //material->SetBinContent(ibin, 4);
    tracking->SetBinContent(ibin, tracking_1trk*2.0*100.0);
    trackcuts->SetBinContent(ibin, 2.5);
    // tofmatching->SetBinContent(ibin, 6);
    // tofmatchingMC->SetBinContent(ibin, 4);
    // pidResponse->SetBinContent(ibin, 4);
    // bgnorm->SetBinContent(ibin, 0.02); //use pt-dep uncert.
    // evtplane->SetBinContent(ibin, 0.3); //use pt-dep uncert.
    // evtmix->SetBinContent(ibin, 0.5); //use pt-dep uncert.  
  }

  //sum all contributions in quadrature
  Double_t syst_sum2 = (tracking_1trk*2.0)*(tracking_1trk*2.0) + trackCuts_rsn*trackCuts_rsn; //tracking^2 + track cuts^2
  Double_t syst_uncorrel_sum2 = (tracking_1trk*tracking_1trk) + (trackCuts_rsn)*(trackCuts_rsn); //tracking_1trk^2 + track cuts_rsn^2, for ratio to K and Pi

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
    
    Double_t contrib_pt2 = 0.0; //default 10% se non c'è
    
    //correlated and uncorrelated uncertainties for particle ratios
    Double_t correlPi2 = 0.0; //uncert. correlated with pi uncert. ^2
    Double_t correlKa2 = 0.0;  //uncert. correlated with K uncert. ^2
    Double_t uncorrelPi2 = 0.0; //uncert. uncorrelated with pi uncert. ^2
    Double_t uncorrelKa2 = 0.0;  //uncert. uncorrelated with K uncert. ^2
  
    //fit range
    if (range_syst>0.0) {
      contrib_pt2+=range_syst*range_syst;
      
      //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=range_syst*range_syst;
      uncorrelKa2+=range_syst*range_syst;
    }
 
    //norm bg
    if (bgnorm_syst>0.0) {
      contrib_pt2+=bgnorm_syst*bgnorm_syst;
 
     //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=bgnorm_syst*bgnorm_syst;
      uncorrelKa2+=bgnorm_syst*bgnorm_syst;
    } 

    //function res bg
    if (func_syst>0.0){
      contrib_pt2+=func_syst*func_syst;      

      //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=func_syst*func_syst;
      uncorrelKa2+=func_syst*func_syst;
    }

    //PID
    if (pid_syst>0.0) {
      contrib_pt2+=pid_syst*pid_syst;
      
      //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=pid_syst*pid_syst;
      uncorrelKa2+=pid_syst*pid_syst;
    }

   //material budget
    if (material_syst>0.0) {
      contrib_pt2+=material_syst*material_syst;

      //removed contribution from correlated uncertainty for particle ratios to pi and K, not for p
      uncorrelPi2+=material_syst_r2pi*material_syst_r2pi;
      uncorrelKa2+=material_syst_r2ka*material_syst_r2ka;
    }

    //hadronic interaction cross section
    if (hadrint_syst>0.0) {
      contrib_pt2+=hadrint_syst*hadrint_syst;

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
      contrib_pt2+=bincount_syst*bincount_syst;

      //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=bincount_syst*bincount_syst;
      uncorrelKa2+=bincount_syst*bincount_syst;
    }
    //event mixing
    if (evmix_syst>0.0) {
      contrib_pt2+=evmix_syst*evmix_syst;
      //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=evmix_syst*evmix_syst;
      uncorrelKa2+=evmix_syst*evmix_syst;
    } 
    //width
    if (width_syst>0.0) {
      contrib_pt2+=width_syst*width_syst;
      //contributes to uncorrelated uncertainty for particle ratios
      uncorrelPi2+=width_syst*width_syst;
      uncorrelKa2+=width_syst*width_syst;
      uncorrelPro2+=width_syst*width_syst;
    } 
#endif


    Double_t totsyst = TMath::Sqrt(syst_sum2+contrib_pt2);
    Double_t totsystUncorrPi = TMath::Sqrt(syst_uncorrel_sum2+uncorrelPi2);
    Double_t totsystUncorrKa = TMath::Sqrt(syst_uncorrel_sum2+uncorrelKa2);
    Double_t totsystUncorrPro = TMath::Sqrt(syst_sum2+contrib_pt2-tracking_1trk*tracking_1trk); //remove contribution from tracking from 1 trk 

    sum2->SetBinContent(ibin, totsyst*100);
    sum2_uncorrPi->SetBinContent(ibin, totsystUncorrPi*100);
    sum2_uncorrKa->SetBinContent(ibin, totsystUncorrKa*100);
    sum2_uncorrPro->SetBinContent(ibin, totsystUncorrPro*100);
    Printf("bin %i tot. perc. %6.4f", ibin, totsyst*100);//TMath::Sqrt(syst_sum2+contrib_pt2));
  }
  
  //assign syst err to data
  TFile * fdata = TFile::Open(Form("%s/%s",fPathCorr.Data(),corrFile.Data()));
  TH1F * data = (TH1F*) fdata->Get(Form("%s",hCorrYieldName.Data()));  
  TH1F * data_Wsyst = (TH1F*) data->Clone(Form("%s%i_syst",hCorrYieldName.Data(),icent));
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
    data_Wsyst->SetBinContent(ibin,yd);
    data_Wsyst->SetBinError(ibin, yd_syst);
    data_Wsyst_Wstat->SetBinContent(ibin,yd);
    data_Wsyst_Wstat->SetBinError(ibin,TMath::Sqrt(yd_syst*yd_syst+yd_stat*yd_stat));
    Printf("bin %i     yield = %e     syst err = %e    stat err = %e", ibin, yd, yd_syst, yd_stat);

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
  sum2->GetYaxis()->SetRangeUser(0.1, 50);
  sum2->GetXaxis()->SetRangeUser(0.0, 15.0);
  
  sum2->Draw();
  range->Draw("same");
  function->Draw("same");
  pid->Draw("same");
  //  bincount->Draw("same");
  //  evtmix->Draw("same");
  //width->Draw("same");
  bgnorm->Draw("same");
  tracking->Draw("same");
  trackcuts->Draw("same");
  material->Draw("same");
  hadrint->Draw("same");
  // if (ipid) {
  //   tofmatching->Draw("same");
  //   tofmatchingMC->Draw("same");
  //   pidResponse->Draw("same");
  // }
  
  TLegend * autolegry = (TLegend*)gPad->BuildLegend(0.25,0.65,0.88,0.88, Form("V0A multiplicity %s",centLabel.Data()));
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

  cs->SaveAs(Form("summaryAllSystUncert_cent%i_%s.png",icent, date.Data()));
  cs->SaveAs(Form("summaryAllSystUncert_cent%i_%s.C",icent,date.Data()));



  //make-up
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

  //save to out file
  TFile * fout = new TFile(Form("finalWsyst_%s_%i.root",date.Data(),icent),"recreate");
  fout->cd();
  material->Write();
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
  hHadrInt4ratio2Pi->Write();
  hHadrInt4ratio2Ka->Write();
  hMaterial4ratio2Pi->Write();
  hMaterial4ratio2Ka->Write();
  data->Write();
  sum2->Write();
  data_Wsyst->Write();
  data_Wsyst_Wstat->Write();
  data_Wsyst_uncorrPi->Write();
  data_Wsyst_Wstat_uncorrPi->Write();
  data_Wsyst_uncorrKa->Write();
  data_Wsyst_Wstat_uncorrKa->Write();
  data_Wsyst_uncorrPro->Write();
  data_Wsyst_Wstat_uncorrPro->Write();
  cs->Write();
  fout->Close(); 
  
  return;
  
}

