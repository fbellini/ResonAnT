void syst_contributions_pA(TString deteco ="TOF", TString date="16set13")
{
  for (Int_t j=0;j<5;j++){
    syst_contributions_pA(j, deteco,date.Data());
  }
}
void syst_contributions_pA(Int_t j=0, TString deteco ="TOF",TString date="16set13")
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  Float_t tofXmin=1.0, tofXmax=8.0;
  Float_t tpcXmin=0.5, tpcXmax=8.0;

  //set input name
  TString fPath = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/systematics/systUncert",deteco.Data());
  //  TString fPathCorr = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/fit_bestRange_fixedW",deteco.Data());
  TString fPathCorr = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s_ANA/ana2s/systematics",deteco.Data());
  TString corrFile = "CORRECTED_best_fit_poly2.root";

  Bool_t isTOF = deteco.Contains("TOF");
  TString hCorrYieldName = Form("h%sCorrected_", deteco.Data());
  
  //PbPb
  // Double_t cent[]={ 0.0, 20.0, 40.0, 60.0, 80.0, 90.0};   
  // Double_t pt[] = { 0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.00 };
  
  //pA analysis
  Double_t cent[]={ 0.0, 20.0, 40.0, 60.0, 80.0, 100.0};   
  Double_t pt[] = { 0.0, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };

  Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;   
  Int_t   ncent  = sizeof(cent) / sizeof(cent[0]) - 1;
  TString centLabel=Form("%i-%i%%", 20*j,20*(j+1));
  
  //cosmetics  
  Color_t color[2][6]={kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack, //tpc
		       kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof
  Int_t marker[2][6]={21, 22, 23, 34, 33, 20, //tpc
		      25, 26, 32, 28, 27, 24}; //tof 
  
  //create axis to reproduce the binning
  TAxis *ptbins = new TAxis(npt, pt);
  TAxis *centbins = new TAxis(ncent, cent);
  
  //********************************************************************/
  /********************************************************************/
 
  //pt independent - material budget
  TH1F * material = new TH1F("material","material budget",npt, pt);
  material->SetLineWidth(4);
  material->SetLineColor(kPink+1);
  material->SetMarkerColor(kPink-1);
  material->SetLineStyle(3);
  material->SetMarkerStyle(0); 

  //hadronic interactions
  TH1F * hadronicint = new TH1F("hadronicint","hadronic int.",npt, pt);
  hadronicint->SetLineWidth(4);
  hadronicint->SetLineColor(kYellow+1);
  hadronicint->SetMarkerColor(kYellow-1);
  hadronicint->SetLineStyle(2);
  hadronicint->SetMarkerStyle(0); 
  
  //global tracking
  TH1F * tracking = new TH1F("tracking","global tracking",npt, pt);
  tracking->SetLineWidth(3);
  tracking->SetLineColor(kCyan+1);
  tracking->SetMarkerColor(kCyan+1);
  tracking->SetLineStyle(4);
  tracking->SetMarkerStyle(0); 

  //tof matching
  TH1F * tofmatching = new TH1F("tofmatching","TOF matching",npt, pt);  
  tofmatching->SetLineWidth(3);
  tofmatching->SetLineColor(kGreen+2);
  tofmatching->SetMarkerColor(kGreen+2);
  tofmatching->SetLineStyle(5);
  tofmatching->SetMarkerStyle(0); 

  //tof matching data vs mc
  TH1F * tofmatchingMC = new TH1F("tofmatchingMC","TOF matching data/MC",npt, pt);  
  tofmatchingMC->SetLineWidth(2);
  tofmatchingMC->SetLineColor(kSpring+5);
  tofmatchingMC->SetMarkerColor(kSpring+5);
  tofmatchingMC->SetLineStyle(0);
  tofmatchingMC->SetMarkerStyle(0); 

  //tof pid response
  TH1F * pidResponse = new TH1F(Form("pidResponse_%i",j),Form("TOF PID response tune"), npt, pt);
  pidResponse->SetLineWidth(3);
  pidResponse->SetLineColor(kOrange+6);
  pidResponse->SetMarkerColor(kOrange+6);
  pidResponse->SetLineStyle(2);
  pidResponse->SetMarkerStyle(0); 

  //pt dependent
  TH1F * bincount = new TH1F(Form("bincount_%i",j),Form("Bin counting"), npt, pt);
  TH1F * evtmix = new TH1F(Form("evtmix_%i",j),Form("Event mixing"), npt, pt);
  TH1F * bgnorm = new TH1F(Form("bgnorm_%i",j),Form("LS bg. normalisation"), npt, pt);
  TH1F * range = new TH1F(Form("range_%i",j),Form("Fit range",j), npt, pt);
  TH1F * function = new TH1F(Form("function_%i",j),Form("Res. bg. fit function",j), npt, pt);
  TH1F * width = new TH1F(Form("width_%i",j),Form("Breit Wigner width & range",j), npt, pt);
  TH1F * pid = new TH1F(Form("pid_%i",j),Form("PID"), npt, pt);

  //backround normalisation
  TFile * fBgNorm = TFile::Open(Form("%s/syst_yields_allCents_LSnorm.root",fPath.Data()));
  TH1F * dummyBN = (TH1F*) fBgNorm->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  bgnorm = (TH1F*) dummyBN->Clone("bgNorm"); 
  bgnorm->SetTitle("LS bg. norm. range");
  bgnorm->SetLineWidth(2);
  bgnorm->SetLineColor(kGreen);
  bgnorm->SetMarkerColor(kGreen);
  bgnorm->SetLineStyle(1);
  bgnorm->SetMarkerStyle(0);
  
  //range of the fit and width
  TFile * fRangeSyst = TFile::Open(Form("%s/syst_yields_range_width_chiMax5.0.root",fPath.Data()));
  TH1F * dummyR = (TH1F*) fRangeSyst->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  range = (TH1F*) dummyR->Clone("range"); 
  range->SetTitle("Fit range and free width");
  range->SetLineWidth(2);
  range->SetLineColor(kPink+8);
  range->SetMarkerColor(kPink+8);
  range->SetLineStyle(1);
  range->SetMarkerStyle(0);

  //res.bg function
  TFile * fFuncSyst = TFile::Open(Form("%s/syst_yields_allCents_ResBg.root",fPath.Data()));
  TH1F * dummyF = (TH1F*) fFuncSyst->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  function = (TH1F*) dummyF->Clone("function"); 
  function->SetTitle("Res. bg. fit function");
  function->SetLineWidth(2);
  function->SetLineColor(kBlack);
  function->SetMarkerColor(kBlack);
  function->SetLineStyle(1);
  function->SetMarkerStyle(0); 

  //pid selection
  TFile * fPidSyst = TFile::Open(Form("%s/syst_yields_allCents_pid25s.root",fPath.Data()));
  TH1F * dummyP = (TH1F*) fPidSyst->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  pid = (TH1F*) dummyP->Clone("function"); 
  pid->SetTitle("PID selection (2#sigma vs 2.5#sigma)");
  pid->SetLineWidth(3);
  pid->SetLineColor(kBlue+1);
  pid->SetMarkerColor(kBlue+1);
  pid->SetLineStyle(2);
  pid->SetMarkerStyle(0); 

 //********************************************************************/
  /********************************************************************/
 
  //sum in quadrature - total uncert
  TH1F * sum2 = new TH1F(Form("sum2_%i",j),Form("Sum pt-independent contr.",j), npt, pt);
  sum2->SetLineWidth(3);
  sum2->SetLineColor(kRed);
  sum2->SetMarkerColor(kRed);
  sum2->SetLineStyle(1);
  sum2->SetMarkerStyle(0); 
  
  //sum in quadrature - total uncorr. uncert.
  TH1F * sum2_uncorr = new TH1F(Form("sum2_uncorr%i",j),Form("Uncorrelated tot. contrib.",j), npt, pt);
  sum2_uncorr->SetLineWidth(3);
  sum2_uncorr->SetLineColor(kBlue);
  sum2_uncorr->SetMarkerColor(kBlue);
  sum2_uncorr->SetLineStyle(1);
  sum2_uncorr->SetMarkerStyle(0); 

  //sum in quadrature - total uncorr. uncert.
  TH1F * sum2_corr = new TH1F(Form("sum2_corr%i",j),Form("Correlated tot. contrib.",j), npt, pt);
  sum2_corr->SetLineWidth(3);
  sum2_corr->SetLineColor(kGray+2);
  sum2_corr->SetMarkerColor(kGray+2);
  sum2_corr->SetLineStyle(1);
  sum2_corr->SetMarkerStyle(0);

 //********************************************************************/
  /********************************************************************/

  for (Int_t ii = 2;ii<npt;ii++){
    Int_t ibin = ii+1;
    material->SetBinContent(ibin, 7.5);
    tracking->SetBinContent(ibin, 8);
    hadronicint->SetBinContent(ibin, 5);
    tofmatching->SetBinContent(ibin, 6);
    tofmatchingMC->SetBinContent(ibin, 4);
    pidResponse->SetBinContent(ibin, 4);
  }
  
  //sum all contributions in quadrature
  Double_t syst_corr2 = 0.0;
  Double_t syst_uncorr2 = 0.0;
  
  for (Int_t ii = 1;ii<npt;ii++){
    //reset values
    syst_corr2 = 0.0;
    syst_uncorr2 = 0.0;
    Int_t ibin = ii+1;
    //pt-dependent
    Double_t material_syst    = material->GetBinContent(ibin)/100.;
    Double_t tracking_syst    = tracking->GetBinContent(ibin)/100.;
    Double_t hadronicint_syst = hadronicint->GetBinContent(ibin)/100.;
    Double_t tofmatch_syst    = tofmatching->GetBinContent(ibin)/100.;
    Double_t tofmatchMC_syst  = tofmatchingMC->GetBinContent(ibin)/100.;
    Double_t pidresp_syst     = pidResponse->GetBinContent(ibin)/100.;
    
    //pt-independent
    //    Double_t width_syst = dummyW->GetBinContent(ibin)/100.;
    Double_t range_syst = dummyR->GetBinContent(ibin)/100.;
    Double_t func_syst = dummyF->GetBinContent(ibin)/100.;
    Double_t pid_syst = dummyP->GetBinContent(ibin)/100.;
    
    //set to 0.0 those for which no syst uncert
    Double_t bincount_syst = 0.0;// dummyBC->GetBinContent(ibin)/100.;
    Double_t evmix_syst = 0.0;//dummyEM->GetBinContent(ibin)/100.;
    Double_t bgnorm_syst = 0.0;//dummyBN->GetBinContent(ibin)/100.;
    
    syst_corr2+= material_syst*material_syst;
    syst_corr2+= tracking_syst*tracking_syst;
    syst_corr2+= hadronicint_syst* hadronicint_syst;

    if (isTOF) {
      syst_uncorr2 += tofmatch_syst*tofmatch_syst;
      syst_uncorr2 += tofmatchMC_syst*tofmatchMC_syst;
      syst_uncorr2 += pidresp_syst*pidresp_syst;
    }
    
    //range and width
    if (range_syst>0.0) {
      syst_uncorr2+=range_syst*range_syst;
    }    
    //bin counting
    if (bincount_syst>0.0) {
      syst_uncorr2+=bincount_syst*bincount_syst;
    }
    //event mixing
    if (evmix_syst>0.0) {
       syst_uncorr2+=evmix_syst*evmix_syst;
     } 
    //function res bg
    if (func_syst>0.0){
      syst_uncorr2+=func_syst*func_syst;      
    }
    //norm bg
    if (bgnorm_syst>0.0) {
      syst_uncorr2+=bgnorm_syst*bgnorm_syst;
    } 
    //PID
    if (pid_syst>0.0) {
      syst_uncorr2+=pid_syst*pid_syst;
    }
    
    Double_t totsyst = TMath::Sqrt(syst_uncorr2+syst_corr2);
    Double_t totcorr = TMath::Sqrt(syst_corr2);
    Double_t totuncorr = TMath::Sqrt(syst_uncorr2);
    sum2->SetBinContent(ibin, totsyst*100);
    sum2_uncorr->SetBinContent(ibin, totuncorr*100);
    sum2_corr->SetBinContent(ibin, totcorr*100);
    Printf("bin %i tot. perc. %6.4f", ibin, totsyst*100);//TMath::Sqrt(syst_sum2+contrib_pt2));
  }
  
  //assign syst err to data
  TFile * fdata = TFile::Open(Form("%s/%s",fPathCorr.Data(),corrFile.Data()));
  TH1F * data = (TH1F*) fdata->Get(Form("%s%i",hCorrYieldName.Data(),j));
  
  TH1F * data_Wsyst = (TH1F*) data->Clone(Form("%s%i_syst",hCorrYieldName.Data(),j));
  for (Int_t ii = 1;ii<npt;ii++){
    Int_t ibin = ii;
    Double_t yd = data->GetBinContent(ibin);
    Double_t perc = sum2->GetBinContent(ibin);
    Double_t errsystem = yd*perc/100.;
    data_Wsyst->SetBinContent(ibin,yd);
    data_Wsyst->SetBinError(ibin,errsystem);
    Printf("bin %i     yield = %e     syst err = %e", ibin,yd, errsystem);
  }
  
  TH1F * data_Wuncorr = (TH1F*) data->Clone(Form("%s%i_uncorr",hCorrYieldName.Data(),j));
  for (Int_t ii = 1;ii<npt;ii++){
    Int_t ibin = ii;
    Double_t yd = data->GetBinContent(ibin);
    Double_t perc = sum2_uncorr->GetBinContent(ibin);
    Double_t errsystem = yd*perc/100.;
    data_Wuncorr->SetBinContent(ibin,yd);
    data_Wuncorr->SetBinError(ibin,errsystem);
    Printf("bin %i     yield = %e     uncorr. syst. err. = %e", ibin,yd, errsystem);
  }
 
  TH1F * data_Wcorr = (TH1F*) data->Clone(Form("%s%i_corr",hCorrYieldName.Data(),j));
  for (Int_t ii = 1;ii<npt;ii++){
    Int_t ibin = ii;
    Double_t yd = data->GetBinContent(ibin);
    Double_t perc = sum2_corr->GetBinContent(ibin);
    Double_t errsystem = yd*perc/100.;
    data_Wcorr->SetBinContent(ibin,yd);
    data_Wcorr->SetBinError(ibin,errsystem);
    Printf("bin %i     yield = %e     corr. syst. err. = %e", ibin,yd, errsystem);
  }
 
 //total systematic uncertainty plot
  TCanvas *cs=new TCanvas("cs","Systematic uncertainty vs p_{t}", 750,600);
  cs->cd();
  sum2->SetTitle("Total syst. uncert.");
  sum2->GetYaxis()->SetTitle("relative uncert. (%)");
  sum2->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  sum2->GetYaxis()->SetRangeUser(0.1, 100);
  if (isTOF) sum2->GetXaxis()->SetRangeUser(tofXmin, tofXmax);
  else sum2->GetXaxis()->SetRangeUser(tpcXmin, tpcXmax);
  
  sum2->Draw();
  //bincount->Draw("same");
  //evtmix->Draw("same");
  //bgnorm->Draw("same");
  //width->Draw("same");
  tracking->Draw("same");
  material->Draw("same");
  hadronicint->Draw("same"); 
  function->Draw("same");
  range->Draw("same");
  pid->Draw("same");
  if (isTOF) {
    tofmatching->Draw("same");
    tofmatchingMC->Draw("same");
    pidResponse->Draw("same");
  }
  
  TLegend * autolegry = (TLegend*)gPad->BuildLegend(0.25,0.65,0.88,0.88, Form("Centrality %s",centLabel.Data()));
  autolegry->SetFillColor(kWhite);
  autolegry->SetLineColor(kWhite);
  autolegry->SetTextFont(42);
  autolegry->SetNColumns(2); 
  //color data
  data->SetMarkerColor(color[isTOF][j]+1);
  data->SetLineColor(color[isTOF][j]+1);
  data->SetMarkerStyle(marker[isTOF][j]);
  data->SetMarkerSize(0.7);
  data->SetLineWidth(2);
  data->SetTitle(Form("%s (stat. uncert.)",centLabel.Data()));
  data_Wsyst->SetMarkerColor(color[isTOF][j]);	       
  data_Wsyst->SetFillStyle(0); //Color(color2[j]);
  data_Wsyst->SetLineWidth(1); //Color(color2[j]);
  data_Wsyst->SetLineColor(color[isTOF][j]);
  data_Wsyst->SetMarkerStyle(0);
  data_Wsyst->SetOption("E2");
  data_Wsyst->SetTitle(Form("%s (syst. uncert.)",centLabel.Data()));

  cs->SaveAs(Form("summaryAllSystUncert_cent%i_%s.png",j, date.Data()));
  cs->SaveAs(Form("summaryAllSystUncert_cent%i_%s.C",j,date.Data()));
  
//uncorrelated systematic uncertainty plot
  TCanvas *cu=new TCanvas("cu","Uncorrelated systematic uncertainty vs p_{T}", 750,600);
  cu->cd();
  sum2_uncorr->SetTitle("uncorrelated syst. uncert.");
  sum2_uncorr->GetYaxis()->SetTitle("relative uncert. (%)");
  sum2_uncorr->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  sum2_uncorr->GetYaxis()->SetRangeUser(0.1, 100);
  if (isTOF) sum2_uncorr->GetXaxis()->SetRangeUser(tofXmin, tofXmax);
  else sum2_uncorr->GetXaxis()->SetRangeUser(tpcXmin, tpcXmax);
  
  sum2_uncorr->Draw();
  //bincount->Draw("same");
  //evtmix->Draw("same");
  //bgnorm->Draw("same");
  //width->Draw("same");
  function->Draw("same");
  range->Draw("same");
  pid->Draw("same");
  if (isTOF) {
    tofmatching->Draw("same");
    tofmatchingMC->Draw("same");
    pidResponse->Draw("same");
  }
  
  TLegend * autolegry_uncorr = (TLegend*)gPad->BuildLegend(0.25,0.65,0.88,0.88, Form("Centrality %s",centLabel.Data()));
  autolegry_uncorr->SetFillColor(kWhite);
  autolegry_uncorr->SetLineColor(kWhite);
  autolegry_uncorr->SetTextFont(42);
  autolegry_uncorr->SetNColumns(2); 
  //color data
  data->SetMarkerColor(color[isTOF][j]+1);
  data->SetLineColor(color[isTOF][j]+1);
  data->SetMarkerStyle(marker[isTOF][j]);
  data->SetMarkerSize(0.7);
  data->SetLineWidth(2);
  data->SetTitle(Form("%s (stat. uncert.)",centLabel.Data()));
  data_Wuncorr->SetMarkerColor(color[isTOF][j]);	       
  data_Wuncorr->SetFillStyle(0); //Color(color2[j]);
  data_Wuncorr->SetLineWidth(1); //Color(color2[j]);
  data_Wuncorr->SetLineColor(color[isTOF][j]);
  data_Wuncorr->SetMarkerStyle(0);
  data_Wuncorr->SetOption("E2");
  data_Wuncorr->SetTitle(Form("%s (syst. uncert.)",centLabel.Data()));

  cu->SaveAs(Form("summaryUncorrSystUncert_cent%i_%s.png", j, date.Data()));
  cu->SaveAs(Form("summaryUncorrSystUncert_cent%i_%s.C", j, date.Data()));

//correlated systematic uncertainty plot
  TCanvas *cc=new TCanvas("cc","Correlated systematic uncertainty vs p_{T}", 750,600);
  cc->cd();
  sum2_corr->SetTitle("correlated syst. uncert.");
  sum2_corr->GetYaxis()->SetTitle("relative uncert. (%)");
  sum2_corr->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  sum2_corr->GetYaxis()->SetRangeUser(0.1, 100);
  if (isTOF) sum2_corr->GetXaxis()->SetRangeUser(tofXmin, tofXmax);
  else sum2_corr->GetXaxis()->SetRangeUser(tpcXmin, tpcXmax);
  
  sum2_corr->Draw();
  tracking->Draw("same");
  material->Draw("same");
  hadronicint->Draw("same"); 
  
  TLegend * autolegry_corr = (TLegend*)gPad->BuildLegend(0.25,0.65,0.88,0.88, Form("Centrality %s",centLabel.Data()));
  autolegry_corr->SetFillColor(kWhite);
  autolegry_corr->SetLineColor(kWhite);
  autolegry_corr->SetTextFont(42);
  autolegry_corr->SetNColumns(2); 
  //color data
  data->SetMarkerColor(color[isTOF][j]+1);
  data->SetLineColor(color[isTOF][j]+1);
  data->SetMarkerStyle(marker[isTOF][j]);
  data->SetMarkerSize(0.7);
  data->SetLineWidth(2);
  data->SetTitle(Form("%s (stat. uncert.)",centLabel.Data()));
  data_Wcorr->SetMarkerColor(color[isTOF][j]);	       
  data_Wcorr->SetFillStyle(0); //Color(color2[j]);
  data_Wcorr->SetLineWidth(1); //Color(color2[j]);
  data_Wcorr->SetLineColor(color[isTOF][j]);
  data_Wcorr->SetMarkerStyle(0);
  data_Wcorr->SetOption("E2");
  data_Wcorr->SetTitle(Form("%s (syst. uncert.)",centLabel.Data()));

  cc->SaveAs(Form("summaryCorrSystUncert_cent%i_%s.png", j, date.Data()));
  cc->SaveAs(Form("summaryCorrSystUncert_cent%i_%s.C", j, date.Data()));



  //save to out file
  TFile * fout = new TFile(Form("finalWsyst_%s_%i.root",date.Data(),j),"recreate");
  fout->cd();
  material->Write();
  tracking->Write();
  hadronicint->Write();
  if (isTOF) {
    tofmatching->Write();
    tofmatchingMC->Write();
    pidResponse->Write();
  }
  pid->Write();
  //bincount->Write();
  //evtmix->Write();
  range->Write();
  bgnorm->Write();
  function->Write();
  data->Write();
  sum2->Write();
  sum2_uncorr->Write();
  sum2_corr->Write();
  data_Wsyst->Write();
  data_Wuncorr->Write();
  data_Wcorr->Write();
  cs->Write();
  cu->Write();
  cc->Write();
  fout->Close(); 
  
  return;
  
}

  //pt depepndent    
  // TFile * fBinCount = TFile::Open(Form("%s/syst_yields_allCents_binCount_E4s_P2s.root",fPath.Data()));
  // TH1F * dummyBC = (TH1F*) fBinCount->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  // bincount = (TH1F*) dummyBC->Clone("binCount"); 
  // bincount->SetTitle("Bin Counting");
  // bincount->SetLineWidth(3);
  // bincount->SetLineColor(kViolet-5);
  // bincount->SetMarkerColor(kViolet-5);
  // bincount->SetLineStyle(3);
  // bincount->SetMarkerStyle(0);
  
  // TFile * fEvMix = TFile::Open(Form("%s/syst_yields_allCents_em.root",fPath.Data()));
  // TH1F * dummyEM = (TH1F*) fEvMix->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  // evtmix = (TH1F*) dummyEM->Clone("evMix"); 
  // evtmix->SetTitle("Event mixing");
  // evtmix->SetLineWidth(3);
  // evtmix->SetLineColor(kBlue);
  // evtmix->SetMarkerColor(kBlue);
  // evtmix->SetLineStyle(2);
  // evtmix->SetMarkerStyle(0);

  // TFile * fWidthSyst = TFile::Open(Form("%s/syst_yields_allCents_freeW08-15.root",fPath.Data()));
  // TH1F * dummyW = (TH1F*) fWidthSyst->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  // width = (TH1F*) dummyW->Clone("width"); 
  // width->SetTitle("Width constraints");
  // width->SetLineWidth(3);
  // width->SetLineColor(kCyan-2);
  // width->SetMarkerColor(kCyan-2);
  // width->SetLineStyle(2);
  // width->SetMarkerStyle(0); 
