void syst_contributions_pA100(Int_t j=0, TString deteco ="TOF",TString date="16set13")
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  //set input name
  TString fPath = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s100/ana2s/systematics/systUncert",deteco.Data());
  TString fPathCorr = Form("/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/%s100/ana2s", deteco.Data());
  TString corrFile = "CORRECTED_best_fit_poly2.root";

  Bool_t isTOF = deteco.Contains("TOF");
  TString hCorrYieldName = Form("correzione%s",deteco.Data());
  
  //PbPb
  // Double_t cent[]={ 0.0, 20.0, 40.0, 60.0, 80.0, 90.0};   
  // Double_t pt[] = { 0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.00 };
  
  //pA analysis
  Double_t cent[]={ 0.0, 100.0};   
  Double_t pt[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 10.00 };
  
  Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;   
  Int_t   ncent  = sizeof(cent) / sizeof(cent[0]) - 1;
  TString centLabel="0-100%";
  
  //cosmetics  
  Color_t color[2][6]={kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack, //tpc
		       kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof
  Int_t marker[2][6]={21, 22, 23, 34, 33, 20, //tpc
		      25, 26, 32, 28, 27, 24}; //tof 
  
  //create axis to reproduce the binning
  TAxis *ptbins = new TAxis(npt, pt);
  TAxis *centbins = new TAxis(ncent, cent);
  
  //pt dependent
  TH1F * bincount = new TH1F(Form("bincount_%i",j),Form("Bin counting"), npt, pt);
  TH1F * evtmix = new TH1F(Form("evtmix_%i",j),Form("Event mixing"), npt, pt);
  TH1F * bgnorm = new TH1F(Form("bgnorm_%i",j),Form("LS bg. normalisation"), npt, pt);
  TH1F * range = new TH1F(Form("range_%i",j),Form("Fit range",j), npt, pt);
  TH1F * function = new TH1F(Form("function_%i",j),Form("Res. bg. fit function",j), npt, pt);
  TH1F * width = new TH1F(Form("width_%i",j),Form("Breit Wigner width & range",j), npt, pt);
  TH1F * pid = new TH1F(Form("pid_%i",j),Form("PID"), npt, pt);

  //pt independent
  TH1F * material = new TH1F("material","material budget",npt, pt);
  material->SetLineWidth(4);
  material->SetLineColor(kPink+1);
  material->SetMarkerColor(kPink-1);
  material->SetLineStyle(3);
  material->SetMarkerStyle(0); 
  TH1F * tracking = new TH1F("tracking","global tracking",npt, pt);
  tracking->SetLineWidth(3);
  tracking->SetLineColor(kCyan+1);
  tracking->SetMarkerColor(kCyan+1);
  tracking->SetLineStyle(4);
  tracking->SetMarkerStyle(0); 
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
  TH1F * pidResponse = new TH1F(Form("pidResponse_%i",j),Form("TOF PID response tune"), npt, pt);
  pidResponse->SetLineWidth(3);
  pidResponse->SetLineColor(kOrange+6);
  pidResponse->SetMarkerColor(kOrange+6);
  pidResponse->SetLineStyle(2);
  pidResponse->SetMarkerStyle(0); 

  TH1F * sum2 = new TH1F(Form("sum2_%i",j),Form("Sum pt-independent contr.",j), npt, pt);
  sum2->SetLineWidth(3);
  sum2->SetLineColor(kRed);
  sum2->SetMarkerColor(kRed);
  sum2->SetLineStyle(1);
  sum2->SetMarkerStyle(0); 

  //pt dependent    
  TFile * fBinCount = TFile::Open(Form("%s/Syst_fromRatio_SystBinCounting_BCtoFIT.root",fPath.Data()));
  TH1F * dummyBC = (TH1F*) fBinCount->Get(Form("hSystVsPtPercentage"));
  bincount = (TH1F*) dummyBC->Clone("binCount"); 
  bincount->SetTitle("Bin Counting");
  bincount->SetLineWidth(3);
  bincount->SetLineColor(kViolet-5);
  bincount->SetMarkerColor(kViolet-5);
  bincount->SetLineStyle(3);
  bincount->SetMarkerStyle(0);
  
  TFile * fEvMix = TFile::Open(Form("%s/Syst_fromRatio_SystBG_EMtoLS.root",fPath.Data()));
  TH1F * dummyEM = (TH1F*) fEvMix->Get(Form("hSystVsPtPercentage"));
  evtmix = (TH1F*) dummyEM->Clone("evMix"); 
  evtmix->SetTitle("Event mixing");
  evtmix->SetLineWidth(3);
  evtmix->SetLineColor(kBlue);
  evtmix->SetMarkerColor(kBlue);
  evtmix->SetLineStyle(2);
  evtmix->SetMarkerStyle(0);

#if 0 //not a systematic effect
  TFile * fBgNorm = TFile::Open(Form("%s/syst_yields_allCents_LSnorm.root",fPath.Data()));
  TH1F * dummyBN = (TH1F*) fBgNorm->Get(Form("hSystVsPtPercentage"));
  bgnorm = (TH1F*) dummyBN->Clone("bgNorm"); 
  bgnorm->SetTitle("LS bg. norm. range");
  bgnorm->SetLineWidth(2);
  bgnorm->SetLineColor(kGreen);
  bgnorm->SetMarkerColor(kGreen);
  bgnorm->SetLineStyle(1);
  bgnorm->SetMarkerStyle(0);
#endif
  
  TFile * fRangeSyst = TFile::Open(Form("%s/syst_yields_SystRMS_WidthRange.root",fPath.Data()));
  TH1F * dummyR = (TH1F*) fRangeSyst->Get(Form("hSystVsPtPercentageOfCentral_0"));
  range = (TH1F*) dummyR->Clone("range"); 
  range->SetTitle("Fit range and free width");
  range->SetLineWidth(2);
  range->SetLineColor(kPink+8);
  range->SetMarkerColor(kPink+8);
  range->SetLineStyle(1);
  range->SetMarkerStyle(0);

  TFile * fFuncSyst = TFile::Open(Form("%s/syst_yields_SystRMSFunctionResBg.root",fPath.Data()));
  TH1F * dummyF = (TH1F*) fFuncSyst->Get(Form("hSystVsPtPercentageOfCentral_0"));
  function = (TH1F*) dummyF->Clone("function"); 
  function->SetTitle("Res. bg. fit function");
  function->SetLineWidth(2);
  function->SetLineColor(kBlack);
  function->SetMarkerColor(kBlack);
  function->SetLineStyle(1);
  function->SetMarkerStyle(0); 

  TFile * fPidSyst = TFile::Open(Form("%s/Syst_fromRatio_SystPID_2sto25s.root",fPath.Data()));
  TH1F * dummyP = (TH1F*) fPidSyst->Get(Form("hSystVsPtPercentage"));
  pid = (TH1F*) dummyP->Clone("function"); 
  pid->SetTitle("PID selection (2#sigma vs 2.5#sigma)");
  pid->SetLineWidth(3);
  pid->SetLineColor(kYellow+1);
  pid->SetMarkerColor(kYellow+1);
  pid->SetLineStyle(2);
  pid->SetMarkerStyle(0); 

  TFile * fWidthSyst = TFile::Open(Form("%s/syst_yields_SystWidth.root",fPath.Data()));
  TH1F * dummyW = (TH1F*) fWidthSyst->Get(Form("hSystVsPtPercentageOfCentral_0"));
  width = (TH1F*) dummyW->Clone("width"); 
  width->SetTitle("Width constraints");
  width->SetLineWidth(3);
  width->SetLineColor(kCyan-2);
  width->SetMarkerColor(kCyan-2);
  width->SetLineStyle(2);
  width->SetMarkerStyle(0); 

  for (Int_t ii = 0;ii<npt;ii++){
    Int_t ibin = ii+1;
    material->SetBinContent(ibin, 6);
    tracking->SetBinContent(ibin, 8);
    tofmatching->SetBinContent(ibin, 6);
    tofmatchingMC->SetBinContent(ibin, 4);
    pidResponse->SetBinContent(ibin, 4);
    // bgnorm->SetBinContent(ibin, 0.02); //use pt-dep uncert.
    // evtplane->SetBinContent(ibin, 0.3); //use pt-dep uncert.
    // evtmix->SetBinContent(ibin, 0.5); //use pt-dep uncert.  
  }

  //sum all contributions in quadrature
  Double_t syst_sum2 = 0.06*0.06 + 0.08*0.08;
  if (isTOF) syst_sum2 +=0.06*0.06 + 0.04*0.04 + 0.04*0.04;

  for (Int_t ii = 0;ii<npt;ii++){
    Int_t ibin = ii+1;
    Double_t bincount_syst = dummyBC->GetBinContent(ibin)/100.;
    Double_t evmix_syst = dummyEM->GetBinContent(ibin)/100.;
    Double_t range_syst = dummyR->GetBinContent(ibin)/100.;
    Double_t func_syst = dummyF->GetBinContent(ibin)/100.;
    Double_t bgnorm_syst;// = dummyBN->GetBinContent(ibin)/100.;
    Double_t width_syst = dummyW->GetBinContent(ibin)/100.;
    Double_t pid_syst = dummyP->GetBinContent(ibin)/100.;
    
    Double_t contrib_pt2 = 0.0; //default 10% se non c'è
   
    //skip those for which no syst uncert
    bincount_syst=0.0;
    //evmix_syst=0.0;
    bgnorm_syst=0.0;
    bgwidth_syst=0.0;

    if (range_syst>0.0) {
      contrib_pt2+=range_syst*range_syst;
    }
 
    //bin counting
    if (bincount_syst>0.0) {
      contrib_pt2+=bincount_syst*bincount_syst;
    }
    //event mixing
    if (evmix_syst>0.0) {
      contrib_pt2+=evmix_syst*evmix_syst;
    } 
    //function res bg
    if (func_syst>0.0){
      contrib_pt2+=func_syst*func_syst;      
    }
    // //norm bg
    if (bgnorm_syst>0.0) {
      contrib_pt2+=bgnorm_syst*bgnorm_syst;
    } 
    //PID
    if (pid_syst>0.0) {
      contrib_pt2+=pid_syst*pid_syst;
    }
    //width
    if (width_syst>0.0) {
      contrib_pt2+=width_syst*width_syst;
    }
    
    Double_t totsyst = TMath::Sqrt(syst_sum2+contrib_pt2);
    sum2->SetBinContent(ibin, totsyst*100);
    Printf("bin %i tot. perc. %6.4f", ibin, totsyst*100);//TMath::Sqrt(syst_sum2+contrib_pt2));
  }
  
  //assign syst err to data
  TFile * fdata = TFile::Open(Form("%s/%s",fPathCorr.Data(),corrFile.Data()));
  TH1F * data = (TH1F*) fdata->Get(Form("%s",hCorrYieldName.Data()));
  
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
  
  //systematic uncertainty plot
  TCanvas *cs=new TCanvas("cs","Systematic uncertainty vs p_{t}", 750,600);
  cs->cd();
  sum2->SetTitle("Total syst. uncert.");
  sum2->GetYaxis()->SetTitle("relative uncert. (%)");
  sum2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  sum2->GetYaxis()->SetRangeUser(0.1, 100);
  if (isTOF) sum2->GetXaxis()->SetRangeUser(1., 8.);
  else sum2->GetXaxis()->SetRangeUser(0.5, 8.);
  
  sum2->Draw();
  range->Draw("same");
  function->Draw("same");
  pid->Draw("same");
  //bincount->Draw("same");
  evtmix->Draw("same");
  //bgnorm->Draw("same");
  //width->Draw("same");
  tracking->Draw("same");
  material->Draw("same");
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
  
  //save to out file
  TFile * fout = new TFile(Form("finalWsyst_%s_%i.root",date.Data(),j),"recreate");
  fout->cd();
  material->Write();
  tracking->Write();
  if (isTOF) {
    tofmatching->Write();
    tofmatchingMC->Write();
    pidResponse->Write();
  }
  pid->Write();
  bincount->Write();
  evtmix->Write();
  range->Write();
  bgnorm->Write();
  function->Write();
  data->Write();
  sum2->Write();
  data_Wsyst->Write();
  cs->Write();
  fout->Close(); 

  return;
  
}

