void syst_contributionMass(Int_t j=0)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  Double_t cent[]={ 0.0, 20.0, 40.0, 60.0, 80.0, 90.0};   
  Double_t pt[] = { 0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.00 };
  Color_t color[]={kRed+1, kOrange, kGreen+2, kBlue+1};
  Color_t color2[]={kRed-10, kYellow-9, kTeal-9, kAzure-9};
  Int_t marker[]={20, 21, 28, 26};

  Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;   
  Int_t   ncent  = sizeof(cent) / sizeof(cent[0]) - 1;
  TString centLabel=Form("%2.0f-%2.0f%%",cent[j],cent[j+1]);
  
  //create axis to reproduce the binning
  TAxis *ptbins = new TAxis(npt, pt);
  TAxis *centbins = new TAxis(ncent, cent);

  TH1F * evtmix = new TH1F(Form("evtmix_%i",j),Form("N. mixed events"), npt, pt);
  TH1F * bgnorm = new TH1F(Form("bgnorm_%i",j),Form("EM bg. normalisation"), npt, pt);
  TH1F * evtplane = new TH1F(Form("evtplane_%i",j),Form("Reaction plane"), npt, pt);

  TH1F * range = new TH1F(Form("range_%i",j),Form("Fit range",j), npt, pt);
  TH1F * function = new TH1F(Form("function_%i",j),Form("Res. bg. fit function",j), npt, pt);
  TH1F * peak = new TH1F(Form("peak_%i",j),Form("Peak fit function",j), npt, pt);

  TH1F * material = new TH1F("material","material budget",npt, pt);
  material->SetLineWidth(3);
  material->SetLineColor(kAzure-5);
  material->SetMarkerColor(kAzure-5);
  material->SetLineStyle(3);
  material->SetMarkerStyle(0); 

  TH1F * tracking = new TH1F("tracking","tracking and track sel.",npt, pt);
  tracking->SetLineWidth(3);
  tracking->SetLineColor(kCyan+1);
  tracking->SetMarkerColor(kCyan+1);
  tracking->SetLineStyle(4);
  tracking->SetMarkerStyle(0); 

  TH1F * tofmatching = new TH1F("tofmatching","TOF matching.",npt, pt);  
  tofmatching->SetLineWidth(3);
  tofmatching->SetLineColor(kGreen+2);
  tofmatching->SetMarkerColor(kGreen+2);
  tofmatching->SetLineStyle(5);
  tofmatching->SetMarkerStyle(0); 

  TH1F * pid = new TH1F(Form("pid_%i",j),Form("PID"), npt, pt);
  pid->SetLineWidth(3);
  pid->SetLineColor(kYellow+1);
  pid->SetMarkerColor(kYellow+1);
  pid->SetLineStyle(2);
  pid->SetMarkerStyle(0); 

  TH1F * sum2 = new TH1F(Form("sum2_%i",j),Form("Sum pt-independent contr.",j), npt, pt);
  sum2->SetLineWidth(2);
  sum2->SetLineColor(kRed);
  sum2->SetMarkerColor(kRed);
  sum2->SetLineStyle(1);
  sum2->SetMarkerStyle(0); 

  for (Int_t ii = 1;ii<npt;ii++){
    Int_t ibin = ii+1;
    material->SetBinContent(ibin, 4);
    tracking->SetBinContent(ibin, 10);
    tofmatching->SetBinContent(ibin, 8);
    pid->SetBinContent(ibin, 0.3);
    bgnorm->SetBinContent(ibin, 0.00); 
    evtplane->SetBinContent(ibin, 0.00); 
    evtmix->SetBinContent(ibin, 0.5); //includes bgnorm & evt plane  
    function->SetBinContent(ibin, 0.2);
  }
  
  
  // TFile * fEvMix = TFile::Open(Form("syst_normYields_cent%i_kstar_3mix_best_poly2.root",j));
  // TH1F * dummyEM = (TH1F*) fEvMix->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  // evtmix = (TH1F*) dummyEM->Clone("evMix"); 
  evtmix->SetTitle("Event mixing");
  evtmix->SetLineWidth(3);
  evtmix->SetLineColor(kBlue);
  evtmix->SetMarkerColor(kBlue);
  evtmix->SetLineStyle(2);
  evtmix->SetMarkerStyle(0);

  /*
    TFile * fBgNorm = TFile::Open(Form("syst_normYields_cent%i_kstar_norm_best_poly2.root",j));
  TH1F * dummyBN = (TH1F*) fBgNorm->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  bgnorm = (TH1F*) dummyBN->Clone("bgNorm"); 
  bgnorm->SetTitle("EM bg. norm. range");
  bgnorm->SetLineWidth(2);
  bgnorm->SetLineColor(kGreen);
  bgnorm->SetMarkerColor(kGreen);
  bgnorm->SetLineStyle(1);
  bgnorm->SetMarkerStyle(0);
 
  TFile * fEvPlane = TFile::Open(Form("syst_normYields_cent%i_kstar_noEP_best_poly2.root",j));
  TH1F * dummyEP = (TH1F*) fEvPlane->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  evtplane = (TH1F*) dummyEP->Clone("evtplane"); 
  evtplane->SetTitle("Evt. plane binning");
  evtplane->SetLineWidth(2);
  evtplane->SetLineColor(kOrange-3);
  evtplane->SetMarkerColor(kOrange-3);
  evtplane->SetLineStyle(9);
  evtplane->SetMarkerStyle(0);
  */
  
  TFile * fRangeSyst = TFile::Open(Form("masssystUncert_cent0%i_poly2_chiMax2.0.root",j));
  TH1F * dummyR = (TH1F*) fRangeSyst->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  range = (TH1F*) dummyR->Clone("range"); 
  range->SetTitle("Fit range");
  range->SetLineWidth(2);
  range->SetLineColor(kMagenta);
  range->SetMarkerColor(kMagenta);
  range->SetLineStyle(5);
  range->SetMarkerStyle(0);

  // TFile * fFuncSyst = TFile::Open(Form("tesimass_cent0%i_chiMax2.0.root",j));
  // TH1F * dummyF = (TH1F*) fFuncSyst->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  // function = (TH1F*) dummyF->Clone("function"); 
  function->SetTitle("Res. bg. fit function");
  function->SetLineWidth(2);
  function->SetLineColor(kBlack);
  function->SetMarkerColor(kBlack);
  function->SetLineStyle(1);
  function->SetMarkerStyle(0); 
  
  TFile * fPeakSyst = TFile::Open(Form("mass_systUncert_peak_cent0%i_chiMax2.0.root",j));
  TH1F * dummyP = (TH1F*) fPeakSyst->Get(Form("hSystVsPtPercentageOfCentral_%i",j));
  peak = (TH1F*) dummyP->Clone("peak"); 
  peak->SetTitle("Peak fit function");
  peak->SetLineWidth(2);
  peak->SetLineColor(kTeal-7);
  peak->SetMarkerColor(kTeal-7);
  peak->SetLineStyle(1);
  peak->SetMarkerStyle(0);
  
  //sum all contributions in quadrature
  Double_t syst_sum2 = 0.003*0.003+0.005*0.005+0.002*0.002;
  for (Int_t ii = 1;ii<npt;ii++){
    Int_t ibin = ii+1;
    Double_t range_syst = dummyR->GetBinContent(ibin)/100.;
    //Double_t func_syst = dummyF->GetBinContent(ibin)/100.;
    Double_t peak_syst = dummyP->GetBinContent(ibin)/100.;
    // Double_t evmix_syst = dummyEM->GetBinContent(ibin)/100.;
    // Double_t bgnorm_syst = dummyBN->GetBinContent(ibin)/100.;
    // Double_t evplane_syst = dummyEP->GetBinContent(ibin)/100.;
    Double_t contrib_pt2 = 0.0; //default 10% se non c'è

    if (range_syst>0.0) {
      contrib_pt2+=(range_syst*range_syst);
    } else {
      contrib_pt2+=(0.003*0.003);
    }

    // if (func_syst>0.0){
    //   contrib_pt2+=(func_syst*func_syst);      
    // } else {
    //   contrib_pt2+=(0.002*0.002);
    // }

    if (peak_syst>0.0){
      contrib_pt2+=(peak_syst*peak_syst);      
    } else {
      contrib_pt2+=(0.003*0.003);
    }
    //comment pt-dependent case
    // if (evmix_syst>0.0) {
    //   contrib_pt2+=evmix_syst*evmix_syst;
    // } else {
    //   contrib_pt2+=0.05*0.05;
    // }

    // if (bgnorm_syst>0.0) {
    //   contrib_pt2+=bgnorm_syst*bgnorm_syst;
    // } else {
    //   contrib_pt2+=0.01*0.01;
    // }

    // if (evplane_syst>0.0) {
    //   contrib_pt2+=evplane_syst*evplane_syst;
    // } else {
    //   if (j==0) contrib_pt2+=0.05*0.05;
    //   else contrib_pt2+=0.03*0.03;
    // }
    
    Double_t totsyst = TMath::Sqrt(syst_sum2+contrib_pt2);
    //if (totsyst>=0.05) 
    sum2->SetBinContent(ibin, totsyst*100);
    Printf("bin %i tot. perc. %6.4f", ibin, totsyst*100);//TMath::Sqrt(syst_sum2+contrib_pt2));
  }
  
  //assign syst err to data
  TFile * fdata = TFile::Open("/Users/bellini/alice/resonances/myKstar/pwglf_train_out/fixedEM_5mix/rawYields/rawYields_tree_best_poly2_aod049_01mar13.root");//Form("final_mass_09mar13.root"));
  TH1F * data = (TH1F*) fdata->Get(Form("hMassVsPt_%i",j));
  
  TH1F * data_Wsyst = (TH1F*) data->Clone(Form("hMassVsPt_%i_syst",j));
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
  sum2->GetYaxis()->SetRangeUser(0.01, 2.);
  sum2->GetXaxis()->SetRangeUser(1., 10.);
  sum2->GetYaxis()->SetTitle("syst. uncert. (% of central value)");
  sum2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  sum2->Draw();
  range->Draw("same");
  function->Draw("same");
  peak->Draw("same");
  evtmix->Draw("same");
  // evtplane->Draw("same");
  // bgnorm->Draw("same");
  // tracking->Draw("same");
  // material->Draw("same");
  // tofmatching->Draw("same");
  pid->Draw("same");
  
  TLegend * autolegry = (TLegend*)gPad->BuildLegend(0.15,0.70,0.6,0.88, Form("Centrality %s",centLabel.Data()));
  autolegry->SetFillColor(kWhite);
  autolegry->SetLineColor(kWhite);
  autolegry->SetTextFont(42);
  autolegry->SetNColumns(2); 
  //color data
  data->SetMarkerColor(color[j]+1);
  data->SetLineColor(color[j]+1);
  data->SetMarkerStyle(marker[j]);
  data->SetMarkerSize(0.7);
  data->SetLineWidth(2);
  data->SetTitle(Form("%s (stat. uncert.)",centLabel.Data()));
  data_Wsyst->SetMarkerColor(color[j]);	       
  data_Wsyst->SetFillStyle(0); //Color(color2[j]);
  data_Wsyst->SetLineWidth(1); //Color(color2[j]);
  data_Wsyst->SetLineColor(color[j]);
  data_Wsyst->SetMarkerStyle(0);
  data_Wsyst->SetOption("E2");
  data_Wsyst->SetTitle(Form("%s (syst. uncert.)",centLabel.Data()));

  cs->SaveAs(Form("summaryAllSystUncert_Mass_cent%i_03mar13.png",j));
  cs->SaveAs(Form("summaryAllSystUncert_Mass_cent%i_03mar13.C",j));

  //save to out file
  TFile * fout = new TFile(Form("finalWsyst_Mass_09mar13_%i.root",j),"recreate");
  fout->cd();
  material->Write();
  tracking->Write();
  tofmatching->Write();
  pid->Write();
  evtmix->Write();
  evtplane->Write();
  range->Write();
  bgnorm->Write();
  function->Write();
  peak->Write();
  data->Write();
  sum2->Write();
  data_Wsyst->Write();
  cs->Write();
  fout->Close(); 

  return;
  
}

