void syst_contributions_pAmulti(Int_t icent=0, 
				TString bg = "MB",
				TString date="10apr18", 
				Int_t barlow = 1, 
				Int_t smooth = -1)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  //set input name
  TString fPath = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiC3_tpc2sPtDep_tof2sveto5smism/syst20180406";
  TString fPathCorr = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiC3_tpc2sPtDep_tof2sveto5smism/norm1.05-1.10/fit_Mixing_VOIGTpoly1_FixRes/fit_r0.994-1.050";
  TString corrFile = "CORRECTED_br_fitResult.root"; //CHANGE-ME
  TString hCorrYieldName = Form("hCorrected_%i",icent); 
  
  //pA analysis
  Double_t cent[]={0.0, 30.0, 60.0, 90.0};   
  Double_t pt[] = {0.0, 0.3, 0.5, 0.7, 0.9, 1.10, 1.30, 1.50, 2.00, 3.00, 4.00, 5.0, 7.0, 10.0};  
  Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;   
  Int_t   ncent  = sizeof(cent) / sizeof(cent[0]) - 1;
  TString centLabel=Form("%3.0f -%3.0f%%",cent[icent], cent[icent+1]);
  
  //cosmetics  
  Color_t color[1][3] = {kRed+1, kSpring+5, kBlue+1};
  Int_t  marker[1][3] = {20, 21, 33}; 
  
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
  TH1F * branching = new TH1F("branching","B.R.",npt, pt);
  branching->SetLineWidth(2);
  branching->SetLineColor(kMagenta);
  branching->SetMarkerColor(kMagenta);
  branching->SetLineStyle(5);
  branching->SetMarkerStyle(0);

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
  TH1F * vertexing = new TH1F("vertexing","Vertexing efficiency",npt, pt);
  vertexing->SetLineWidth(2);
  vertexing->SetLineColor(kAzure+3);
  vertexing->SetMarkerColor(kAzure+3);
  vertexing->SetLineStyle(3);
  vertexing->SetMarkerStyle(0); 
 
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
  
  TFile * fHadrSyst = TFile::Open(Form("%s/systHadrInt.root",fPath.Data()));
  TH1F * dummyHI = (TH1F*) fHadrSyst->Get(Form("hSystVsPtPercentageOfCentral"));
  hadrint = (TH1F*) dummyHI->Clone("hadrint"); 
  hadrint->SetTitle("Hadronic int.");
  hadrint->SetLineWidth(2);
  hadrint->SetLineColor(kOrange+1);
  hadrint->SetMarkerColor(kOrange+1);
  hadrint->SetLineStyle(1);
  hadrint->SetMarkerStyle(0); 
  
  TString normfile = Form("systematicsB1_Background_cent%i.root", icent);
  TString sysHistoName = "hSystVsPtPercentageOfCentral_rms";
  TFile * fBgNorm = TFile::Open(Form("%s/%s",fPath.Data(),normfile.Data()));
  TH1F * dummyBN = (TH1F*) fBgNorm->Get(sysHistoName.Data());
  bgnorm = (TH1F*) dummyBN->Clone("bgNorm");
  bgnorm->Scale(100.);
  bgnorm->SetTitle("MEB normalisation");
  bgnorm->SetLineWidth(2);
  bgnorm->SetLineColor(kPink+8);
  bgnorm->SetMarkerColor(kPink+8);
  bgnorm->SetLineStyle(8);
  bgnorm->SetMarkerStyle(0);
  
  TString rangefile = Form("systematicsB1_Fit_range_cent%i.root", icent);
  sysHistoName = "hSystVsPtPercentageOfCentral_rms";
  TFile * fRangeSyst = TFile::Open(Form("%s/%s",fPath.Data(),rangefile.Data()));
  TH1F * dummyR = 0x0;
  if (fRangeSyst) {
    dummyR = (TH1F*) fRangeSyst->Get(sysHistoName.Data());
    range = (TH1F*) dummyR->Clone("range");
    range->Scale(100.);
    range->SetTitle("Fit range");
    range->SetLineWidth(3);
    range->SetLineColor(kTeal+2);
    range->SetMarkerColor(kTeal+2);
    range->SetLineStyle(1);
    range->SetMarkerStyle(0);
  }
  
  // TString fcnfile = Form("systematicsB1_Fit_params_cent%i.root", icent);
  // sysHistoName = "hSystVsPtPercentageOfCentral_rms";
  // TFile * fFuncSyst = TFile::Open(Form("%s/%s",fPath.Data(),fcnfile.Data()));
  // TH1F * dummyF = 0x0;
  // if (fFuncSyst) {
  //   dummyF = (TH1F*) fFuncSyst->Get(sysHistoName.Data());
  //   function = (TH1F*) dummyF->Clone("function"); 
  //   function->SetTitle("Res. bg. fit function");
  //   function->SetLineWidth(2);
  //   function->SetLineColor(kBlue+2);
  //   function->SetMarkerColor(kBlue+2);
  //   function->SetLineStyle(1);
  //   function->SetMarkerStyle(0); 
  // }
  
  // TString PIDfile = Form("systematicsB1_PID_cent%i.root", icent);
  // sysHistoName = "hSystVsPtPercentageOfCentral_rms";
  // TFile * fPidSyst = TFile::Open(Form("%s/%s",fPath.Data(),PIDfile.Data()));
  // TH1F * dummyP = 0x0;
  // if (fPidSyst) {
  //   dummyP = (TH1F*) fPidSyst->Get(sysHistoName.Data());
  //   pid = (TH1F*) dummyP->Clone("PID"); 
  //   pid->SetTitle("PID");
  //   pid->SetLineWidth(3);
  //   pid->SetLineColor(kBlack);
  //   pid->SetMarkerColor(kBlack);
  //   pid->SetLineStyle(2);
  //   pid->SetMarkerStyle(0); 
  // }
  
//ITS-TPC matching efficiency sys uncertainty in pp 5 TeV LHC15n
  //https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsTrackSystematicUncertaintyBookkeping
  //https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsTrackSystematicUncertainty
  //Summary of recommended values of the systematic error for the single track:
  //Summary of recommended values of the systematic error: 2.5-4% for 0.5<pT<5 GeV /c, 2.5-1% for 5<pT<10, 1% for 10<pT<15, values for pT>15 should be studied for each specific analysis.
  //When considering a resonance in two-body decay, multiply the uncert. on the single track times 2
  //as the two daughters are considered independently tracked thus independently affected by the systemtics
  Double_t tracking_1trkBelow2GeV = 0.03;
  Double_t tracking_1trk2to5GeV = 0.04;
  Double_t tracking_1trk5to8GeV = 0.02;
  Double_t tracking_1trkAbove8GeV = 0.015;

  //track cuts systematics uncertainty inherited from the analysis of pi,K,p and pt-independent
  //When considering a resonance in two-body decay, multiply the uncert. on the single track times 2
  //as the two daughters are considered independently tracked thus independently affected by the systemtics
  Double_t trackCuts_rsn = 0.025;

  for (Int_t ii = 0;ii<npt;ii++){
    //uncert. on phi-->KK BR equal to 1%    
    Int_t ibin = ii+1;
    branching->SetBinContent(ibin, 1.0);
    
    if (pt[ibin] < 2.2) {
      tracking->SetBinContent(ibin, 5.0);
    } else if (pt[ibin]>2.2 && pt[ibin]< 4.5) {
      tracking->SetBinContent(ibin, 5.0);
    } else if (pt[ibin]>4.0 && pt[ibin]< 6.5) {
      tracking->SetBinContent(ibin, 8.0);
    } else if (pt[ibin]>6.0) {
      tracking->SetBinContent(ibin, 10.0);
    }  //FIXME
    trackcuts->SetBinContent(ibin, trackCuts_rsn*2.0*100.0);
  }

  //sum all contributions in quadrature
  Double_t syst_ptFullyCorr_sum2 = 0.01; //(BR uncert = 1%)^2 
  
  //estimate uncertainty per each pt bin
  for (Int_t ii = 0;ii<npt;ii++){
    Int_t ibin = ii+1;
    
#if 0 //not a systematic effect    
    Double_t bincount_syst = dummyBC->GetBinContent(ibin)/100.;
    Double_t evmix_syst = dummyEM->GetBinContent(ibin)/100.;
    Double_t width_syst = dummyW->GetBinContent(ibin)/100.;
    Double_t pid_syst = dummyP->GetBinContent(ibin)/100.;
    Double_t func_syst = dummyF->GetBinContent(ibin)/100.;
#endif
    Double_t tracking_sys = tracking->GetBinContent(ibin)/100.;
    Double_t material_syst = dummyMT->GetBinContent(ibin)/100.;
    Double_t hadrint_syst = dummyHI->GetBinContent(ibin)/100.;
    Double_t range_syst = dummyR->GetBinContent(ibin)/100.;
    Double_t bgnorm_syst = dummyBN->GetBinContent(ibin)/100.;
    
    Double_t ptuncorr2 = tracking_sys*tracking_sys;
    
    //track cuts
    //ptuncorr2+=(trackCuts_rsn*trackCuts_rsn);
    
    //fit range
    if (range_syst>0.0) {
      ptuncorr2+=range_syst*range_syst;
    }
    
    //norm bg
    if (bgnorm_syst>0.0) {
      ptuncorr2+=bgnorm_syst*bgnorm_syst;
    } 

    // //function res bg
    // if (func_syst>0.0){
    //   ptuncorr2+=func_syst*func_syst;      
    // }

    // //PID
    // if (pid_syst>0.0) {
    //   ptuncorr2+=pid_syst*pid_syst;
    // }

   //material budget
    if (material_syst>0.0) {
      ptuncorr2+=material_syst*material_syst;
    }

    //hadronic interaction cross section
    if (hadrint_syst>0.0) {
      ptuncorr2+=hadrint_syst*hadrint_syst;
    }

    
    Double_t totsyst = TMath::Sqrt(syst_ptFullyCorr_sum2+ptuncorr2);
    sum2->SetBinContent(ibin, totsyst*100.);
    Printf("bin %i tot. perc. %6.4f", ibin, totsyst*100.);
    
    Double_t totsystUncorr = TMath::Sqrt(ptuncorr2);
    sum2_uncorr->SetBinContent(ibin, totsystUncorr*100.);
  }
  
  //assign syst err to data
  TFile * fdata = TFile::Open(Form("%s/%s",fPathCorr.Data(),corrFile.Data()));
  TH1F * data = (TH1F*) fdata->Get(Form("%s",hCorrYieldName.Data()));  
  TH1F * data_Wsyst = (TH1F*) data->Clone(Form("%s%i_syst",hCorrYieldName.Data(),icent));
  TH1F * data_Wsyst_uncorr = (TH1F*) data->Clone(Form("%s%i_syst_uncorr",hCorrYieldName.Data(),icent));
  TH1F * data_Wsyst_Wstat = (TH1F*) data->Clone(Form("%s%i_syst_stat",hCorrYieldName.Data(),icent));
  
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
  }
  
  //systematic uncertainty plot
  TCanvas *cs=new TCanvas("cs","Systematic uncertainty vs p_{t}", 750,600);
  cs->cd();
  sum2->SetTitle("Total syst. uncert.");
  sum2->GetYaxis()->SetTitle("relative uncert. (%)");
  sum2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  sum2->GetYaxis()->SetRangeUser(0.1, 35);
  sum2->GetXaxis()->SetRangeUser(0.0, 9.9);

  sum2->Draw();
  range->Draw("same");
  // function->Draw("same");
  // pid->Draw("same");
  //  bincount->Draw("same");
  //  evtmix->Draw("same");
  //width->Draw("same");
  branching->Draw("same");
  bgnorm->Draw("same");
  tracking->Draw("same");
  //trackcuts->Draw("same");
  material->Draw("same");
  hadrint->Draw("same");
  //sum2_uncorr->Draw("same");
  
  TLegend * autolegry = (TLegend*)gPad->BuildLegend(0.25,0.60,0.88,0.88, Form("V0M %s",centLabel.Data()));
  autolegry->SetFillColor(kWhite);
  autolegry->SetLineColor(kWhite);
  autolegry->SetTextFont(42);
  autolegry->SetNColumns(2); 
  //color data
  data->SetMarkerColor(color[1][icent]+1);
  data->SetLineColor(color[1][icent]+1);
  data->SetMarkerStyle(marker[1][icent]);
  data->SetMarkerSize(0.7);
  data->SetLineWidth(2);
  data->SetTitle(Form("%s (stat. uncert.)",centLabel.Data()));
  data_Wsyst->SetMarkerColor(color[1][icent]);	       
  data_Wsyst->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst->SetLineColor(color[1][icent]);
  data_Wsyst->SetMarkerStyle(0);
  data_Wsyst->SetOption("E2");
  data_Wsyst->SetTitle(Form("%s (syst. uncert.)",centLabel.Data()));
  // TPaveText * barlowLabel = new TPaveText(0.25,0.57,0.88,0.65,"NDC");
  // barlowLabel->SetFillColor(kWhite);
  // barlowLabel->SetBorderSize(0);
  // barlowLabel->SetTextFont(42);
  //if (barlow>0) barlowLabel->InsertText(Form("Stat. compatibility test for N_{b} = #Delta_{i}/#sigma_{i} = %i",barlow));
  //else barlowLabel->InsertText("No stat. compatibility test");
  cs->cd();
  //barlowLabel->Draw();

  TString imagefilename;
  if (barlow>0) imagefilename = Form("summaryAllSystUncert_%s_B%i_cent%i_%s", bg.Data(), barlow, icent, date.Data()) ;
  else imagefilename = Form("summaryAllSystUncert_%s_cent%i_%s", bg.Data(), icent, date.Data());
  if (smooth>0) imagefilename.ReplaceAll("summaryAllSystUncert",Form("summaryAllSystUncert_smooth%i",smooth));

  cs->Print(Form("%s.png", imagefilename.Data()));
  cs->Print(Form("%s.C", imagefilename.Data()));

  //make-up
  data_Wsyst_uncorr->SetMarkerColor(color[1][icent]);	       
  data_Wsyst_uncorr->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst_uncorr->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst_uncorr->SetLineColor(color[1][icent]);
  data_Wsyst_uncorr->SetMarkerStyle(0);
  data_Wsyst_uncorr->SetOption("E2");
  data_Wsyst_uncorr->SetTitle(Form("%s (syst. uncert., p_{T}-uncorr.)",centLabel.Data()));
  
  data_Wsyst_Wstat->SetMarkerColor(color[1][icent]);	       
  data_Wsyst_Wstat->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst_Wstat->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst_Wstat->SetLineColor(color[1][icent]);
  data_Wsyst_Wstat->SetMarkerStyle(0);
  data_Wsyst_Wstat->SetTitle(Form("%s (#sqrt{syst^{2}+stat^{2}})",centLabel.Data()));

  TCanvas *cunc=new TCanvas("cunc","Summary of uncertainty vs p_{t}", 750,600);
  cunc->cd();
  sum2->Draw();
  sum2_uncorr->Draw("HIST same");
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
  //  vertexing->Write();
  branching->Write();
  tracking->Write();
  trackcuts->Write();
  // if (ipid) {
  //   tofmatching->Write();
  //   tofmatchingMC->Write();
  //   pidResponse->Write();
  // }
  //  bincount->Write();
  // evtmix->Write();
  // pid->Write();
  range->Write();
  bgnorm->Write();
  // function->Write();
  hadrint->Write();
  statunc->Write();
  data->Write();
  sum2->Write();
  sum2_uncorr->Write();
  data_Wsyst->Write();
  data_Wsyst_Wstat->Write();
  data_Wsyst_uncorr->Write();  
  cunc->Write();
  cs->Write();
  fout->Close(); 
  
  return;
  
}

