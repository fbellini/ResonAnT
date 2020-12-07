#include "/Users/fbellini/alice/macros/cosmetics/MakeUp.C"
void SmoothenSysPtRange(TH1F * hist, Float_t ptmin, Float_t ptmax);

void syst_contributionsXeXecent(Int_t icent=0, 
				TString date="04dec20", 
				Bool_t smooth = 1,
				Bool_t useReweigthed = 0)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  myOptions(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.17);
  //set input name
  //PREL TString fPath = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiA3_tpc2sPtDep_tof2sveto5smism/systematics";
  //PREL TString fPathCorr = "~/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiA3_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070";
  TString fPath = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/final/analysis/systematics";
  TString fPathCorr = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/final/analysis/phifinal_default_LowBdca/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070";
  TString corrFile = "CORRECTED_br_fitResult.root"; 
  TString hCorrYieldName = Form("hCorrected_%i",icent); 

  if (useReweigthed){
    fPathCorr = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/final/analysis/phifinal_default_LowBdca/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070";
    //PREL "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiA3_tpc2sPtDep_tof2sveto5smism/spectra";
    corrFile = "REWEIGHT_CORRECTED_br_fitResult.root";
  }
  
  //Binning for XeXe analysis - preliminary
  Double_t cent[] = {0.0, 10.0, 30.0, 50.0, 70.0, 90.0};   
  Double_t pt[]   = {0.0, 0.3, 0.5, 0.7, 0.9, 1.10, 1.30, 1.50, 2.00, 3.00, 4.00, 5.0, 7.0, 10.0};  
  Int_t    npt    = sizeof(pt) / sizeof(pt[0]) - 1;   
  Int_t   ncent   = sizeof(cent) / sizeof(cent[0]) - 1;
  TString centLabel = Form("%3.0f-%3.0f%%",cent[icent], cent[icent+1]);
  
  //cosmetics  
  Color_t color[5] = {kRed+1, kOrange, kSpring+5, kBlue+1, kMagenta+2};
  Int_t  marker[5] = {20, 21, 34, 33, 22}; 
  
  //create axis to reproduce the binning
  TAxis *ptbins = new TAxis(npt, pt);
  TAxis *centbins = new TAxis(ncent, cent);
  
  //statistical 
  TH1F * statunc = new TH1F(Form("statunc_%i",icent),Form("Statistical uncertainty"), npt, pt);
  statunc->SetLineWidth(3);
  statunc->SetLineColor(kGreen+2);
  statunc->SetMarkerColor(kGreen+2);
  statunc->SetLineStyle(1);
  statunc->SetMarkerStyle(0);

  //total systematics
  TH1F * sum2 = new TH1F(Form("sum2_%i",icent),Form("Sum"), npt, pt);
  sum2->SetLineWidth(3);
  sum2->SetLineColor(kRed);
  sum2->SetMarkerColor(kRed);
  sum2->SetLineStyle(1);
  sum2->SetMarkerStyle(0); 

  //total uncorrelated with pi
  TH1F * sum2_uncorr = new TH1F(Form("sum2_uncorr_%i",icent),Form("Uncorrelated sys. uncert."), npt, pt);
  sum2_uncorr->SetLineWidth(3);
  sum2_uncorr->SetLineColor(kGray+2);
  sum2_uncorr->SetMarkerColor(kGray+2);
  sum2_uncorr->SetLineStyle(1);
  sum2_uncorr->SetMarkerStyle(0); 
  
  //pt dependent
  TH1F * tracking = new TH1F(Form("tracking_%i",icent), Form("ITS-TPC matching"), npt, pt);
  TH1F * trackcuts = new TH1F(Form("trackcuts_%i",icent), Form("Track cuts"), npt, pt);
  TH1F * material = new TH1F(Form("material_%i",icent), Form("Material Budget"), npt, pt);
  TH1F * hadrint = new TH1F(Form("hadrint_%i",icent),Form("Hadronic inter."), npt, pt);
  TH1F * bgnorm = new TH1F(Form("bgnorm_%i",icent),Form("Bg. norm."), npt, pt);
  TH1F * range = new TH1F(Form("range_%i",icent),Form("Fit range"), npt, pt);
  TH1F * function = new TH1F(Form("bgfunc_%i",icent),Form("Res. bg. fit function"), npt, pt);
  TH1F * parameters = new TH1F(Form("params_%i",icent),Form("Fit parameters"), npt, pt);
  TH1F * pid = new TH1F(Form("pid_%i",icent), Form("PID"), npt, pt);
  TH1F * branching = new TH1F("branching","B.R.",npt, pt);
  
  branching->SetLineWidth(2);
  branching->SetLineColor(kMagenta);
  branching->SetMarkerColor(kMagenta);
  branching->SetLineStyle(5);
  branching->SetMarkerStyle(0);
  for (int j=1; j<npt+1; j++){
    branching->SetBinContent(j, 0.01);
    branching->SetBinError(j, 0.0);
  }

  //for (int j=1; j<npt+1; j++){
  //  pid->SetBinContent(j, 0.08);
  //  pid->SetBinError(j, 0.0);
  //}
  // TString PIDfile = Form("systematicsB1_PID_cent%i.root", icent);
  // sysHistoName = "hSystVsPtPercentageOfCentral_rms";
  // TFile * fPidSyst = TFile::Open(Form("%s/%s",fPath.Data(),PIDfile.Data()));
  // TH1F * dummyP = 0x0;
  // if (fPidSyst) {
  //   dummyP = (TH1F*) fPidSyst->Get(sysHistoName.Data());
  //   pid = (TH1F*) dummyP->Clone("PID");
  // }
  TString pidfile = Form("systematicsB1_PID_cent%i.root", icent);
  TFile * fPidSyst = TFile::Open(Form("%s/%s",fPath.Data(), pidfile.Data()));
  //if (fPidSyst) {
  TH1F * dummyPid = (TH1F*) fPidSyst->Get("hSystVsPtPercentageOfCentral_rms");
  pid = (TH1F*) dummyPid->Clone("pid");
  pid->SetTitle("PID");
  pid->SetLineWidth(4);
  pid->SetLineColor(kBlack);
  pid->SetMarkerColor(kBlack);
  pid->SetLineStyle(3);
  pid->SetMarkerStyle(0);
  
  float smallq = 0.001;

  if (smooth) {
    SmoothenSysPtRange(pid, 0.5+smallq,1.0-smallq);
    SmoothenSysPtRange(pid, 1.0+smallq,1.5-smallq);
    SmoothenSysPtRange(pid, 1.5+smallq,3.0-smallq);
    if (icent==4) SmoothenSysPtRange(pid, 3.0+smallq, 7.-smallq);
    else SmoothenSysPtRange(pid, 3.0+smallq, 7.-smallq);
  }
  //ITS TPC matching from
  //https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsTrackSystematicUncertaintyBookkeping
  //https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsTrackSystematicUncertainty  
  TFile * fITSTPCmatching = TFile::Open(Form("%s/systITSTPCmatch.root",fPath.Data()));
  TH1F * dummyTtracking = (TH1F*) fITSTPCmatching->Get(Form("hSystVsPtPercentageOfCentral"));
  tracking = (TH1F*) dummyTtracking->Clone("ITSTPCmatching");
  tracking->SetTitle("ITS-TPC matching");
  tracking->Scale(0.01);
  tracking->SetLineWidth(2);
  tracking->SetLineColor(kOrange);
  tracking->SetMarkerColor(kOrange);
  tracking->SetLineStyle(9);
  tracking->SetMarkerStyle(0); 

  TFile * fMaterial = TFile::Open(Form("%s/systMaterial.root",fPath.Data()));
  TH1F * dummyMT = (TH1F*) fMaterial->Get(Form("hSystVsPtPercentageOfCentral"));
  material = (TH1F*) dummyMT->Clone("material"); 
  material->SetTitle("Material budget");
  material->Scale(0.01);
  material->SetLineWidth(4);
  material->SetLineColor(kPink+2);
  material->SetMarkerColor(kPink+2);
  material->SetLineStyle(3);
  material->SetMarkerStyle(0); 
  
  TFile * fHadrSyst = TFile::Open(Form("%s/systHadrInt.root",fPath.Data()));
  TH1F * dummyHI = (TH1F*) fHadrSyst->Get(Form("hSystVsPtPercentageOfCentral"));
  hadrint = (TH1F*) dummyHI->Clone("hadrint"); 
  hadrint->SetTitle("Hadronic int.");
  hadrint->Scale(0.01);
  hadrint->SetLineWidth(2);
  hadrint->SetLineColor(kOrange+1);
  hadrint->SetMarkerColor(kOrange+1);
  hadrint->SetLineStyle(1);
  hadrint->SetMarkerStyle(0);

  TFile * fTrackcuts = TFile::Open(Form("%s/systTrackCuts%i.root", fPath.Data(),icent));
  TH1F * dummyTrackcuts = (TH1F*) fTrackcuts->Get(Form("hSystVsPtPercentageOfCentral"));
  trackcuts = (TH1F*) dummyTrackcuts->Clone("trackCuts");
  trackcuts->SetTitle("Track cuts");
  trackcuts->Scale(0.01);
  trackcuts->SetLineWidth(2);
  trackcuts->SetLineColor(kAzure+10);
  trackcuts->SetMarkerColor(kAzure+10);
  trackcuts->SetLineStyle(2);
  trackcuts->SetMarkerStyle(0); 
  
  TString normfile = Form("systematicsB1_Bg_norm_cent%i.root", icent);
  TFile * fBgNorm = TFile::Open(Form("%s/%s",fPath.Data(),normfile.Data()));
  TH1F * dummyBN = (TH1F*) fBgNorm->Get("hSystVsPtPercentageOfCentral_rms");
  bgnorm = (TH1F*) dummyBN->Clone("bgNorm");
  bgnorm->SetTitle("MEB normalisation");
  bgnorm->SetLineWidth(2);
  bgnorm->SetLineColor(kPink+8);
  bgnorm->SetMarkerColor(kPink+8);
  bgnorm->SetLineStyle(8);
  bgnorm->SetMarkerStyle(0);
  if (smooth) {
    SmoothenSysPtRange(bgnorm, 0.5+smallq,1.3-smallq);
    SmoothenSysPtRange(bgnorm, 1.3+smallq, 3.0-smallq);
    if (icent==4) SmoothenSysPtRange(bgnorm, 3.0+smallq, 7.0-smallq);
    else SmoothenSysPtRange(bgnorm, 3.0+smallq, 10.0-smallq);
  }

  TString rangefile = Form("systematicsB1_Fit_range_cent%i.root", icent);
  TFile * fRangeSyst = TFile::Open(Form("%s/%s",fPath.Data(),rangefile.Data()));
  TH1F * dummyR = 0x0;
  if (fRangeSyst) {
    dummyR = (TH1F*) fRangeSyst->Get("hSystVsPtPercentageOfCentral_rms");
    range = (TH1F*) dummyR->Clone("range");
    range->SetTitle("Fit range");
    range->SetLineWidth(3);
    range->SetLineColor(kTeal+2);
    range->SetMarkerColor(kTeal+2);
    range->SetLineStyle(1);
    range->SetMarkerStyle(0);
  }
  if (smooth) {
    SmoothenSysPtRange(range, 1.5+smallq,3.0-smallq);
    if (icent==4) {
      SmoothenSysPtRange(range, 0.7+smallq,1.5-smallq);
      SmoothenSysPtRange(range, 3.0+smallq,7.-smallq);
    } else { 
      SmoothenSysPtRange(range, 0.5+smallq,1.5-smallq);
      SmoothenSysPtRange(range, 3.0+smallq,10.-smallq);
    }
  }

  TFile * fFuncSyst = TFile::Open(Form("%s/systematicsB1_Bg_fit_cent%i.root", fPath.Data(), icent));
  TH1F * dummyF = 0x0;
  if (fFuncSyst) {
    dummyF = (TH1F*) fFuncSyst->Get("hSystVsPtPercentageOfCentral_rms");
    function = (TH1F*) dummyF->Clone("function"); 
    function->SetTitle("Res. bg. fit function");
    function->SetLineWidth(2);
    function->SetLineColor(kBlue);
    function->SetMarkerColor(kBlue);
    function->SetLineStyle(7);
    function->SetMarkerStyle(0); 
  }
  if (smooth) { 
    SmoothenSysPtRange(function, 1.0+smallq,1.5-smallq);
    SmoothenSysPtRange(function, 1.5+smallq,3.0-smallq);
    if (icent==4) SmoothenSysPtRange(function, 3.0+smallq,7.-smallq);
    else SmoothenSysPtRange(function, 3.0+smallq,10.-smallq);
  }

  TFile * fParamSyst = TFile::Open(Form("%s/systematicsB1_Fit_params_cent%i.root", fPath.Data(), icent));
  TH1F * dummyPar = 0x0;
  if (fParamSyst) {
    dummyPar = (TH1F*) fParamSyst->Get("hSystVsPtPercentageOfCentral_rms");
    parameters = (TH1F*) dummyPar->Clone("FitParams"); 
    parameters->SetTitle("Fit parameters");
    parameters->SetLineWidth(2);
    parameters->SetLineColor(kBlue+2);
    parameters->SetMarkerColor(kBlue+2);
    parameters->SetLineStyle(1);
    parameters->SetMarkerStyle(0); 
  }

  if (smooth) {
    if (icent==4) SmoothenSysPtRange(parameters, 0.7+smallq, 1.5-smallq);
    else if (icent==1) SmoothenSysPtRange(parameters, 0.5+smallq, 3.0-smallq);
    else SmoothenSysPtRange(parameters, 0.5+smallq, 1.5-smallq);
    SmoothenSysPtRange(parameters, 3.0+smallq, 7.0-smallq);
  }
  
  //sum all contributions in quadrature
  Double_t syst_sum2 = 0.; 
  
  //estimate uncertainty per each pt bin
  for (Int_t ii = 0; ii<npt; ii++){
    Int_t ibin = ii+1;
    syst_sum2 = 0.0;
    Double_t branching_sys = branching->GetBinContent(ibin);
    Double_t tracking_sys  = tracking->GetBinContent(ibin);
    Double_t trackcuts_sys = trackcuts->GetBinContent(ibin);
    Double_t material_syst = material->GetBinContent(ibin);
    Double_t hadrint_syst  = hadrint->GetBinContent(ibin);
    Double_t range_syst    = range->GetBinContent(ibin);
    Double_t bgnorm_syst   = bgnorm->GetBinContent(ibin);
    Double_t parameters_syst = parameters->GetBinContent(ibin);
    Double_t function_syst = function->GetBinContent(ibin);
    Double_t pid_syst = pid->GetBinContent(ibin);

    syst_sum2 += (branching_sys*branching_sys);
    syst_sum2 += (tracking_sys*tracking_sys);
    if (trackcuts_sys>0) syst_sum2 += (trackcuts_sys*trackcuts_sys);
    if (material_syst>0) syst_sum2 += (material_syst*material_syst);
    if (hadrint_syst>0)  syst_sum2 += (hadrint_syst*hadrint_syst);
    if (range_syst>0) syst_sum2 += (range_syst*range_syst);
    if (bgnorm_syst>0) syst_sum2 += (bgnorm_syst*bgnorm_syst);
    if (parameters_syst>0) syst_sum2 += (parameters_syst*parameters_syst);
    if (function_syst>0) syst_sum2 += (function_syst*function_syst);      
    if (pid_syst>0) syst_sum2 += (pid_syst*pid_syst);
 
    Double_t totsyst = TMath::Sqrt(syst_sum2);
    sum2->SetBinContent(ibin, totsyst);
    Printf("bin %i tot. perc. %6.4f", ibin, totsyst);
    
    Double_t totsystUncorr = TMath::Sqrt(syst_sum2);
    sum2_uncorr->SetBinContent(ibin, totsystUncorr);
  }
  
  //assign syst err to data
  TFile * fdata = TFile::Open(Form("%s/%s",fPathCorr.Data(),corrFile.Data()));
  TH1F * data = (TH1F*) fdata->Get(Form("%s",hCorrYieldName.Data()));  
  TH1F * data_Wsyst = (TH1F*) data->Clone(Form("%s%i_syst",hCorrYieldName.Data(),icent));
  TH1F * data_Wsyst_uncorr = (TH1F*) data->Clone(Form("%s%i_syst_uncorr",hCorrYieldName.Data(),icent));
  TH1F * data_Wsyst_Wstat = (TH1F*) data->Clone(Form("%s%i_syst_stat",hCorrYieldName.Data(),icent));
  
  for (Int_t ii = 1;ii<npt+1;ii++){   
    Int_t ibin = ii;
    Double_t yd = data->GetBinContent(ibin);
    Double_t yd_stat = data->GetBinError(ibin);
    Double_t perc = sum2->GetBinContent(ibin);
    Double_t yd_syst = yd*perc;
    statunc->SetBinContent(ibin, (yd>0? (yd_stat/yd) : 0.0));
    data_Wsyst->SetBinContent(ibin,yd);
    data_Wsyst->SetBinError(ibin, yd_syst);
    data_Wsyst_Wstat->SetBinContent(ibin, yd);
    data_Wsyst_Wstat->SetBinError(ibin,TMath::Sqrt(yd_syst*yd_syst+yd_stat*yd_stat));
    Printf("bin %i     yield = %e     syst err = %e    stat err = %e", ibin, yd, yd_syst, yd_stat);
    
    Double_t perc_uncorr = sum2_uncorr->GetBinContent(ibin);
    Double_t yd_syst_uncorr = yd*perc_uncorr;
    data_Wsyst_uncorr->SetBinContent(ibin,yd);
    data_Wsyst_uncorr->SetBinError(ibin, yd_syst_uncorr);
  }

  
  //systematic uncertainty plot
  TCanvas *cs=new TCanvas("cs","Systematic uncertainty vs p_{t}", 750,600);
  cs->cd();
  sum2->SetTitle("Total syst. uncert.");
  sum2->GetYaxis()->SetTitle("relative uncert.");
  sum2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  sum2->GetYaxis()->SetRangeUser(0.001, 0.6);
  sum2->GetXaxis()->SetRangeUser(0.5, 9.9);
  gPad->SetTicky();
  gPad->SetTickx();
  sum2->Draw();
  branching->Draw("same");
  tracking->Draw("same");
  trackcuts->Draw("same");
  material->Draw("same");
  hadrint->Draw("same");
  bgnorm->Draw("same");
  range->Draw("same");
  function->Draw("same");
  parameters->Draw("same");
  pid->Draw("same");
  TLegend * autolegry = (TLegend*)gPad->BuildLegend(0.18,0.62,0.88,0.88, Form("V0M %s",centLabel.Data()));
  autolegry->SetFillColor(kWhite);
  autolegry->SetLineColor(kWhite);
  autolegry->SetTextFont(42);
  autolegry->SetTextSize(0.03);
  autolegry->SetNColumns(3); 
  autolegry->Draw();

  TCanvas *cspectra = new TCanvas("cspectra","Systematic uncertainty vs p_{t}", 750,600);
  Beautify(data, color[icent], 1, 2, marker[icent], 1.3);
  Beautify(data_Wsyst, color[icent], 1, 2, 1, 1.3);
  Beautify(data_Wsyst_Wstat, color[icent], 1, 2, marker[icent], 1.3);
  data_Wsyst->SetTitle(Form("%s (syst. uncert.)",centLabel.Data()));
  data_Wsyst_Wstat->SetTitle(Form("%s (#sqrt{syst^{2}+stat^{2}})",centLabel.Data()));
  cspectra->cd();  data_Wsyst->Draw("E2"); data->Draw("same");


  TString imagefilename = Form("summaryAllSystUncert_cent%i_%s", icent, date.Data());
  if (smooth>0) imagefilename.ReplaceAll("summaryAllSystUncert",Form("summaryAllSystUncert_smooth%i",smooth));
  cs->Print(Form("%s.pdf", imagefilename.Data()));
  cs->Print(Form("%s.png", imagefilename.Data()));

  
  TCanvas *cunc=new TCanvas("cunc","Summary of uncertainty vs p_{t}", 750,600);
  cunc->cd();
  sum2->Draw();
  statunc->Draw("HIST same");
  
  TLegend * autolegry2 = (TLegend*)gPad->BuildLegend(0.25,0.65,0.88,0.88);
  autolegry2->SetFillColor(kWhite);
  autolegry2->SetLineColor(kWhite);
  autolegry2->SetTextFont(42);
  autolegry2->SetNColumns(1); 
  autolegry2->Draw();
  //save to out file
  TString outfilename = Form("finalWsyst_%s_%i.root", date.Data(),icent);
  if (smooth>0) outfilename.ReplaceAll("finalWsyst", Form("finalWsyst_smooth%i",smooth));

  TString imagefilename2 = Form("summaryUncert_cent%i_%s", icent, date.Data());
  if (smooth>0) imagefilename2.ReplaceAll("summaryUncert", Form("summaryUncert_smooth%i",smooth));
  //cunc->SaveAs(Form("%s.png", imagefilename2.Data()));
  cunc->SaveAs(Form("%s.eps", imagefilename2.Data()));

  TFile * fout = new TFile(outfilename.Data(),"recreate");
  fout->cd();
  branching->Write();
  tracking->Write();
  trackcuts->Write();
  material->Write();
  hadrint->Write();
  bgnorm->Write();
  range->Write();
  function->Write();
  parameters->Write();
  pid->Write();
  statunc->Write();
  data->Write();
  sum2->Write();
  sum2_uncorr->Write();
  data->Write();
  data_Wsyst->Write();
  data_Wsyst_Wstat->Write();
  data_Wsyst_uncorr->Write();  
  cunc->Write();
  cs->Write();
  fout->Close(); 
  
  return;
  
}

void SmoothenSysPtRange(TH1F * hist, Float_t ptmin, Float_t ptmax)
{
  if (!hist) return;
  if (ptmax<=ptmin) return;
  Int_t nmin = hist->GetXaxis()->FindBin(ptmin);
  Int_t nmax = hist->GetXaxis()->FindBin(ptmax);
  //calculate average over bins
  Float_t avg = 0.0;
  for (int j = nmin; j < nmax; j++) {
    Printf("(bin=%i - y=%f)", hist->GetBinContent(j));
    avg += hist->GetBinContent(j); 
    }
  avg = avg/(nmax-nmin);
  Printf("Smoothing %i bins of %s, from %f to %f GeV/c: result avg %5.2f", nmax-nmin, hist->GetName(), ptmin, ptmax, avg);
  //TF1 * f1 = new TF1("f1","pol0", ptmin, ptmax);
  //TFitResultPtr fitresult = hist->Fit(f1,"R0SQWL","", ptmin, ptmax);
  //avg = f1->GetParameter(0);
  //set
  for (int j = nmin; j <= nmax; j++) {
    hist->SetBinContent(j, avg);
    hist->SetBinError(j, 0.0);
  }
  return;
}
