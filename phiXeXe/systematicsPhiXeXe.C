// #include "/Users/fbellini/alice/macros/MakeUp.C"
// #include "/Users/fbellini/alice/macros/SetStyle.C"
#include "/Users/fbellini/alice/macros/GetPlotRatio.C"

enum ESysSource_t {kPID = 0,
		   kBgNorm, 
		   kBgFit,
		   kLSB,
		   kBreit,
		   kParams,
		   kBinCount,
		   kFitRange,
		   kTrackCuts,
		   kNtypes};


TList * GetListOfAlternativeHistograms(ESysSource_t source = ESysSource_t::kNtypes, Short_t cent = 0, TString histoPrefix = "hRawYieldVsPt_", TString binning = "A3");
TH1D  * GetAlternativeHist(Int_t altID = 0, TString  alternativeFile = "", Short_t cent = 0, TString histoPrefix = "hRawYieldVsPt_");
TH1D  * GetCentralValueHistRaw(Short_t cent = 0, TString histoPrefix = "hRawYieldVsPt_", TString binning = "A3");
TH1D  * GetCentralValueHistCorr(Short_t cent = 0, TString histoPrefix = "hCorrected_");
TString GetSourceName(ESysSource_t source = ESysSource_t::kNtypes);
TCanvas * createCanvas(TString name = "c1", Short_t nrows = 1, Short_t ncols = 1);

Float_t GetSysMaxDeviationPt(TH1D * hDevPt = 0x0);
Float_t GetSysUnifRangePt(TH1D * hDevPt = 0x0);

void runSystematics();


void systematicsPhiXeXe(Short_t cent = 0, 
			TString BarlowOpt = "", // empty for no Barlow, B1 for Barlow 1sigma, B2 for Barlow 2sigma
			ESysSource_t source = ESysSource_t::kFitRange,
			TString binning = "A3")
{
  //style for display
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetFrameBorderSize(0);

  Int_t centrality[5] = {0, 10, 30, 60, 90};
  //Get histo with central value
  TH1D * hCentral = 0X0;
  if (source == ESysSource_t::kPID) hCentral = (TH1D *) GetCentralValueHistCorr(cent, "hCorrected_");
  else hCentral = (TH1D *) GetCentralValueHistRaw(cent, "hRawYieldVsPt_", binning.Data());

  const Int_t nBins = hCentral->GetNbinsX();
  Printf("::::: N. of pT Bins = %i", nBins);

  //Get list for histos with alternative values for given source
  TList * alternativeList = (TList *) GetListOfAlternativeHistograms(source, cent,  "hRawYieldVsPt_", binning.Data());
  if (!alternativeList) return;
  else Printf("::::: We are ready to extract systematics!");

  //Get list for ratio of variation/default vs pt
  TList * ratioList = new TList();
  ratioList->SetName(Form("lRatios%s", GetSourceName(source).Data()));

  //Create histograms with variations and add them to a new list
  TList * varList = new TList();
  varList->SetName(Form("lVar%s", GetSourceName(source).Data()));

  TH1D * hVar[nBins];
  TH1D * hSigma2[nBins];
  TH1D * hNbarlow[nBins];
  TH1D * hDevRatio[nBins];
  TH1D * hDevRatioB1[nBins];
  TH1D * hDevRatioB2[nBins];
  
  for (Int_t ipt = 1; ipt<nBins+1; ipt++){
    if (hCentral->GetBinContent(ipt)<=0) continue;    
    Double_t varMax = hCentral->GetBinContent(ipt)*0.50;
    Double_t maxSigma = TMath::Power(TMath::Sqrt(varMax), 2.0);
    
    hVar[ipt] = new TH1D(Form("hVar_%s_pt%i", GetSourceName(source).Data(), ipt), Form("hVar_%s_pt%i; #Deltax_{i} = x_{i} - x_{default}", GetSourceName(source).Data(), ipt), 100, -varMax, varMax);
    varList->Add(hVar[ipt]);

    hSigma2[ipt] = new TH1D(Form("hSigma2_%s_pt%i", GetSourceName(source).Data(), ipt), Form("hSigma2_%s_pt%i; #Delta#sigma_{i} = #sqrt{|#sigma_{i}^{2} - #sigma_{default}^{2}|}", GetSourceName(source).Data(), ipt), 100, -maxSigma, maxSigma);
    varList->Add(hSigma2[ipt]);

    hNbarlow[ipt] = new TH1D(Form("hNbarlow_%s_pt%i", GetSourceName(source).Data(), ipt), Form("hNbarlow_%s_pt%i; N_{B} = #Deltax_{i} / #Delta#sigma_{i}", GetSourceName(source).Data(), ipt), 200, -10., 10.);
    varList->Add(hNbarlow[ipt]);

    hDevRatio[ipt] = new TH1D(Form("hDevRatio_%s_pt%i", GetSourceName(source).Data(), ipt), Form("hDevRatio_%s_pt%i; 1 - x_{i}/x_{default}", GetSourceName(source).Data(), ipt), 100, -0.5, 0.5);
    varList->Add(hDevRatio[ipt]);

     hDevRatioB1[ipt] = new TH1D(Form("hDevRatioB1_%s_pt%i", GetSourceName(source).Data(), ipt), Form("hDevRatioB1_%s_pt%i; 1 - x_{i}/x_{default}", GetSourceName(source).Data(), ipt), 100, -0.5, 0.5);
    varList->Add(hDevRatioB1[ipt]);

    hDevRatioB2[ipt] = new TH1D(Form("hDevRatioB2_%s_pt%i", GetSourceName(source).Data(), ipt), Form("hDevRatioB2_%s_pt%i; 1 - x_{i}/x_{default}", GetSourceName(source).Data(), ipt), 100, -0.5, 0.5);
    varList->Add(hDevRatioB2[ipt]);
  }

  //create histogram for sys uncert
  TH1D * hsys_rms = (TH1D*) hCentral->Clone("hSystVsPtPercentageOfCentral_rms");
  hsys_rms->Reset("ICES");
  hsys_rms->SetTitle("Relative sys. uncertainty; #it{p}_{T} (GeV/#it{c});  relative uncertainty;");
  hsys_rms->GetYaxis()->SetRangeUser(0., 1.);
  TH1D * hsys_unif = (TH1D*) hsys_rms->Clone("hSystVsPtPercentageOfCentral_unif");
  TH1D * hsys_max = (TH1D*) hsys_rms->Clone("hSystVsPtPercentageOfCentral_max");

  //Number of alternative histograms
  Short_t nAltern = alternativeList->GetEntries(); 

  //prep for display
  TCanvas * c1 = new TCanvas("c1","c1", 800, 600);
  TString centLabel = Form("%i-%i %%",centrality[cent],centrality[cent+1]);
  TLegend * leg = new TLegend(0.2, 0.85-0.05*nAltern/3, 0.75, 0.9, Form("%s", centLabel.Data()));
  leg->SetNColumns(3);
  myLegendSetUp(leg, 0.035);
  
  //loop over list and over pt to get variations in each pt bin
  for (Int_t ih = 0; ih<alternativeList->GetEntries(); ih++){
    
    TH1D * htmp = (TH1D *) alternativeList->At(ih);
    if (!htmp) return;

    TH1D * centralClone = (TH1D*) hCentral->Clone();
    TH1D * tmpClone = (TH1D*) htmp->Clone();
    TH1D * htmpratio = (TH1D*) GetPlotRatio(tmpClone, centralClone, 0, "", htmp->GetName(), hCentral->GetName(), "alternative / default", 0, 0, 0.0, -1.0, "B"); 
    htmpratio->SetName(Form("Ratio_%itoDefault", ih));
    htmpratio->SetLineColor(kAzure+9-ih);
    htmpratio->SetMarkerColor(kAzure+9-ih);
    htmpratio->SetMarkerStyle(20+ih);
    htmpratio->GetYaxis()->SetRangeUser(0.5, 1.5);
    if (source == ESysSource_t::kParams) htmpratio->GetYaxis()->SetRangeUser(0.2, 1.8);
    
    htmpratio->GetYaxis()->SetTitle("alternative / default");
    ratioList->Add(htmpratio);

    //plot ratio to default
    c1->cd();
    if (ih==0) htmpratio->Draw();
    else htmpratio->Draw("same");
    leg->AddEntry(htmpratio, Form("%s%i", GetSourceName(source).Data(), ih), "p");

    for (Int_t ipt = 1; ipt<nBins+1; ipt++){
      if (hCentral->GetBinContent(ipt)<=0) continue;
      
      Double_t var = htmp->GetBinContent(ipt) - hCentral->GetBinContent(ipt);
      Double_t sigma2diff = TMath::Sqrt( TMath::Abs(TMath::Power(htmp->GetBinError(ipt), 2.0) - TMath::Power(hCentral->GetBinError(ipt), 2.0)));
      Double_t nbarlow = 0.0;
      if (sigma2diff!=0) nbarlow = var / sigma2diff;
      Double_t devRatio = 1 - htmp->GetBinContent(ipt) / hCentral->GetBinContent(ipt);
      
      hVar[ipt]->Fill(var);
      hVar[ipt]->SetLineColor(kOrange+7-ipt);
      hVar[ipt]->SetLineWidth(2);
      hSigma2[ipt]->Fill(sigma2diff);
      hNbarlow[ipt]->Fill(nbarlow);
      hDevRatio[ipt]->Fill(devRatio);

      //fill deviation plot only if Nbarlow >1 or >2, otherwise it is considered as statistical fluctuation
      if (TMath::Abs(nbarlow)>1.) hDevRatioB1[ipt]->Fill(devRatio);
      if (TMath::Abs(nbarlow)>2.) hDevRatioB2[ipt]->Fill(devRatio);
    } 
  }

  Printf("::::: N. examined alternative settings for sys. due to %s = %i",  GetSourceName(source).Data(), alternativeList->GetEntries());

  //RMS or MAX deviation to be implemented to get fractional sys
  for (Int_t ipt = 1; ipt<nBins+1; ipt++){
    if (hCentral->GetBinContent(ipt)<=0) continue;
    TH1D * hDevPt = (TH1D*) varList->FindObject(Form("hDevRatio%s_%s_pt%i", BarlowOpt.Data(), GetSourceName(source).Data(), ipt));
    if (!hDevPt) {Printf("Invalid histogram with deviations for pt bin %i", ipt); continue;}
    
    Float_t relSysPt = 0.0;
    Float_t maxDevPt = 0.0;
    Float_t unifDevPt = 0.0;
    
    Int_t nUsed = hDevPt->GetEntries();
    // compute sys - zero alternative after Barlow
    if (nUsed<1) {
      Printf("\n:::: Pt bin %i --> %i alternative measurements pass Barlow. Rel. sys. = 0.", ipt, nUsed);
      continue;
    }

    // compute sys - 1 alternative after Barlow - take it as sys
    if (nUsed==1) {
      relSysPt = TMath::Abs(hDevPt->GetMean()); //absolute value symmetrises the uncertainty
      Printf("\n:::: Pt bin %i --> %i alternative measurements passed Barlow. Rel. sys. = %4.3f",ipt, nUsed, relSysPt);
    }

    // compute sys - >1 alternative after Barlow - RMS
    if (nUsed>1) {
      relSysPt = hDevPt->GetRMS();
      Printf("\n:::: Pt bin %i --> %i alternative measurements passed Barlow. \nRel. sys. RMS = %4.3f", ipt, nUsed, relSysPt);
      
      // compute sys - >1 alternative after Barlow - max deviation
      maxDevPt = GetSysMaxDeviationPt(hDevPt);
      Printf("Rel. sys. MAX = %4.3f", maxDevPt);

      // compute sys - >1 alternative after Barlow - (max-min)/√12 --> FIXME: should take the range in yield and not in variation
      unifDevPt = GetSysUnifRangePt(hDevPt);
      Printf("Rel. sys. UNIF = %4.3f", unifDevPt);
    }
    
    hsys_max->SetBinContent(ipt, maxDevPt);  hsys_max->SetBinError(ipt, 0.0);
    hsys_rms->SetBinContent(ipt, relSysPt); hsys_rms->SetBinError(ipt, 0.0);
    hsys_unif->SetBinContent(ipt, unifDevPt);  hsys_unif->SetBinError(ipt, 0.0);
    
  }
  c1->cd();
  leg->Draw();

  //save summary plots to file
  TFile * fout = new TFile(Form("systematics%s_%s_cent%i.root", BarlowOpt.Data(), GetSourceName(source).Data(), cent), "recreate");
  fout->cd();
  hCentral->Write();
  hsys_rms->Write();
  hsys_max->Write();
  hsys_unif->Write();
  alternativeList->Write();
  ratioList->Write();
  varList->Write();
  c1->Write();
  //fout->Close();
  c1->Print(Form("systematics%s_%s_cent%i.pdf", BarlowOpt.Data(), GetSourceName(source).Data(), cent));
  return;
}

void runSystematics()
{
  // run all systematics for all centralities
  for (Int_t is = 1; is<ESysSource_t::kNtypes; is ++){
    for (Int_t ic = 0; ic<4; ic ++){
      systematicsPhiXeXe(ic, "", (ESysSource_t) is, "A3");
    }
  }
  return; 
}

TString GetSourceName(ESysSource_t source)
{  
  switch (source) {
  case ESysSource_t::kPID:
    return "PID";
  case ESysSource_t::kBgNorm:
    return "Bg_norm";
  case ESysSource_t::kBgFit:
    return "Bg_fit";
  case ESysSource_t::kLSB:
    return "LSB";
  case ESysSource_t::kBreit:
    return "Breit-Wigner";
  case ESysSource_t::kParams:
    return "Fit_params";
  case ESysSource_t::kBinCount:
    return "Bin_count";
  case ESysSource_t::kFitRange:
    return "Fit_range";
  case ESysSource_t::kTrackCuts:
    return "Track_cuts";
  default: 
    return "none";
  }
  
}


TH1D * GetAlternativeHist(Int_t altID, TString alternativeFile, Short_t cent, TString histoPrefix)
{
  //Get central value from file
  if (!histoPrefix || !alternativeFile) return 0x0;
 
  TFile * fin = TFile::Open(alternativeFile.Data());
  if (!fin) return 0x0;
  TH1D * h = (TH1D *) fin->Get(Form("%s%i", histoPrefix.Data(), cent));
  if (!h) return 0;
  h->SetName(Form("%s%i_alt%i", histoPrefix.Data(), cent, altID));
  h->SetTitle(Form("%s%i_alt%i", histoPrefix.Data(), cent, altID));
  Printf("::::: Alternative measurement loaded \n----- id %i -- file: %s", altID, alternativeFile.Data());
  return h;
}


TH1D * GetCentralValueHistRaw(Short_t cent, TString histoPrefix, TString binning)
{
  //Get central value from file
  
  if (!histoPrefix) return 0x0;
  
  TString centralValueFile = Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070/RAW_fitResult.root", binning.Data());
  TFile * fin = TFile::Open(centralValueFile.Data());
  if (!fin) return 0x0;
  TH1D * h = (TH1D *) fin->Get(Form("%s%i", histoPrefix.Data(), cent));
  if (!h) return 0;
  Printf("::::: Central measurement loaded \n----- file: %s", centralValueFile.Data());
  return h;
}

TH1D * GetCentralValueHistCorr(Short_t cent, TString histoPrefix)
{
  //Get central value from file
  
  if (!histoPrefix) return 0x0;
  TString centralValueFile = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiC3_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_Wfix/fit_r0.994-1.070/CORRECTED_br_fitResult.root";
  
  TFile * fin = TFile::Open(centralValueFile.Data());
  if (!fin) return 0x0;
  TH1D * h = (TH1D *) fin->Get(Form("%s%i", histoPrefix.Data(), cent));
  if (!h) return 0;
  Printf("::::: Central measurement loaded \n----- file: %s", centralValueFile.Data());
  return h;
}


TCanvas * createCanvas(TString name, Short_t nrows, Short_t ncols)
{
  TCanvas * c = new TCanvas(name.Data(), name.Data(), 800*ncols, 600*nrows);
  c->Divide(nrows, ncols);
  return c;
}

TList * GetListOfAlternativeHistograms(ESysSource_t source, Short_t cent, TString histoPrefix, TString binning)
{

  if (source<0) return 0x0;
  TList * list = new TList();
  list->SetName(GetSourceName(source));
  Int_t i = -1;

  switch (source) {
  case ESysSource_t::kPID:
    list->AddLast((TH1D*)GetAlternativeHist(++i, Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0414pidSys/phi%s_tpc2sPtDep_tof4sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070/CORRECTED_br_fitResult.root", binning.Data()), cent, "hCorrected_"));
    list->AddLast((TH1D*)GetAlternativeHist(++i, "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0414pidSys/phiC3_tpc3sPtDep_tof3sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070/CORRECTED_br_fitResult.root", cent, "hCorrected_"));
    break;

  case ESysSource_t::kBgNorm:
    list->AddLast((TH1D*)GetAlternativeHist(++i, Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.06-1.08/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    list->AddLast((TH1D*)GetAlternativeHist(++i, Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.05-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    list->AddLast((TH1D*)GetAlternativeHist(++i, Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.10-1.20/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    break;
    
  case ESysSource_t::kBgFit:
    list->AddLast((TH1D*)GetAlternativeHist(++i, Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly2_fixW/fit_r0.992-1.100/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    break;
    
  case ESysSource_t::kLSB:
    list->AddLast((TH1D*)GetAlternativeHist(++i, Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Like_VOIGTpoly2_fixW/fit_r0.992-1.100/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    break;

  case ESysSource_t::kBreit:
    list->AddLast((TH1D*)GetAlternativeHist(++i, Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_BREITpoly1_allFree/fit_r0.994-1.070/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    break;
    
  case ESysSource_t::kParams:
    list->AddLast((TH1D*) GetAlternativeHist(++i, Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_allFixed/fit_r0.994-1.070/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    list->AddLast((TH1D*) GetAlternativeHist(++i,  Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_ResRMS1_fixW/fit_r0.994-1.070/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data())); 
    list->AddLast((TH1D*) GetAlternativeHist(++i, Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_ResGaus0_fixW/fit_r0.994-1.070/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    list->AddLast((TH1D*) GetAlternativeHist(++i,  Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_ResRMS3_fixW/fit_r0.994-1.070/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    list->AddLast((TH1D*) GetAlternativeHist(++i,  Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_ResRMS1/fit_r0.994-1.070/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data())); 
    list->AddLast((TH1D*) GetAlternativeHist(++i,  Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_ResGaus0/fit_r0.994-1.070/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    list->AddLast((TH1D*) GetAlternativeHist(++i,  Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_ResRMS3/fit_r0.994-1.070/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    break;

    
  case ESysSource_t::kBinCount:
    list->AddLast((TH1D*)GetAlternativeHist(++i, Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.994-1.070/RAW_fitResult.root", binning.Data()), cent, "hRawYieldBCVsPt_"));
    break;
    
  case ESysSource_t::kFitRange:
    list->AddLast((TH1D*)GetAlternativeHist(++i,  Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.992-1.100/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    list->AddLast((TH1D*)GetAlternativeHist(++i,  Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.996-1.070/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    list->AddLast((TH1D*)GetAlternativeHist(++i,  Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phi%s_tpc2sPtDep_tof2sveto5smism/norm1.07-1.10/fit_Mixing_VOIGTpoly1_fixW/fit_r0.998-1.060/RAW_fitResult.root", binning.Data()), cent, histoPrefix.Data()));
    break;
    
  case ESysSource_t::kTrackCuts:
    break;
    
  default: 
    Printf(":::: No systematic source selected. Doing nothing.");
  }
  
  if (i<0) {
    Printf("::::: WARNING: your list of alternative histograms seems to be empty!");  
    list = 0x0;
  }
  return list;
}


Float_t GetSysMaxDeviationPt(TH1D * hDevPt)
{
  if (!hDevPt) return 0.0;
		 
  Float_t maxDevPt = 0.0;
  for (int i = 1; i<hDevPt->GetNbinsX(); i++) {
    Float_t xvar = hDevPt->GetXaxis()->GetBinCenter(i);
    if ( (hDevPt->GetBinContent(i)>0) && (TMath::Abs(xvar) > maxDevPt) )
      maxDevPt = TMath::Abs(xvar);
  }
  return maxDevPt;
}

Float_t GetSysUnifRangePt(TH1D * hDevPt)
{
  if (!hDevPt) return 0.0;
		 
  Float_t maxDev = 0.0;
  Float_t minDev = 999.0;
  for (int i = 1; i<hDevPt->GetNbinsX(); i++) {
    Float_t xvar = hDevPt->GetXaxis()->GetBinCenter(i);
    if (hDevPt->GetBinContent(i)>0) {
      if (xvar < minDev) minDev = xvar;
      if (TMath::Abs(xvar) > maxDev) maxDev = xvar;
    }
  }
  return (maxDev - minDev)/TMath::Sqrt(12.);
}
