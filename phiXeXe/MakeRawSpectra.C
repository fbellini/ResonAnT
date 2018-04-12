/* fbellini@cern.ch, 05.03.2018 */

#include "/Users/fbellini/alice/macros/MakeUp.C"

void MakeRawSpectra(Float_t fitLow = 0.994, Float_t fitUp = 1.050,
		    TString anaPath = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710",
		    TString binning = "C3",
                    TString pid = "tpc2sPtDep_tof2sveto5smism",
		    TString bgType = "Mixing",
		    TString fitFcn = "VOIGTpoly1_fixRes",
		    Float_t normLow = 1.050, Float_t normUp = 1.100,
                    const Float_t ptmin = 0.0, const Float_t ptmax = 9.999,
		    const Float_t cutChi2 = 4.0, Bool_t skipLastBin = 0,
		    Short_t legendEntryStyle = -1)
{
  //graphics
  gStyle->SetTextFont(42);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  myOptions(0);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  TGaxis::SetMaxDigits(2);
  Int_t legendEntryCounter = 1;

  Color_t color[1][3] = {kRed+2, kSpring-2, kBlue+2};
  Int_t  marker[1][3] = {20, 21, 33}; 

  Color_t colorFunc = kBlue+2;
  Color_t colorPID = kTeal+5;
  Color_t colorBg = kRed+2;
  Color_t colorNorm = kAzure+7;

  TString legendEntry = ""; //default, centrality will be appended
  switch (legendEntryStyle) {
  case 0:
    legendEntry = Form("Fit: %4.3f-%4.3f GeV/#it{c}", fitLow, fitUp);
    break;
  case 1 : 
    legendEntry = Form("%s ", pid.Data());
    break;
  case 2 :
    legendEntry = Form("%s ", bgType.Data());
    break;
  case 3 :
    legendEntry = Form("Norm: %4.3f-%4.3f GeV/#it{c} ", normLow, normUp); 
    break;
  case 4:
    legendEntry = Form("%s ", fitFcn.Data());
    break;
  default :
    break;
  }
  
  TString projPath = Form("%s/phi%s_%s/proj_%s.root", anaPath.Data(), binning.Data(), pid.Data(), binning.Data());
  TString fitResultPath = Form("%s/phi%s_%s/norm%3.2f-%3.2f/fit_%s_%s/fit_r%4.3f-%4.3f", anaPath.Data(), binning.Data(), pid.Data(), normLow, normUp, bgType.Data(), fitFcn.Data(), fitLow, fitUp);
  TFile * fin[6] = {0x0,0x0,0x0,0x0,0x0,0x0};

  //get bins
  TFile * f= TFile::Open(projPath.Data());
  if (!f) return;
  //get bins
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  Int_t npt = ptbins->GetNbins();
  const Int_t dimpt = npt+1;
  Double_t pt[dimpt];
  for (Int_t k=0; k<dimpt;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
  }
  TAxis *centbins = (TAxis*)f->Get("centbins");
  Int_t ncent = centbins->GetNbins();
  const Int_t dimcent = ncent+1;
  Double_t cent[dimcent]; 
  for (Int_t k=0; k<dimcent;k++){
    cent[k] = centbins->GetBinLowEdge(k+1);
  }

  // Define output
  TH1D*frameMassVsPt=new TH1D("MassVsPt","Mass vs. #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); #it{M} (GeV/#it{c}^{2})", npt, pt);
  Beautify(frameMassVsPt, kBlack, 1, 2, 1, 1.0);
  frameMassVsPt->GetYaxis()->SetRangeUser(1.016, 1.024);
  frameMassVsPt->GetYaxis()->SetTitleOffset(1.3);
  frameMassVsPt->GetXaxis()->SetRangeUser(ptmin, ptmax);

  TH1D*frameWidthVsPt=new TH1D("WidthVsPt","Width vs. #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); #Gamma (GeV/#it{c}^{2})", npt, pt);
  Beautify(frameWidthVsPt, kBlack, 1, 2, 1, 1.0);
  frameWidthVsPt->GetYaxis()->SetRangeUser(0.0, 0.012);
  frameWidthVsPt->GetXaxis()->SetRangeUser(ptmin, ptmax);
   
  TH1D*frameSigmaVsPt=new TH1D("SigmaVsPt","Sigma vs. #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); #sigma_{Voigt} (GeV/#it{c}^{2})", npt, pt);
  Beautify(frameSigmaVsPt, kBlack, 1, 2, 1, 1.0);
  frameSigmaVsPt->GetYaxis()->SetRangeUser(0.0, 2.5e-3);
  frameSigmaVsPt->GetXaxis()->SetRangeUser(ptmin, ptmax);

  TH1D*frameYieldVsPt=new TH1D("YieldVsPt","Raw yield vs. #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); d#it{N}/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}", npt, pt);
  Beautify(frameYieldVsPt, kBlack, 1, 2, 1, 1.0);
  frameYieldVsPt->GetYaxis()->SetRangeUser(1., 3.0e5);
  frameYieldVsPt->GetYaxis()->SetTitleOffset(1.3);
  frameYieldVsPt->GetXaxis()->SetRangeUser(ptmin, ptmax);

  TH1D*frameYieldBCVsPt=new TH1D("YieldBCVsPt","Raw yield (BC) vs. #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); d#it{N}/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}", npt, pt);
  Beautify(frameYieldBCVsPt, kBlack, 1, 2, 1, 1.0);
  frameYieldBCVsPt->GetYaxis()->SetRangeUser(1., 3.0e5);
  frameYieldBCVsPt->GetYaxis()->SetTitleOffset(1.3);
  frameYieldBCVsPt->GetXaxis()->SetRangeUser(ptmin, ptmax);

  TH1D*frameTailFractionBCVsPt=new TH1D("YieldBCVsPt","Tails / raw yield (BC) vs. #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); fraction of yield in tails", npt, pt);
  Beautify(frameTailFractionBCVsPt, kBlack, 1, 2, 1, 1.0);
  frameTailFractionBCVsPt->GetYaxis()->SetRangeUser(0., 0.15);
  frameTailFractionBCVsPt->GetYaxis()->SetTitleOffset(1.3);
  frameTailFractionBCVsPt->GetXaxis()->SetRangeUser(ptmin, ptmax);
  
  TH1D*frameChi2VsPt=new TH1D("Chi2VsPt","#chi^{2} vs. #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); #chi^{2}/NDF", npt, pt);
  Beautify(frameChi2VsPt, kBlack, 1, 2, 1, 1.0);
  frameChi2VsPt->GetYaxis()->SetRangeUser(0.0, cutChi2);
  frameChi2VsPt->GetXaxis()->SetRangeUser(ptmin, ptmax);

  TH1D*frameSoverBVsPt=new TH1D("SoverBVsPt","S/B; #it{p}_{T} (GeV/#it{c}); #it{S/B}", npt, pt);
  Beautify(frameSoverBVsPt, kBlack, 1, 2, 1, 1.0);
  frameSoverBVsPt->GetYaxis()->SetRangeUser(0.0, 300.0);
  frameSoverBVsPt->GetXaxis()->SetRangeUser(ptmin, ptmax);

  TH1D*frameSignificanceVsPt=new TH1D("SignificanceVsPt","S/#sqrt{S+B}; #it{p}_{T} (GeV/#it{c}); #it{S/#sqrt{S+B}}", npt, pt);
  Beautify(frameSignificanceVsPt, kBlack, 1, 2, 1, 1.0);
  frameSignificanceVsPt->GetYaxis()->SetRangeUser(0.0, 15.0);
  frameSignificanceVsPt->GetXaxis()->SetRangeUser(ptmin, ptmax);

  TLine * pdgmass = new TLine(0.,1.01995,10.,1.01995);
  pdgmass->SetLineStyle(2);
  pdgmass->SetLineWidth(3);
  pdgmass->SetLineColor(kBlack);
  
  TLine * pdgwidth = new TLine(0., 0.00487, 10.,0.00487);
  pdgwidth->SetLineStyle(2);
  pdgwidth->SetLineWidth(3);
  pdgwidth->SetLineColor(kBlack);

  TCanvas *csummary = new TCanvas("csummary","summary of fit", 1200, 900);
  csummary->Divide(3,2);
    
  TCanvas *cm = new TCanvas("MassVsPt","Mass vs #it{p}_{T}", 700,600);
  cm->cd();
  frameMassVsPt->Draw();

  TCanvas *cwi=new TCanvas("WidthVsPt","Width vs #it{p}_{T}", 700,600);
  cwi->cd();
  frameWidthVsPt->Draw();

  TCanvas *cre=new TCanvas("SigmaVsPt","Sigma vs #it{p}_{T}", 700,600);
  cre->cd();
  frameSigmaVsPt->Draw();

  TCanvas *cry=new TCanvas("cry","Raw yield vs #it{p}_{T}", 600,700);
  cry->cd(); gPad->SetLogy();
  frameYieldVsPt->Draw();

  TCanvas *cbc=new TCanvas("cbc","Raw yield (BC) vs #it{p}_{T}", 600,700);
  cbc->cd(); gPad->SetLogy();
  frameYieldBCVsPt->Draw();

  TCanvas *ctail=new TCanvas("ctail","Tail fraction vs #it{p}_{T}", 600,700);
  ctail->cd(); 
  frameTailFractionBCVsPt->Draw();
  
  TCanvas *cc2=new TCanvas("cc2","#Chi^{2} vs #it{p}_{T}", 700,600);
  cc2->cd();
  frameChi2VsPt->Draw();

  TCanvas *cSoB=new TCanvas("cSoB","S/B vs #it{p}_{T}", 700,600);
  cSoB->cd();
  frameSoverBVsPt->Draw();
  
  TCanvas *csign=new TCanvas("csign","Significance vs #it{p}_{T}", 700,600);
  csign->cd();
  frameSignificanceVsPt->Draw();

  
  TH1D *hMassVsPt[6]  =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1D *hWidthVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1D *hSigmaVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1D *hRawYieldVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1D *hRawYieldBCVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1D *hTailFractionBCVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1D *hChi2VsPt[6]     =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1D *hSoverBVsPt[6]   =  {0x0,0x0,0x0,0x0,0x0,0x0};
  TH1D *hSignificanceVsPt[6] =  {0x0,0x0,0x0,0x0,0x0,0x0};


  TLegend *mleg = new TLegend(0.75,0.65,0.95,0.65+0.05*ncent);
  myLegendSetUp(mleg, 0.05);

  TFile * fout = new TFile("RAW_fitResult.root","recreate");
  TString fitResultFile;
  /*loop on centrality bins*/
  for (Int_t icentbin=0;icentbin<ncent;icentbin++){

    //FIXME fitResultFile = Form("%s/result_c%i.root", fitResultPath.Data(), icentbin);
    fitResultFile = Form("result_c%i.root", icentbin);
    fin[icentbin] = TFile::Open(fitResultFile.Data());
    if (!fin[icentbin]) { Printf("Error: cannot find file with fit result. CHECK!"); return; }
    else Printf(">>>>>>> opened file %s", fin[icentbin]->GetName());
    
    TString centLabel=Form("%3.0f-%3.0f %",cent[icentbin],cent[icentbin+1]);
    Printf("\n\n\n*******************\n %s \n*******************", centLabel.Data());
    TString drawOpt = "same";
    
    hMassVsPt[icentbin] = (TH1D*) fin[icentbin]->Get("mass");
    hWidthVsPt[icentbin] = (TH1D*) fin[icentbin]->Get("gamma");
    hSigmaVsPt[icentbin] = (TH1D*) fin[icentbin]->Get("sigma");
    hRawYieldVsPt[icentbin] = (TH1D*) fin[icentbin]->Get("rawIntegral");
    hRawYieldBCVsPt[icentbin] = (TH1D*) fin[icentbin]->Get("rawBC"); hRawYieldBCVsPt[icentbin]->SetMarkerStyle(24);
    hTailFractionBCVsPt[icentbin] = (TH1D*) fin[icentbin]->Get("tailBCfrac");
    hChi2VsPt[icentbin] = (TH1D*) fin[icentbin]->Get("chi2oNDF");
    hSoverBVsPt[icentbin] = (TH1D*) fin[icentbin]->Get("hSoverB");
    hSignificanceVsPt[icentbin] = (TH1D*) fin[icentbin]->Get("hSignif");

    cry->cd();
    hRawYieldVsPt[icentbin]->Draw(drawOpt.Data());
    cbc->cd();
    hRawYieldBCVsPt[icentbin]->Draw(drawOpt.Data());
    hRawYieldVsPt[icentbin]->Draw(drawOpt.Data());
    
    ctail->cd();
    hTailFractionBCVsPt[icentbin]->Draw(drawOpt.Data());
    cm->cd();
    hMassVsPt[icentbin]->Draw(drawOpt.Data()); pdgmass->Draw();
    cwi->cd();
    hWidthVsPt[icentbin]->Draw(drawOpt.Data()); pdgwidth->Draw();
    cre->cd();
    hSigmaVsPt[icentbin]->Draw(drawOpt.Data());  
    cc2->cd();
    hChi2VsPt[icentbin]->Draw(drawOpt.Data());	
    cSoB->cd();
    hSoverBVsPt[icentbin]->Draw(drawOpt.Data());
    csign->cd();
    hSignificanceVsPt[icentbin]->Draw(drawOpt.Data());
    //plot summary
    csummary->cd(1); gPad->SetLogy();
    hRawYieldVsPt[icentbin]->Draw(drawOpt.Data());
    hRawYieldBCVsPt[icentbin]->Draw(drawOpt.Data());
    csummary->cd(2);
    hTailFractionBCVsPt[icentbin]->Draw(drawOpt.Data());
    csummary->cd(3);
    hChi2VsPt[icentbin]->Draw(drawOpt.Data());
    csummary->cd(4);
    hMassVsPt[icentbin]->Draw(drawOpt.Data());  pdgmass->Draw();
    csummary->cd(5);
    hWidthVsPt[icentbin]->Draw(drawOpt.Data()); pdgwidth->Draw();
    csummary->cd(6);
    hSigmaVsPt[icentbin]->Draw(drawOpt.Data()); 
    
    fout->cd();
    hRawYieldVsPt[icentbin]->Write(Form("hRawYieldVsPt_%i",icentbin));
    hRawYieldBCVsPt[icentbin]->Write(Form("hRawYieldBCVsPt_%i",icentbin));
    hTailFractionBCVsPt[icentbin]->Write(Form("hTailFractionBCVsPt_%i",icentbin));
    hMassVsPt[icentbin]->Write(Form("hMassVsPt_%i",icentbin));
    hSigmaVsPt[icentbin]->Write(Form("hSigmaVsPt_%i",icentbin));
    hChi2VsPt[icentbin]->Write(Form("hChi2VsPt_%i",icentbin));
    hSoverBVsPt[icentbin]->Write(Form("hSoverBVsPt_%i",icentbin));
    hSignificanceVsPt[icentbin]->Write(Form("hSignificanceVsPt_%i",icentbin));
  }
  
  TString fimgname("result.eps");
  cry->Print(Form("YIELD_%s",fimgname.Data()));
  cbc->Print(Form("YIELDBC_%s",fimgname.Data()));
  ctail->Print(Form("TAILS_%s",fimgname.Data()));
  cm->Print(Form("MASS_%s",fimgname.Data()));
  cre->Print(Form("RES_%s",fimgname.Data()));
  cwi->Print(Form("WIDTH_%s",fimgname.Data()));
  cc2->Print(Form("CHI2_%s",fimgname.Data()));
  cSoB->Print(Form("SoB_%s",fimgname.Data()));
  csign->Print(Form("Signif_%s",fimgname.Data()));

  fout->cd();
  cry->Write();
  cbc->Write();
  ctail->Write();
  cm->Write();
  cre->Write();
  cwi->Write();
  cc2->Write();
  cSoB->Write();
  csign->Write();
  fout->Close();

  
  return;
}
