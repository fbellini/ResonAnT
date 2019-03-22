#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "RooPlot.h"
#include "RooProduct.h"
using namespace RooFit ;

Int_t PDG = 313;
TDatabasePDG *pdg     = TDatabasePDG::Instance();
TParticlePDG *part    = pdg->GetParticle(PDG);
Double_t      pdgMass    = part->Mass(); // const Double_t pdgMass = 0.89594;
Double_t      pdgWidth   = part->Width(); // 0.0487;// const Double_t pdgWidth = 0.0487

void getIntegralCorrectionFactorHisto(TString strategy = "histo",TString infilename="proj/projMC_eff_kstar_train42.root")
{
  Color_t color[4]={kRed,kYellow,kGreen,kBlue};
  TGaxis::SetMaxDigits(3);
  // //open input file
  if (!infilename){
    printf("ERROR: invalid file name provided. \n");
    return;
  }
  TFile* fin=TFile::Open(infilename.Data());
  if (!fin) {
    printf("ERROR: invalid or corrupted input file. Exiting.\n");
    return;
  }

  if (strategy.Contains("fit")) gSystem->Load("libRooFit");
  
  //get bins
  TAxis *ptbins = (TAxis*)fin->Get("ptbins");
  Int_t npt = ptbins->GetNbins();
  TAxis *centbins = (TAxis*)fin->Get("centbins");
  Int_t ncent = centbins->GetNbins();

  const Int_t dimpt = npt+1;
  Double_t pt[dimpt];
  for (Int_t k=0; k<dimpt;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
    //Printf("%5.2f",pt[k]);
  }
  const Int_t dimcent = ncent+1;
  Double_t cent[dimcent]; 
  for (Int_t k=0; k<dimcent;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
    //Printf("%5.2f",cent[k]);
  }

  //output file
  TString fileout(Form("integralCorrectionMC_%s.root",strategy.Data()));
  TFile * fout=new TFile(fileout.Data(),"recreate");
  TH1F*hCorrVsRange[4][10];
  TH1F*hCorrFactorVsRange[4][10];
  for (Int_t icentbin=0; icentbin<ncent;icentbin++){   
    for (Int_t iptbin=0; iptbin<npt;iptbin++){         
      TH1D * h = (TH1D*)fin->Get(Form("TruesPM_ptBin0%i_centBin0%i", iptbin,icentbin));
      TAxis * xa = (TAxis*) h->GetXaxis();
      
      Int_t nx = xa->GetNbins();
      Double_t binw = xa->GetBinWidth(1);
      
      hCorrVsRange[icentbin][iptbin]=new TH1F(Form("hCorr_cent%i_pt%i",icentbin,iptbin),Form("pt bin %i (%2.0f-%2.0f %); M_{low};fraction of signal in M_{low} #leq M #leq infinity", iptbin, cent[icentbin],cent[icentbin+1]), 0.9/binw, 0.6, 1.5);
      hCorrFactorVsRange[icentbin][iptbin]=new TH1F(Form("hCorrFactor_cent%i_pt%i",icentbin,iptbin),Form("pt bin %i (%2.0f-%2.0f %); M_{low};correction factor",iptbin, cent[icentbin],cent[icentbin+1]), 0.9/binw, 0.6, 1.5);
      Double_t totalErr = 0.0;
      Double_t total = h->IntegralAndError(xa->FindBin(0.6),xa->FindBin(1.5), totalErr);

      for (Int_t irange=1; irange<nx; irange++){
	Double_t tailErr, tail;
	if (strategy.Contains("fit")){//estimate tail from integral of fit function
	  Double_t fitParams[8];
	  //	RooPlot *plotframe = (RooPlot*)fitKStar(hnsig, fitParams, hnsig->GetBinLowEdge(irange), 1.5);
	  if (getIntegralFromFitKStar(h, fitParams, h->GetBinLowEdge(irange), 1.5)>0.0){
	    tail = fitParams[4]-fitParams[7];
	    total = fitParams[4];
	  }
	} else {
	  tail = h->IntegralAndError(1,irange,tailErr);
	}
	Double_t frac = tail*100./total;	
	hCorrVsRange[icentbin][iptbin]->SetBinContent(irange,frac);
	//hCorrVsRange[icentbin][iptbin]->SetBinError(irange,tailErr);
	hCorrVsRange[icentbin][iptbin]->SetLineColor(color[icentbin]+5-iptbin);
	hCorrVsRange[icentbin][iptbin]->SetMarkerColor(color[icentbin]+5-iptbin);
	hCorrFactorVsRange[icentbin][iptbin]->SetBinContent(irange, 1.0-frac/100.);
	//hCorrFactorVsRange[icentbin][iptbin]->SetBinError(irange,tailErr);
	hCorrFactorVsRange[icentbin][iptbin]->SetLineColor(color[icentbin]+5-iptbin);
	hCorrFactorVsRange[icentbin][iptbin]->SetMarkerColor(color[icentbin]+5-iptbin);
      }
      fout->cd();
      hCorrVsRange[icentbin][iptbin]->GetYaxis()->SetRangeUser(0.1,100.);
      hCorrVsRange[icentbin][iptbin]->GetYaxis()->SetNdivisions(510);
      hCorrVsRange[icentbin][iptbin]->GetXaxis()->SetRangeUser(0.65,1.15);
      hCorrVsRange[icentbin][iptbin]->Write();
      hCorrFactorVsRange[icentbin][iptbin]->Write();
    }
  }
  return;
}

//------------------------------------------------------------------------------------------------
//RooPlot*fitKStar(TH1D * hnsig, Double_t * fitParams, Double_t xlow = 0.6, Double_t xup=1.5)
Double_t getIntegralFromFitKStar(TH1D * hnsig, Double_t * fitParams, Double_t xlow = 0.6, Double_t xup=1.5)
{
  /*
    performs fit on K* plot histo
    where it is assumed that already (norm+)subtracted distribution is given as input
  */
  //  gROOT->LoadMacro("/Users/bellini/alice/macro/SetGraphicStyle.C");
  //  SetGraphicStyle(0,0,0);
// #ifdef __CINT__
//   gROOT->ProcessLine(".x $ASD/kstar/roofit/RooRelBW.cxx+") ;
//   gROOT->ProcessLine(".x $ASD/kstar/roofit/RooRelBWPS2.cxx+") ;
// #endif
  
  if (!hnsig){
    printf("fitKstar - ERROR: input histogram %s not found in file . returning...\n", hSignalName.Data() );
    return;
  }
  
  Double_t histo_integral = hnsig->Integral();
  Printf("Histo entries = %e", histo_integral);
  RooRealVar x("m","M (GeV/c^{2})", 0.6, 1.5);
  x.setRange("tail", xlow, xup);
  x.setRange("full", 0.6, 1.5);
  
  RooRealVar width("width","width",pdgWidth,0.9*pdgWidth,1.5*pdgWidth);
  width.setConstant(kTRUE);
  RooRealVar mean("mean","mean",pdgMass, 0.5*pdgMass, 1.5*pdgMass);
  
  RooDataHist data("data","data", RooArgList(x), hnsig);
  RooBreitWigner breit("breit","signal Breit-Wigner",x, mean, width);
  RooRealVar nsig("nsig","signal fraction",histo_integral*0.5, 0., histo_integral) ;
  // RooRelBW relbw("relbw","relativistic BW", x, mean, width);
  // RooRelBWPS2 relbwps("relbwps","BW (X) PSfactor", x, pt, mean, width) ;
  RooAddPdf modelBW("modelBW","modelBW",RooArgList(breit),RooArgList(nsig)) ;
  
  //Perform fit 
  RooFitResult * rf = modelBW.fitTo(data, Extended(kTRUE),  SumW2Error(kFALSE), Save());
  RooArgList pars(* modelBW.getParameters(RooArgSet(x) ) );
  
  //make fitted function (need to add normalized term to pdf)
  RooArgSet prodSet(modelBW); //prodSet.add(nsig);
  RooProduct unNormPdf("fitted Function", "fitted Function", prodSet);
  TF1 * f2 = unNormPdf.asTF(RooArgList(x), pars);
  Double_t signalMass=((RooRealVar*) pars.find("mean"))->getVal();
  Double_t  signalMassErr=((RooRealVar*) pars.find("mean"))->getError();
  Double_t signalWidth=((RooRealVar*) pars.find("width"))->getVal();
  Double_t  signalWidthErr=((RooRealVar*) pars.find("width"))->getError();
  Double_t nSignal=((RooRealVar*) pars.find("nsig"))->getVal();
  Double_t nSignalErr=((RooRealVar*) pars.find("nsig"))->getError();
  Double_t chi2=0.0;//=xframe->chiSquare("modelBW", "data", 3);
  
  //get integrals 
  Double_t integ2_full = f2->Integral(0.6, 1.5);
  Double_t integ2 = nSignal*f2->Integral(xlow, xup, 0)/integ2_full;
  Double_t dinteg2 = nSignal*f2->IntegralError(xlow, xup, 0, rf->covarianceMatrix().GetMatrixArray())/integ2_full;
  // LM: for drawing need to clone 
  // new TCanvas(); f2->DrawClone();
  //draw in frame
  // RooPlot * xframe = x.frame(0.6, 1.5);
  // data.plotOn(xframe, Name("data"), LineColor(kBlack), MarkerColor(kBlack), MarkerSize(0.2), LineWidth(2), DrawOption("E1X0")) ;//nota: add Name("data") to get chi2
  // modelBW.plotOn(xframe, Name("modelBW"), LineColor(kRed),LineWidth(2),Range(0.6, 1.5));//nota: add Name("model") to get chi2
  // modelBW.paramOn(xframe,Layout(0.55));
  
  if (fitParams){
    fitParams[0]=signalMass;
    fitParams[1]=signalMassErr;
    fitParams[2]=signalWidth;
    fitParams[3]=signalWidthErr;
    fitParams[4]=nSignal;
    fitParams[5]=nSignalErr;
    fitParams[6]=chi2;
    fitParams[7]=integ2;
    fitParams[8]=dinteg2;
  }
  
  //xframe->Draw("E");
  //return xframe;
  return integ2;
}
