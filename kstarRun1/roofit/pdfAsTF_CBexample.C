#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
using namespace RooFit ;

// build fit function from RooFit pdf
struct FitFunc { 
   FitFunc(TF1 * f, double binWidth) : _f(f), _binW(binWidth) {}
   double operator() (const double *x, const double *p) { 
      // need to unnormalized function
      double ntot = (p[0]+p[7]+p[8]);
      double fval = _f->EvalPar(x,p)*ntot*_binW;
      return fval;
   }
   TF1 * _f;
   double _binW;
};

void pdfAsTF_CBexample(bool doFit = true){
			       
  gStyle->SetOptTitle(0) ;
  gStyle->SetOptStat(0) ;
  gStyle->SetPalette(1) ;
  gStyle->SetCanvasColor(10) ;
  gStyle->SetFrameFillColor(10) ;
  
  TFile *f = new TFile("Histo.root");
  TH1D *histo = (TH1D*) f->Get("histo"); 
    
  //------------------------------------
  // define x range and frame
  //------------------------------------
  
  RooRealVar x("DimuonMass","Mass",2.0,5.,"GeV/c^{2}") ;
  //RooRealVar x("DimuonMass","Mass",1.5,5.,"GeV/c^{2}") ;
  RooPlot* frame = x.frame(Title("Test Roofit")) ;
  TCanvas *c1 = new TCanvas("c1","c1",20,20,600,600);
  gPad->SetLogy(1);
  
  //------------------------------------
  // import histo 
  //------------------------------------
  
  RooDataSet *datat; RooDataHist *datah;
  datah = new RooDataHist("datah","datah",x,Import(*histo));
  datah->plotOn(frame,Name("datah"),DataError(RooAbsData::SumW2),MarkerColor(kRed),MarkerSize(0.8),XErrorSize(0.));   
  datah->statOn(frame);

  //------------------------------------
  // define Crystal Ball
  //------------------------------------
  RooRealVar cbmean("cb mean","cb_mean",3.096,2.9,3.4);
  RooRealVar cbsigma("cb sigma","cb_sigma",0.08,0.05,0.16);
  RooRealVar cbn("cb n","cb_n",3.6,0.1,100.);
  RooRealVar cbalpha("cb alpha","cb_alpha",1,0.1,10.);
  RooCBShape cb("cb","cb",x,cbmean,cbsigma,cbalpha,cbn) ;
  
  cbalpha.setVal(0.98);  
  cbalpha.setConstant(kTRUE);
  cbn.setVal(5.2); 
  cbn.setConstant(kTRUE);
  //----------------------------------------------------
  // fit expo1 
  //-----------------------------------------------------
   
  RooRealVar alpha1("bck1 alpha","bck1_alpha",-2.282635,-5.,0.) ;
  RooRealVar nbck1("nbck1","nbck1",500,0.,1e10);
  RooExponential exp1("exp1","exp1",x,alpha1) ;
  x.setRange("rangeBck1",1.5,2.7);
  RooExtendPdf ebck1_step1("ebck1_step1","ebck1_step1",exp1,nbck1);
  RooAddPdf bck1_step1("bck1_step1","exp bck1",RooArgList(ebck1_step1));
  
  //-----------------------------------------------------
  // fit expo2 
  // ----------------------------------------------------- 
  
  RooRealVar alpha2("bck2 alpha","bck2_alpha",-0.559789,-10.,0.) ;
  RooRealVar nbck2("nbck2","nbck2",500,0.,1e10);
  RooExponential exp2("exp2","exp2",x,alpha2) ;
  RooExtendPdf ebck2_step1("ebck2_step1","ebck2_step1",exp2,nbck2);
  RooAddPdf bck2_step1("bck2_step1","exp bck2",RooArgList(ebck2_step1));
  x.setRange("rangeBck2",3.5,5.);
  
  //------------------------------------ ----------------
  // define expo1+expo2
  //----------------------------------------------------- 
  
  RooAddPdf bck("bck","exp1+exp2",RooArgList(ebck1_step1,ebck2_step1));  
  
  //------------------------------------
  // fit with signal+bck
  //------------------------------------
  
  RooRealVar nsig("N.JPsi","nsignal",500,0.,1000000);
  RooExtendPdf esig("esig","esig",cb,nsig);
  RooAddPdf sum("sum","bck1+bck2+cb",RooArgList(esig,ebck1_step1,ebck2_step1));


  // convert RooAbsPdf to a TF1
  //------------------------------------
 
  // TF1 want to have addition instead of a pdf 
  RooProduct b1("b1","b1",RooArgList(exp1,nbck1));
  RooProduct b2("b2","b2",RooArgList(exp2,nbck2));
  RooProduct sig("sig","sig",RooArgList(esig,nsig));

  // un-normalized version
  //RooAddition usum("usum","usum", RooArgList(sig,b1,b2));
  //normalized version
  RooAddPdf sum("sum","bck1+bck2+cb",RooArgList(cb,exp1,exp2),RooArgList(nsig,nbck1,nbck2));

  RooArgSet * obs  = sum.getObservables(*datah);
  RooArgSet * par  = sum.getParameters(*datah);
  par->Print("V");
  RooArgList p (*par);

  // take normalized pdf 
  TF1 * fnorm = sum.asTF(*obs,p, *obs );

  Double_t binWidth= histo->GetBinWidth(1);
  FitFunc wf( fnorm, binWidth );
  TF1 * functot = new TF1("ftot", &wf, 2.,5.,p.getSize(), "FitFunc");

  
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1);
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
  

  for (int i = 0; i < p.getSize(); ++i) { 
     RooRealVar & v = (RooRealVar&) p[i];
     functot->SetParLimits(i,v.getMin(), v.getMax() );
     if (v.isConstant() )   functot->FixParameter(i,v.getVal()); 
     v.Print();
  }
  // parameters from roofit fit (variable must be defined in 2.,5 range
  functot->SetParameter(0,2084.13);
   functot->SetParameter(1,-2.34744);
   functot->SetParameter(2,-0.941247);
  // // 3 is constant
   functot->SetParameter(4,3.12931);
   functot->SetParameter(6,0.0715441);
  functot->SetParameter(7,159566.);
  functot->SetParameter(8,13491.7);
  
  histo->Draw("e");
  TF1 * functot0 = (TF1*) functot->Clone();
  //functot0->SetLineColor(kBlue);
  functot0->Draw("SAME");
  double xmin,xmax;
  functot->SetRange(2.,5.);
     
  //functot->FixParameter(functot->GetParNumber("cb n"),5.2); 

  if (!doFit) 
     ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1);


  TFitResultPtr r = histo->Fit(functot,"RSI N");
  
  
  histo->GetXaxis()->SetRangeUser(2.,5.);
  histo->SetMinimum(1.);
  gPad->SetLogy(1);
  histo->Draw("e");
  histo->SetMarkerStyle(20);
  histo->SetMarkerColor(2);
  histo->SetMarkerSize(0.7);
  
  functot->DrawClone("same");

// compute integrals 
TF1 * psifix = sig.asTF(*obs);

// TF1 *psifix = new TF1("psifix",FuncJpsi,0.,5.,9); 
// for(int i=0;i<9;i++) psifix->SetParameter(i,functot->GetParameter(i));
 psifix->SetLineColor(kGreen);
 psifix->DrawClone("same");
 
//Double_t NPsi=psifix->Integral(0.,5.);
 cout << "NPSI = " << functot->GetParameter(0) << " +/- " << functot->GetParError(0) << endl;

  }
