
#include <stdio.h>
#include "TCanvas.h"
//#include "RooPlot.h"

// retrieve reference values in the database PDG
Int_t PDG = 313;
TDatabasePDG *pdg     = TDatabasePDG::Instance();
TParticlePDG *part    = pdg->GetParticle(PDG);
Double_t      pdgMass    = part->Mass(); // const Double_t pdgMass = 0.89594;
Double_t      pdgWidth   = part->Width(); // 0.0487;// const Double_t pdgWidth = 0.0487;


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void levyFitCents( Double_t rangefitMin = 0.5,
                   Double_t rangefitMax = 8.0)
{
  for (Int_t ic=0;ic<5;ic++){
    HistogramFit("TPC_CORRECTED_best_fit_poly2.root", 
    		 //"CORRECTED_NOBR__best_fit_poly2noBR.root",
		 Form("hTPCCorrected_%i",ic), 
    		 rangefitMin, 
    		 rangefitMax, 
    		 "$HOME/alice/resonances/kstar_pA5.02TeV/ratio2K/levy", 
    		 ic, 
    		 kTRUE /*BR*/);
    HistogramFit("TOF_CORRECTED_best_fit_poly2.root", 
    		 Form("hTOFCorrected_%i",ic), 
    		 rangefitMin, 
    		 rangefitMax, 
    		 "$HOME/alice/resonances/kstar_pA5.02TeV/ratio2K/levy", 
    		 ic, 
    		 kTRUE /*BR*/);
  }
  return;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
void HistogramFit (TString inputfile = "TPC_CORRECTED_best_fit_poly2.root", 
                   TString histogram = "hTPCCorrected_0",  
                   Double_t rangefitMin = 1.,
                   Double_t rangefitMax = 8.0,
		   TString percorsoHist = "$HOME/alice/resonances/kstar_pA5.02TeV/ratio2K/levy",
		   Int_t ic = 1,
		   Bool_t correctForBR = 0)
{
  gStyle->SetOptStat(00001);
  gStyle->SetOptFit(1);
  gStyle->SetTextFont(42);

  //open file with histogram
  TFile *file = TFile::Open(Form("%s/%s", percorsoHist.Data(),inputfile.Data()));
  TH1D *hist = (TH1D*) file->Get(histogram.Data())->Clone("histo");
  if (!hist) return;
  hist->SetTitle(Form("K* %i-%i",ic*20, (ic+1)*20));

  //patch to exclude bin from 0.3-0.5
  if (hist->GetBinContent(2)>0 && hist->GetBinLowEdge(2)==0.3) {
    hist->SetBinContent(2,0.0);
    hist->SetBinError(2,0.0);
  }
  TH1D *fit2hist = (TH1D*) file->Get(histogram.Data())->Clone("ratio");
  fit2hist->SetTitle("fit/hist");
  fit2hist->Reset("ICES");
  fit2hist->SetLineColor(kBlack);
  fit2hist->SetMarkerColor(kBlack);
  fit2hist->SetMarkerStyle(1);
  fit2hist->SetLineWidth(2);
  fit2hist->GetYaxis()->SetRangeUser(0.7,1.3);
  
  TLine * l11 = new TLine(rangefitMin, 1.1, rangefitMax, 1.1);
  l11->SetLineWidth(1); l11->SetLineStyle(7); l11->SetLineColor(kRed);
  TLine * l09 = new TLine(rangefitMin, 0.9, rangefitMax, 0.9); 
  l09->SetLineWidth(1); l09->SetLineStyle(7); l09->SetLineColor(kRed);
  
  Double_t valueMin = rangefitMin;
  Double_t valueMax = rangefitMax;
  if (ic==4 && valueMax>6.0) valueMax = 6.0; 
  //cout << "Insert minimun value for bin: " << endl;
  //cin >> valueMin;
  //cout << "Insert maximum value for bin: " << endl;
  //cin >> valueMax;
  
  Int_t ibinmin = hist->GetXaxis()->FindBin(valueMin);
  Int_t ibinmax = hist->GetXaxis()->FindBin(valueMax);
  //Debug
  Printf("###### Debug   ibinmin = %i ---->  lowEdge = %f", ibinmin,  hist->GetXaxis()->GetBinLowEdge(ibinmin));
  Printf("###### Debug   ibinmax = %i ---->  upEdge = %f", ibinmax,  hist->GetXaxis()->GetBinLowEdge(ibinmax));
  //se fit va da ibinmin low edge a ibinmax low edge => ok
  //se fit va da ibinmin low edge a ibinmax up edge => ibinmax --> ibinmax-1
  
  //draw
  TCanvas *c1 = new TCanvas(Form("c_%i",ic), "Histogram to fit",700,600);
  TCanvas *c2 = new TCanvas(Form("cr_%i",ic), "Ratio fit to histogram",700,600);
  //c1->Divide(1,2);
  hist->GetXaxis()->SetTitle(" p_{T} (GeV/c)");
  hist->GetYaxis()->SetTitle(" d^{2}N/dydp_{T}");   	
  c1->cd();
  gPad->SetLogy();
  // hist->SetLineColor(kRed);
  // hist->SetMarkerColor(kRed);
  // hist->SetMarkerStyle(22);
  // hist->SetMarkerSize(0.8);
  hist->GetXaxis()->SetRangeUser(0.0,10.0);
  hist->GetYaxis()->SetRangeUser(6e-5,5.);
  hist->Draw();

  //*********************************************************** FIT **********************************************************
  Double_t Levy(Double_t *pt, Double_t *par)
  {
    Double_t ldNdy  = par[0];
    Double_t lTemp  = par[1];
    Double_t lPower = par[2];
    Double_t lMass  = par[3];
	     
    Double_t lBigCoef = ((lPower - 1) * (lPower - 2)) / (lPower * lTemp * (lPower * lTemp + lMass * (lPower - 2)));
    Double_t lInPower = 1 + (TMath::Sqrt(pt[0] * pt[0] + lMass * lMass) - lMass) / (lPower * lTemp);
	     
    return ldNdy * pt[0] * lBigCoef * TMath::Power(lInPower, (-1) * lPower);
  }
  TF1 *fitLevy = new TF1("fitLevy",Levy,0.0,11,4);
  fitLevy->SetLineWidth(2);
  fitLevy->SetLineStyle(7);
  fitLevy->SetLineColor(kBlue+2);
  
  Double_t    mass      = pdgMass;
  Double_t    initYield = 0.0001;
  Double_t    initTemp  = 0.3;
  Double_t    initPower = 2.0;
	
  // first try without starting values for the parameters
  // This defaults to 1 for each param.
  // this results in an ok fit for the polynomial function
  // however the non-linear part (lorenzian) does not
  // respond well.
  fitLevy->SetParName(0,"Yield");
  //fitLevy->SetParLimits(0, 0.0, 10000.0);
  fitLevy->SetParameter(0, initYield);
    
  fitLevy->SetParName(1,"Temp");
  //fitLevy->SetParLimits(1, 0.0, 10000.0);
  fitLevy->SetParameter(1, initTemp);
	
  fitLevy->SetParName(2,"n");
  fitLevy->SetParLimits(2, 0.0, 10000.0);
  fitLevy->SetParameter(2, initPower);
	
  fitLevy->SetParName(3,"Mass");
  //fitLevy->SetParLimits(3, pdgMass*0.5, 1.5*pdgMass);
  fitLevy->FixParameter(3, mass);
	
  if (correctForBR) hist->Scale(1./0.666); //commented if using BR-corrected spectra
  hist->Fit("fitLevy","VM","VM",rangefitMin,rangefitMax);
	
  TF1 *dlevy = new TF1("dlevy",Levy,0.0,11,4);

  //Double_t Integral(Double_t a, Double_t b, const Double_t* params = 0, Double_t epsilon = 1e-12)
  Double_t IntegralHist = hist->Integral(ibinmin,ibinmax,"width"); //extremes are included in integral
  cout << "The integral of the histogram between bin " << ibinmin << " and bin " << ibinmax << " is: " << IntegralHist << endl;
  cout << "The integral of the histogram between x =" << hist->GetBinLowEdge(ibinmin) << " and x =" << hist->GetBinLowEdge(ibinmax) << " is: " << IntegralHist << endl;
  Double_t IntegralFit = fitLevy->Integral(rangefitMin,rangefitMax);
  Double_t IntegralFitExtrapolated = fitLevy->Integral(0.0,10.0);
  cout << "The integral of the fit function between " << rangefitMin << " and " << rangefitMax << " is: " << IntegralFit << endl;
  cout << "The integral of the fit function between 0 and 10.0 is: " << IntegralFitExtrapolated << endl;
  	
  for (Int_t i= ibinmin; i<ibinmax+1; i++){
    Float_t ratio = 0.0;
    if (hist->GetBinContent(i)==0.0) {
      Printf("oooops");
      continue;
    }
    ratio = (fitLevy->Eval(hist->GetBinCenter(i)) ) / hist->GetBinContent(i);
    fit2hist->SetBinContent(i, ratio);
  }
  
  Double_t levyPar[4], levyErr[3];
  fitLevy->GetParameters(levyPar);
  levyErr[0] = fitLevy->GetParError(0);
  levyErr[1] = fitLevy->GetParError(1);
  levyErr[2] = fitLevy->GetParError(2);

  TPaveText *result = new TPaveText(0.45, 0.65, 0.88, 0.88,"NDC"); 
  result->SetLineWidth(0);
  result->SetBorderSize(0);
  result->SetFillColor(kWhite);
  result->AddText(Form("Fit range: %4.2f-%4.2f GeV/c",rangefitMin,rangefitMax));
  result->AddText(Form("dN/dy = %10.4f #pm %10.4f", levyPar[0],levyErr[0]));
  result->AddText(Form("T = %10.4f #pm %10.4f", levyPar[1],levyErr[1]));
  result->AddText(Form("n = %10.4f #pm %10.4f", levyPar[2],levyErr[2]));
  result->AddText(Form("#chi^{2}/NDF = %10.4f / %i", fitLevy->GetChisquare(),fitLevy->GetNDF()));
  result->AddText(Form("<pt> = %10.6f ", dlevy->Moment(1,0.0,8.0, fitLevy->GetParameters())));
  
  c1->cd();
  fitLevy->Draw("same");
  fitLevy->Print();
  result->Draw();
  
  c2->cd();
  fit2hist->SetMarkerStyle(20);  //  c1->cd(2);
  fit2hist->Draw("histp");
  l11->Draw("same");
  l09->Draw("same");       
  
  //save
  TString nomeOutput = Form("cent%i_levy%3.1f-%3.1f", ic,rangefitMin,rangefitMax);
  TString OutputFit = Form("%s_%s_inNOBR.root", (inputfile.Contains("TOF")? "TOF":"TPC") ,nomeOutput.Data());
  TFile *FitOut = TFile::Open(OutputFit.Data(),"RECREATE");
  FitOut->cd();
  hist->Write();
  fitLevy->Write();
  if (fit2hist) fit2hist->Write();
  c1->Write();

  c1->SaveAs(Form("FIT_%s_%s.png", (inputfile.Contains("TOF")? "TOF":"TPC") ,nomeOutput.Data()));
  c2->SaveAs(Form("RATIOFit2Hist_%s_%s.png", (inputfile.Contains("TOF")? "TOF":"TPC") ,nomeOutput.Data()));
  

  return;
}
