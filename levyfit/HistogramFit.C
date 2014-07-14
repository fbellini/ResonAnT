
#include <stdio.h>
#include "TCanvas.h"
//#include "RooPlot.h"

// retrieve reference values in the database PDG
Int_t PDG = 313;
TDatabasePDG *pdg     = TDatabasePDG::Instance();
TParticlePDG *part    = pdg->GetParticle(PDG);
Double_t      pdgMass    = part->Mass(); // const Double_t pdgMass = 0.89594;
Double_t      pdgWidth   = part->Width(); // 0.0487;// const Double_t pdgWidth = 0.0487;



void HistogramFit (TString nomeOutput = "levy_05set13", //"prova"
		   TString inputfile = "Correction3AnalysisOutput_yCMpos_definitive_BR.root",     //"rawYields_fitEM_LXG_0.65-1.20_EMnorm1.30-1.50_combined0-100_AnalysisTrain_13b.root",
                   TString histogram = "correzioneTOF",   //"hRawYieldsVsPt",
                   TString fitFunction = "Levy",                //"gauss, pol1, pol0, Levy",
                   Double_t rangefitMin = 1.0,
                   Double_t rangefitMax = 10.0,
                   TString percorsoHist = "$HOME/Dropbox/LL/BestFIt_forApproval/LevyFit")

{
    
	//open file with histogram
    TFile *file = TFile::Open(Form("%s/%s", percorsoHist.Data(),inputfile.Data()));

    TH1D *histTofit = (TH1D*) file->Get(histogram.Data())->Clone("histo");
    histTofit->SetTitle("K* 0-100%");
	
    TH1D *fit2hist = (TH1D*) file->Get(histogram.Data())->Clone("ratio");
    fit2hist->SetTitle("fit/hist");
    fit2hist->Reset("ICES");
    fit2hist->SetLineColor(kBlack);
    fit2hist->SetMarkerColor(kBlack);
    fit2hist->SetMarkerStyle(1);
    fit2hist->SetLineWidth(2);
    TLine * l11 = new TLine(rangefitMin, 1.1, rangefitMax, 1.1); l11->SetLineWidth(1); l11->SetLineStyle(7); l11->SetLineColor(kRed);
    TLine * l09 = new TLine(rangefitMin, 0.9, rangefitMax, 0.9); l09->SetLineWidth(1); l09->SetLineStyle(7); l09->SetLineColor(kRed);

    TCanvas *c1 = new TCanvas("c1", "Histogram to fit",14,35,700,600);
    c1->Divide(1,2);
    //histTofit->GetXaxis()->SetTitle(" p_{t} (GeV/c)");
    //histTofit->GetYaxis()->SetTitle(" dN/dp_{t} ");
   	
    c1->cd(1);
    gPad->SetLogy();
    histTofit->SetLineColor(kRed);
    histTofit->SetMarkerColor(kRed);
    histTofit->SetMarkerStyle(22);
    histTofit->SetMarkerSize(0.8);
    histTofit->Draw();
    
    Double_t valueMin = rangefitMin;
    Double_t valueMax = rangefitMax;
    //cout << "Insert minimun value for bin: " << endl;
    //cin >> valueMin;
    //cout << "Insert maximum value for bin: " << endl;
    //cin >> valueMax;
    
    //    Int_t nbins = histTofit->GetXaxis()->GetNbins();
    // Double_t binlow = histTofit->GetXaxis()->GetBinLowEdge(1);
    // Double_t binhigh = histTofit->GetXaxis()->GetBinLowEdge(2);
    // Double_t width = binhigh - binlow;    
    //    cout << "Nbins hist1: " << nbins << ". Lowest bin hist1: " << binlow << ". Highest bin hist1: " << binhigh << ". Bin width hist1: " << width <<  endl;
    //Int_t ibinmin = 1 + ( valueMin-(histTofit->GetXaxis()->GetBinLowEdge(1)) )/width;
    //Int_t ibinmax = nbins - ( (histTofit->GetXaxis()->GetBinLowEdge(nbins)) - valueMax )/width;

    //replaces lines from 51 to 57
    Int_t ibinmin = histTofit->GetXaxis()->FindBin(valueMin);
    Int_t ibinmax = histTofit->GetXaxis()->FindBin(valueMax);
    //Debug
    Printf("###### Debug   ibinmin = %i ---->  lowEdge = %f", ibinmin,  histTofit->GetXaxis()->GetBinLowEdge(ibinmin));
    Printf("###### Debug   ibinmax = %i ---->  upEdge = %f", ibinmax,  histTofit->GetXaxis()->GetBinLowEdge(ibinmax));
    //se fit va da ibinmin low edge a ibinmax low edge => ok
    //se fit va da ibinmin low edge a ibinmax up edge => ibinmax --> ibinmax-1
 
    //*********************************************************** FIT **********************************************************
    if(fitFunction.Contains("Levy"))
      {
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
	fitLevy->SetLineColor(kCyan-3);
	
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
	
	//histTofit->Scale(1./0.66); //commented if using BR-corrected spectra
	histTofit->Fit("fitLevy","VM","VM",rangefitMin,rangefitMax);
	
	//Double_t Integral(Double_t a, Double_t b, const Double_t* params = 0, Double_t epsilon = 1e-12)
	Double_t IntegralHist = histTofit->Integral(ibinmin,ibinmax,"width"); //extremes are included in integral
	cout << "The integral of the histogram between bin " << ibinmin << " and bin " << ibinmax << " is: " << IntegralHist << endl;
	Double_t IntegralFit = fitLevy->Integral(rangefitMin,rangefitMax);
	Double_t IntegralFitExtrapolated = fitLevy->Integral(0.0,10.0);
	cout << "The integral of the fit function between " << rangefitMin << " and " << rangefitMax << " is: " << IntegralFit << endl;
	cout << "The integral of the fit function between 0 and 10.0 is: " << IntegralFitExtrapolated << endl;
    
	// gStyle->SetOptStat(000000010);
	// gStyle->SetOptFit(1);    
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	for (Int_t i= ibinmin; i<ibinmax+1; i++){
	  Float_t ratio = fitLevy->Eval(histTofit->GetBinCenter(i)) / histTofit->GetBinContent(i);
	  fit2hist->SetBinContent(i, ratio);
	}
	fit2hist->GetXaxis()->SetRangeUser(rangefitMin, rangefitMax);
	fit2hist->GetYaxis()->SetRangeUser(0.7, 1.3);
	fit2hist->GetYaxis()->SetNdivisions(512);
	c1->cd(2);
	gPad->SetGridy();
	fit2hist->Draw();
	
	// l11->Draw("same");
	// l09->Draw("same");
	//salvo output in un file
	/*TString OutputRawFit = Form("OutputRawFit_%s_Levy.root",nomeOutput.Data());
	
	TFile *RawFitOut = TFile::Open(OutputRawFit.Data(),"RECREATE");
	     RawFitOut->cd();
	     
	     
	     histTofit->Write();
	     hEfficiency->Write();
	     correzione->Write();
	     
	   RawFitOut->Close();*/
       } //fine Levy
       
       else
       {
           TF1 *fit = new TF1(fitFunction.Data(),fitFunction.Data(),0.0,11);
           fit->SetLineWidth(2);
           fit->SetLineColor(kCyan-3);
           histTofit->Fit(fitFunction.Data(),"VM","VM",rangefitMin,rangefitMax);
           Double_t IntegralHist = histTofit->Integral(ibinmin,ibinmax);
           cout << "The integral of the histogram between bin " << ibinmin << " and bin " << ibinmax << " is: " << IntegralHist << endl;
           Double_t IntegralFit = fit->Integral(rangefitMin,rangefitMax);
           cout << "The integral of the fit function between " << rangefitMin << " and " << rangefitMax << " is: " << IntegralFit << endl;
           
           gStyle->SetOptStat(000000010);
           gStyle->SetOptFit(1);

       }
       

    //salvo output in un file
    TString OutputFit = Form("OutputFit_%s_%s.root",nomeOutput.Data(),fitFunction.Data());
    
    TFile *FitOut = TFile::Open(OutputFit.Data(),"RECREATE");
    FitOut->cd();
    histTofit->Write();
    if (fit2hist) fit2hist->Write();
    FitOut->Close();

    return;
}
