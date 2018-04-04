//fbellini@cern.ch
//for XeXe spectra - 04/03/2018
//
#include "/Users/fbellini/alice/macros/MakeUp.C"

TPaveText * AddPaveTextXeXe();
TPaveText * AddPaveTextStatOnly();

void MakeNormCorrSpectra(TString spectraFileName = "RAW_ana0221_default.root",                    
			    Int_t tpcNs = 2,
                            Int_t tofNsveto = 3,
                            TString suffix = "",
                            Bool_t correctBR = 1)
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

  Color_t color[1][3] = {kRed+1, kSpring+5, kBlue+1};
  Int_t  marker[1][3] = {20, 21, 33}; 

  const TString effHistName = "hEffVsPt";
  Float_t trgAndVtxEff_pPbMinBias = 1.0; //trigger and vertex reco efficiency for pPb 0-100%
  Float_t trgAndVtxEff_80to100 = 1.0;
  Float_t branchingRatio = 0.489; // BR = 0.489 ± 0.005;
    
  Float_t nEvents_ana0221[3] = {2.866550e+05, 2.870050e+05, 2.872150e+05}; //ana0221
  Float_t nEvents_ana0301[3] = { 2.361130e+05, 2.354870e+05, 2.359310e+05};
  Int_t centBinning[4] = {0, 30, 60, 90}; 
    
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TCanvas * c1 = new TCanvas("c1","raw norm. spectra",600,700);
  TCanvas * ceff = new TCanvas("eff","efficiency correction",600,500);
  TCanvas * ccorr = new TCanvas("corr","corrected spectra",600,700);
  TString corrspectraFileName = spectraFileName;
  corrspectraFileName.ReplaceAll("RAW",Form("CORRECTED_%s",(correctBR?"br":"NOBR")));
  corrspectraFileName.ReplaceAll(".root",Form("%s.root",suffix.Data()));
  TFile * fout = new TFile(corrspectraFileName.Data(),"recreate");
    
  for (Int_t ic = 0; ic<3;ic++){

    TString centLabel = Form("%i-%i%%",centBinning[ic], centBinning[ic+1]);
    Printf("\n\n\n*************************\n Spectra %s \n*************************",centLabel.Data());
    
    TString effFileName = "";
    TString effFilePath = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/sim";
    if (spectraFileName.Contains("ana0221")) {
      effFilePath = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/sim/ana0221mc"; // no trailing /
      effFileName = "eff_C3_tpc2s_tof3sveto.root";
    } else
      if (spectraFileName.Contains("ana0301")) {
	effFilePath = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/sim/ana03011mc"; // no trailing /
	effFileName = "eff_C3_tpc2sPtDep_tof3sveto.root";
      }
    
    TString spectraHistName(Form("hRawYieldVsPt_%i",ic));
    TFile *fraw = TFile::Open(spectraFileName.Data());
    if (!fraw) return;
    TH1F * hraw = (TH1F*) fraw->Get(spectraHistName.Data());
    if (!hraw) return;    
    hraw->SetTitle(Form("Raw yields %s", centLabel.Data()));
    hraw->GetYaxis()->SetTitle("1/N_{evt} d#it{N_{raw}}/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");
    hraw->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hraw->GetXaxis()->SetRangeUser(0.0,10.5);
    hraw->GetYaxis()->SetRangeUser(1.e-5,10.);
    hraw->SetLineColor(color[0][ic]);
    hraw->SetMarkerColor(color[0][ic]);
    hraw->SetMarkerStyle(marker[0][ic]);
    hraw->SetLineWidth(1);
        
    //normalise by event number
    if (spectraFileName.Contains("ana0221")){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEvents_ana0221[ic]);
      Printf(":::: Normalization by accepted event number for ana0221");
    } else if (spectraFileName.Contains("ana0301")){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEvents_ana0301[ic]);
      Printf(":::: Normalization by accepted event number for ana0301");
    }
    
    
    TH1F * hcorr = (TH1F*) hraw->Clone(Form("hCorrected_%i",ic));
    hcorr->GetYaxis()->SetTitle("1/N_{evt} d#it{N}/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");
    hcorr->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hraw->GetXaxis()->SetRangeUser(0.0,10.5);
    hraw->GetYaxis()->SetRangeUser(1.e-7,10.);
    hcorr->SetTitle(Form("%s", centLabel.Data()));
    hcorr->SetLineColor(color[0][ic]);
    hcorr->SetMarkerColor(color[0][ic]);
    hcorr->SetMarkerStyle(marker[0][ic]);
    hcorr->SetMarkerSize(1.);
    hcorr->SetLineWidth(1);
          
    //get efficiency 
    TFile *feff = TFile::Open(Form("%s/%s", effFilePath.Data(), effFileName.Data()));
    if (!feff) {
      Printf("CANNOT FIND EFF FILE!!!");
      return;
    }
    TH1F * heff = (TH1F*) feff->Get(Form("%s%i", effHistName.Data(), ic));
    if (!heff) { Printf("Invalid efficiency histogram name."); return;}
    heff->SetTitle(Form("Efficiency (%s)", centLabel.Data()));
    heff->SetName(Form("hEff_%i", ic));
    heff->GetXaxis()->SetRangeUser(0.0,10.5);
    heff->GetYaxis()->SetRangeUser(0.0, 1.0);
    heff->GetYaxis()->SetTitle("efficiency");
    heff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    heff->SetLineColor(color[0][ic]);
    heff->SetMarkerColor(color[0][ic]);
    heff->SetMarkerStyle(marker[0][ic]);
    heff->SetLineWidth(1);

    //---------------------------------------------
    // CORRECTIONS
    //---------------------------------------------

	   
    //apply correction for vertex and trigger efficiency for p-Pb min bias (0-100%) //nb: this refers to the number of events, thus it multiplies the spectrum
    // if (train%100 == 0) {
    //   hraw->Scale(trgAndVtxEff_pPbMinBias);
    //   Printf(":::: Correction for trigger and vtx reco efficiency for pPb min bias (0-100%) => scaled by %4.3f", trgAndVtxEff_pPbMinBias);
    // }
  
    //correct for efficiency
    if (hcorr->Divide(heff))  Printf(":::: Efficiency correction from file: %s/%s", effFilePath.Data(), effFileName.Data());
        
    //correct for BR
    if (correctBR) {
      hcorr->Scale(1./(branchingRatio));
      Printf(":::: Yields corrected for BR => scaled by 1/%5.3f", branchingRatio);
    }
    
    //correct for dy
    hcorr->Scale(1.);
    Printf(":::: Yields scaled for DeltaY = 1.0 to get d^2N/dydp_t");
                
    //DRAW
    c1->cd();
    gPad->SetLogy();
    if (ic==0)  hraw->Draw();
    else hraw->Draw("same");
        
    ccorr->cd();
    gPad->SetLogy();
    if (ic==0)  hcorr->Draw();
    else hcorr->Draw("same");
        
    ceff->cd();
    if (ic==0)  heff->Draw();
    else heff->Draw("same");
        
    fout->cd();
    heff->Write();
    hraw->Write();
    hcorr->Write();
  }
    
  TPaveText * paveSTAT = new TPaveText(0.25, 0.25, 0.6, 0.30);
  BeautifyPave(paveSTAT, 0.04);
  paveSTAT->AddText("Stat. uncertainty only");
  
  ceff->cd();
  TLegend * effleg = (TLegend*) gPad->BuildLegend(0.6,0.65,0.85,0.8);
  myLegendSetUp(effleg, 0.04);
    
  ccorr->cd();
  TLegend * corrleg = (TLegend*) gPad->BuildLegend(0.6,0.70,0.85,0.85);
  myLegendSetUp(effleg, 0.04);

  TPaveText * pave = new TPaveText(0.6, 0.8, 0.89, 0.9);
  BeautifyPave(pave, 0.04);
  pave->AddText("Xe-Xe, #sqrt{#it{s}_{NN}} = 5.44 TeV");
  c1->cd();
  pave->Draw();
  paveSTAT->Draw();
  
  TLegend * rawleg = (TLegend*) gPad->BuildLegend(0.6,0.65,0.85,0.8);
  myLegendSetUp(rawleg, 0.04);
  paveSTAT->Draw();
  pave->Draw();
  
  fout->cd();
  ceff->Write();
  ccorr->Write();
  c1->Write();
  fout->Close();   
    
  Printf("DONE");
  return;
}

  
