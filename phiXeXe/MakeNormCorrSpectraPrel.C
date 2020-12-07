//fbellini@cern.ch
//for XeXe spectra - 04/03/2018
//
#include "/Users/fbellini/alice/macros/cosmetics/MakeUp.C"

TPaveText * AddPaveTextXeXe();
TPaveText * AddPaveTextStatOnly();

void MakeNormCorrSpectra(TString spectraFileName = "RAW_fitResult.root",                    
        TString effFilePath =  "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/sim/final",
			 //TString effFilePath =  "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0503ec/simulation/LHC17j7_3",//"/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0423/simulation", 
			 TString pid = "default_LowBdca", //"tpc2sPtDep_tof3sveto5smism",
			 TString binning = "C3",
			 Bool_t reweightEff = 0,
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

  Color_t color[1][5] = {kOrange, kSpring+5, kTeal+5, kBlue+1, kMagenta+3};
  Int_t  marker[1][5] = {20, 21, 33, 34, 45}; 

  const TString effHistName = "hEffVsPt";
  Float_t trgAndVtxEff_pPbMinBias = 1.0; //trigger and vertex reco efficiency for pPb 0-100%
  Float_t trgAndVtxEff_80to100 = 1.0;
  Float_t branchingRatio = 0.489; // BR = 0.489 ± 0.005;
    
  Float_t nEvents_ana0221[3] = { 2.866550e+05, 2.870050e+05, 2.872150e+05}; //ana0221
  Float_t nEvents_ana0301[3] = { 2.361130e+05, 2.354870e+05, 2.359310e+05};
  Float_t nEvents_ana0406[3] = { 4.367050e+05, 4.367330e+05, 4.367250e+05};
  Float_t nEvents_ana0406A3[4] = {1.456080e+05, 2.910970e+05, 4.367330e+05, 4.367250e+05};
  Float_t nEvents_ana0414A[5] = {1.392300e+05, 2.784580e+05, 2.785530e+05, 2.784210e+05, 2.785030e+05};
  Float_t nEvents_ana0414B[5] = {1.462360e+05, 2.924550e+05, 2.924940e+05, 2.924800e+05, 2.924840e+05};
//  Int_t centBinning[5] = {0, 30, 60, 90, 100}; 
  Int_t centBinning[7] = {0, 10, 30, 50, 70, 90, 100}; 
  Float_t nEvents_ana0423[5] = {1.441320e+05, 2.882710e+05, 2.881750e+05, 2.883070e+05, 2.883840e+05};
  Float_t nEvents_ana0503[5] ={ 1.445520e+05, 2.890270e+05, 2.889130e+05, 2.885410e+05, 2.858640e+05};
  if (binning.Contains("A3")) {
    centBinning[1] = 10;
    centBinning[2] = 30;
    centBinning[3] = 60;
    centBinning[4] = 90;
    centBinning[5] = 90;
  }

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TCanvas * c1 = new TCanvas("c1","raw norm. spectra",600,700);
  TCanvas * ceff = new TCanvas("eff","efficiency correction",600,500);
  TCanvas * ccorr = new TCanvas("corr","corrected spectra",600,700);
  TString corrspectraFileName = spectraFileName;
  corrspectraFileName.ReplaceAll("RAW",Form("CORRECTED_%s",(correctBR?"br":"NOBR")));
  if (reweightEff)  corrspectraFileName.Prepend("REWEIGHT_");
  corrspectraFileName.ReplaceAll(".root",Form("%s.root",suffix.Data()));
  TFile * fout = new TFile(corrspectraFileName.Data(),"recreate");
  Int_t maxc = 5; if (binning.Contains("A3")) maxc = 4;
		    
  for (Int_t ic = 0; ic<maxc;ic++){

    TString centLabel = Form("%i-%i%%",centBinning[ic], centBinning[ic+1]);
    Printf("\n\n\n*************************\n Spectra %s \n*************************",centLabel.Data());
    
    TString effFileName = Form("eff_o_%s.root", pid.Data()); 
    if (binning.Contains("A3")) effFileName.ReplaceAll("C3", "A3");
    
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
    if (binning.Contains("C3")) {
      if  (effFilePath.Contains("ana0503")){
        hraw->SetName(Form("hCombinedNorm_%i",ic));
        hraw->Scale(1./nEvents_ana0503[ic]);
        Printf(":::: Normalization by accepted event number for ana0503");
      } else if (effFilePath.Contains("ana0423")){
        hraw->SetName(Form("hCombinedNorm_%i",ic));
        hraw->Scale(1./nEvents_ana0423[ic]);
        Printf(":::: Normalization by accepted event number for ana0423");
      } else if (effFilePath.Contains("ana0414") && (pid.Contains("tpc3sPtDep_tof3sveto5smism") || pid.Contains("tpc2sPtDep_tof4sveto5smism"))){
        hraw->SetName(Form("hCombinedNorm_%i",ic));
        hraw->Scale(1./nEvents_ana0414B[ic]);
        Printf(":::: Normalization by accepted event number for ana0414B - pid %s", pid.Data());
      } else if (effFilePath.Contains("ana0414") && (pid.Contains("tpc2sPtDep_tof3sveto"))){      
        hraw->SetName(Form("hCombinedNorm_%i",ic));
        hraw->Scale(1./nEvents_ana0414A[ic]);
        Printf(":::: Normalization by accepted event number for ana0414A - pid %s", pid.Data());
      }
    } else {
      if (binning.Contains("A3") && (effFilePath.Contains("ana0406"))){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEvents_ana0406A3[ic]);
      Printf(":::: Normalization by accepted event number for ana0406");
    } else if (effFilePath.Contains("ana0221")){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEvents_ana0221[ic]);
      Printf(":::: Normalization by accepted event number for ana0221");
    } else if (effFilePath.Contains("ana0301")){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEvents_ana0301[ic]);
      Printf(":::: Normalization by accepted event number for ana0301");
    } else if (effFilePath.Contains("ana0406")){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEvents_ana0406[ic]);
      Printf(":::: Normalization by accepted event number for ana0406");
      } 
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
    TFile * feff = 0x0;
    if (reweightEff) feff = TFile::Open(Form("%s/reweight_%s", effFilePath.Data(), effFileName.Data()));
    else feff = TFile::Open(Form("%s/%s", effFilePath.Data(), effFileName.Data()));
    if (!feff) {
      Printf("CANNOT FIND EFF FILE!!!");
      return;
    }
    TH1F * heff = (TH1F*) feff->Get(Form("%s%i", effHistName.Data(), ic));
    if (!heff) { Printf("Invalid efficiency histogram name."); return;}
    heff->SetTitle(Form("Efficiency (%s)", centLabel.Data()));
    heff->SetName(Form("hEff_%i", ic));
    heff->GetXaxis()->SetRangeUser(0.0, 10.5);
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
    if (hcorr->Divide(heff))  Printf(":::: Efficiency correction from file: %s/%s", effFilePath.Data(), feff->GetName());
        
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

  
