//for pA spectra - 16/12/2013
//
//isTOF == 0 (TPC only), 1 (TOF only), 2 (combined with tof veto)
//

void MakeNormCorrSpectraPPB(TString spectraFileName = "RAW_best_fit_poly2.root", 
			    Int_t train=89,
			    Int_t tpcNs = 2,
			    Int_t tofNsveto = 3,
			    Int_t isTOF = 2, 
			    TString suffix = "", 
			    Bool_t correctBR = 1,
                Int_t check = 1 //use 0 for normal, use 1 for check with fit of correction, use 2 for fit func in full range
			    )
{
  const TString effHistName = "hEffVsPt";
  Float_t trgAndVtxEff_pPbMinBias = 0.978; //trigger and vertex reco efficiency for pPb 0-100%
  Float_t trgAndVtxEff_80to100 = (100-5.*2.2)/100.;
  Float_t branchingRatio = 0.666;
  Float_t pid_eff_tpc2s = (0.954*0.954);
  Float_t pid_eff_tpc25s = (0.9876*0.9876);
  Float_t pid_eff_tpc3s = (0.999*0.999);
  
  Float_t nEventsCombined89[5] = {1.157537e+07, 1.158668e+07, 1.158418e+07, 1.194203e+07, 1.085574e+07}; //aod139 LHC13b,c train 8+9
  Float_t nEventsCombined89100[5] = {57544008., -1.0, -1.0, -1.0, -1.0}; //aod139 LHC13b,c train 8+9   //0-100%
  
  Color_t color[3][5]={ kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, //tpc
			kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, //tof
			kRed+1, kPink+6, kGreen+1, kAzure+1, kBlue+2 }; //combined tpc3s_tof3sveto
  //			kOrange+7, kPink+6, kGreen+1, kAzure+1, kBlue+4 }; //combined tpc3s_tof3sveto
  Int_t marker[3][5]={21, 22, 23, 34, 33, //tpc
		      25, 26, 32, 28, 27, //tof
		      21, 22, 32, 28, 24}; //combined tpc3s_tof3sveto
  
  //TPC_TOF matching eff
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TCanvas * c1 = new TCanvas("c1","raw norm. spectra",600,700);
  TCanvas * ceff = new TCanvas("eff","efficiency correction",600,500);
  TCanvas * ccorr = new TCanvas("corr","corrected spectra",600,700);
  TString corrspectraFileName = spectraFileName;
  corrspectraFileName.ReplaceAll("RAW",Form("CORRCHECK_%s",(correctBR?"br":"NOBR")));
  corrspectraFileName.ReplaceAll(".root",Form("%s.root",suffix.Data()));
  TFile * fout = new TFile(corrspectraFileName.Data(),"recreate");
  
  //K*+antK* matching eff
  for (Int_t ic = 0; ic<5;ic++){
    TString centLabel = Form("%i-%i%%",ic*20,(ic+1)*20);
    Printf("*************************\n*************************\n Spectra %s \n*************************",centLabel.Data());
    TString effFileName;
    if (check==1) //USE MIN BIAS EFFICIENCY BUT CORRECT FOR MULTI DEP FOR 60-80% and 80-100%
        effFileName = Form("$HOME/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/eff_train12-20/binB/efficiency_RsnOut_tpc%is_tof%isveto_cent%03i-%03i.root", tpcNs, tofNsveto, 0, 100);
    else
        effFileName = Form("$HOME/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/eff_train12-20/multi-binB/efficiency_RsnOut_tpc%is_tof%isveto_cent%03i-%03i.root", tpcNs, tofNsveto, ic*20, 20*(ic+1));
      
    TString spectraHistName(Form("hRawYieldVsPt_%i",ic));
    TFile *fraw = TFile::Open(spectraFileName.Data());
    if (!fraw) return;
    TH1F * hraw = (TH1F*) fraw->Get(spectraHistName.Data());
    if (!hraw) return;
    
    if (train==89100){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEventsCombined89100[ic]);
      Printf(":::: Normalization by accepted event number from train 8+9  - 0-100%");
    } 
    
    if (train==89){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEventsCombined89[ic]);
      Printf(":::: Normalization by accepted event number from train 8+9");
    } 

    if (train%100 == 0) { //apply correction for vertex and trigger efficiency for p-Pb min bias (0-100%) //nb: this refers to the number of events, thus it multiplies the spectrum 
      hraw->Scale(trgAndVtxEff_pPbMinBias);
      Printf(":::: Correction for trigger and vtx reco efficiency for pPb min bias (0-100%) => scaled by %4.3f", trgAndVtxEff_pPbMinBias);
    }
   
    hraw->SetTitle(Form("Raw yields %s", centLabel.Data()));
    hraw->GetYaxis()->SetTitle("1/N_{evt}* dN/dp_{T} (0<y_{CMS}<0.5)");
    //hraw->GetYaxis()->SetTitleOffset(1.4);
    hraw->GetXaxis()->SetRangeUser(0.0,10.);
    hraw->GetYaxis()->SetRangeUser(1.e-6,2.);
    hraw->SetLineColor(color[isTOF][ic]);
    hraw->SetMarkerColor(color[isTOF][ic]);
    hraw->SetMarkerStyle(marker[isTOF][ic]);
    hraw->SetLineWidth(1);
    
    //get efficiency
    TFile *feff = TFile::Open(effFileName.Data());
    if (!feff) {
      Printf("TROUBLES WITH EFF FILE!!!"); 
      return;
    }
    TH1F * heff = (TH1F*) feff->Get(effHistName.Data());
    if (!heff) return;
    heff->SetTitle(Form("Efficiency (%s)", centLabel.Data()));
    heff->SetName(Form("hCombEff_%i", ic));
    heff->GetXaxis()->SetRangeUser(0.0,16.0);
    heff->GetYaxis()->SetRangeUser(0.0, 1.0);
    heff->GetYaxis()->SetTitle("efficiency");
    heff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    heff->SetLineColor(color[isTOF][ic]);
    heff->SetMarkerColor(color[isTOF][ic]);
    heff->SetMarkerStyle(marker[isTOF][ic]);
    heff->SetLineWidth(1);
    
    //correct for efficiency
    TH1F * hcorr = (TH1F*) hraw->Clone(Form("hCorrected_%i",ic));
      
    if (check==1 && ic>2) {
          Printf(":::: Min bias efficiency used - now correcting for multiplicity dependence of efficiency");
          TF1 * func2 = new TF1(Form("func2_%i",ic),"[0]+exp(-x+[1])",0.0, 15.);
          func2->SetLineColor(kBlue);
          if (ic == 3) {
              //Par [0] = 1.005 +/- 0.003, Par [1] = -3.638 +/- 0.242
              func2->SetParameter(0, 1.005); func2->SetParError(0, 0.003);
              func2->SetParameter(1, -3.638); func2->SetParError(1, 0.242);
          }
          if ( ic == 4) {
              //Par [0] = 1.014 +/- 0.005, Par [1] = -3.072 +/- 0.218
              func2->SetParameter(0, 1.014); func2->SetParError(0, 0.005);
              func2->SetParameter(1, -3.072); func2->SetParError(1, 0.218);
          }
        TCanvas * ccorreff = new TCanvas("ccorreff","ccorreff", 800, 600);
        ccorreff->cd();
        func2->Draw();
        if (heff->Multiply(func2))  Printf("Correction function applied to mb eff: f(x) = [0]+exp(-x+[1])");
        fout->cd();
        func2->Write();
    }
    if (hcorr->Divide(heff))  Printf(":::: Efficiency correction from file: %s", effFileName.Data());
   
      //correct for BR
    if (correctBR) {
      hcorr->Scale(1./(branchingRatio));
      Printf(":::: Yields corrected for BR => scaled by 1/%5.3f", branchingRatio);
    }
    //correct for dy
    hcorr->Scale(2.); 
    Printf(":::: Yields scaled for dy=0.5 to get d^2N/dydp_t");
    
    //correct for .p.le antip.le
    hcorr->Scale(0.5); 
    Printf(":::: Yields scaled for 1/2 to get 1/2 d^2N/dydp_t (p.le/antip.le)");
    
    hcorr->SetTitle(Form("%s", centLabel.Data()));
    hcorr->GetYaxis()->SetTitle(Form("d^{2}#it{N}/(#it{p}_{T}d#it{y}) (GeV/#it{c})^{-1}"));//, correctBR?"* 1/B.R.":""));
    hcorr->GetYaxis()->SetRangeUser(1e-6,2.);
    //hcorr->GetYaxis()->SetTitleOffSet(1.4);
    hcorr->GetXaxis()->SetRangeUser(0.0,14.9);
    hcorr->SetLineColor(color[isTOF][ic]);
    hcorr->SetMarkerColor(color[isTOF][ic]);
    hcorr->SetMarkerStyle(marker[isTOF][ic]);
    hcorr->SetMarkerSize(0.7);
    hcorr->SetLineWidth(1);
    
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
  
  gROOT->LoadMacro("$ASD/AddPaveText.C");
  ceff->cd();
  TLegend * effleg = (TLegend*) gPad->BuildLegend(0.12,0.65,0.55,0.88);
  effleg->SetBorderSize(0);
  effleg->SetFillColor(kWhite);
  
  ccorr->cd();
  TLegend * corrleg = (TLegend*) gPad->BuildLegend(0.55,0.6,0.99,0.88,"V0A multiplicity (Pb side)");
  corrleg->SetBorderSize(0);
  corrleg->SetFillColor(kWhite);
  corrleg->SetFillStyle(0);
  AddPaveText_KStar_pPb("bl");
  AddPaveTextStatOnly();


  c1->cd();
  TLegend * rawleg = (TLegend*) gPad->BuildLegend(0.4,0.7,0.88,0.88);
  rawleg->SetBorderSize(0);
  rawleg->SetFillColor(kWhite);
  AddPaveTextStatOnly();

    fout->cd();
    ceff->Write();
    ccorr->Write();
    c1->Write();
    fout->Close();
  
  Printf("DONE");
  return;
}


