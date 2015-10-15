//for pA spectra - 16/12/2013
//
//isTOF == 0 (TPC only), 1 (TOF only), 2 (combined with tof veto)
//

void MakeNormCorrSpectraPPB(TString spectraFileName = "RAW_best_fit_poly2.root",
                            Int_t train=5253,
                            Int_t tpcNs = 2,
                            Int_t tofNsveto = 3,
                            Int_t isTOF = 2,
                            TString suffix = "",
                            Bool_t correctBR = 1
                            )
{
  const TString effHistName = "hEffVsPt";
  Float_t trgAndVtxEff_pPbMinBias = 0.964; //trigger and vertex reco efficiency for pPb 0-100%
  Float_t trgAndVtxEff_80to100 = (100-5.*2.2)/100.;
  Float_t branchingRatio = 0.666;
  Float_t pid_eff_tpc2s = (0.954*0.954);
  Float_t pid_eff_tpc25s = (0.9876*0.9876);
  Float_t pid_eff_tpc3s = (0.999*0.999);
    
  Float_t nEventsTpc215216[5] = {1.9624771e+07, 1.9589295e+07, 1.956609e+07, 2.0104316e+07, 1.7385471e+07}; //aod LHC13b,c train 215,216
  Float_t nEventsTpc215216100[5] = {9.90721600000000000e+07 ,9.62699430000000000e+07, -1.0, -1.0, -1.0}; //aod LHC13b,c train 215,216
  Float_t nEventsTpc231232[5] = {2.016376e+07, 2.012588e+07, 2.012634e+07, 2.089858e+07, 2.008804e+07}; //aod LHC13b,c train 231,232
  Float_t nEventsTof231232[5], nEventsTof215216[5];
  for (int i=0; i<5;i++) {
    //4.606e-03 = fraction of events from run 195390 wrt the full mb stat (good runs) LHC13b+c
    nEventsTof215216[i] = nEventsTpc215216[i]*(1-4.606e-03);
    nEventsTof231232[i] = nEventsTpc231232[i]*(1-4.606e-03);
    //Printf("cent %i: nTpcEvents = %e     nTofEvents = %e", i,nEventsTpc215216[i],nEventsTof215216[i]);
  }
    
  Float_t nEventsCombined12[5] = {1.5955585e+07, 1.5979663e+07, 1.5965594e+07, 1.6360261e+07, 1.3933877e+07}; //aod LHC13b,c train LF_pPb_AOD 1+2
  Float_t nEventsCombined12100[5] = {7.8194980e+07, -1.,-1.,-1.,-1.};// //aod LHC13b,c train LF_pPb_AOD 1+2 0-100%
  Float_t nEventsCombined67[5] = {1.436007e+07, 1.437908e+07, 1.437454e+07, 1.482758e+07, 1.347172e+07}; //aod139 LHC13b,c train 6+7
  Float_t nEventsCombined67100[5] = {7.03489090000000000e+07, -1.0, -1.0, -1.0, -1.0}; //aod139 LHC13b,c train 6+7   //0-100%
  Float_t nEventsCombined89[5] = {1.157537e+07, 1.158668e+07, 1.158418e+07, 1.194203e+07, 1.085574e+07}; //aod139 LHC13b,c train 8+9
  Float_t nEventsCombined89100[5] = {57544008., -1.0, -1.0, -1.0, -1.0}; //aod139 LHC13b,c train 8+9   //0-100%
    
  Float_t nEventsCombined5253[5] = {2.005112e+07, 2.007708e+07, 2.007164e+07, 2.069941e+07, 1.872713e+07}; //aod LHC13b,c train LF_pPb_AOD 52+53
  Float_t nEventsCombined5253100[5] = {9.962638e+07, -1., -1., -1., -1.}; //aod LHC13b,c train LF_pPb_AOD 52+53



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
  corrspectraFileName.ReplaceAll("RAW",Form("CORRECTED_%s",(correctBR?"br":"NOBR")));
  corrspectraFileName.ReplaceAll(".root",Form("%s.root",suffix.Data()));
  TFile * fout = new TFile(corrspectraFileName.Data(),"recreate");
    
  //K*+antK* matching eff
  for (Int_t ic = 0; ic<5;ic++){
    TString centLabel = Form("%i-%i%%",ic*20,(ic+1)*20);
    Printf("*************************\n*************************\n Spectra %s \n*************************",centLabel.Data());
    TString effFileName;
        
    switch (isTOF) {
    case 0 :
      if (train==215216) effFileName = Form("$HOME/alice/resonances/kstar_pA5.02TeV/MC_CF_pPb82/eff/efficiency_RsnOut_tof2s_centBin0%i.root",ic);
      if (train==231232) effFileName = Form("$HOME/alice/resonances/kstar_pA5.02TeV/MC_CF_pPb99/eff/efficiency_RsnOut_tof25s_centBin0%i.root",ic);
      break;
                
    case 1 :
      if (train==215216) effFileName = Form("$HOME/alice/resonances/kstar_pA5.02TeV/MC_CF_pPb82/eff/efficiency_RsnOut_quality2011_centBin0%i.root",ic);
      // if (train==12) effFileName = Form("$HOME/alice/resonances/kstar_pA5.02TeV/LF_pPb_1-2/LSvsEMpng_binningOld/eff/efficiency_RsnOut_comb_cent000-100.root  ",ic);
      //effFileName = Form("$HOME/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/efficiency_RsnOut_comb_centBin0%i.root",ic);//combined
      break;
                
    case 2:
      if (train==12100) effFileName = Form("$HOME/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/efficiency_RsnOut_comb_cent000-100.root  ");
      if (train==12) effFileName    = Form("$HOME/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/efficiency_RsnOut_comb_cent%03i-%03i.root",20*ic, 20*(ic+1));
      if (train==67100) effFileName = Form("$HOME/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/efficiency_RsnOut_tpc3s_tof3sveto_cent000-100.root");
      if (train==67) effFileName    = Form("$HOME/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/efficiency_RsnOut_tpc3s_tof3sveto_cent%03i-%03i.root",20*ic, 20*(ic+1));
      if (train==89100 || train==5253100) {
	effFileName = Form("$HOME/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/eff_train12-20/binB/efficiency_RsnOut_tpc%is_tof%isveto_cent000-100.root", tpcNs, tofNsveto);
      }
      if (train==89 || train==5253) {
	effFileName = Form("$HOME/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/eff_train12-20/multi-binB/efficiency_RsnOut_tpc%is_tof%isveto_cent%03i-%03i.root", tpcNs, tofNsveto, ic*20, 20*(ic+1));
      }
      break;
                
    default :
      Printf("No efficiency file assigned to unknown analysis strategy");
    }
        
    TString spectraHistName(Form("hRawYieldVsPt_%i",ic));
    //0-100% TString spectraHistName(Form("hRawYieldVsPtTPC"));
        
    TFile *fraw = TFile::Open(spectraFileName.Data());
    if (!fraw) return;
    TH1F * hraw = (TH1F*) fraw->Get(spectraHistName.Data());
    if (!hraw) return;
        
    //normalise by event number
    if (train==215216){
      hraw->SetName(Form("h%sNorm_%i",(isTOF==1)?"TOF":"TPC",ic));
      if (isTOF==1) hraw->Scale(1./nEventsTof215216[ic]);
      else hraw->Scale(1./nEventsTpc215216[ic]);
      Printf(":::: Normalisation by accepted event number from train 215+216");
    }
        
    if (train==231232){
      hraw->SetName(Form("h%sNorm_%i",(isTOF==1)?"TOF":"TPC",ic));
      if (isTOF==1) hraw->Scale(1./nEventsTof231232[ic]);
      else hraw->Scale(1./nEventsTpc231232[ic]);
      Printf(":::: Normalisation by accepted event number from train 231+232");
    }
        
    if (train==12){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEventsCombined12[ic]);
      Printf(":::: Normalization by accepted event number from train 1+2");
    }
    if (train==12100){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEventsCombined12100[ic]);
      Printf(":::: Normalization by accepted event number from train 1+2 - 0-100%");
    }
    if (train==67){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEventsCombined67[ic]);
      Printf(":::: Normalization by accepted event number from train 6+7");
    }
        
    if (train==67100){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEventsCombined67100[ic]);
      Printf(":::: Normalization by accepted event number from train 6+7  - 0-100%");
    }
        
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

    if (train==5253100){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEventsCombined5253100[ic]);
      Printf(":::: Normalization by accepted event number from train 52+53 - 0-100%");
    }
        
    if (train==5253){
      hraw->SetName(Form("hCombinedNorm_%i",ic));
      hraw->Scale(1./nEventsCombined5253[ic]);
      Printf(":::: Normalization by accepted event number from train 52+53");
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
    if (train==12 || train==67){
      heff->SetTitle(Form("Efficiency (%s)", centLabel.Data()));
      heff->SetName(Form("hCombEff_%i", ic));
    } else {
      heff->SetTitle(Form("%s efficiency (%s)", isTOF?"TOF":"TPC",centLabel.Data()));
      heff->SetName(Form("h%sEff_%i", isTOF?"TOF":"TPC",ic));
    }
    heff->GetXaxis()->SetRangeUser(0.0,16.0);
    heff->GetYaxis()->SetRangeUser(0.0, 1.0);
    heff->GetYaxis()->SetTitle("efficiency");
    heff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    heff->SetLineColor(color[isTOF][ic]);
    heff->SetMarkerColor(color[isTOF][ic]);
    heff->SetMarkerStyle(marker[isTOF][ic]);
    heff->SetLineWidth(1);
        
    //if tpc, multiply eff quality by 2sigma PID factor
    if (isTOF==0) {
      if (train==231232) {
	heff->Scale(pid_eff_tpc25s);
	Printf(":::: Scaled TPC efficiency quality by PID gaussian factor (2.5sigma)= %f", pid_eff_tpc25s);
      } else {
	if (train==215216) {
	  heff->Scale(pid_eff_tpc2s);
	  Printf(":::: Scaled TPC efficiency quality by PID gaussian factor (2.5sigma)= %f", pid_eff_tpc2s);
	} else {
	  Printf(":::: Efficiency not scaled for any gaussian factor: ok for combined ");
	}
      }
    }
    //correct for efficiency
    TH1F * hcorr = (TH1F*) hraw->Clone(Form("hCorrected_%i",ic));
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
        
    // if (isRoofit) {
    //   hraw->SetLineColor(kGreen+2);
    //   hraw->SetMarkerColor(kGreen+2);
    // }
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

