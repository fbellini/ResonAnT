/***** files with spectra with stat. uncert. only

       TString fileTOF="./pAexpress215-216/TOF_ANA/ana2s/fit_bestRange_fixedW/CORRECTED_best_fit_poly2.root",
       TString fileTPC="./pAexpress215-216/TPC_ANA/ana2s/fit_bestRange_fixedW/CORRECTED_best_fit_poly2.root",

******/

void compareTOFTPCcorrected_pA(TString fileTOF="./TOF_ANA/ana2s/FINAL/finalWsyst_24set13_0.root",
			       TString fileTPC="./TPC_ANA/ana2s/FINAL/finalWsyst_24set13_0.root",
			       TString suffix = "")
{
  // TFile * fout = new TFile(Form("compareTPCTOF_pA%s.root",suffix.Data()),"recreate");
  for (Int_t j=0;j<5;j++){
    fileTOF.ReplaceAll("_0.root", Form("_%i.root",j));
    fileTPC.ReplaceAll("_0.root", Form("_%i.root",j));
    compareTOFTPCcorrected_cent(j,fileTOF.Data(),fileTPC.Data());
    
  }
  return;
}


void compareTOFTPCcorrected_cent(Int_t centBin=0,
				 TString fileTOF="./TOF_ANA/ana2s/FINAL/finalWsyst_24set13_0.root",
				 TString fileTPC="./TPC_ANA/ana2s/FINAL/finalWsyst_24set13_0.root"
				 ){
  fileTOF.ReplaceAll("_0.root", Form("_%i.root",centBin));
  fileTPC.ReplaceAll("_0.root", Form("_%i.root",centBin));
  
  TFile* ftpc = TFile::Open(fileTPC.Data(),"read");
  TFile* ftof = TFile::Open(fileTOF.Data(),"read");

  //color TOF spectra
  Color_t color[2][6]={kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, kBlack, //tpc
		       kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof
  Int_t marker[2][6]={21, 22, 23, 34, 33, 20, //tpc
		      25, 26, 32, 28, 27, 24}; //tof

  Double_t pt[] = { 0.0, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };
  TString tofhname = Form("hTOFCorrected_%i", centBin); 
  TString tpchname = Form("hTPCCorrected_%i",centBin);
  
  TH1F* htof_stat = (TH1F*)ftof->Get(tofhname.Data());
  TH1F* htpc_stat =(TH1F*)ftpc->Get(tpchname.Data());
  TH1F* htof_syst = (TH1F*)ftof->Get(Form("%s_syst",tofhname.Data()));
  TH1F* htpc_syst =(TH1F*)ftpc->Get(Form("%s_syst",tpchname.Data()));
  TH1F* htof_corr = (TH1F*)ftof->Get(Form("%s_corr",tofhname.Data()));
  TH1F* htpc_corr =(TH1F*)ftpc->Get(Form("%s_corr",tpchname.Data()));
  TH1F* htof_uncorr = (TH1F*)ftof->Get(Form("%s_uncorr",tofhname.Data()));
  TH1F* htpc_uncorr =(TH1F*)ftpc->Get(Form("%s_uncorr",tpchname.Data()));

  if (!htof_stat || !htof_syst  ||  !htof_uncorr || ! htof_corr){
    Printf("::::: Tof data unaccessible");
    return;
  }
  if (!htpc_stat || !htpc_syst  ||  !htpc_uncorr || ! htpc_corr){
    Printf("::::: Tpc data unaccessible");
    return;
  }

  if (!(htof_stat->GetXaxis()->GetNbins()==htpc_stat->GetXaxis()->GetNbins())) 
    Printf("::::: Warning: TPC and TOF have different number of bins");
 
  TAxis * ptAxis = htof_stat->GetXaxis();
  Int_t nPtBins = ptAxis->GetNbins();
  
  //compare with stat uncert.
  TH1D * copyTOF = new TH1D("TOF","TOF",nPtBins, ptAxis->GetXbins()->GetArray());
  copyTOF->Sumw2();
  TH1D * copyTPC = new TH1D("TPC","TPC",nPtBins, ptAxis->GetXbins()->GetArray());
  copyTPC->Sumw2();

  //compare with syst uncert.
  TH1D * copyTOF_uncorr = new TH1D("TOF_uncorr","TOF_uncorr",nPtBins, ptAxis->GetXbins()->GetArray());
  copyTOF_uncorr->Sumw2();
  TH1D * copyTPC_uncorr = new TH1D("TPC_uncorr","TPC_uncorr",nPtBins, ptAxis->GetXbins()->GetArray());
  copyTPC_uncorr->Sumw2();

//copy content in common range
  for (Int_t ibin=0; ibin< nPtBins+1; ibin++){
    if (htof_stat->GetBinContent(ibin)==0) continue;
    if (htpc_stat->GetBinContent(ibin)==0) continue;
    copyTOF->SetBinContent(ibin, htof_stat->GetBinContent(ibin));
    copyTPC->SetBinContent(ibin, htpc_stat->GetBinContent(ibin));
    copyTOF->SetBinError(ibin, htof_stat->GetBinError(ibin));
    copyTPC->SetBinError(ibin, htpc_stat->GetBinError(ibin));
    
    copyTOF_uncorr->SetBinContent(ibin, htof_uncorr->GetBinContent(ibin));
    copyTPC_uncorr->SetBinContent(ibin, htpc_uncorr->GetBinContent(ibin));
    copyTOF_uncorr->SetBinError(ibin, htof_uncorr->GetBinError(ibin));
    copyTPC_uncorr->SetBinError(ibin, htpc_uncorr->GetBinError(ibin));
  } 

  TString centLabel=Form(" (%i-%i%%)",centBin*20,(centBin+1)*20);
  TH1D * ratioS = new TH1D(Form("ratioS_%i",centBin),Form("stat. uncert. sumw2 %s; p_{T}(GeV/#it{c}); ratio TOF/TPC",centLabel.Data()), nPtBins, ptAxis->GetXbins()->GetArray());
  TH1D * ratioB = new TH1D(Form("ratioB_%i",centBin),Form("stat. uncert binomial %s; p_{T}(GeV/#it{c}); ratio TOF/TPC",centLabel.Data()), nPtBins, ptAxis->GetXbins()->GetArray());
  TH1D * ratio_uncorr = new TH1D(Form("ratio_uncorr_%i",centBin),Form("syst. uncert. sumw2 %s; p_{T}(GeV/#it{c}); ratio TOF/TPC",centLabel.Data()), nPtBins, ptAxis->GetXbins()->GetArray());
  
  //sum in quadrature of the weights = sumw2
  ratioS->Divide(copyTOF, copyTPC, 1.0, 1.0); //summed like uncorr errors
  
  //binomial errors
  ratioB->Divide(copyTOF, copyTPC, 1.0, 1.0, "B"); //summed like uncorr errors
  
  //uncorrelated syst
  ratio_uncorr->Divide(copyTOF_uncorr, copyTPC_uncorr, 1.0, 1.0); //summed like uncorr errors
  
  ratioS->SetLineColor(kBlack+1);
  ratioS->SetMarkerColor(kBlack+1);
  ratioS->SetMarkerStyle(24);

  ratioB->SetLineColor(kRed);
  ratioB->SetMarkerColor(kRed);
  ratioB->SetMarkerStyle(20);
  
  ratio_uncorr->SetFillColor(kGray+1);
  ratio_uncorr->SetFillStyle(3001);
  ratio_uncorr->SetLineColor(kGray+1);
  ratio_uncorr->SetMarkerColor(kRed+2);
  ratio_uncorr->SetMarkerStyle(24);
  
  TLine * ln[3];
  for (Int_t i=0;i<3;i++){
    ln[i] = new TLine(0.0, 0.9+0.1*i, 10.0, 0.9+0.1*i);
    ln[i]->SetLineStyle(7);
    ln[i]->SetLineColor(kBlack);
    ln[i]->SetLineWidth(1);
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  ratio_uncorr->GetYaxis()->SetRangeUser(0.0,2.);
  ratioB->GetYaxis()->SetRangeUser(0.0,2.);

  TCanvas * ccomp = new TCanvas(Form("compare_%i",centBin),Form("compare_%i",centBin), 900,500); 
  //ccomp->Divide(1,3);
  // ratioS->Draw();
  // ccomp->cd(2);
  //ccomp->cd(3);
  ccomp->cd(1);
  ratio_uncorr->Draw("E2");
  ratioB->Draw("same");
  leg = (TLegend*) gPad->BuildLegend(0.5,0.7,0.89,0.89,centLabel.Data());
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  ln[0]->Draw("same");
  ln[1]->Draw("same");
  ln[2]->Draw("same");
  gSystem->Exec("mkdir -p compareTOF2TPC");
  ccomp->SaveAs(Form("compareTOF2TPC/compareTOF2TPC_correctedY_%i.png",centBin));
  ccomp->SaveAs(Form("compareTOF2TPC/compareTOF2TPC_correctedY_%i.C",centBin));
  return;
}

Double_t diff2(Double_t x1, Double_t x2)
{
  Double_t dummy=0.0;
  if (x1 >= x2){
    dummy = x1*x1-x2*x2;
  } else {
    dummy = x2*x2-x1*x1;
  }
  return TMath::Sqrt(dummy);	
}
