void GetEfficiencyRatio2minb()
{
  Color_t color[3][5]={ kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-6, //tpc
			kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, //tof
			kBlack, kPink+6, kGreen+1, kAzure+1, kBlue+3 }; //combined tpc3s_tof3sveto
  //			kOrange+7, kPink+6, kGreen+1, kAzure+1, kBlue+4 }; //combined tpc3s_tof3sveto
  Int_t marker[3][5]={21, 22, 23, 34, 33, //tpc
		      25, 26, 32, 28, 27, //tof
		      21, 22, 32, 28, 24}; //combined tpc3s_tof3sveto
  
  TFile * fmb = TFile::Open("$HOME/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/eff_train12-20/binB/efficiency_RsnOut_tpc2s_tof3sveto_cent000-100.root");
  TH1D * hmb = (TH1D*) fmb->Get("hEffVsPt");
  hmb->SetLineColor(kBlack);
  hmb->SetMarkerColor(kBlack);
  TH1D * hmb_ks = (TH1D*) fmb->Get("hEffVsPtKstar");
  TH1D * hmb_aks = (TH1D*) fmb->Get("hEffVsPtAntiKstar");
    
  TH1D * hmbRelUncert = (TH1D*) hmb->Clone("hEffRelUncert");
  hmbRelUncert->Reset();
  for (int j=1;j<hmb->GetNbinsX()+1;j++) {
        Float_t binw = (hmb->GetBinError(j)*100.)/hmb->GetBinContent(j);
        hmbRelUncert->SetBinContent(j, binw);
  }
    
  TH1D * heff[5];
  TH1D * hEffRelUncert[5];
  TH1D * heff_ks[5];
  TH1D * heff_aks[5];
  TH1D * heff2mb[5];
  TCanvas * c1 = new TCanvas("c1","c1", 800, 600);
  
  for (Int_t ic=0;ic<5;ic++){
    TFile * fin = TFile::Open(Form("efficiency_RsnOut_tpc2s_tof3sveto_cent%03i-%03i.root", ic*20, (ic+1)*20));
    heff[ic] =  (TH1D*)fin->Get("hEffVsPt");
    heff[ic]->SetName(Form("hEffVsPt_%i",ic));
    heff2mb[ic] = (TH1D*) heff[ic]->Clone(Form("heff2mb_%i",ic));
    heff2mb[ic]->Divide(hmb);
    heff2mb[ic]->GetYaxis()->SetRangeUser(0.7,1.3);
    heff2mb[ic]->SetLineColor(color[2][ic]);
    heff2mb[ic]->SetMarkerColor(color[2][ic]);
    heff2mb[ic]->SetMarkerStyle(marker[2][ic]);
    heff2mb[ic]->GetYaxis()->SetTitle("ratio");
    heff2mb[ic]->SetTitle(Form("Ratio efficiency %i-%i%% / min.bias", ic*20, (ic+1)*20));
    heff2mb[ic]->SetLineWidth(1);
    c1->cd();
    if (ic==0) heff2mb[ic]->Draw();
    else heff2mb[ic]->Draw("same");
      
    hEffRelUncert[ic] =  (TH1D*)fin->Get("hEffVsPt");
    hEffRelUncert[ic]->SetName(Form("hEffRelUncert_%i",ic));
    hEffRelUncert[ic]->Reset();
    for (int j=1;j<heff[ic]->GetNbinsX()+1;j++) {
        Float_t binw = (heff[ic]->GetBinError(j)*100.)/heff[ic]->GetBinContent(j);
        hEffRelUncert[ic]->SetBinContent(j, binw);
    }
    hEffRelUncert[ic]->SetLineColor(color[2][ic]);
    hEffRelUncert[ic]->SetMarkerColor(color[2][ic]);
    hEffRelUncert[ic]->SetMarkerStyle(marker[2][ic]);
    hEffRelUncert[ic]->GetYaxis()->SetTitle("rel. stat. uncert. (%)");
    hEffRelUncert[ic]->SetTitle(Form("Relative stat. uncertainty of efficiency", ic*20, (ic+1)*20));
    hEffRelUncert[ic]->SetLineWidth(2);
    hEffRelUncert[ic]->GetYaxis()->SetRangeUser(0.,25.);
      
  }
    
  c1->cd();
  TLegend * leg = (TLegend*) c1->BuildLegend();
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  gStyle->SetOptTitle(0);
  TFile * fout = new TFile("ratios2mb_tpc2s_tof3sveto.root","recreate");
  fout->cd();
  hmb->Write();
  hmbRelUncert->Write();
  for (Int_t ic=0;ic<5;ic++){
        heff[ic]->Write();
        heff2mb[ic]->Write();
        hEffRelUncert[ic]->Write();
  }
  return;
}
