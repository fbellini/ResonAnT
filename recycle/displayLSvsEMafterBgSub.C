//particle id = 0 (sum), 1 (kstar), 2 (antikstar)


void displayLSvsEMafterBgSub(TString filename = "_kstar_EMnorm1.30-1.50_cut1818_train215-216.root", TString pngsuffix = "tof3s",  Bool_t is0100 = kTRUE, Int_t colorcombin = 5)
{
  TGaxis::SetMaxDigits(2);  
  gStyle->SetOptTitle(1);
  Color_t colors[2][7] = { kRed+1, kGreen+2, kAzure+1, kOrange+1, kViolet+6, kCyan+1, kBlack, 
			   kRed-5, kGreen-6, kAzure+3, kOrange+4, kViolet-6, kCyan-1, kRed};
  TLine * zeroline = new TLine (0.7, 0.0, 1.3, 0.0); 
  zeroline->SetLineStyle(3); 
  zeroline->SetLineWidth(1);
  //TOF and TPC standalone analyses
  //  Double_t pt[] = { 0.0, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };
  //define binning in pT  - 300MeV bins - binning A 
  //Double_t pt[] = {0.0, 0.15, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };
  //define binning in pT  - 200MeV bins - binning B
  Double_t pt[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00, 12.0, 15.0};
 //define binning in pT  - 200MeV bins - binning C
  //Double_t pt[] = {0.0, 0.1, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 10.00 };
  //binning D
  //Double_t pt[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00, 12.0, 15.0 };


  Char_t particleName[3][12]={"Sum","KStar","AntiKStar"};
  Int_t particle;
  if (filename.Contains("_antikstar_")) particle = 2;
  else if (filename.Contains("_kstar_")) particle = 1;
  else particle = 0;
  
  TFile * fks = TFile::Open(filename.Data());
  if (!fks) return;
  
  TCanvas * canv[5];
  TCanvas * tmp = new TCanvas("tmp","tmp",800,600);
  Int_t ncent = 5;
  Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;   
    
  for (Int_t c=0;c<ncent;c++){    
    canv[c] = new TCanvas(Form("c%i_%s",c,pngsuffix.Data()),Form("c%i_%s",c,pngsuffix.Data()),1200,1000);
    canv[c]->Divide(6,4);
  
    for (Int_t p=0;p<npt;p++){
      TLegend * cleg = new TLegend(0.65,0.65,0.89,0.89);
      cleg->SetFillColor(kWhite);
      cleg->SetBorderSize(0);
      Double_t ymin=1e10, ymax=-1e10;
      TH1D * hks[2];
      for (Int_t isLS = 0; isLS<2; isLS++){
	hks[isLS] = (TH1D*) fks->Get(Form("sub_norm_%s_ptBin%02i_centBin%02i",(isLS?"Like":"Mixing"),p,c));
	TString titleks(Form("%s",hks[isLS]->GetTitle())); 
	hks[isLS]->SetTitle(Form("%2.1f<p_{T}<%2.1f GeV/c (%i-%i%%)",pt[p],pt[p+1], (is0100?0:c*20), (is0100?100:(c+1)*20)));
	hks[isLS]->GetXaxis()->SetTitleSize(0.05);
	hks[isLS]->GetXaxis()->SetLabelSize(0.06);
	hks[isLS]->GetYaxis()->SetTitleSize(0.05);
	hks[isLS]->GetYaxis()->SetLabelSize(0.06);
	hks[isLS]->SetLineWidth(1);
	hks[isLS]->SetMarkerColor(colors[isLS][colorcombin]);
	hks[isLS]->SetLineColor(colors[isLS][colorcombin]);
	hks[isLS]->SetMarkerStyle((isLS?20:28));
	hks[isLS]->Rebin(2);
	if (p==0) {
	  cleg->AddEntry(hks[isLS],(isLS?"LS":"EM"),"lpf");
	}
	hks[isLS]->GetXaxis()->SetTitle("M_{K#pi} (GeV/c^{2})");
	hks[isLS]->GetYaxis()->SetTitle("");//counts / (0.01 GeV/c^{2})");
	hks[isLS]->GetXaxis()->SetRangeUser(0.7, 1.3);
	canv[c]->cd(p+1);
	hks[isLS]->Draw((isLS?"same":""));
	hks[isLS]->GetXaxis()->SetTitleSize(0.04);
	hks[isLS]->GetXaxis()->SetLabelSize(0.04);
	hks[isLS]->GetYaxis()->SetTitleSize(0.04);
	hks[isLS]->GetYaxis()->SetLabelSize(0.04);
	if (isLS) zeroline->Draw("same");
  	tmp->cd();
	hks[isLS]->Draw((isLS?"same":""));
	if (isLS) zeroline->Draw("same");
	if (hks[isLS]->GetBinContent(hks[isLS]->GetMinimumBin())<ymin) ymin = hks[isLS]->GetBinContent(hks[isLS]->GetMinimumBin());
	if (hks[isLS]->GetBinContent(hks[isLS]->GetMaximumBin())>ymax) ymax = hks[isLS]->GetBinContent(hks[isLS]->GetMaximumBin());
      }
      hks[0]->GetYaxis()->SetRangeUser((ymin<0? ymin*1.2:0.0), 1.2*ymax);
      hks[1]->GetYaxis()->SetRangeUser((ymin<0? ymin*1.2:0.0), 1.2*ymax);
      tmp->SaveAs(Form("./LSvsEMpng_%s_c%i_p%i.png",particleName[particle], c,p));
    }
    canv[c]->SaveAs(Form("./LSvsEM_AfterBgSub_%s_%s_c%i.png", particleName[particle], pngsuffix.Data(),c));
  }

  return;
}
