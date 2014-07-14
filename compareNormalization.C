void compareNormalization(Int_t ipt, Int_t icent=0)
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  Char_t files[7][50] = {
    "_sum_EMnorm0.00-0.00_cut2424_train67.root",
    "_sum_EMnorm0.70-0.75_cut2424_train67.root",
    "_sum_EMnorm1.10-1.20_cut2424_train67.root",
    "_sum_EMnorm1.20-1.30_cut2424_train67.root",
    "_sum_EMnorm1.20-1.40_cut2424_train67.root",
    "_sum_EMnorm1.30-1.40_cut2424_train67.root",
    "_sum_EMnorm1.30-1.50_cut2424_train67.root" };

  Char_t labels[7][20] = {
    "EM norm positive",
    "EM norm 0.70-0.75",
    "EM norm 1.10-1.20",
    "EM norm 1.20-1.30",
    "EM norm 1.20-1.40",
    "EM norm 1.30-1.40",
    "EM norm 1.30-1.50"
  };

  Color_t colors[7] = {kRed+2, kMagenta, kBlue, kGreen+1, kAzure+1, kYellow+2, kRed};
  Int_t markers[7] = { 25, 20, 24, 23, 26, 27, 28};
  const Int_t nfiles = 7;
  TFile * fin[7];
  TH1D * hsig;
  TH1D * hls;
  TH1D * hem[7];
  TH1D * hsubem[7];
  TH1D * hsubls;
  TString hsigName = Form("Signal_ptBin%02i_centBin%02i",ipt, icent);
  TString hemName = Form("norm_Mixing_ptBin%02i_centBin%02i", ipt, icent);
  TString hlsName = Form("norm_Like_ptBin%02i_centBin%02i", ipt, icent);
  TString hsubemName = Form("sub_norm_Mixing_ptBin%02i_centBin%02i",ipt, icent);
  TString hsublsName = Form("sub_norm_Like_ptBin%02i_centBin%02i",ipt, icent);

  for (Int_t j=0;j<nfiles;j++){
    fin[j] = TFile::Open(files[j]);
    if (j==0) {
      hsig = (TH1D*) fin[j]->Get(hsigName.Data());
      hsubls = (TH1D*) fin[j]->Get(hsublsName.Data());
      hls = (TH1D*) fin[j]->Get(hlsName.Data());
    }
    hem[j] = (TH1D*) fin[j]->Get(hemName.Data());
    hem[j]->SetLineColor(colors[j]);
    hem[j]->SetMarkerColor(colors[j]);
    hem[j]->SetMarkerStyle(markers[j]);
    hem[j]->SetMarkerSize(0.7);

    hsubem[j] = (TH1D*) fin[j]->Get(hsubemName.Data());
    hsubem[j]->SetLineColor(colors[j]);
    hsubem[j]->SetMarkerColor(colors[j]);
    hsubem[j]->SetMarkerStyle(markers[j]);
    hsubem[j]->SetMarkerSize(0.7);
  }
  hsig->SetLineColor(kBlack);
  hsig->SetMarkerColor(kBlack);
  hsig->SetMarkerStyle(20);
  hsig->SetMarkerSize(0.7);

  hsubls->SetLineColor(kBlack);
  hsubls->SetMarkerColor(kBlack);
  hsubls->SetMarkerStyle(20);
  hsubls->SetMarkerSize(0.7);

  hls->SetLineColor(kGray+2);
  hls->SetFillColor(kGray+1);
  hls->SetFillStyle(3001);
  hls->SetMarkerColor(kGray+2);
  hls->SetMarkerStyle(0);
  hls->SetMarkerSize(0.7);
  
  Double_t emDisplayYrange[2] = {0.1,1e6};
  GetHistArrayYrange(hem, nfiles, emDisplayYrange);
  Double_t subDisplayYrange[2] = {0.1,1e6};
  GetHistArrayYrange(hsubem, nfiles, subDisplayYrange);

  //display distributions before background subtraction
  TCanvas * cem = new TCanvas("cem", "cem", 800, 600);
  TLegend *lem = new TLegend(0.11,0.57,0.43,0.89);
  lem->SetFillColor(kWhite);
  lem->SetLineColor(kWhite);
  hls->GetYaxis()->SetRangeUser((emDisplayYrange[0]<0? emDisplayYrange[0]*1.2 : 0.0),1.2*emDisplayYrange[1]);
  cem->cd();
  hls->Draw("hist");
  lem->AddEntry(hls,"Like-sign not norm","lpf");
  //  h->GetYaxis()->SetRangeUser((ymin<0? ymin*1.2:0.0), 1.2*ymax);
  
  for (Int_t j=0;j<nfiles;j++){
    hem[j]->Draw("same");
    lem->AddEntry(hem[j],labels[j],"lpf");
  }
  hsig->Draw("same");
  lem->AddEntry(hsig,"Unlike-sign","lpf");
  lem->Draw();

  //display distributions after background subtraction
  TCanvas * csub = new TCanvas("csub", "csub", 800, 600);
  TLegend *lsub = new TLegend(0.65,0.63,0.89,0.89);
  lsub->SetFillColor(kWhite);
  lsub->SetLineColor(kWhite);
  hsubls->GetYaxis()->SetRangeUser((subDisplayYrange[0]<0? subDisplayYrange[0]*1.2 : 0.0),1.2*subDisplayYrange[1]);
  csub->cd();
  hsubls->Draw();
  lsub->AddEntry(hsubls,"Like-sign not norm","lpf");
  for (Int_t j=0;j<nfiles;j++){
    hsubem[j]->Draw("same");
    lsub->AddEntry(hsubem[j],labels[j],"lpf");
  }
  lsub->Draw();
  TLine * zeroline = new TLine (0.6, 0.0, 1.5, 0.0); 
  zeroline->SetLineStyle(3); 
  zeroline->SetLineWidth(1);
  zeroline->Draw("same");


  //plot histo of norm factors
  TH1D * hnormem[7]; 
  Double_t pt[] = { 0.0, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };
  Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;   
  Double_t normfactor_em=0.0; 
  Int_t treept = 0, treecent = 0;

  //prepare for display
  TCanvas * cnorm= new TCanvas("cnorm", "cnorm", 800, 600);
  TLegend *lnorm = new TLegend(0.65,0.63,0.89,0.89);
  lnorm->SetFillColor(kWhite);
  lnorm->SetLineColor(kWhite);
  
  //make norm factors plots
  for (Int_t j=0;j<nfiles;j++){
    hnormem[j] = new TH1D(Form("hnorm_%i",j),Form("Norm. factors - %s", labels[j]), npt, pt);
    TTree * tree = (TTree*) fin[j]->Get("ntree");
    if (!tree) Printf("Invalid tree");
    tree->SetBranchAddress("factor_em", &normfactor_em);
    tree->SetBranchAddress("pt_bin", &treept);
    tree->SetBranchAddress("cent_bin", &treecent);
    for (Int_t ie=0;ie<tree->GetEntries(); ie++) {
      tree->GetEntry(ie);
      //tree->Show(ie);
      if (treecent!=icent) continue;
      hnormem[j]->SetBinContent(treept+1, normfactor_em);
      hnormem[j]->SetBinError(treept+1, normfactor_em*1e-4);
    }
    cnorm->cd();
    hnormem[j]->Draw((j==0?"*","same*"));
    lnorm->AddEntry(hnormem[j], labels[j],"lpf");
    hnormem[j]->SetLineColor(colors[j]);
    hnormem[j]->SetMarkerColor(colors[j]);
    hnormem[j]->SetMarkerStyle(markers[j]);
    hnormem[j]->SetMarkerSize(0.7);
    hnormem[j]->GetYaxis()->SetRangeUser(0.0,1.0);
    hnormem[j]->GetYaxis()->SetTitle("EM normalization factor");
    hnormem[j]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    
  }
  cnorm->cd();
  lnorm->Draw();
  return;	
}

void GetHistArrayYrange(TH1D ** histArray, Int_t dim, Double_t * yrange)
{
  Double_t ymin=1e10, ymax=-1e10; 
  if (!histArray) return;
  if (dim<1) return;
  for (Int_t i=0;i<dim;i++){
    if (!histArray[i]) { Printf("invalid histogram %i", i); return;}
    Double_t mintmp = histArray[i]->GetBinContent(histArray[i]->GetMinimumBin());
    Double_t maxtmp = histArray[i]->GetBinContent(histArray[i]->GetMaximumBin());
    if (mintmp<ymin) ymin = mintmp;
    if (maxtmp>ymax) ymax = maxtmp;
  }
  yrange[0] = ymin;
  yrange[1] = ymax;
}
