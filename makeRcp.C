void makeRcp()
{
  TFile * fout = new TFile("rcp_31may13.root","recreate");
  makeRcp("tof", fout);
  makeRcp("tpc", fout);
  fout->Close();
  return;
}

void makeRcp(TString det = "tof", TFile * fout)
{
  if (!fout) return;
  Double_t nEvts81[4]={3.090159e+06, 3.077256e+06, 3.081290e+06, 3.095443e+06};
  
  //various
  Color_t color[4] = {kRed+1, kOrange+2, kGreen+1, kAzure+3};
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //read Kstar histos
  TH1D* hks[4];
  TH1D* haks[4];
  TH1D* hrawks[4];
  TH1D* hrawaks[4];
  if (det.Contains("tof")){
    TFile * f1 = TFile::Open("/Users/bellini/alice/resonances/myKstar/pwglf_train_out/train81/analysisBWPS3_31may13.root"); //kstar
    for (Int_t ic=0;ic<4;ic++){
      hks[ic] = (TH1D*) f1->Get(Form("corrected_kstar_%i", ic))->Clone();
      haks[ic] = (TH1D*) f1->Get(Form("corrected_antikstar_%i", ic))->Clone(); 
      hrawks[ic] = (TH1D*) f1->Get(Form("raw_ks_cent%i", ic))->Clone(Form("norm_raw_ks_cent%i", ic));
      hrawaks[ic] = (TH1D*) f1->Get(Form("raw_aks_cent%i", ic))->Clone(Form("norm_raw_aks_cent%i", ic)); 
      hrawks[ic]->Scale(1./nEvts81[ic]);
      hrawaks[ic]->Scale(1./nEvts81[ic]);
    }
  } else {
    TFile * f2 = TFile::Open("/Users/bellini/alice/resonances/myKstar/pwglf_train_out/resultTPC_SS/10mar2013/spec_0to20_stat.root");
    hks[0] = (TH1D*)f2->Get("hout");
    TFile * f3 = TFile::Open("/Users/bellini/alice/resonances/myKstar/pwglf_train_out/resultTPC_SS/10mar2013/spec_20to40_stat.root");
    hks[1] = (TH1D*)f3->Get("hout");
    TFile * f4 = TFile::Open("/Users/bellini/alice/resonances/myKstar/pwglf_train_out/resultTPC_SS/10mar2013/spec_40to60_stat.root");
    hks[2] = (TH1D*)f4->Get("hout");
    TFile * f5 = TFile::Open("/Users/bellini/alice/resonances/myKstar/pwglf_train_out/resultTPC_SS/10mar2013/spec_60to80_stat.root");
    hks[3] = (TH1D*)f5->Get("hout");      
  }
  
  //Get ratios for each bin
  gROOT->LoadMacro("$ASD/GetPlotRatio.C");
  TH1D* ratioks[4];
  TH1D* ratioaks[4];
  TH1D* ratiorawks[4];
  TH1D* ratiorawaks[4];
  TCanvas * cus[4];
  
  for (Int_t ic=0;ic<4;ic++){
    cus[ic]= new TCanvas(Form("cus%i",ic),Form("cent. %i",ic), 600,500);
    ratioks[ic] = (TH1D*) GetPlotRatio(hks[ic],hks[0],0, Form("K* cent. %02i-%02i%%",ic*20,(ic+1)*20),"K* cent. 0-20%%","ratio");
    ratioks[ic]->SetNameTitle(Form("ks_%s_ratioC2x_cent%i",det.Data(),ic),Form("%s K*: (%02i-%02i%%)/(0-20%%)",det.Data(), ic*20,(ic+1)*20));
    ratioks[ic]->GetYaxis()->SetRangeUser(0.001,1.2);
    ratioks[ic]->GetYaxis()->SetTitle("ratio");
    ratioks[ic]->SetFillColor(kWhite);
    if (det.Contains("tof")) {
      ratioks[ic]->SetLineColor(color[ic]);
      ratioks[ic]->SetMarkerColor(color[ic]);
    } else {
      ratioks[ic]->SetLineColor(kBlack);
      ratioks[ic]->SetMarkerColor(kBlack);
    }
    ratioks[ic]->SetMarkerStyle(20+ic);
    ratioks[ic]->SetMarkerSize(0.8);
  
    if (det.Contains("tof")) {
      ratioaks[ic] = (TH1D*) GetPlotRatio(haks[ic],haks[0],0, Form("TOF #bar{K*} cent. %02i-%02i%%",ic*20,(ic+1)*20),"#bar{K*} cent. 0-20%%");
      ratioaks[ic]->SetNameTitle(Form("aks_%s_ratioC2x_cent%i",det.Data(),ic),Form("%s #bar{K*}: (%02i-%02i%%)/(0-20%%)",det.Data(),ic*20,(ic+1)*20));
      ratioaks[ic]->GetYaxis()->SetTitle("ratio");
      ratioaks[ic]->GetYaxis()->SetRangeUser(0.001,1.2);
      ratioaks[ic]->SetLineColor(color[ic]);
      ratioaks[ic]->SetFillColor(kWhite);
      ratioaks[ic]->SetMarkerColor(color[ic]);
      ratioaks[ic]->SetMarkerStyle(24+ic);
      ratioaks[ic]->SetMarkerSize(0.8);

      ratiorawks[ic] = (TH1D*) GetPlotRatio(hrawks[ic],hrawks[0],0, Form("raw/N_{evts} K* cent. %02i-%02i%%",ic*20,(ic+1)*20),"raw/N_{evts} K* cent. 0-20%%","ratio");
      ratiorawks[ic]->SetNameTitle(Form("rawks_%s_ratioC2x_cent%i",det.Data(),ic),Form("%s raw/N_{evts} K*: (%02i-%02i%%)/(0-20%%)",det.Data(),ic*20,(ic+1)*20));
      ratiorawks[ic]->GetYaxis()->SetTitle("ratio");
      ratiorawks[ic]->GetYaxis()->SetRangeUser(0.001,1.2);
      ratiorawks[ic]->SetLineColor(color[ic]-2);
      ratiorawks[ic]->SetLineWidth(2);
      ratiorawks[ic]->SetLineStyle(7);
      ratiorawks[ic]->SetFillColor(kWhite);
      ratiorawks[ic]->SetMarkerColor(color[ic]-2);
      ratiorawks[ic]->SetMarkerStyle(20+ic);
      ratiorawks[ic]->SetMarkerSize(0.8);
      
      ratiorawaks[ic] = (TH1D*) GetPlotRatio(hrawaks[ic],hrawaks[0],0, Form("raw/N_{evts} K* cent. %02i-%02i%%",ic*20,(ic+1)*20),"raw/N_{evts} K* cent. 0-20%%","ratio");
      ratiorawaks[ic]->SetNameTitle(Form("rawaks_%s_ratioC2x_cent%i",det.Data(),ic),Form("%s raw/N_{evts} #bar{K*}: (%02i-%02i%%)/(0-20%%)",det.Data(),ic*20,(ic+1)*20));
      ratiorawaks[ic]->GetYaxis()->SetTitle("ratio");
      ratiorawaks[ic]->GetYaxis()->SetRangeUser(0.001,1.2);
      ratiorawaks[ic]->SetLineColor(color[ic]-2);
      ratiorawaks[ic]->SetLineWidth(2);
      ratiorawaks[ic]->SetLineStyle(7);
      ratiorawaks[ic]->SetFillColor(kWhite);
      ratiorawaks[ic]->SetMarkerColor(color[ic]-2);
      ratiorawaks[ic]->SetMarkerStyle(24+ic);
      ratiorawaks[ic]->SetMarkerSize(0.8);
    }
    cus[ic]->cd();
    ratioks[ic]->Draw();
    if (det.Contains("tof")) ratioaks[ic]->Draw("same");
    lus = (TLegend*) cus[ic]->BuildLegend(0.58, 0.7,0.89,0.89,"");
    lus->SetFillColor(kWhite);
    lus->SetLineColor(kWhite);
    lus->SetTextFont(42);
    fout->cd(); 
    ratioks[ic]->Write();
    if (det.Contains("tof")) {
      ratioaks[ic]->Write();
      ratiorawks[ic]->Write();
      ratiorawaks[ic]->Write();
    }
    cus[ic]->Write();    
  }
  return;
}
