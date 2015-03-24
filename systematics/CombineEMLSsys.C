void CombineEMLSsys( TString centLabel = "0-100%", TString suffix = "28apr14_0")
{

  TString em_name=Form("finalWsyst_smooth2_EM_%s.root", suffix.Data());
  TString ls_name=Form("finalWsyst_smooth2_LS_%s.root", suffix.Data());
  TString hmaterial_name("material");
  TString htracking_name("tracking");
  TString htrackcuts_name("trackcuts");
  TString hPID_name("PID");
  TString hrange_name("range");
  TString hbgNorm_name("bgNorm");
  TString hfunction_name("function");
  TString hhadrint_name("hadrint");
  TString hsum2_0_name("sum2_0");
  TString hsum2_uncorr_0_name("sum2_uncorr_0");

  TFile * fls = TFile::Open(ls_name.Data(),"read"); if (!fls) { Printf("No LS file found"); return;}
  TFile * fem = TFile::Open(em_name.Data(),"read"); if (!fem) { Printf("No EM file found"); return;}
  
  TH1D *  hmaterial[3]; 
  TH1D *  htracking[3];
  TH1D *  htrackcuts[3];
  TH1D *  hPID[3];
  TH1D *  hrange[3];
  TH1D *  hbgNorm[3];
  TH1D *  hfunction[3];
  TH1D *  hhadrint[3];
  TH1D *  hsum2_0[3];
  TH1D *  hsum2_uncorr_0[3];

  TFile * fin = 0x0;
  for (Int_t j=0;j<2;j++){
    if (j==0) fin=fls;
    else fin=fem;
    if (!fin) return;
    hmaterial[j] = (TH1D*) fin->Get(hmaterial_name.Data()); 
    htracking[j] = (TH1D*) fin->Get(htracking_name.Data());
    htrackcuts[j] = (TH1D*) fin->Get(htrackcuts_name.Data());
    hPID[j] = (TH1D*) fin->Get(hPID_name.Data());
    hrange[j] = (TH1D*) fin->Get(hrange_name.Data());
    hbgNorm[j] = (TH1D*) fin->Get(hbgNorm_name.Data());
    hfunction[j] = (TH1D*) fin->Get(hfunction_name.Data());
    hhadrint[j] = (TH1D*) fin->Get(hhadrint_name.Data());
    hsum2_0[j] = (TH1D*) fin->Get(hsum2_0_name.Data());
    hsum2_uncorr_0[j] = (TH1D*) fin->Get(hsum2_uncorr_0_name.Data());
    
  }

  gROOT->LoadMacro("$HOME/alice/macro/combinePlotsFromBinN.C");
  Int_t combineFirstNbins = 2;
  hmaterial[2]  = ( TH1D*) combinePlotsFromBinN(hmaterial[0], hmaterial[1], combineFirstNbins); 
  htracking[2]  = ( TH1D*) combinePlotsFromBinN(htracking[0], htracking[1], combineFirstNbins);
  htrackcuts[2] = ( TH1D*) combinePlotsFromBinN(htrackcuts[0], htrackcuts[1], combineFirstNbins);
  hPID[2]       = ( TH1D*) combinePlotsFromBinN(hPID[0], hPID[1], combineFirstNbins);
  hrange[2]     = ( TH1D*) combinePlotsFromBinN(hrange[0],hrange[1], combineFirstNbins);
  hbgNorm[2]    = ( TH1D*) combinePlotsFromBinN(hbgNorm[0], hbgNorm[1], combineFirstNbins);
  hfunction[2]  = ( TH1D*) combinePlotsFromBinN(hfunction[0], hfunction[1], combineFirstNbins);
  hhadrint[2]   = ( TH1D*) combinePlotsFromBinN(hhadrint[0], hhadrint[1], combineFirstNbins);
  hsum2_0[2]    = ( TH1D*) combinePlotsFromBinN(hsum2_0[0], hsum2_0[1], combineFirstNbins);
  hsum2_uncorr_0[2] = ( TH1D*) combinePlotsFromBinN(hsum2_uncorr_0[0], hsum2_uncorr_0[1], combineFirstNbins);

 //systematic uncertainty plot
  TCanvas *cs=new TCanvas("cs","Systematic uncertainty vs p_{t}", 750,600);
  gStyle->SetOptTitle(0);
  cs->cd();
  hsum2_0[2]->SetTitle("Total systematic uncertainty");
  hsum2_0[2]->GetYaxis()->SetTitle("relative uncert. (%)");
  hsum2_0[2]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hsum2_0[2]->GetYaxis()->SetRangeUser(0.1, 20);
  hsum2_0[2]->GetXaxis()->SetRangeUser(0.0, 14.9);
  
  hsum2_0[2]->Draw("hist");
  hsum2_uncorr_0[2]->Draw("hist same");
  htracking[2]->Draw("hist same");
  hrange[2]->Draw("hist same");
  htrackcuts[2]->Draw("hist same");
  hfunction[2]->Draw("hist same");
  hmaterial[2]->Draw("hist same");
  hPID[2]->Draw("hist same");
  hhadrint[2]->Draw("hist same");
  hbgNorm[2]->Draw("hist same");

  TLegend * autolegry = (TLegend*)gPad->BuildLegend(0.25,0.65,0.88,0.88, Form("V0A multiplicity class %s",centLabel.Data()));
  autolegry->SetFillColor(kWhite);
  autolegry->SetLineColor(kWhite);
  autolegry->SetTextFont(42);
  autolegry->SetNColumns(2); 
  cs->SaveAs(Form("summaryCombinedEMLSsys_%s.png", suffix.Data()));
  cs->SaveAs(Form("summaryCombinedEMLSsys_%s.C", suffix.Data()));

  TFile * fout = new TFile(Form("combinedEMLSsystematicsPlots_%s.root", suffix.Data()),"recreate");
  fout->cd();
  hmaterial[2]->Write();
  htracking[2]->Write();
  htrackcuts[2]->Write();
  hPID[2]->Write();
  hrange[2]->Write();
  hbgNorm[2]->Write();
  hfunction[2]->Write();
  hhadrint[2]->Write();
  hsum2_0[2]->Write();
  hsum2_uncorr_0[2]->Write();
  cs->Write();
  fout->Close(); 
  
  Printf("::::::::::::::: Average hmaterial uncert: %4.2f", GetAverageUncert(hmaterial[2]));
  Printf("::::::::::::::: Average htracking uncert: %4.2f", GetAverageUncert(htracking[2]));
  Printf("::::::::::::::: Average htrackcuts uncert: %4.2f", GetAverageUncert(htrackcuts[2]));
  Printf("::::::::::::::: Average hPID uncert: %4.2f", GetAverageUncert(hPID[2]));
  Printf("::::::::::::::: Average hrange uncert: %4.2f", GetAverageUncert(hrange[2]));
  Printf("::::::::::::::: Average htracking uncert: %4.2f", GetAverageUncert(htracking[2]));
  Printf("::::::::::::::: Average hbgNorm uncert: %4.2f", GetAverageUncert(hbgNorm[2]));
  Printf("::::::::::::::: Average hfunction uncert: %4.2f", GetAverageUncert(hfunction[2]));
  Printf("::::::::::::::: Average hhadrint uncert: %4.2f", GetAverageUncert(hhadrint[2]));
  Printf("::::::::::::::: Average hsum2_0 uncert: %4.2f", GetAverageUncert(hsum2_0[2]));
  Printf("::::::::::::::: Average sum2_uncorr_0 uncert: %4.2f", GetAverageUncert(hsum2_uncorr_0[2]));

  return;
}


Float_t GetAverageUncert(TH1D* h)
{
  TH1F * tmp = new TH1F("tmp","tmp", 100., 0., 100.);
  for (Int_t i=1; i<h->GetNbinsX()+1;i++){
    tmp->Fill(h->GetBinContent(i));
  }
  Float_t mean = tmp->GetMean();
  delete tmp;
  return mean;
}


