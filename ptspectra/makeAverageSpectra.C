Int_t makeAverageSpectra(char * list = "likeOut.lst", Float_t chi2cut=2.5)
{
  TString infile ; 
  Int_t filesCounter=0;
  Int_t bincounter[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  FILE * files = fopen(list, "r") ; 
  
  TH1D * raws = new TH1D();
  TCanvas * c1 = new TCanvas("superp","superp",600,600);
  Bool_t isEM=0;
  while ( infile.Gets(files) ){
    infile.Prepend("$HOME/alice/resonances/myKstar/ESD_pp2011/data/");
    
    TFile * ifile = TFile::Open(infile.Data());
    if (!ifile) break;
    filesCounter++;
    if (infile.Contains("Mixing")) isEM=1;
    TH1D* dummyraw = (TH1D*) ifile->Get("raw");
    TH1D* dummybg = (TH1D*) ifile->Get("bg");
    TH1D* dummychi = (TH1D*) ifile->Get("chi2");
    
    dummyraw->SetLineColor(29+filesCounter);
    dummyraw->SetLineWidth(1);
    dummyraw->SetMarkerColor(29+filesCounter);
  
    //copies only axis ranges
    if (filesCounter==1){
      raws = (TH1D*) dummyraw->Clone();
      raws->Reset("ICES");
      c1->cd();
      gPad->SetLogy();
      dummyraw->SetTitle(Form("Raw spectra by fit range and function variation - %s bg",(isEM? "EM":"LS")));
      dummyraw->GetXaxis()->SetTitle("p_{t} (GeV/c)");
      dummyraw->GetYaxis()->SetTitle("raw dN/dp_{t}");
      dummyraw->Draw();
      
    }else{
      dummyraw->Draw("same");
    }
    c1->SaveAs(Form("superposition_%s.png",(isEM? "EM":"LS")));
    Printf("================== File %i",filesCounter);
    for (Int_t ibin = 1; ibin<dummyraw->GetXaxis()->GetNbins(); ibin++ ){
      Float_t chi2 = dummychi->GetBinContent(ibin);
      
      if (chi2>0 && chi2<chi2cut) {
	bincounter[ibin]++;
	Double_t rawnew = dummyraw->GetBinContent(ibin);
	Double_t rawnewerr = dummyraw->GetBinError(ibin);
	Double_t currentRaw= raws->GetBinContent(ibin);
	Double_t currentRawerr= raws->GetBinError(ibin);
	Printf("bin %i: before = %e -> after = %e", ibin, currentRaw, rawnew+currentRaw );
	currentRaw+=rawnew;
	raws->SetBinContent(ibin,currentRaw);
	raws->SetBinError(ibin,currentRawerr+rawnewerr);
      }//check chi2      
      else {
	Printf("---- Chi2 for bin %i too high = %f : skipping", ibin, chi2);
      }
    }//loop on pt bins
  }
  for (Int_t ibin = 1; ibin<raws->GetXaxis()->GetNbins(); ibin++ ){
    Double_t content = raws->GetBinContent(ibin);
    Double_t contenterr = raws->GetBinError(ibin);
    raws->SetBinContent(ibin, content/bincounter[ibin]);
    raws->SetBinError(ibin, contenterr/bincounter[ibin]);
  }
  raws->SetName("hRawYieldVsPt_0");
  raws->SetLineColor(kBlack);
  raws->SetLineWidth(2);
  raws->SetMarkerStyle(20);
  raws->SetMarkerSize(1.);
  raws->SetMarkerColor(kBlack);
  c1->cd();
  //raws->Draw("same");
  TFile * fout = new TFile(Form("uncorrectedYields_%s.root",(isEM? "EM":"LS")),"recreate");
  raws->Write();
  fout->Close();
  Printf("processed files n. %i",filesCounter);
  Printf("Output saved in %s/%s",gSystem->pwd(),fout->GetName());
  return 0;
}


