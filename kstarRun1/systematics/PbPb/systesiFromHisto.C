enum EQuantity{ kYields,
		  kMass,
		  kWidth};

void systesiFromHisto( TString fPath = "normYields", TString fCentral="cent0_central_spectra.root",  TString fSyst="cent0_syst_spectra.root", Int_t selCent = 0, Int_t quantityID = EQuantity::kYields, Float_t chi2cut=2.0, TString filebins="/Users/bellini/alice/resonances/myKstar/pwglf_train_out/fixedEM_5mix/systematics/sub_EMnorm1.30-1.50_aod049_kstar.root", Bool_t display=1)
{
  enum EFitFunction{ kPOLY2,
		     kPOLY3,
		     kLandau,
		     kPOLY1,
		     kEXP,
		     kVoigt,
		     kRELBW,
		     kBW,
		     kData};
  
  Color_t color[]={kRed, kOrange, kGreen+2, kBlue, kMagenta, kBlack};
  Int_t marker[]={20, 21, 28, 22, 23};

  TString hRawYieldsName("hRawYieldVsPt_");
  TString hMassName("hMassVsPt_");
  TString hWidthName("hWidthVsPt_");

  //get bins
  Int_t npt_axis = 0, ncent_axis=0; 
  TFile *f=TFile::Open(filebins.Data());
  if (!ptbins || !centbins) return;
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  npt_axis = ptbins->GetNbins();  
  TAxis *centbins = (TAxis*)f->Get("centbins");
  ncent_axis = centbins->GetNbins();
  f->Close();
  const Int_t npt = npt_axis+1;
  const Int_t ncent = ncent_axis+1;
  Double_t pt[npt];
  for (Int_t k=0; k<npt;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
  }
  Double_t cent[ncent]; 
  for (Int_t k=0; k<ncent;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
  }
  
  TString centLabel=Form("%2.0f-%2.0f%%",cent[selCent],cent[selCent+1]);
  
  //define output
  TString outname(fSyst.Data());
  outname.Prepend("systUncert/syst_");
  if (quantityID == EQuantity::kMass) outname.ReplaceAll("syst_","syst_mass_");
  TFile * outfile = new TFile(outname.Data(),"recreate"); 
  
  
  //Get spectra from files 
  TFile * finCentral= TFile::Open(Form("%s/%s", fPath.Data(),fCentral.Data()),"read");
  if (!finCentral) {
    Printf("Something wrong in output file 1. Check if they're there!");
    return;
  }
  
  //Get pot name according to quantity to be displayed
  TString histoName;
  if (quantityID == EQuantity::kYields) histoName = hRawYieldsName;
  if (quantityID == EQuantity::kMass) histoName = hMassName;
  if (quantityID == EQuantity::kWidth) histoName = hWidthName;
  
  TH1F* hCentral = (TH1F*) finCentral->Get(histoName.Append(Form("%i",selCent))); 
  TString centralTitle(finCentral->GetName()); centralTitle.ReplaceAll(".root", ""); //centralTitle.ReplaceAll("_", " ");
  hCentral->SetName(Form("central%i",selCent));
  hCentral->SetTitle(centralTitle.Data());
  hCentral->SetMarkerStyle(20);
  hCentral->SetLineWidth(1);

  TFile * finSyst= TFile::Open(Form("%s/%s", fPath.Data(),fSyst.Data()),"read");
  if (!finSyst) {
    Printf("Something wrong in output file 2. Check if they're there!");
    return;
  }
  TH1F* spectrum2 = (TH1F*) finSyst->Get(histoName.Data());
  TString secondTitle(finSyst->GetName()); secondTitle.ReplaceAll(".root", ""); //secondTitle.ReplaceAll("_", " ");
  spectrum2->SetName(Form("secondary%i",selCent));
  spectrum2->SetTitle(secondTitle.Data());
  spectrum2->SetMarkerStyle(25);
  spectrum2->SetLineWidth(1);
  
  TH1F* ratio = (TH1F*) hCentral->Clone("ratioCoS");
  ratio->Divide(spectrum2);
  ratio->SetTitle("central/secondary");
  ratio->SetMarkerStyle(28);
  ratio->SetLineWidth(2);
  //save original plots in file
  outfile->cd();
  hCentral->Write();
  spectrum2->Write();
  ratio->Write();
  
  //define output
  TString titleQvsPt;
  TString titleQvsPt_yaxis;
  if (quantityID==EQuantity::kMass) {
    titleQvsPt = Form("Mass %s", centLabel.Data());
    titleQvsPt_yaxis = Form("M (GeV/c^{2})");
  }
  if (quantityID==EQuantity::kYields) {
    titleQvsPt = Form("Raw yields %s",centLabel.Data());
    titleQvsPt_yaxis = Form("1/N_{ev}*d^{2}N/dp_{t}dy (|y|<0.5)");
  }

  TH1F*frameYieldVsPt=new TH1F("frameYieldVsPt","; p_{t} (GeV/c); ", npt_axis, pt);
  TH1F*frameSystVsPt=new TH1F("frameSystVsPt","; p_{t} (GeV/c); syst. error", npt_axis, pt);
  TH1F* hYieldVsPt = new TH1F(Form("hYieldVsPt_%i",selCent),"stat. uncert. only; p_{t} (GeV/c);", npt_axis, pt);
  hYieldVsPt->SetTitle(titleQvsPt.Data());
  hYieldVsPt->GetYaxis()->SetTitle(titleQvsPt_yaxis.Data());
  hYieldVsPt->SetMarkerColor(color[selCent]+1);
  hYieldVsPt->SetLineColor(color[selCent]+1);
  hYieldVsPt->SetMarkerStyle(1);
  hYieldVsPt->SetLineWidth(2);
  
  //  hYieldVsPt->SetFillStyle(0);
  // hYieldVsPt->SetLineColor(kBlack/*color[selCent]*/);	       
  // hYieldVsPt->SetMarkerColor(kBlack/*color[selCent]*/);	       
  // hYieldVsPt->SetMarkerStyle(marker[selCent]);	       
  
  TH1F* hYieldVsPtWsyst = new TH1F(Form("hYieldVsPtWsyst_%i",selCent),"#sqrt{#sigma_{stat}^{2}+#sigma_{syst}^{2}}; p_{t} (GeV/c);", npt_axis, pt);
  hYieldVsPtWsyst->GetYaxis()->SetTitle(titleQvsPt_yaxis.Data());
 // hYieldVsPtWsyst->SetLineColor(color[selCent]+2);	       
  hYieldVsPtWsyst->SetMarkerColor(color[selCent]-9);	       
  // hYieldVsPtWsyst->SetMarkerStyle(1);
  hYieldVsPtWsyst->SetFillColor(color[selCent]-9);
  hYieldVsPtWsyst->SetLineColor(color[selCent]-9);
  hYieldVsPtWsyst->SetMarkerStyle(0);
  hYieldVsPtWsyst->SetOption("E2");

  TH1F* hSystVsPt = new TH1F(Form("hSystVsPt_%i",selCent),"; p_{t} (GeV/c); syst. uncert.", npt_axis, pt);
  hSystVsPt->SetLineColor(color[selCent]);	       
  hSystVsPt->SetMarkerColor(color[selCent]);	       
  hSystVsPt->SetMarkerStyle(marker[selCent]);	       
  hSystVsPt->SetLineWidth(1);
  
  TH1F* hSystVsPtPercentageOfCentral = new TH1F(Form("hSystVsPtPercentageOfCentral_%i",selCent),"; p_{t} (GeV/c); syst. uncert. (% of central value)", npt_axis, pt);
  hSystVsPtPercentageOfCentral->SetTitle(Form("centrality %s", centLabel.Data()));
  hSystVsPtPercentageOfCentral->SetLineColor(color[selCent]);	       
  hSystVsPtPercentageOfCentral->SetMarkerColor(color[selCent]);	       
  hSystVsPtPercentageOfCentral->SetMarkerStyle(marker[selCent]);
  hSystVsPtPercentageOfCentral->SetLineWidth(1);

  //read input 
  Double_t yield1, yield2, erry1, erry2;//fitParams1[9],fitParams2[9];
  Int_t centBinID1,ptBinID1,
    centBinID2,ptBinID2;
  Int_t nentries = hCentral->GetNbinsX()+1;
  for (Int_t ientry = 2;ientry<nentries;ientry++){
    yield1=hCentral->GetBinContent(ientry);
    erry1=hCentral->GetBinError(ientry);
    yield2=spectrum2->GetBinContent(ientry);
    erry2=spectrum2->GetBinError(ientry);
    
    Double_t diffy=0.0, delta2=0.0, syst=0.0, nsigma=0.0, sumerr2=0.0; 
    diffy = TMath::Abs(yield1-yield2);
    delta2 = TMath::Abs(erry2*erry2-erry1*erry1);
    sumerr2 = erry2*erry2+erry1*erry1;
    if (delta2>0.0) nsigma = diffy/TMath::Sqrt(delta2);
    //syst= (nsigma-1)*TMath::Sqrt(delta2);
    Printf("Delta=%e , diffy=%e , nsigma=%5.2f", TMath::Sqrt(delta2), diffy, nsigma);
    if (diffy*diffy-sumerr2>0) syst = TMath::Sqrt(diffy*diffy-sumerr2)/2.0;
    else syst = 0.0;
    if (syst<=0) {
      Printf("Error compatible with stat - syst error negligible");
      syst=1e-5;
    }
    //when fill histo add +2 because definition of bins start
      Float_t dpt = ptbins->GetBinWidth(ientry);
      hYieldVsPt->SetBinContent(ientry, yield1/dpt);
      hYieldVsPt->SetBinError(ientry, erry1/dpt);
      hYieldVsPtWsyst->SetBinContent(ientry, yield1/dpt);
      hYieldVsPtWsyst->SetBinError(ientry, TMath::Sqrt(erry1*erry1/(dpt*dpt)+syst*syst));
      hSystVsPt->SetBinContent(ientry, syst/dpt);
      hSystVsPtPercentageOfCentral->SetBinContent(ientry, (syst*100)/yield1);
  }//loop on bins
  
  if (display){
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    
    //superposition of the two versions of the spectra
    TCanvas *cspectra=new TCanvas("cspectra","Raw yield vs p_{T}", 750,600);
    cspectra->cd();
    hCentral->GetYaxis()->SetRangeUser(1e-5, 15.);
    if (quantityID==EQuantity::kMass)     
      hCentral->GetYaxis()->SetRangeUser(0.86, 0.94);
    if (quantityID==EQuantity::kYields) 
      hCentral->GetYaxis()->SetRangeUser(1e-5, 15.); 

    hCentral->GetYaxis()->SetTitleOffset(1.2);
    hCentral->Draw();
    spectrum2->Draw("same");
    if (quantityID==EQuantity::kYields) gPad->SetLogy();
    TLegend * autolegspectra = (TLegend*)gPad->BuildLegend(0.3,0.75,0.88,0.88, Form("Centrality %s",centLabel.Data()));
    autolegspectra->SetFillColor(kWhite);
    autolegspectra->SetLineColor(kWhite);
    autolegspectra->SetTextFont(42);

    TPaveText * labelkind = new TPaveText(0.20,0.91,0.88,0.98,"NDC");
    labelkind->InsertText(secondTitle.Data());

    //central spectrum with errors
    TCanvas *cry=new TCanvas("cry","Raw yield vs p_{T} with systematic errors", 750,600);
    cry->cd(1);
   if (quantityID==EQuantity::kYields) 
     hYieldVsPtWsyst->GetYaxis()->SetRangeUser(1e-5, 15.); 
   if (quantityID==EQuantity::kMass)     
     hYieldVsPtWsyst->GetYaxis()->SetRangeUser(0.86, 0.94);
   hYieldVsPtWsyst->GetYaxis()->SetTitleOffset(1.2);
   hYieldVsPtWsyst->Draw();
   hYieldVsPt->Draw("same");
   labelkind->Draw("same");
   gPad->SetLogy();
   TLegend * autolegry = (TLegend*)gPad->BuildLegend(0.6,0.7,0.88,0.88, Form("Centrality %s",centLabel.Data()));
    autolegry->SetFillColor(kWhite);
    autolegry->SetLineColor(kWhite);
    autolegry->SetTextFont(42);
    
    //systematic uncertainty plot
    TCanvas *cs=new TCanvas("cs","Systematic error vs p_{T}", 750,600);
    // cs->Divide(1,2);
    // cs->cd(1);
    // frameSystVsPt->GetYaxis()->SetRangeUser(0.1, 1e6);
    // gPad->SetLogy();
    // frameSystVsPt->Draw();
    // hSystVsPt->Draw("same");
    // cs->cd(2);
    cs->cd();
    hSystVsPtPercentageOfCentral->SetLineWidth(2);
    hSystVsPtPercentageOfCentral->GetYaxis()->SetRangeUser(0.1, 50);
    hSystVsPtPercentageOfCentral->Draw();
    labelkind->Draw("same");
  }
  outfile->cd();
  if (display) {
    gStyle->SetOptTitle(0);
    cry->Write();
    cspectra->Write();
    cs->Write();
    fSyst.ReplaceAll(".root",".png");
    if (quantityID==EQuantity::kYields) {
      cry->SaveAs(Form("systUncert/YIELDSwSyst_%s",fSyst.Data()));
      cspectra->SaveAs(Form("systUncert/SPECTRA_%s",fSyst.Data()));
      cs->SaveAs(Form("systUncert/SYSTUNC_%s",fSyst.Data()));
    }
    if (quantityID==EQuantity::kMass) {
      cry->SaveAs(Form("systUncert/MASSwSyst_%s",fSyst.Data()));
      cspectra->SaveAs(Form("systUncert/MASSsyst_%s",fSyst.Data()));
      cs->SaveAs(Form("systUncert/SYSTUNCMASS_%s",fSyst.Data()));
    }
  }
  frameYieldVsPt->Write();
  hYieldVsPt->Write();
  frameSystVsPt->Write();
  hYieldVsPtWsyst->Write();
  hSystVsPt->Write();
  hSystVsPtPercentageOfCentral->Write();
  return;  
  
}
