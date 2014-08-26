
Float_t fitEfficiency2minbias(Float_t * pars, Float_t * parErrs,
                   Int_t icent=0, Float_t fitXmin = 1.0, Float_t fitXmax = 15.0){
    
    Color_t color[6]={ kRed+1, kPink+6, kGreen+1, kAzure+1, kBlue+2, kBlack}; //combined tpc3s_tof3sveto
    Int_t marker[6]= { 21, 22, 32, 28, 24, 20};

    TFile * fin = TFile::Open(Form("/Users/fbellini/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/eff_train12-20/multi-binB/ratios2mb_tpc2s_tof3sveto.root"));
    
    if (!fin) return;
    else Printf("Opening file %s", fin->GetName());
 
    TH1F * heff;
    if (icent<0) heff = (TH1F*) fin->Get("heff2mb");
    else heff = (TH1F*) fin->Get(Form("heff2mb_%i",icent));
    if (!heff) return;
    
    //fitting function
    TF1 * func2 = new TF1("func2","[0]+exp(-x+[1])",0.0, 15.);//[0]+[1]*x+[2]*x^2
    Printf("Fitting function defined f(x) = [0]+exp(-x+[1])");
    func2->SetLineColor(kBlue);

    TFitResultPtr result = heff->Fit(func2, "SEM0RQ", "", fitXmin, fitXmax);
    if (result)
        Printf("::::: Fitting histogram %s in pT range [%3.1f,%3.1f] GeV/c", heff->GetName(), fitXmin, fitXmax);
    else return;
    
    Float_t chi2ndf = func2->GetChisquare()/func2->GetNDF();
    //Float_t pars[2]; Float_t parErrs[2];

    for (int j=0;j<2;j++) {
        pars[j]    = result->Parameter(j);
        parErrs[j] = result->ParError(j);
        Printf("Par %i = %5.3f +/- %5.3f", j, pars[j], parErrs[j]);
    }
    
    Printf("Chi^2/NDF = %4.2f", chi2ndf);
    //result->Print("V");
    TCanvas *c1 = new TCanvas("c1","fit correction",0,0,700,500);
    heff->Draw();
    func2->Draw("SAME");
    c1->Update();
    
    TH1F * hratio = (TH1F*) heff->Clone(Form("hist2fit_%i", icent<0? "mb" : icent));
    hratio->SetLineColor((icent<0? color[5] : color[icent]));
    hratio->SetMarkerColor((icent<0? color[5] : color[icent]));
    hratio->SetMarkerStyle((icent<0? marker[5] : marker[icent]));
    hratio->GetYaxis()->SetTitle("ratio hist./fit");
    hratio->SetTitle(Form("ratio hist./fit - multi bin %i", icent));

    hratio->Divide(func2);
    hratio->GetYaxis()->SetRangeUser(0.7, 1.3);
    TLine * l0 = new TLine(fitXmin, 0.9, fitXmin, 1.1); l0->SetLineColor(kBlue); l0->SetLineStyle(2); l0->SetLineWidth(2);
    TLine * l1 = new TLine(fitXmax, 0.9, fitXmax, 1.1); l1->SetLineColor(kBlue); l1->SetLineStyle(2); l1->SetLineWidth(2);
    TCanvas *c2 = new TCanvas("c2","hist./fit",0,0,700,500);
    hratio->Draw();
    gPad->SetGridy();
    l0->Draw("same"); l1->Draw("same");
    
    c1->Print(Form("fitCorr2mb_%i_%3.1f-%3.1f.png",icent, fitXmin, fitXmax));
    c2->Print(Form("fitCorr2mbRatio_%i_%3.1f-%3.1f.png",icent, fitXmin, fitXmax));
    
    return chi2ndf;
}

void fitEfficiencyCent(Float_t fitXmin = 0.0, Float_t fitXmax = 15.0){
    
    gStyle->SetOptStat(0);
    //gPad->SetTickx(); gPad->SetTicky();
    
    TH1F * h0 = new TH1F("h0","parameter [0]; multiplicity bin index; ", 6, -1., 5.);
    h0->GetYaxis()->SetRangeUser(0.8, 1.2);
    h0->GetYaxis()->SetLabelSize(0.05);
    h0->GetXaxis()->SetLabelSize(0.05); h0->GetXaxis()->SetTitleSize(0.05);
    h0->SetLineWidth(2);
    
    TH1F * h1 = new TH1F("h1","parameter [1]; multiplicity bin index; ", 6, -1., 5.);
    h1->GetYaxis()->SetRangeUser(-10., 10.);
    h1->SetLineWidth(2);
    h1->GetYaxis()->SetLabelSize(0.05);
    h1->GetXaxis()->SetLabelSize(0.05); h1->GetXaxis()->SetTitleSize(0.05);
    h1->SetLineWidth(2);
    
    TH1F * hchi2 = new TH1F("hchi2","chi2/Ndof; multiplicity bin index; ", 6, -1., 5.);
    hchi2->SetLineWidth(2);
    hchi2->GetYaxis()->SetRangeUser(0., 10.);
    hchi2->GetYaxis()->SetLabelSize(0.05);
    hchi2->GetXaxis()->SetLabelSize(0.05); hchi2->GetXaxis()->SetTitleSize(0.05);
    hchi2->SetLineWidth(2);
    
    h0->SetBinContent(1, 1.0); h0->SetBinError(1, 0.0);
    h1->SetBinContent(1, -1e4); h1->SetBinError(1, 0.0);
    for (Int_t icent = 0; icent<5; icent++){
        Float_t pars[2]; Float_t parErrs[2];
        Float_t chi = fitEfficiency2minbias(pars, parErrs, icent, fitXmin, fitXmax);
        int ibin = h0->GetXaxis()->FindBin(icent);
        h0->SetBinContent(ibin, pars[0]); h0->SetBinError(ibin, parErrs[0]);
        h1->SetBinContent(ibin, pars[1]); h1->SetBinError(ibin, parErrs[1]);
        hchi2->SetBinContent(ibin, chi);
    }
    TPaveText * ptrange = new TPaveText(0.15, 0.72, 0.89, 0.88, "NDC");
    ptrange->InsertText("Fit range:");
    ptrange->InsertText(Form("%3.1f < #it{p}_{T} < %3.1f GeV/#it{c}", fitXmin, fitXmax));
    ptrange->SetFillColor(kWhite); ptrange->SetBorderSize(0);
    
    TCanvas *cs = new TCanvas("cs","fit correction",0,0,900,400);
    cs->Divide(3,1);
    cs->cd(1); h0->Draw(); ptrange->Draw("same");
    cs->cd(2); h1->Draw();
    cs->cd(3); hchi2->Draw();
    cs->Print(Form("fitCorrParams_%3.2f-%3.1f.png", fitXmin,fitXmax));
    return;
}
