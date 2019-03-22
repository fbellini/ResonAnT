void fitEfficiencyCent(Float_t fitXmin = 1.0, Float_t fitXmax = 15.0){
    
    gStyle->SetOptStat(0);
    //gPad->SetTickx(); gPad->SetTicky();
    
    TH1F * h0 = new TH1F("h0","parameter [0]; multiplicity bin index; ", 6, -1., 5.);
    h0->GetYaxis()->SetRangeUser(0.5, 0.7);
    h0->GetYaxis()->SetLabelSize(0.05);
    h0->GetXaxis()->SetLabelSize(0.05); h0->GetXaxis()->SetTitleSize(0.05);
    h0->SetLineWidth(2);
    
    TH1F * h1 = new TH1F("h1","parameter [1]; multiplicity bin index; ", 6, -1., 5.);
    h1->GetYaxis()->SetRangeUser(2., 3.);
    h1->SetLineWidth(2);
    h1->GetYaxis()->SetLabelSize(0.05);
    h1->GetXaxis()->SetLabelSize(0.05); h1->GetXaxis()->SetTitleSize(0.05);
    h1->SetLineWidth(2);
    
    TH1F * h2 = new TH1F("h2","parameter [2]; multiplicity bin index; ", 6, -1., 5.);
    h2->GetYaxis()->SetRangeUser(-2.5, -1.5);
    h2->SetLineWidth(2);
    h2->GetYaxis()->SetLabelSize(0.05);
    h2->GetXaxis()->SetLabelSize(0.05); h2->GetXaxis()->SetTitleSize(0.05);
    h2->SetLineWidth(2);
    
    TH1F * hchi2 = new TH1F("hchi2","chi2/Ndof; multiplicity bin index; #Chi^{2}/N_{dof}", 6, -1., 5.);
    hchi2->SetLineWidth(2);
    hchi2->GetYaxis()->SetRangeUser(0., 10.);
    hchi2->GetYaxis()->SetLabelSize(0.05);
    hchi2->GetXaxis()->SetLabelSize(0.05); hchi2->GetXaxis()->SetTitleSize(0.05);
    hchi2->SetLineWidth(2);
    
    for (Int_t icent = -1; icent<5; icent++){
        Float_t pars[3]; Float_t parErrs[3];
        Float_t chi = fitEfficiency(pars, parErrs, icent, fitXmin, fitXmax);
        int ibin = h0->GetXaxis()->FindBin(icent);
        h0->SetBinContent(ibin, pars[0]); h0->SetBinError(ibin, parErrs[0]);
        h1->SetBinContent(ibin, pars[1]); h1->SetBinError(ibin, parErrs[1]);
        h2->SetBinContent(ibin, pars[2]); h2->SetBinError(ibin, parErrs[2]);
        hchi2->SetBinContent(ibin, chi);
    }
    TPaveText * ptrange = new TPaveText(0.15, 0.72, 0.89, 0.88, "NDC");
    ptrange->InsertText("Fit range:");
    ptrange->InsertText(Form("%3.1f < #it{p}_{T} < %3.1f GeV/#it{c}", fitXmin, fitXmax));
    ptrange->SetFillColor(kWhite); ptrange->SetBorderSize(0);
    
    TCanvas *cs = new TCanvas("cs","fit efficiency",0,0,1300,400);
    cs->Divide(4,1);
    cs->cd(1); h0->Draw(); ptrange->Draw("same");
    cs->cd(2); h1->Draw();
    cs->cd(3); h2->Draw();
    cs->cd(4); hchi2->Draw();
    cs->Print(Form("fitParams_%3.2f-%3.1f.png", fitXmin,fitXmax));
    return;
}

Float_t fitEfficiency(Float_t * pars, Float_t * parErrs,
                   Int_t icent=0, Float_t fitXmin = 1.0, Float_t fitXmax = 15.0) {
    
    Color_t color[6]={ kRed+1, kPink+6, kGreen+1, kAzure+1, kBlue+2, kBlack}; //combined tpc3s_tof3sveto
    Int_t marker[6]= { 21, 22, 32, 28, 24, 20};

    TFile * fin;
    if (icent>=0)
        fin = TFile::Open(Form("/Users/fbellini/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/eff_train12-20/multi-binB/efficiency_RsnOut_tpc2s_tof3sveto_cent%03i-%03i.root", icent*20, (icent+1)*20));
    else
        fin = TFile::Open("/Users/fbellini/alice/resonances/kstar_pA5.02TeV/MC_LF_pPb/eff_train12-20/binB/efficiency_RsnOut_tpc2s_tof3sveto_cent000-100.root");
    
    if (!fin) return;
    else Printf("Opening file %s", fin->GetName());
 
    TH1F * heff = (TH1F*) fin->Get("hEffVsPt");
    if (!heff) return;
    
    for (int j=1;j<heff->GetNbinsX()+1;j++) {
        Float_t binw = heff->GetXaxis()->GetBinWidth(j);
        //Printf("Bin %i [%3.1f,%3.1f], w = %3.1f) Rel uncert = %4.2f %%", j, heff->GetXaxis()->GetBinLowEdge(j), heff->GetXaxis()->GetBinUpEdge(j), binw, heff->GetBinError(j)*100./heff->GetBinContent(j));
    }
    //fitting function
    TF1 * func = new TF1("func","[0]-1/(x^[1]-[2])",0.0, 15.);
    //TF1 * func = new TF1("func","[0]-1/(x^[1]-[2])+1/(x^[1]-[3])",0.0, 15.);
    Printf("Fitting function defined f(x) = [0]-1/(x^[1]-[2])");
    func->SetLineColor(kBlue);
    //func->SetParLimits(0, 0.55, 0.65);
    //Printf("Parameter [0] limited in range [0.55, 0.65]");
    
    TFitResultPtr result = heff->Fit(func, "SEMRQ", "", fitXmin, fitXmax);
    if (result)
        Printf("Fitting histogram %s in pT range [%3.1f,%3.1f] GeV/c", heff->GetName(), fitXmin, fitXmax);
    else return;
    
    Float_t chi2ndf = func->GetChisquare()/func->GetNDF();
    
    for (int j=0;j<3;j++) {
        pars[j]    = result->Parameter(j);
        parErrs[j] = result->ParError(j);
        Printf("Par %i = %5.3f +/- %5.3f", j, pars[j], parErrs[j]);
    }
    
    Printf("Chi^2/NDF = %4.2f", chi2ndf);
    //result->Print("V");
    TCanvas *c1 = new TCanvas("c1","fit efficiency",0,0,700,500);
    heff->Draw();
    c1->Update();
    
    TH1F * hratio = (TH1F*) heff->Clone(Form("hist2fit_%i", icent<0? "mb" : icent));
    hratio->SetLineColor((icent<0? color[5] : color[icent]));
    hratio->SetMarkerColor((icent<0? color[5] : color[icent]));
    hratio->SetMarkerStyle((icent<0? marker[5] : marker[icent]));
    hratio->GetYaxis()->SetTitle("ratio data/fit");
    hratio->SetTitle(Form("ratio data/fit - multi bin %i", icent));

    hratio->Divide(func);
    hratio->GetYaxis()->SetRangeUser(0.7, 1.3);
    TLine * l0 = new TLine(fitXmin, 0.9, fitXmin, 1.1); l0->SetLineColor(kBlue); l0->SetLineStyle(2); l0->SetLineWidth(2);
    TLine * l1 = new TLine(fitXmax, 0.9, fitXmax, 1.1); l1->SetLineColor(kBlue); l1->SetLineStyle(2); l1->SetLineWidth(2);
    TCanvas *c2 = new TCanvas("c2","fit/efficiency",0,0,700,500);
    hratio->Draw();
    gPad->SetGridy();
    l0->Draw("same"); l1->Draw("same");
    
    c1->Print(Form("fitEff_%i_%3.1f-%3.1f.png",icent, fitXmin, fitXmax));
    c2->Print(Form("fitEffRatio_%i_%3.1f-%3.1f.png",icent, fitXmin, fitXmax));
    
    return chi2ndf;
}