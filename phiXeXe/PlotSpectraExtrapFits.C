
#include "/Users/fbellini/alice/macros/cosmetics/SetStyle.C"
#include "/Users/fbellini/alice/macros/cosmetics/MakeUp.C"

void PlotSpectraExtrapFits()
{
    TFile fbgbw("FITSPECTRUM_bgbw_0.5-10.0_09may18.root", "read");
    TFile flevy("FITSPECTRUM_levy_0.5-10.0_09may18.root", "read");
    TFile ffermi("FITSPECTRUM_fermi_0.5-10.0_09may18.root", "read");
    TFile fboltz("FITSPECTRUM_boltz_0.5-10.0_09may18.root", "read");
    TFile fmtexp("FITSPECTRUM_mtexp_0.5-10.0_09may18.root", "read");

    TH1F * hphi_stat[4];
    TH1F * hphi_syst[4];
    TF1 * bgbw[4];
    TF1 * levy[4];
    TF1 * fermi[4];
    TF1 * boltz[4]; 
    TF1 * mtexp[4];
    TCanvas * cs[4];

    Color_t color[4] = {kRed, kOrange, kGreen+1, kBlue};
    int Marker_Style[4] = {20, 21, 34, 24};
    float Marker_Size = 1.2;
    int fLine_width = 2;
    int cent[] = {0, 10, 30, 60, 90};
    for (int j = 0; j<4; j++) {
    

        hphi_stat[j] = (TH1F*) fbgbw.Get(Form("hPhiXeXe_cent%i_stat",j));
        hphi_syst[j] = (TH1F*) fbgbw.Get(Form("hPhiXeXe_cent%i_syst",j));
        bgbw[j] = (TF1*)fbgbw.Get(Form("bgbw%i",j)); 
        levy[j] = (TF1*)flevy.Get(Form("levy%i",j));
        fermi[j] = (TF1*)ffermi.Get(Form("fermi%i",j));
        boltz[j] = (TF1*)fboltz.Get(Form("boltz%i",j));
        mtexp[j] = (TF1*)fmtexp.Get(Form("mtexp%i",j));

        bgbw[j]->SetLineColor(kRed); bgbw[j]->SetLineStyle(1); bgbw[j]->SetLineWidth(fLine_width);
        levy[j]->SetLineColor(kBlue); levy[j]->SetLineStyle(3); levy[j]->SetLineWidth(fLine_width);
        fermi[j]->SetLineColor(kSpring-5); fermi[j]->SetLineStyle(9); fermi[j]->SetLineWidth(fLine_width);
        boltz[j]->SetLineColor(kMagenta); boltz[j]->SetLineStyle(6); boltz[j]->SetLineWidth(fLine_width);
        mtexp[j]->SetLineColor(kViolet-6); mtexp[j]->SetLineStyle(10); mtexp[j]->SetLineWidth(fLine_width);

        hphi_stat[j]->SetTitle("stat. uncert.");        
        hphi_syst[j]->SetTitle("syst. uncert.");        
        Beautify(hphi_stat[j], color[j], 1, 2, Marker_Style[j], Marker_Size);
        Beautify(hphi_syst[j], color[j], 1, 2, Marker_Style[j], Marker_Size);
        hphi_syst[j]->SetFillColorAlpha(color[j], 0.2);
        hphi_syst[j]->SetFillStyle(1001);

        SetStyle();
        cs[j] = new TCanvas(Form("cs%i",j),Form("cs%i",j), 900, 1200);
        cs[j]->cd(); 
        gPad->SetMargin(0.2, 0.05, 0.2, 0.05);
        gPad->SetLogy(); 
        gPad->SetTicky();
        gPad->SetTicky();
        hphi_syst[j]->Draw("E2"); hphi_stat[j]->GetXaxis()->SetRangeUser(0.01, 9.99);
        hphi_stat[j]->Draw("same");
        bgbw[j]->Draw("same");
        levy[j]->Draw("same");
        fermi[j]->Draw("same");
        boltz[j]->Draw("same");
        mtexp[j]->Draw("same");

        TLegend * leg = new TLegend(0.6, 0.6, 0.9, 0.9, Form("#phi, %i-%i %%", cent[j],cent[j+1]));
        myLegendSetUp(leg, 0.04);
        leg->AddEntry(hphi_syst[j], "syst. uncert.", "f");
        leg->AddEntry(hphi_stat[j], "stat. uncert.", "lp");
        leg->AddEntry(bgbw[j], "bgbw", "l");
        leg->AddEntry(levy[j], "levy", "l");
        leg->AddEntry(fermi[j], "Fermi-Dirac", "l");
        leg->AddEntry(boltz[j], "Boltz", "l");
        leg->AddEntry(mtexp[j], "m_{T}-exp", "lp");
        leg->Draw();
        cs[j]->Print(Form("yield_syst_functions_%i.png",j));
        cs[j]->Print(Form("yield_syst_functions_%i.pdf",j));
    }    
}