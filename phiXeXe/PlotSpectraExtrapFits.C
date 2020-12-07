
#include "/Users/fbellini/alice/macros/cosmetics/SetStyle.C"
#include "/Users/fbellini/alice/macros/cosmetics/MakeUp.C"

void PlotSpectraExtrapFits()
{
    /*TFile fbgbw("FITSPECTRUM_bgbw_0.5-10.0_09may18.root", "read");
    TFile flevy("FITSPECTRUM_levy_0.5-10.0_09may18.root", "read");
    TFile fbose("FITSPECTRUM_bose_0.5-10.0_09may18.root", "read");
    TFile fboltz("FITSPECTRUM_boltz_0.5-10.0_09may18.root", "read");
    TFile fmtexp("FITSPECTRUM_mtexp_0.5-10.0_09may18.root", "read");
*/
    TFile fbgbw("FITSPECTRUM_bgbw_0.5-10.0_30oct20.root", "read");
    TFile flevy("FITSPECTRUM_levy_0.5-10.0_30oct20.root", "read");
    TFile fbose("FITSPECTRUM_bose_0.5-10.0_30oct20.root", "read");
    TFile fboltz("FITSPECTRUM_boltz_0.5-10.0_30oct20.root", "read");
    TFile fmtexp("FITSPECTRUM_mtexp_0.5-10.0_30oct20.root", "read");

    TH1F * hphi_stat[5];
    TH1F * hphi_syst[5];
    TF1 * bgbw[5];
    TF1 * levy[5];
    TF1 * bose[5];
    TF1 * boltz[5]; 
    TF1 * mtexp[5];
    TCanvas * cs[5];

    Color_t color[5] = {kOrange, kSpring+5, kTeal+5, kBlue+1, kMagenta+3};
    int Marker_Style[5] = {20, 21, 34, 24, 22};
    float Marker_Size = 1.2;
    int fLine_width = 2;
    int cent[] = {0, 10, 30, 50, 70, 90};
    for (int j = 0; j<5; j++) {
        hphi_stat[j] = (TH1F*) fbgbw.Get(Form("hPhiXeXe_cent%i_stat",j));
        hphi_syst[j] = (TH1F*) fbgbw.Get(Form("hPhiXeXe_cent%i_syst",j));
        bgbw[j] = (TF1*)fbgbw.Get(Form("bgbw%i",j)); 
        levy[j] = (TF1*)flevy.Get(Form("levy%i",j));
        bose[j] = (TF1*)fbose.Get(Form("bose%i",j));
        boltz[j] = (TF1*)fboltz.Get(Form("boltz%i",j));
        mtexp[j] = (TF1*)fmtexp.Get(Form("mtexp%i",j));

        bgbw[j]->SetLineColor(kRed); bgbw[j]->SetLineStyle(1); bgbw[j]->SetLineWidth(fLine_width);
        levy[j]->SetLineColor(kBlue); levy[j]->SetLineStyle(3); levy[j]->SetLineWidth(fLine_width);
        bose[j]->SetLineColor(kSpring-5); bose[j]->SetLineStyle(9); bose[j]->SetLineWidth(fLine_width);
        boltz[j]->SetLineColor(kMagenta); boltz[j]->SetLineStyle(6); boltz[j]->SetLineWidth(fLine_width);
        mtexp[j]->SetLineColor(kViolet-6); mtexp[j]->SetLineStyle(10); mtexp[j]->SetLineWidth(fLine_width);

        hphi_stat[j]->SetTitle("stat. uncert.");        
        hphi_syst[j]->SetTitle("syst. uncert.");        
        Beautify(hphi_stat[j], color[j], 1, 2, Marker_Style[j], Marker_Size);
        Beautify(hphi_syst[j], color[j], 1, 2, Marker_Style[j], Marker_Size);
        hphi_syst[j]->SetFillColorAlpha(color[j], 0.4);
        hphi_syst[j]->SetFillStyle(1001);

        SetStyle();
        cs[j] = new TCanvas(Form("cs%i",j),Form("cs%i",j), 900, 1200);
        cs[j]->cd(); 
        gPad->SetMargin(0.2, 0.05, 0.2, 0.05);
        gPad->SetLogy(); 
        gPad->SetTicky();
        gPad->SetTicky();
        hphi_syst[j]->Draw("E2"); 
        hphi_stat[j]->GetXaxis()->SetRangeUser(0.01, 9.99);
        hphi_stat[j]->Draw("same");
        bgbw[j]->Draw("same");
        levy[j]->Draw("same");
        bose[j]->Draw("same");
        boltz[j]->Draw("same");
        mtexp[j]->Draw("same");

        TLegend * leg = new TLegend(0.6, 0.6, 0.9, 0.9, Form("#phi, %i-%i %%", cent[j],cent[j+1]));
        myLegendSetUp(leg, 0.04);
        leg->AddEntry(hphi_syst[j], "syst. uncert.", "f");
        leg->AddEntry(hphi_stat[j], "stat. uncert.", "lp");
        leg->AddEntry(bgbw[j], "BG-bw", "l");
        leg->AddEntry(levy[j], "Levy", "l");
        leg->AddEntry(bose[j], "Bose", "l");
        leg->AddEntry(boltz[j], "Boltz", "l");
        leg->AddEntry(mtexp[j], "m_{T}-exp", "lp");
        leg->Draw();
        cs[j]->Print(Form("yield_syst_functions_%i.png",j));
        cs[j]->Print(Form("yield_syst_functions_%i.pdf",j));
    }    
}