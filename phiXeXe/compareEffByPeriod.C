#include "/Users/fbellini/alice/macros/MakeUp.C"
TH1F * GetPlotRatio(TH1* num,
		    TH1* den,
		    Bool_t display,
		    TString savepngName,
		    TString numName,
		    TString denName,
		    TString titleY,
		   TString errorType);

Int_t compareEffByPeriod()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleX(0.6);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatH(0.3);
  gStyle->SetStatW(0.25);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetHistLineWidth(0);
  gStyle->SetLineWidth(1);
  
  Char_t fileName[4][200] = {"~/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/simulation/LHC17j7_ZDCfix_extra/eff_C3_tpc2sPtDep_tof2sveto5smism.root",
			     "~/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/simulation/LHC17j7_ZDCfix/eff_C3_tpc2sPtDep_tof2sveto5smism.root",
			     "~/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/simulation/LHC17j7/eff_C3_tpc2sPtDep_tof2sveto5smism.root",
			     "~/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/simulation/eff_C3_tpc2sPtDep_tof2sveto5smism.root"};
  
  TFile * fin[4];
  TH1F * eff[3][4];
  TH1F * ratio[3][3];
  TString histPref = "hEffVsPt";
  Color_t color[4] = {kAzure-7, kOrange-3, kGreen+1, kBlack};
  Int_t marker[4] = {24, 25, 26, 20};
  Char_t * prodName[4][40] = {"LHC17j7_ZDCfix_extra", "LHC17j7_ZDCfix", "LHC17j7", "all productions"};
  
  for (int iprod = 0; iprod<4; iprod++) {
    fin[iprod] = new TFile(fileName[iprod]);
    if (!fin[iprod]) return iprod;
    for (int icent = 0; icent <3; icent++){
      eff[icent][iprod] = (TH1F*) fin[iprod]->Get(Form("%s%i",histPref.Data(), icent));
    } 
  }

  //prepare display
  TCanvas * ceff = new TCanvas("ceff", "ceff", 1600, 800);
  ceff->Divide(3,2);

  TLegend * leg[3];
  //get eff ratios
  for (int icent = 0; icent<3; icent++){
    leg[icent] = new TLegend(0.2, 0.65,0.7,0.95, Form("A#times#epsilon (%i-%i%%)", 30*icent, 30*(icent+1)));
    myLegendSetUp(leg[icent], 0.05);

    for (int iprod = 0; iprod<4; iprod++) {
      ceff->cd(icent+1);
      
      Beautify(eff[icent][iprod], color[iprod], 1, 2, marker[iprod], 1.3);
      if (iprod==0) leg[icent]->AddEntry(eff[icent][iprod], "LHC17j_ZDCfix_extra", "p");
      if (iprod==1) leg[icent]->AddEntry(eff[icent][iprod], "LHC17j_ZDCfix", "p");
      if (iprod==2) leg[icent]->AddEntry(eff[icent][iprod], "LHC17j", "p");
      if (iprod==3) leg[icent]->AddEntry(eff[icent][iprod], "all", "p");
      
      if (iprod==0) eff[icent][iprod]->Draw();
      else eff[icent][iprod]->Draw("same");

      if (iprod<3){
	ratio[icent][iprod] = GetPlotRatio(eff[icent][iprod], eff[icent][3], 0, "", Form("%s", prodName[iprod]), "all productions", "efficiency", "sumw2");
	ratio[icent][iprod]->GetYaxis()->SetRangeUser(0.6, 1.4);
	ratio[icent][iprod]->GetYaxis()->SetTitle("ratio to integrated sample");
	
	ceff->cd(4+icent);
	if (iprod==3) ratio[icent][iprod]->Draw();
	else ratio[icent][iprod]->Draw("same");
	
	Beautify(ratio[icent][iprod], color[iprod], 1, 2, marker[iprod], 1.3);
      }
      
    }
    ceff->cd(icent+1);
    leg[icent]->Draw();
  }
  return -1;
}


TH1F * GetPlotRatio(TH1* num,
		    TH1* den,
		    Bool_t display,
		    TString savepngName,
		    TString numName,
		    TString denName,
		    TString titleY,
		    TString errorType)
{
  Int_t rebinFactorNum = 0;
  Int_t rebinFactorDen = 0;
  Float_t showMin=0.0;
  Float_t showMax=-1.e10;
  /* returns ratio of two plots passed as arguments and displays if requested */
  if (!num) {
    Printf("invalid numerator passed as argument");
    return 0;
  }
  if (!den) {
    Printf("invalid numerator passed as argument");
    return 0;
  }

  TCanvas * cr = 0x0;
  TPad * pad1 = 0x0;
  TPad * pad2 = 0x0;
  if (display){
    cr = new TCanvas("cr","compare",800, 800);
    cr->cd();
    pad1 = new TPad("pad1","This is pad1",0.001,0.5,0.999,0.999);
    pad1->SetFillColor(0);
    pad1->SetBorderMode(0);
    pad1->SetBorderSize(0);
    pad1->SetMargin(0.15,0.05,0.01,0.02);

    pad2 = new TPad("pad2","This is pad2",0.001,0.001, 0.999,0.5);
    pad2->SetFillColor(0);
    pad2->SetBorderMode(0);
    pad2->SetBorderSize(0);
    pad2->SetMargin(0.15,0.05,0.2,0.01);
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    pad1->SetLogy();
    if (numName.CompareTo("",TString::kExact)==0) numName = num->GetName();
    if (numName.CompareTo("",TString::kExact)==0) denName = den->GetName();
    num->SetTitle(numName.Data());

    den->SetTitle(denName.Data());
    num->SetFillStyle(0);
    den->SetFillStyle(0);
    den->GetYaxis()->SetTitle(titleY.Data());
  
    if (showMax>showMin) den->GetXaxis()->SetRangeUser(showMin,showMax);
    den->Draw("");
    num->Draw("same");
    den->GetYaxis()->SetTitle(titleY.Data());

    TLegend * leg = (TLegend*) pad1->BuildLegend(0.7,0.75,0.93,0.93);
    myLegendSetUp(leg, 0.06);
    leg->SetEntryOption("pl");
  }

  //get ratio
  // num->Sumw2();
  // den->Sumw2();
  if (rebinFactorNum>0) num->Rebin(rebinFactorNum);
  if (rebinFactorDen>0) den->Rebin(rebinFactorDen);
  num->SetTitle(numName.Data());
  den->SetTitle(denName.Data());

  TH1F * ratio = (TH1F*) num->Clone("LSB");
  ratio->SetTitle(Form("%s/%s",num->GetTitle(), den->GetTitle()));
  ratio->SetMarkerStyle(20);
  if (errorType.Contains("num")) {
    ratio->Reset();
    for (Int_t j=1; j<num->GetNbinsX()+1;j++){
      if (den->GetBinContent(j)>0) {
	ratio->SetBinContent(j, num->GetBinContent(j)/den->GetBinContent(j));
	ratio->SetBinError(j, num->GetBinError(j)/num->GetBinContent(j));
      } else {
	ratio->SetBinContent(j,0);
      	ratio->SetBinError(j,0);
      }
    }
  } else {
    ratio->Divide(num, den, 1., 1., errorType.Data());
  }
  
  if (display) {
    Float_t xmax = ratio->GetXaxis()->GetBinUpEdge(ratio->GetNbinsX());
    TLine * at1 = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1., xmax, 1.);
    at1->SetLineStyle(7);
    at1->SetLineWidth(1);

    TLine * at120 = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1.2, xmax, 1.2);
    at120->SetLineStyle(2);
    at120->SetLineWidth(1);
    
    TLine * at080 = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 0.8, xmax, 0.8);
    at080->SetLineStyle(2);
    at080->SetLineWidth(1);

    ratio->GetYaxis()->SetTitle(Form("%s / %s",num->GetTitle(), den->GetTitle()));
    ratio->GetXaxis()->SetTitle("#it{p}_{T}(GeV/#it{c})");
    if (showMax>showMin) ratio->GetXaxis()->SetRangeUser(showMin,showMax);
    pad2->cd();
    ratio->Draw("");
    TPaveText *pave = new TPaveText(0.6, 0.25, 0.8, 0.32, "NDC");
    if (errorType.Contains("num"))
      pave->InsertText(Form("Relative uncert. of numerator."));
    else
      pave->InsertText(Form("%s uncert.", errorType.Data()));
    pave->SetBorderSize(0);
    pave->SetTextFont(42);
    pave->SetFillStyle(0);
    pave->Draw();
    at1->Draw();
    at120->Draw();
    at080->Draw();
    // TLegend * leg2 = (TLegend*) pad2->BuildLegend(0.6,0.88, 0.88,0.94);
    // leg2->SetEntryOption("p");
    // myLegendSetUp(leg2, 0.06);
  }
  
  if (savepngName.Contains(".png") || savepngName.Contains(".eps")) cr->Print(savepngName.Data());
  return ratio;
}
