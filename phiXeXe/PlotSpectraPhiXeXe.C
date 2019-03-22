#include "/Users/fbellini/alice/macros/MakeUp.C"
#include "/Users/fbellini/alice/macros/SetStyle.C"
#include "/Users/fbellini/alice/macros/GetGraphFromHisto.C"

TH1F * GetPhiXeXeSpectrumCent(Int_t cent = 0, Bool_t sys = 0);

void PlotSpectraPhiXeXe(Bool_t isPreliminary = 1)
{

  //plotting macro with good cosmetics for spectra
  SetStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  
  //cosmetics  
  Color_t color[4] = {kRed+1, kOrange, kSpring+5, kAzure+2};
  Int_t  Marker_Style[4] = {20, 21, 34, 33}; 
  Int_t cent[] = {0, 10, 30, 60, 90}; 
  TCanvas *c1 = new TCanvas("c1","spectra", 800, 950);
  c1->SetLogy();
  c1->SetTickx();
  c1->SetTicky();

  TH1F *frame = new TH1F("frame","frame", 105, 0.0, 10.5);
  frame->SetMinimum(2e-5);
  frame->SetMaximum(50.0);
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(1.6);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetRangeUser(9.e-6, 1e3);
  frame->GetXaxis()->SetTitleSize(0.05);
  frame->SetTitle("; #it{p}_{T} (GeV/#it{c}); 1/#it{N}_{ev} d^{2}#it{N}/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");

  TPaveText *phitext = new TPaveText(0.35,0.85,0.45,0.90,"brNDC");
  phitext->SetBorderSize(0);
  phitext->SetFillColor(0);
  phitext->SetFillStyle(0);
  phitext->SetTextAlign(32);
  phitext->SetTextSize(0.05);
  phitext->SetTextFont(42);
  phitext->InsertText("#bf{#phi}(1020)");
  
  TPaveText *titletext = new TPaveText(0.52,0.67,0.93,0.77,"brNDC");
  titletext->SetBorderSize(0);
  titletext->SetFillColor(0);
  titletext->SetFillStyle(0);
  titletext->SetTextAlign(12);
  titletext->SetTextSize(0.04);
  titletext->SetTextFont(42);
  titletext->AddText("#bf{ALICE Preliminary}");
  titletext->AddText("Xe-Xe, #sqrt{#it{s}_{NN}} = 5.44 TeV");
  //titletext->InsertText("|#it{y} | < 0.5");
  titletext->Draw();

  TPaveText *titletextunc = new TPaveText(0.22,0.18,0.27,0.22,"brNDC");
  titletextunc->SetBorderSize(0);
  titletextunc->SetFillColor(0);
  titletextunc->SetFillStyle(0);
  titletextunc->SetTextAlign(12);
  titletextunc->SetTextSize(0.035);
  titletextunc->SetTextFont(42);
  titletextunc->InsertText("Uncertainties: stat.(bar), syst.(box)");  
  
  TLegend * leg = new TLegend(0.52, 0.8, 0.93, 0.90);
  myLegendSetUp(leg, 0.035);
  leg->SetTextAlign(12);
  leg->SetNColumns(2);
  c1->cd();
  frame->Draw();
  titletext->Draw();
  titletextunc->Draw();
  phitext->Draw();

  Int_t Line_Style = 1;
  Int_t Line_Width = 2;
  Int_t Fill_Style = 0;
  
  TH1F * hStat[4];
  TH1F * hSyst[4];
  TGraphErrors * gStat[4];
  TGraphErrors * gSyst[4];
  Float_t multiplyingScalingFactor[4] = {2., 1., 1., 1.0};

  TFile * fout = new TFile("Preliminary_Spectra_phi_XeXe544TeV.root", "recreate");
  for (int ic = 0; ic < 4; ic++){
    hStat[ic] = (TH1F*) GetPhiXeXeSpectrumCent(ic, 0);
    hSyst[ic] = (TH1F*) GetPhiXeXeSpectrumCent(ic, 1);
    Beautify(hStat[ic], color[ic], Line_Style, Line_Width, Marker_Style[ic], 1.5);
    Beautify(hSyst[ic], color[ic], Line_Style, Line_Width, Marker_Style[ic], 0.1);
    hSyst[ic]->SetFillStyle(0);

    fout->cd();
    hStat[ic]->Write();
    hSyst[ic]->Write();

    gStat[ic] = (TGraphErrors *) GetGraphFromHisto(hStat[ic], kTRUE, multiplyingScalingFactor[ic]);
    gSyst[ic] = (TGraphErrors *) GetGraphFromHisto(hSyst[ic], kFALSE, multiplyingScalingFactor[ic]);
    BeautifyGraph(gStat[ic], color[ic], color[ic], Fill_Style, Line_Style, Line_Width, Marker_Style[ic], 1.5);
    BeautifyGraph(gSyst[ic], color[ic], color[ic], Fill_Style, Line_Style, Line_Width, Marker_Style[ic], 0.1);
    
    TString legLabel;
    if (multiplyingScalingFactor[ic]!=1.0) legLabel = Form("%i-%i%% #times%2.0f", cent[ic], cent[ic+1], multiplyingScalingFactor[ic]);
    else legLabel = Form("%i-%i%%", cent[ic], cent[ic+1] );
    leg->AddEntry(gStat[ic], legLabel.Data(), "pf");
    
    c1->cd();
    gSyst[ic]->Draw("E2 same");
    gStat[ic]->Draw("PZ same");    
  }
  
  c1->cd();
  leg->Draw();

  fout->cd();
  c1->Draw();
  
  TString nameimg = "PhiXeXe_Spectra";
  if (isPreliminary) nameimg.Prepend("Preliminary_");
  c1->Print(Form("%s.eps",nameimg.Data()));
  c1->Print(Form("%s.pdf",nameimg.Data()));
	   
  return;

}

//-------------------------------------------------------
TH1F * GetPhiXeXeSpectrumCent(Int_t cent, Bool_t sys)
{
  //
  // get phi spectrum from file
  ///Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiA3_tpc2sPtDep_tof2sveto5smism/spectra/finalWsyst_smooth1_01may18_*.root
  //
  TH1F * hist = 0x0;
  TFile * fin = TFile::Open(Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/ana0406esd710/phiA3_tpc2sPtDep_tof2sveto5smism/spectra/finalWsyst_smooth1_09may18_%i.root", cent));
  if (!fin) return hist;
  if (sys) hist = (TH1F*) fin->Get(Form("hCorrected_%i%i_syst",cent, cent));
  else hist = (TH1F*) fin->Get(Form("hCorrected_%i",cent));
  return hist;
}
