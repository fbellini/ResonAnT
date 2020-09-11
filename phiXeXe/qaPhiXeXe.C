#include "/Users/fbellini/alice/macros/cosmetics/SetStyle.C"
#include "/Users/fbellini/alice/macros/cosmetics/MakeUp.C"
#include "/Users/fbellini/alice/macros/utils/fitSlices.C"

void compareLowBdca(TString finname = "RsnOut.root");

void qaPhiXeXe(Int_t nCuts = 0, TString finname = "RsnOut.root")
{
  gStyle->SetOptStat(10);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatH(0.3);
  gStyle->SetStatW(0.25);
  
  
  TString listName = "RsnOut_tpc2sPtDep_tof3s";
  switch (nCuts){
  case 1:
    listName = "RsnOut_tpc2s_tof3sveto";
    break;
  case 2:
    listName = "RsnOut_tpc3s_tof3sveto";
    break;
  case 3:
    listName = "RsnOut_tpc2s_tof4sveto";
    break;
  case 4:
    listName = "RsnOut_tpc3s";
    break;
  case 5:
    listName = "RsnOut_tof3s";
    break;
  case 6:
    listName = "RsnOut_tpc2sPtDep";
    break;
  case 7:
    listName = "RsnOut_tpc2sPtDep_tof3sveto";
    break;
  case 8:
    listName = "RsnOut_tpc2sPtDep_tof4sveto5smism";
    break;
  case 9:
    listName = "RsnOut_tpc3sPtDep_tof3sveto5smism";
    break;
  case 10:
    listName = "RsnOut_tpc2sPtDep_tof3sveto5smism";//buggy label -->it's a 3 sigma veto instead
    break;
  case 11:
    listName = "RsnOut_tpc2sPtDep_tof3sveto_elRej";
    break;
  case 101:
    listName = "RsnOut_default_LowBdca";
    break;
  default:
    listName = "RsnOut_tpc2sPtDep_tof3s";
  }

  TFile * fin = new TFile(finname.Data());
  if (!fin) return;
  TList * lin = (TList*) fin->Get(listName.Data());
  if (!lin) return;

  TH2F * hTPCpidQual = (TH2F*) lin->FindObject("cutQ_bit5.TPC_nsigmaK_VsPtpc_pTPC_K");
  hTPCpidQual->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV#it{c})");
  hTPCpidQual->GetYaxis()->SetTitle("N#sigma_{K}^{TPC}");
  hTPCpidQual->GetYaxis()->SetTitleSize(0.06);
  hTPCpidQual->GetYaxis()->SetLabelSize(0.06);
  hTPCpidQual->GetXaxis()->SetTitleSize(0.06);
  hTPCpidQual->GetXaxis()->SetLabelSize(0.06);
  hTPCpidQual->GetXaxis()->SetRangeUser(0.1, 11.0);
  
  TH2F * hTPCpidCut = (TH2F*) lin->FindObject("cutKa.TPC_nsigmaK_VsPtpc_pTPC_K");
  hTPCpidCut->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV#it{c})");
  hTPCpidCut->GetYaxis()->SetTitle("N#sigma_{K}^{TPC}");
  hTPCpidCut->SetTitle(listName.Data());
  hTPCpidCut->GetYaxis()->SetTitleSize(0.06);
  hTPCpidCut->GetYaxis()->SetLabelSize(0.06);
  hTPCpidCut->GetXaxis()->SetTitleSize(0.06);
  hTPCpidCut->GetXaxis()->SetLabelSize(0.06);
  hTPCpidCut->GetXaxis()->SetRangeUser(0.1, 11.0);

  TH2F * hTOFpidQual = (TH2F*) lin->FindObject("cutQ_bit5.TOF_nsigmaK_vsP_p_K");
  hTOFpidQual->GetXaxis()->SetTitle("#it{p} (GeV#it{c})");
  hTOFpidQual->GetYaxis()->SetTitle("N#sigma_{K}^{TOF}");
  hTOFpidQual->GetYaxis()->SetTitleSize(0.06);
  hTOFpidQual->GetYaxis()->SetLabelSize(0.06);
  hTOFpidQual->GetXaxis()->SetTitleSize(0.06);
  hTOFpidQual->GetXaxis()->SetLabelSize(0.06);
  hTOFpidQual->GetXaxis()->SetRangeUser(0.1, 11.0);
  fitSlices(hTPCpidQual);
  
  TH2F * hTOFpidCut = (TH2F*) lin->FindObject("cutKa.TOF_nsigmaK_vsP_p_K");
  hTOFpidCut->GetXaxis()->SetTitle("#it{p} (GeV#it{c})");
  hTOFpidCut->GetYaxis()->SetTitle("N#sigma_{K}^{TOF}");
  hTOFpidCut->SetTitle(listName.Data());
  hTOFpidCut->GetYaxis()->SetTitleSize(0.06);
  hTOFpidCut->GetYaxis()->SetLabelSize(0.06);
  hTOFpidCut->GetXaxis()->SetTitleSize(0.06);
  hTOFpidCut->GetXaxis()->SetLabelSize(0.06);
  hTOFpidCut->GetXaxis()->SetRangeUser(0.1, 11.0);

  TF1 fg("fg","gaus",-2.,2.); // fit range +- 6 sigma
  TLine l;
  TObjArray arrTPC;
  TObjArray arrTOF;
  fg.SetParameters(1,0,1);
  hTPCpidQual->FitSlicesY(&fg,0,-1,0,"NQR",&arrTPC);

  TH1 *hMtpc=(TH1*)arrTPC.At(1);
  hMtpc->SetMarkerStyle(20);
  hMtpc->SetMarkerSize(.5);

  TH1 *hStpc=(TH1*)arrTPC.At(2);
  hStpc->SetMarkerStyle(20);
  hStpc->SetMarkerSize(.5);
  hStpc->SetMarkerColor(kRed);
  hStpc->SetLineColor(kRed);
  

  hTOFpidQual->FitSlicesY(&fg,0,-1,0,"NQR",&arrTOF);
  TH1 *hMtof=(TH1*)arrTOF.At(1);
  hMtof->SetMarkerStyle(20);
  hMtof->SetMarkerSize(.5);
  
  TH1 *hStof=(TH1*)arrTOF.At(2);
  hStof->SetMarkerStyle(20);
  hStof->SetMarkerSize(.5);
  hStof->SetMarkerColor(kRed);
  hStof->SetLineColor(kRed);
  hMtof->Draw("sames");
  hStof->Draw("sames");
  
  
  TCanvas * c = new TCanvas("c", "c", 1200, 800);
  c->Divide(2,2);
  c->cd(1); hTPCpidQual->Draw("colz");
  hMtpc->DrawClone("sames");
  hStpc->DrawClone("sames");
  l.SetLineColor(kBlack);
  l.DrawLine(.1,0,10,0);
  l.SetLineColor(kRed);
  l.DrawLine(.1,1,10,1);
  gPad->SetLogz();gPad->SetLogx();
  
  c->cd(2); hTPCpidCut->Draw("colz");
  gPad->SetLogz();gPad->SetLogx();

  c->cd(3); hTOFpidQual->Draw("colz");
  hMtof->DrawClone("sames");
  hStof->DrawClone("sames");
  l.SetLineColor(kBlack);
  l.DrawLine(.1,0,10,0);
  l.SetLineColor(kRed);
  l.DrawLine(.1,1,10,1);
  gPad->SetLogz();gPad->SetLogx();
 
  
  c->cd(4); hTOFpidCut->Draw("colz"); 
  gPad->SetLogz();gPad->SetLogx();

  return;
}



void compareLowBdca(TString finname){
  SetStyle();
  TFile * fin = new TFile(finname.Data());
  if (!fin) return;

  TList * lin1 = (TList*) fin->Get("RsnOut_default_LowBdca");
  TList * lin2 = (TList*) fin->Get("RsnOut_tpc2sPtDep_tof3sveto5smism");
  if (!lin1 || !lin2) return;

  TH2F * hLowB = (TH2F*) lin1->FindObject("cutKa.DCAxyVsPt_pt_DCAxy");
  TH2F * hNominalB = (TH2F*) lin2->FindObject("cutKa.DCAxyVsPt_pt_DCAxy");

  if (!hLowB || !hNominalB) return;
  int istart = hLowB->GetXaxis()->FindBin(0.2);
  int istop = hLowB->GetXaxis()->FindBin(1.5);
  TH1F * hLowB_py = (TH1F *) hLowB->ProjectionY("hLowB_py", istart, istop);
  Beautify(hLowB_py, kRed, 1, 2, 20, 1.3);
  hLowB_py->SetTitle("low B");
  TH1F * hNominalB_py = (TH1F *) hNominalB->ProjectionY("hNominalB_py", istart, istop);
  Beautify(hNominalB_py, kBlue, 1, 2, 21, 1.3);
  hNominalB_py->SetTitle("standard");

  TCanvas *cdca = new TCanvas("cdca", "DCAxy", 800, 900);
  cdca->cd();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  gPad->SetMargin(0.2, 0.05, 0.2, 0.07);
  hLowB_py->Draw("hist");
  hLowB_py->GetXaxis()->SetRangeUser(-0.5, 0.5);
  hLowB_py->GetXaxis()->SetTitle("DCA_{xy} (cm)");
  hLowB_py->GetYaxis()->SetTitle("tracks");
  hNominalB_py->Draw("hist same");
  TLegend * leg = (TLegend*) gPad->BuildLegend(0.7, 0.7, 0.9, 0.9, "DCAxy cut");
  myLegendSetUp(leg, 0.04);
  
  return;
}