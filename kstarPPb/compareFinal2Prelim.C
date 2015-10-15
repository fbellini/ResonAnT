void compareFinal2Prelim(bool useSys=0,
			 TString final = "$HOME/alice/resonances/kstar_pA5.02TeV/output_LF5455/multi/finalCand_kstar_pPb_smoothSys.root",
			 TString prelim = "$HOME/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/preliminary/prelim_kstar_pPb_smoothSys.root")
{

  TString uncert="";

  if (useSys==1) uncert="_sys"; 
  gROOT->LoadMacro("$ASD/GetPlotRatio.C");
  TFile *_file0 = TFile::Open(final.Data());
  if (!_file0) {Printf("Problem with file %s", final.Data()); return;}

  TString prefixfinal = "hKstar"; 
  TString prefixprelim = "hKstar"; 
  
  TH1F * h100p = (TH1F*) _file0->Get(Form("%s_100%s", prefixfinal.Data(), uncert.Data()));
  TH1F * h0p = (TH1F*) _file0->Get(Form("%s_0%s", prefixfinal.Data(), uncert.Data()));
  TH1F * h1p = (TH1F*) _file0->Get(Form("%s_1%s", prefixfinal.Data(), uncert.Data()));
  TH1F * h2p = (TH1F*) _file0->Get(Form("%s_2%s", prefixfinal.Data(), uncert.Data()));
  TH1F * h3p = (TH1F*) _file0->Get(Form("%s_3%s", prefixfinal.Data(), uncert.Data()));
  TH1F * h4p = (TH1F*) _file0->Get(Form("%s_4%s", prefixfinal.Data(), uncert.Data()));

  TFile *_file1 = TFile::Open(prelim.Data());//
  if (!_file1) {Printf("Problem with file %s", prelim.Data()); return;}

  TH1F * h100old = (TH1F*) _file1->Get(Form("hKstar_100%s",uncert.Data()));
  TH1F * h0old = (TH1F*) _file1->Get(Form("%s_0%s",prefixprelim.Data(),uncert.Data()));
  TH1F * h1old = (TH1F*) _file1->Get(Form("%s_1%s",prefixprelim.Data(),uncert.Data()));
  TH1F * h2old = (TH1F*) _file1->Get(Form("%s_2%s",prefixprelim.Data(),uncert.Data()));
  TH1F * h3old = (TH1F*) _file1->Get(Form("%s_3%s",prefixprelim.Data(),uncert.Data()));
  TH1F * h4old = (TH1F*) _file1->Get(Form("%s_4%s",prefixprelim.Data(),uncert.Data()));
  
  useSys=0;
  TH1D * r100 = GetPlotRatio(h100p, h100old, 1, "final","preliminary","yields",0,0,Form("final2prel_100_%sUncert_12oct15.png",(useSys?"":"Stat")),0.0,14.9,(useSys?"num":"B"));
  TH1D * r0 = GetPlotRatio(h0p, h0old, 1, "final","preliminary","yields",0,0,Form("final2prel_0_%sUncert_12oct15.png",(useSys?"":"Stat")),0.0,14.9,(useSys?"num":"B")); 
  TH1D * r1 = GetPlotRatio(h1p, h1old, 1, "final","preliminary","yields",0,0,Form("final2prel_1_%sUncert_12oct15.png",(useSys?"":"Stat")),0.0,14.9,(useSys?"num":"B"));
  TH1D * r2 = GetPlotRatio(h2p, h2old, 1, "final","preliminary","yields",0,0,Form("final2prel_2_%sUncert_12oct15.png",(useSys?"":"Stat")),0.0,14.9,(useSys?"num":"B")) ;
  TH1D * r3 = GetPlotRatio(h3p, h3old, 1, "final","preliminary","yields",0,0,Form("final2prel_3_%sUncert_12oct15.png",(useSys?"":"Stat")),0.0,14.9,(useSys?"num":"B")) ;
  TH1D * r4 = GetPlotRatio(h4p, h4old, 1, "final","preliminary","yields",0,0,Form("final2prel_4_%sUncert_12oct15.png",(useSys?"":"Stat")),0.0,5.9,(useSys?"num":"B")) ;
  r100->SetName(Form("final2prel_100_%sUncert",(useSys?"":"Stat")));
  r0->SetName(Form("final2prel_0_%sUncert",(useSys?"":"Stat")));
  r1->SetName(Form("final2prel_1_%sUncert",(useSys?"":"Stat")));
  r2->SetName(Form("final2prel_2_%sUncert",(useSys?"":"Stat")));
  r3->SetName(Form("final2prel_3_%sUncert",(useSys?"":"Stat")));
  r4->SetName(Form("final2prel_4_%sUncert",(useSys?"":"Stat")));
  
  useSys=1;
  TH1D * r100s = GetPlotRatio(h100p, h100old, 1, "final","preliminary","yields",0,0,Form("final2prel_100_num%sUncert_12oct15.png",(useSys?"":"Stat")),0.0,14.9,(useSys?"num":"B")); 
  TH1D * r0s = GetPlotRatio(h0p, h0old, 1, "final","preliminary","yields",0,0,Form("final2prel_0_%sUncert_12oct15.png",(useSys?"":"Stat")),0.0,14.9,(useSys?"num":"B")); 
  TH1D * r1s = GetPlotRatio(h1p, h1old, 1, "final","preliminary","yields",0,0,Form("final2prel_1_%sUncert_12oct15.png",(useSys?"":"Stat")),0.0,14.9,(useSys?"num":"B"));
  TH1D * r2s = GetPlotRatio(h2p, h2old, 1, "final","preliminary","yields",0,0,Form("final2prel_2_%sUncert_12oct15.png",(useSys?"":"Stat")),0.0,14.9,(useSys?"num":"B")) ;
  TH1D * r3s = GetPlotRatio(h3p, h3old, 1, "final","preliminary","yields",0,0,Form("final2prel_3_%sUncert_12oct15.png",(useSys?"":"Stat")),0.0,14.9,(useSys?"num":"B")) ;
  TH1D * r4s = GetPlotRatio(h4p, h4old, 1, "final","preliminary","yields",0,0,Form("final2prel_4_%sUncert_12oct15.png",(useSys?"":"Stat")),0.0,5.9,(useSys?"num":"B")) ;
  r100s->SetName(Form("final2prel_100_%sUncert",(useSys?"":"Stat")));
  r0s->SetName(Form("final2prel_0_%sUncert",(useSys?"":"Stat")));
  r1s->SetName(Form("final2prel_1_%sUncert",(useSys?"":"Stat")));
  r2s->SetName(Form("final2prel_2_%sUncert",(useSys?"":"Stat")));
  r3s->SetName(Form("final2prel_3_%sUncert",(useSys?"":"Stat")));
  r4s->SetName(Form("final2prel_4_%sUncert",(useSys?"":"Stat")));
  
  Color_t color[6] = {kOrange+1, kSpring+5, kTeal-5, kAzure+1, kBlue+1, kGray};
  r100s->SetFillStyle(3001);  r100s->SetFillColor(color[5]);
  r0s->SetFillStyle(3001);  r0s->SetFillColor(color[0]);
  r4s->SetFillStyle(3001);  r4s->SetFillColor(color[4]);



  TLine * l0 = new TLine(0.0, 1., 14.9, 1.); 
  l0->SetLineStyle(3);
  l0->SetLineColor(kBlack);
  TLine * l1 = new TLine(0.0, 0.9, 14.9, 0.9); 
  l1->SetLineStyle(2);
  l1->SetLineColor(kGreen+1);
  TLine * l2 = new TLine(0.0, 1.1, 14.9, 1.1); 
  l2->SetLineStyle(2);
  l2->SetLineColor(kGreen+1);

  TLine * ls0 = new TLine(0.0, 1., 5.9, 1.); 
  ls0->SetLineStyle(3);
  ls0->SetLineColor(kBlack);
  TLine * ls1 = new TLine(0.0, 0.9, 5.9, 0.9); 
  ls1->SetLineStyle(2);
  ls1->SetLineColor(kGreen+1);
  TLine * ls2 = new TLine(0.0, 1.1, 5.9, 1.1); 
  ls2->SetLineStyle(2);
  ls2->SetLineColor(kGreen+1);

  gROOT->LoadMacro("$ASD/fig_template.C");
  SetStyle();
  TCanvas * summary = new TCanvas("summary","summary", 600, 600);
  summary->Divide(1,3);
  summary->cd(1);
  r100s->Draw("E2");
  r100->Draw("same");
  l1->Draw(); l0->Draw(); l2->Draw();
  summary->cd(2);
  r0s->Draw("E2");
  r0->Draw("same");
  l1->Draw(); l0->Draw(); l2->Draw();

  // summary->cd(3);
  // r1s->Draw("E2");
  // r1->Draw("same");
  // summary->cd(4);
  // r2s->Draw("E2");
  // r2->Draw("same");
  // summary->cd(5);
  // r3s->Draw("E2");
  // r3->Draw("same");
  summary->cd(3);
  r4s->Draw("E2");
  r4->Draw("same");
  ls1->Draw(); ls0->Draw(); ls2->Draw();

  TFile * fout = new TFile("final2prel_12oct15.root","recreate");
  fout->cd();
  r100->Write();
  r100s->Write();
  r0->Write();
  r0s->Write();
  r1->Write();
  r1s->Write();
  r2->Write();
  r2s->Write();
  r3->Write();
  r3s->Write();
  r4->Write();
  r4s->Write();
  
  return;
}
