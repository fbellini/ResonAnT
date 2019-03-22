void makeKstarVsAntiKStarRatios(TString kstarfilename="/Users/bellini/alice/resonances/myKstar/pwglf_train_out/kstar_kppim/sub_EMnorm1.30-1.50_aod049_kstar.root", Int_t startBin = -1, Int_t stopBin = -1, Int_t startCBin = -1, Int_t stopCBin = -1, TString suffix = "")
{
  //"/Users/bellini/alice/resonances/myKstar/pwglf_train_out/kstar_kppim/sub_EMnorm1.30-1.50_aod049_kstar.root"
  //"/Users/bellini/alice/resonances/myKstar/pwglf_train_out/antikstar_kmpip/sub_EMnorm1.30-1.50_aod049_antikstar.root"
  Double_t pt[] = { 0.0, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.00 };

  if (startBin<0) startBin=0;
  if (stopBin<0) stopBin=11; 
  if (startCBin<0) startBin=0;
  if (stopCBin<0) stopBin=6;
 
  //various
  Color_t color[6] = {kRed+3, kBlue+3, kGreen+1, kMagenta+3, kOrange, kCyan+3};
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //read Kstar histos
  TFile * fileSub1 = TFile::Open(kstarfilename.Data()); //kstar
  TH1D* hUS1[55];
  TH1D* hEM1[55];
  TH1D* hSub1[55];
  TH1D* hLS1[55];
  TH1D* hSubLS1[55];

  for (Int_t ic=startCBin;ic<stopCBin;ic++){
    for (Int_t ip=startBin;ip<stopBin;ip++){
      Int_t index = stopBin*ic + ip;
      if (ip<10){
	hUS1[index] = (TH1D*) fileSub1->Get(Form("Signal_ptBin0%i_centBin0%i", ip, ic))->Clone();
	hEM1[index] = (TH1D*) fileSub1->Get(Form("norm_Mixing_ptBin0%i_centBin0%i", ip, ic))->Clone();
	hSub1[index] = (TH1D*) fileSub1->Get(Form("sub_norm_Mixing_ptBin0%i_centBin0%i", ip, ic))->Clone();
	hLS1[index] = (TH1D*) fileSub1->Get(Form("norm_Like_ptBin0%i_centBin0%i", ip, ic))->Clone();
	hSubLS1[index] = (TH1D*) fileSub1->Get(Form("sub_norm_Like_ptBin0%i_centBin0%i", ip, ic))->Clone();
      } else {
	hUS1[index] = (TH1D*) fileSub1->Get(Form("Signal_ptBin%i_centBin0%i", ip, ic))->Clone();
	hEM1[index] = (TH1D*) fileSub1->Get(Form("norm_Mixing_ptBin%i_centBin0%i", ip, ic))->Clone();
	hSub1[index] = (TH1D*) fileSub1->Get(Form("sub_norm_Mixing_ptBin%i_centBin0%i", ip, ic))->Clone();
	hLS1[index] = (TH1D*) fileSub1->Get(Form("norm_Like_ptBin%i_centBin0%i", ip, ic))->Clone();
	hSubLS1[index] = (TH1D*) fileSub1->Get(Form("sub_norm_Like_ptBin%i_centBin0%i", ip, ic))->Clone();
      }
    }
  }
  
  //read antiKstar histos
  TFile * fileSub2 = TFile::Open(kstarfilename.ReplaceAll("kstar_","antikstar_")); //antikstar
  TH1D* hUS2[55];
  TH1D* hEM2[55];
  TH1D* hSub2[55];
  TH1D* hLS2[55];
  TH1D* hSubLS2[55];

  for (Int_t ic=startCBin;ic<stopCBin;ic++){
    for (Int_t ip=startBin;ip<stopBin;ip++){
      Int_t index = stopBin*ic + ip;
       if (ip<10){
	hUS2[index] = (TH1D*) fileSub2->Get(Form("Signal_ptBin0%i_centBin0%i", ip, ic))->Clone();
	hEM2[index] = (TH1D*) fileSub2->Get(Form("norm_Mixing_ptBin0%i_centBin0%i", ip, ic))->Clone();
	hSub2[index] = (TH1D*) fileSub2->Get(Form("sub_norm_Mixing_ptBin0%i_centBin0%i", ip, ic))->Clone();
      	hLS2[index] = (TH1D*) fileSub2->Get(Form("norm_Like_ptBin0%i_centBin0%i", ip, ic))->Clone();
	hSubLS2[index] = (TH1D*) fileSub2->Get(Form("sub_norm_Like_ptBin0%i_centBin0%i", ip, ic))->Clone();
} else {
	hUS2[index] = (TH1D*) fileSub2->Get(Form("Signal_ptBin%i_centBin0%i", ip, ic))->Clone();
	hEM2[index] = (TH1D*) fileSub2->Get(Form("norm_Mixing_ptBin%i_centBin0%i", ip, ic))->Clone();
	hSub2[index] = (TH1D*) fileSub2->Get(Form("sub_norm_Mixing_ptBin%i_centBin0%i", ip, ic))->Clone();
      	hLS2[index] = (TH1D*) fileSub2->Get(Form("norm_Like_ptBin%i_centBin0%i", ip, ic))->Clone();
	hSubLS2[index] = (TH1D*) fileSub2->Get(Form("sub_norm_Like_ptBin%i_centBin0%i", ip, ic))->Clone();
       }    
    }
  }
  
  //Get ratios for each bins
  gROOT->LoadMacro("$ASD/GetPlotRatio.C");
  TH1D* ratioUS[55];
  TH1D* ratioEM[55];
  TH1D* ratioSub[55];
  TH1D* ratioLS[55];
  TH1D* ratioSubLS[55];

  TFile *file_rUS = new TFile(Form("ratioUS_%s.root",suffix.Data()),"recreate");
  TFile *file_rEM = new TFile(Form("ratioEM_%s.root",suffix.Data()),"recreate");
  TFile *file_rSub = new TFile(Form("ratioSub_%s.root",suffix.Data()),"recreate");
  TFile *file_rLS = new TFile(Form("ratioLS_%s.root",suffix.Data()),"recreate");
  TFile *file_rSubLS = new TFile(Form("ratioSubLS_%s.root",suffix.Data()),"recreate");

  TCanvas * cus[4];
  TCanvas * cem[4];
  TCanvas * cls[4];
  
  for (Int_t ic=startCBin;ic<stopCBin;ic++){
    cus[ic]= new TCanvas(Form("cus%i",ic),Form("cent. %i",ic), 900,500);
    cem[ic]= new TCanvas(Form("cem%i",ic),Form("cent. %i",ic), 900,500);
    cls[ic]= new TCanvas(Form("cls%i",ic),Form("cent. %i",ic), 900,500);

    for (Int_t ip=startBin;ip<stopBin;ip++){
      Int_t index = stopBin*ic + ip;
      ratioUS[index] = (TH1D*) GetPlotRatio(hUS2[index],hUS1[index],0, "K^{-}#pi^{+}","K^{+}#pi^{-}");
      ratioUS[index]->SetNameTitle(Form("ratioUS_cent%i_pt%i",ic,ip),Form("%2.1f<p_{T}<%3.1f (%i-%i%%)",pt[ip], pt[ip+1],ic*20,(ic+1)*20));
      ratioUS[index]->SetLineColor(color[0]-ip);
      ratioUS[index]->SetFillColor(kWhite);
      ratioUS[index]->SetMarkerColor(color[0]-ip);
      ratioUS[index]->SetMarkerStyle(20+ip);
      ratioUS[index]->SetMarkerSize(0.8);
      ratioUS[index]->GetYaxis()->SetRangeUser(0.5,1.2);

      ratioEM[index] = (TH1D*) GetPlotRatio(hEM2[index],hEM1[index],0, "K-#pi+ EM","K+#pi- EM");
      ratioEM[index]->SetNameTitle(Form("ratioEM_cent%i_pt%i",ic,ip),Form("%2.1f<p_{T}<%3.1f (%i-%i%%)",pt[ip], pt[ip+1],ic*20,(ic+1)*20));
      ratioEM[index]->SetLineColor(color[1]-ip);
      ratioEM[index]->SetFillColor(kWhite);
      ratioEM[index]->SetMarkerColor(color[1]-ip);
      ratioEM[index]->SetMarkerStyle(20+ip);
      ratioEM[index]->SetMarkerSize(0.8);
      ratioEM[index]->GetYaxis()->SetRangeUser(0.5,1.2);
      
      ratioSub[index] = (TH1D*) GetPlotRatio(hSub2[index],hSub1[index],0, "K-#pi+ Sub","K+#pi- Sub");
      ratioSub[index]->SetNameTitle(Form("ratioSub_cent%i_pt%i",ic,ip),Form("%2.1f<p_{T}<%3.1f (%i-%i%%)",pt[ip], pt[ip+1],ic*20,(ic+1)*20));
      ratioSub[index]->SetLineColor(color[2]-ip);
      ratioSub[index]->SetFillColor(kWhite);
      ratioSub[index]->SetMarkerColor(color[2]-ip);
      ratioSub[index]->SetMarkerStyle(20+ip);
      ratioSub[index]->SetMarkerSize(0.8);

      ratioLS[index] = (TH1D*) GetPlotRatio(hLS2[index],hLS1[index],0, "K-#pi- LS","K+#pi+ LS");
      ratioLS[index]->SetNameTitle(Form("ratioLS_cent%i_pt%i",ic,ip),Form("%2.1f<p_{T}<%3.1f (%i-%i%%)",pt[ip], pt[ip+1],ic*20,(ic+1)*20));
      ratioLS[index]->SetLineColor(color[3]-ip);
      ratioLS[index]->SetFillColor(kWhite);
      ratioLS[index]->SetMarkerColor(color[3]-ip);
      ratioLS[index]->SetMarkerStyle(20+ip);
      ratioLS[index]->SetMarkerSize(0.8);
      ratioLS[index]->GetYaxis()->SetRangeUser(0.5,1.2);

      ratioSubLS[index] = (TH1D*) GetPlotRatio(hSubLS2[index],hSubLS1[index],0, "K-#pi- Sub","K+#pi+ Sub");
      ratioSubLS[index]->SetNameTitle(Form("ratioSubLS_cent%i_pt%i",ic,ip),Form("%2.1f<p_{T}<%3.1f (%i-%i%%)",pt[ip], pt[ip+1],ic*20,(ic+1)*20));
      ratioSubLS[index]->SetLineColor(color[4]-ip);
      ratioSubLS[index]->SetFillColor(kWhite);
      ratioSubLS[index]->SetMarkerColor(color[4]-ip);
      ratioSubLS[index]->SetMarkerStyle(20+ip);
      ratioSubLS[index]->SetMarkerSize(0.8);

      file_rUS->cd(); ratioUS[index]->Write();
      file_rEM->cd(); ratioEM[index]->Write();
      file_rSub->cd(); ratioSub[index]->Write();
      file_rLS->cd(); ratioLS[index]->Write();
      file_rSubLS->cd(); ratioSubLS[index]->Write();

      cus[ic]->cd();
      ratioUS[index]->Draw(((ip==1)?"":"same"));
      cem[ic]->cd();
      ratioEM[index]->Draw(((ip==1)?"":"same"));
      cls[ic]->cd();
      ratioLS[index]->Draw(((ip==1)?"":"same"));
    }
    file_rUS->cd();
    cus[ic]->Write();
    file_rEM->cd();
    cem[ic]->Write();
    file_rLS->cd();
    cls[ic]->Write();
    lus = (TLegend*) cus[ic]->BuildLegend(0.78, 0.2,0.99,0.89,"US: K^{-}#pi^{+} / K^{+}#pi^{-}");
    lem = (TLegend*) cem[ic]->BuildLegend(0.78, 0.2,0.99,0.89,"EM: K^{-}#pi^{+} / K^{+}#pi^{-}");
    lls = (TLegend*) cls[ic]->BuildLegend(0.78, 0.2,0.99,0.89,"LS: K^{-}#pi^{+} / K^{+}#pi^{-}");
    lus->SetFillColor(kWhite);
    lus->SetLineColor(kWhite);
    lem->SetFillColor(kWhite);
    lem->SetLineColor(kWhite);
    lls->SetLineColor(kWhite);
    lls->SetFillColor(kWhite);
    lus->SetTextFont(42);
    lem->SetTextFont(42);
    lls->SetTextFont(42);
    
    cus[ic]->SaveAs(Form("ratioUS_c%i_pt%i-%i.png",ic,startBin,stopBin));
    cem[ic]->SaveAs(Form("ratioEM_c%i_pt%i-%i.png",ic,startBin,stopBin));
    cls[ic]->SaveAs(Form("ratioLS_c%i_pt%i-%i.png",ic,startBin,stopBin));
  }
  return;

}

//--------------------------------------------------
void makeNormRawSpectraRatio()
{

  Color_t color[]={kRed+1, kOrange+1, kGreen+2, kBlue+2};
  
  TFile * file1 = TFile::Open("/Users/bellini/alice/resonances/myKstar/pwglf_train_out/kstar_kppim/normYields/normYields_treebest_poly2.root"); //kstar
  TH1D* hUS1[4];

  TFile * file2 = TFile::Open("/Users/bellini/alice/resonances/myKstar/pwglf_train_out/antikstar_kmpip/normYields/normYields_treetreeBest.root"); //antikstar
  TH1D* hUS2[4];

 gROOT->LoadMacro("$ASD/GetPlotRatio.C");
 TH1D* ratioUS[4];
 TFile *file_rUS = new TFile("ratioUS.root","recreate");
 TCanvas * cus = new TCanvas(Form("cus"),Form("ratios"), 900,500);

 for (Int_t ic=0;ic<4;ic++){
   Int_t index = ic;
   hUS1[index] = (TH1D*) file1->Get(Form("hRawYieldVsPt_%i", ic))->Clone();
   hUS2[index] = (TH1D*) file2->Get(Form("hRawYieldVsPt_%i", ic))->Clone();

   ratioUS[index] = (TH1D*) GetPlotRatio(hUS2[index],hUS1[index],1, "K-#pi+","K+#pi-");
   ratioUS[index]->SetNameTitle(Form("ratioNormRawY_cent%i",ic),Form("ratioNormRawY_cent%i",ic));
   ratioUS[index]->SetLineColor(color[ic]);
   ratioUS[index]->SetFillColor(kWhite);
   ratioUS[index]->SetMarkerColor(color[ic]);
   ratioUS[index]->SetMarkerStyle(20+ic);
   ratioUS[index]->SetMarkerSize(0.8);
   ratioUS[index]->GetYaxis()->SetRangeUser(0.5,1.5);
   file_rUS->cd(); ratioUS[index]->Write();
   cus->cd();
   if (ic==0) ratioUS[index]->Draw();
   else ratioUS[index]->Draw("same");
   cus->BuildLegend();  
 }
 
 file_rUS->cd(); cus->Write();
 return;
}
