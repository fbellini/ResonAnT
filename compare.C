/*
fbellini@cern.ch - 29/06/2012

*/
TString macroDir = "/Users/bellini/alice/macro/kstar";

void compare(Char_t*file_ss="ss_kstar",Char_t*file_fb="fb_kstar", Int_t ipt =-1, Int_t icent=-1, Bool_t saveImg=0){
  
  gStyle->SetTextFont(42);

  TFile * fss= TFile::Open(Form("%s.root",file_ss));
  if (!fss) return;
  TFile * ffb= TFile::Open(Form("%s.root",file_fb));
  if (!ffb) return;

  TFile * fssE= TFile::Open("analysis.root");
  if (!fssE) return;
  TFile * ffbE= TFile::Open("analysisAOD_0-80.root");
  if (!ffbE) return;

  TH1I *hEss = (TH1I*) fssE->Get("RsnOut")->FindObject("hEventStat")->Clone();
  Double_t nEvents_ss = 1./(hEss->GetBinContent(4));
  Printf("SS entries = %e", hEss->GetBinContent(4)) ;
  TH1I *hEfb = (TH1*) ffbE->Get("RsnOut_Tof")->FindObject("hEventStat");
  Double_t nEvents_fb = 1./(hEfb->GetBinContent(4));
  Printf("FB entries = %e", hEfb->GetBinContent(4));
  Double_t ratioScaleF = hEfb->GetBinContent(4)/hEss->GetBinContent(4);

  fssE->Close();
  ffbE->Close();

  TAxis *ptbins = (TAxis*)fss->Get("ptbins");
  Int_t nPtBins = ptbins->GetNbins();
  TAxis *centbins = (TAxis*)fss->Get("centbins");
  Int_t nCentBins = centbins->GetNbins();
  
  TString fout_name = "compare_ss_fb.root";
  TFile * fout = new TFile(fout_name.Data(),"recreate");
  fout->cd();
  ptbins->Write("ptbins");
  centbins->Write("centbins");
  
  Bool_t nextBin;
  Char_t prefix[10] = "compare";
  TCanvas * cdisplay = new TCanvas(Form("display_%s",prefix), Form("%s bg - pt bin %i - centrality bin %i", prefix, ipt,icent),600,600);
  cdisplay->Divide(2,1);

  for (Int_t icentbin=0;icentbin<nCentBins;icentbin++){
    //if only one bin selected skip the others
    if ((icent>=0) && (icentbin!=icent)) continue;
    
    for (Int_t iptbin=1;iptbin<nPtBins;iptbin++){
      //if only one bin selected skip the others      
      if ((ipt>0) && (iptbin!=ipt)) continue;
     
      Double_t lowPt=ptbins->GetBinLowEdge(iptbin+1);
      Double_t upPt=ptbins->GetBinLowEdge(iptbin+2);
      Double_t lowC=centbins->GetBinLowEdge(icentbin+1);
      Double_t upC=centbins->GetBinLowEdge(icentbin+2);
      
      TString hs1_name=Form("Data_UnlikePM_ptBin%02i_centBin%02i",iptbin,icentbin);
      TString hs2_name=Form("Data_UnlikeMP_ptBin%02i_centBin%02i",iptbin,icentbin);
      TString hb1em_name=Form("Data_MixingPM_ptBin%02i_centBin%02i",iptbin,icentbin);
      TString hb2em_name=Form("Data_MixingMP_ptBin%02i_centBin%02i",iptbin,icentbin);
      TString hb1ls_name=Form("Data_LikePP_ptBin%02i_centBin%02i",iptbin,icentbin);
      TString hb2ls_name=Form("Data_LikeMM_ptBin%02i_centBin%02i",iptbin,icentbin);
      
      TH1D * hs1_ss = (TH1D*) fss->Get(hs1_name.Data())->Clone();
      TH1D * hs1_fb = (TH1D*) ffb->Get(hs1_name.Data())->Clone();
      Printf("Signal+background SE kstar: %s",hs1_name.Data());

      TH1D * hs2_ss = (TH1D*) fss->Get(hs2_name.Data())->Clone();
      TH1D * hs2_fb = (TH1D*) ffb->Get(hs2_name.Data())->Clone();
      Printf("Signal+background SE antikstar: %s",hs2_name.Data());
      hs1_ss->Add(hs2_ss);
      hs1_fb->Add(hs2_fb);
      
      TH1D * hs_ratio;
      hs_ratio = (TH1D*) hs1_ss->Clone();
      hs_ratio->Divide(hs1_fb);
      hs_ratio->SetName(Form("ratio_Signal_ptBin%02i_centBin%02i",iptbin,icentbin));
      hs_ratio->SetTitle (Form("ratio of signal - ptBin %02i, centBin%02i",iptbin,icentbin));
      hs_ratio->Scale(ratioScaleF);
      HistoMakeUp(hs_ratio);  
      fout->cd();
      hs_ratio->Write();
     
      
      /***************************************************************************
	 EVENT MIXING BACKGROUND
	 background subtraction formula: sub = s-b
	 s = ls_pm + ls_mp       (only if K*+anti-K* is enabled, otherwise only one of the two ls)
	 b = (em_pm + em_mp)*0.5      
      ***************************************************************************/
      
      TH1D * hb1em_ss = (TH1D*) fss->Get(hb1em_name.Data())->Clone();
      TH1D * hb2em_ss = (TH1D*) fss->Get(hb2em_name.Data())->Clone();
      Printf("Background em 1: %s",hb1em_name.Data());

      TH1D * hb1em_fb = (TH1D*) ffb->Get(hb1em_name.Data())->Clone();
      TH1D * hb2em_fb = (TH1D*) ffb->Get(hb2em_name.Data())->Clone();
      Printf("Background em 2: %s",hb2em_name.Data());

      TH1D * hbem_ss =  (TH1D*) BgMixing(hb2em_ss,hb2em_ss);
      hbem_ss->SetName(Form("Mixing_ptBin%02i_centBin%02i",iptbin,icentbin));
      TH1D * hbem_fb =  (TH1D*) BgMixing(hb2em_fb,hb2em_fb);
      hbem_fb->SetName(Form("Mixing_ptBin%02i_centBin%02i",iptbin,icentbin));

      HistoMakeUp(hbem_ss);
      HistoMakeUp(hbem_fb);

      TH1D * hbem_ratio;
      hbem_ratio = (TH1D*) hbem_ss->Clone();
      hbem_ratio->Divide(hbem_fb);
      hbem_ratio->SetName(Form("ratio_bgEM_ptBin%02i_centBin%02i",iptbin,icentbin));
      hbem_ratio->SetTitle (Form("ratio of bg EM - ptBin %02i, centBin%02i",iptbin,icentbin));
      hbem_ratio->Scale(ratioScaleF);
    //normalize
      // TH1D * hbem_norm = normValuesInterval(hbem, hs, 1.3, 1.5, 1.);
      // hbem_norm->SetName(Form("norm_%s",hbem->GetName()));
      // HistoMakeUp(hbem_norm);  
      // //subtract
      // TH1D * sub_em =(TH1D*) subtractBackgnd(hs, hbem_norm);
      // sub_em->SetName(Form("sub_%s",hbem_norm->GetName()));
      // HistoMakeUp(sub_em);  
      //rebin 
      // if (isRebin) {
      // 	hs_tris->Rebin(2);
      // 	hbem_norm->Rebin(2);
      // 	sub_em->Rebin(2);
      // }
      // hs_tris->SetTitle(Form("EM bg: %4.2f<p_{t}<%4.2f GeV/c (%3.1f-%3.1f central)",lowPt,upPt,lowC,upC));
      // sub_em->SetTitle(Form("EM bg: %4.2f<p_{t}<%4.2f GeV/c (%3.1f-%3.1f central)",lowPt,upPt,lowC,upC));
      //save to file
      fout->cd();
      hbem_ratio->Write();
      // hbem_norm->Write();
      // sub_em->Write();
      //draw
      /***************************************************************************
	LIKE SIGN BACKGROUND
	background subtraction formula: sub = s-b
	s = ls_pm + ls_mp       (only if K*+anti-K* is enabled, otherwise only one of the two ls)
	b = Sqrt (ls_pp * ls_mm)
      ***************************************************************************/
           
      TH1D * hb1ls_ss = (TH1D*) fss->Get(hb1ls_name.Data())->Clone();
      TH1D * hb2ls_ss = (TH1D*) fss->Get(hb2ls_name.Data())->Clone();
      Printf("Background ls 1: %s",hb1ls_name.Data());
      
      TH1D * hb1ls_fb = (TH1D*) ffb->Get(hb2ls_name.Data())->Clone();
      TH1D * hb2ls_fb = (TH1D*) ffb->Get(hb2ls_name.Data())->Clone();
      Printf("Background ls 2: %s",hb2ls_name.Data());

      TH1D * hbls_ss = (TH1D*) BgLike(hb1ls_ss, hb2ls_ss);
      hbls_ss->SetName(Form("Like_ptBin%02i_centBin%02i",iptbin,icentbin));
      
      TH1D * hbls_fb = (TH1D*) BgLike(hb1ls_fb, hb2ls_fb);
      hbls_fb->SetName(Form("Like_ptBin%02i_centBin%02i",iptbin,icentbin));
      HistoMakeUp(hbls_ss);
      HistoMakeUp(hbls_fb);
      //TH1D * hbls_norm = (TH1D*) hbls->Clone();
      //normalize hbls
      // normValuesInterval(hbls, hs_bis, 0.7, 0.8, 1.);//hbls->Clone();
      // BestNormalization(hs_bis, hbls_norm, 0.010, 5.0, 0.01, hs_bis->GetXaxis()->FindBin(0.6), hs_bis->GetXaxis()->FindBin(1.5) );
      // hbls_norm->SetName(Form("norm_%s",hbls->GetName()));
      // HistoMakeUp(hbls);
      // HistoMakeUp(hbls_norm);
      //subtract
      // TH1D * sub_ls =(TH1D*)subtractBackgnd(hs_bis,hbls_norm);
      // sub_ls->SetName(Form("sub_%s",hbls_norm->GetName()));
      // HistoMakeUp(sub_ls);
      // //rebin 
      // if (isRebin) {
      // 	hbls->Rebin(2);
      // 	hbls_norm->Rebin(2);
      // 	sub_ls->Rebin(2);
      // }
      // hs_bis->SetTitle(Form("LS bg: %4.2f<p_{t}<%4.2f GeV/c (%3.1f-%3.1f central)",lowPt,upPt,lowC,upC));
      // sub_ls->SetTitle(Form("LS bg: %4.2f<p_{t}<%4.2f GeV/c (%3.1f-%3.1f central)",lowPt,upPt,lowC,upC));

      TH1D * hbls_ratio;
      hbls_ratio = (TH1D*) hbls_ss->Clone();
      hbls_ratio->Divide(hbls_fb);
      hbls_ratio->SetName(Form("ratio_bgLS_ptBin%02i_centBin%02i",iptbin,icentbin));
      hbls_ratio->SetTitle (Form("ratio of bg LS - ptBin %02i, centBin%02i",iptbin,icentbin));
      hbls_ratio->Scale(ratioScaleF);

      //save to file
      fout->cd();
      hbls_ratio->Write();
      // hbls_norm->Write();
      // sub_ls->Write();
      //draw
    //   switch (display){
    //   case 1:
    // 	cdisplay[1]->cd();
    // 	hs_tris->Draw();
    // 	hbls->Draw("same");
    // 	hbls_norm->Draw("same");
    // 	sub_ls->Draw("same");
    // 	break;
    //   case 2:
    // 	cdisplay[1]->cd(iptbin+1);
    // 	hs_tris->Draw();
    // 	hbls->Draw("same");
    // 	hbls_norm->Draw("same");
    // 	sub_ls->Draw("same");
    // 	break;
    //   case 3:
    //   cdisplay[1]->cd(icentbin+1);
    // 	hs_tris->Draw();
    // 	hbls->Draw("same");
    // 	hbls_norm->Draw("same");
    // 	sub_ls->Draw("same");
    // 	break;
    //   default:
    // 	break;
    //   }
    }
    
    // fout->cd();
    // if (cdisplay[0]){
    // cdisplay[0]->Write();
    // if (saveImg) cdisplay[0]->SaveAs(Form("%s_canvas.png",cdisplay[0]->GetName()));
    // }
    // if (cdisplay[1]){
    // cdisplay[1]->Write();
    // if (saveImg) cdisplay[1]->SaveAs(Form("%s_canvas.png",cdisplay[1]->GetName()));
    // }
    // cdisplay->Update();
    // cout << "Continue? (1/0)" << endl;
    // cin >> nextBin;
    //  switch (nextBin) 
    //    {
    //    case 1:
    //  	 continue;
    //    default:
    //  	 return;
    //    }
    
  }
  Printf("========================================= Succesfully saved output file %s",fout->GetName());
  return;
}

/*****************************************/
TH1D* BgLike(TH1D *hpp, TH1D *hmm)
{
   Int_t i, nbins = hpp->GetNbinsX();
   Double_t y1, y2, e1, e2, bg, err;
   
   TH1D *out = (TH1D*)hpp->Clone();
   
   out->SetName("tmp");
   out->Clear();
   out->SetEntries(hpp->GetEntries() + hmm->GetEntries());
   
   for (i = 1; i <= nbins; i++) {
      y1 = hpp->GetBinContent(i);
      y2 = hmm->GetBinContent(i);
      e1 = hpp->GetBinError(i);
      e2 = hmm->GetBinError(i);
      
      if (y1 < 1.0 || y2 < 1.0) {
         out->SetBinContent(i, 0.0);
         out->SetBinError(i, 0.0);
      } else {
         bg = TMath::Sqrt(y1*y2);
         err = y1*y1*e2*e2 + y2*y2*e1*e1;
         err /= y1*y2;
         out->SetBinContent(i, bg);
         out->SetBinError(i, TMath::Sqrt(err));
      }
   }
   
   return out;
}

/*****************************************/
TH1D* SumSignal(TH1D *hpm, TH1D *hmp)
{
   Int_t i, nbins = hpm->GetNbinsX();
   Double_t y1, y2, e1, e2, sig, err;   
   TH1D *out = (TH1D*)hpm->Clone();  
   out->SetName("tmp");
   out->Clear();
   out->SetEntries(hpm->GetEntries() + hmp->GetEntries());   
   for (i = 1; i <= nbins; i++) {
      y1 = hpm->GetBinContent(i);
      y2 = hmp->GetBinContent(i);
      e1 = hpm->GetBinError(i);
      e2 = hmp->GetBinError(i);
      if (y1 < 1.0 || y2 < 1.0) {
         out->SetBinContent(i, 0.0);
         out->SetBinError(i, 0.0);
      } else {
         sig = (y1+y2)*0.5;
         err = e2*e2 + e1*e1;
         out->SetBinContent(i, sig);
         out->SetBinError(i, TMath::Sqrt(err));
      }
   }   
   return out;
}
/*****************************************/
TH1D* BgMixing(TH1D *hpm, TH1D *hmp)
{
   Int_t i, nbins = hpm->GetNbinsX();
   Double_t y1, y2, e1, e2, bg, err;   
   TH1D *out = (TH1D*)hpm->Clone();  
   out->SetName("tmp");
   out->Clear();
   out->SetEntries(hpm->GetEntries() + hmp->GetEntries());   
   for (i = 1; i <= nbins; i++) {
      y1 = hpm->GetBinContent(i);
      y2 = hmp->GetBinContent(i);
      e1 = hpm->GetBinError(i);
      e2 = hmp->GetBinError(i);
      if (y1 < 1.0 || y2 < 1.0) {
         out->SetBinContent(i, 0.0);
         out->SetBinError(i, 0.0);
      } else {
         bg = (y1+y2)*0.5;
         err = e2*e2 + e1*e1;
         out->SetBinContent(i, bg);
         out->SetBinError(i, TMath::Sqrt(err));
      }
   }   
   return out;
}

/*****************************************/
void HistoMakeUp(TH1D*histo){
  if (!histo) return;
  TString hname = Form("%s",histo->GetName());

  // histo->Rebin(2);

  TAxis * xaxis = (TAxis*)histo->GetXaxis();
  xaxis->SetTitle("M_{inv} (GeV/c^{2})");

  TAxis * yaxis = (TAxis*)histo->GetYaxis();
  yaxis->SetTitle("pairs");

  Color_t color;
  if (hname.Contains("Signal") 
      || hname.Contains("UnlikeMP")
      || hname.Contains("UnlikePM") )
    color = kBlack;
  if (hname.Contains("Mixing"))
    color = kRed;
  if (hname.Contains("Like"))
    color = kGreen+1;
  if (hname.Contains("Like") && hname.Contains("norm"))
    color = kAzure-4;
  if (hname.Contains("Like") && hname.Contains("sub"))
    color = kBlue;
  if (hname.Contains("Mixing") && hname.Contains("sub"))
    color = kMagenta;
  
  histo->SetLineColor(color);
  histo->SetLineWidth(2);
  histo->SetMarkerColor(color);
  return;
}
