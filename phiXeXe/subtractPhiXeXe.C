/*
fbellini@cern.ch - 25/10/2012

Starting from projected histogram per centrality bins and pt bins:
- normalize background
- compute signal (K*+anti-k*), LS and EM background
- subtract bg
- display

*/

#include "/Users/fbellini/alice/macros/SetStyle.C"
#include "/Users/fbellini/alice/macros/MakeUp.C"

const Double_t kBigNumber=1E10;
const Double_t kSmallNumber=1E-10;

TH1D*    BgLike(TH1D* hpp = NULL, TH1D* hmm = NULL);
TH1D*    SumSignal(TH1D* hpm = NULL, TH1D* hmp = NULL);
TH1D*    BgMixing(TH1D *hpm = NULL, TH1D* hmp = NULL);
void     HistoMakeUp(TH1D* histo = NULL, Int_t markerStyle=-1);
Double_t BestNormalization(TH1D *hSig, TH1D *hBg, Double_t normMin, Double_t normMax, Double_t normStep, Int_t firstBin, Int_t lastBin);
TH1D*    normalize2Integral(TH1D* hist=NULL, TH1D* histRef=NULL, Double_t factor=1.);
Float_t  GetRangeValuesNormalizationFactor(TH1D* hist=NULL, TH1D* histRef=NULL, Double_t valueMin=1.3, Double_t valueMax=1.5);
TH1D*    normValuesInterval(TH1D* hist=NULL, TH1D* histRef=NULL, Double_t valueMin=1.3, Double_t valueMax=1.5, Double_t factor=1., Double_t *outputFactor = 0);
TH1D*    subtractBackgnd(TH1D* hist=NULL, TH1D* hist2=NULL);


enum EHistStyle {kSig = 0,
		 kLSBPP,
		 kLSBMM,
		 kLSB,
		 kLSBnorm,
		 kLSBsub,
		 kMEB,
		 kMEBnorm,
		 kMEBsub};

Int_t subtractPhiXeXe(	 TString projectionFile="proj_A3.root", 
			 Float_t emNormInf = 1.05,
			 Float_t emNormSup = 1.15, 
			 Double_t desiredIMbinWidth = 0.002, //in GeV/cˆ2
			 Int_t ipt = -1,
			 Int_t istartpt = 1,
			 Int_t icent = -1,
			 Float_t scaleEMNormFactor = 1.0,
		     	 Bool_t doNormLS = 0, 
			 Bool_t enableCanvas = 1,
			 Bool_t isDisplayDist = 1,
			 Bool_t displayEMonly = 0)
{
  //general - set style
   // initial setup
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTextFont(42);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(1);
  gStyle->SetLabelOffset(0.005,"yx");
  gStyle->SetLabelSize(0.09,"xyz");
  gStyle->SetTitleSize(0.09,"xyz");
  gStyle->SetTitleOffset(1.1,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetEndErrorSize(0); //sets in #of pixels the lenght of the tick at the end of the error bar
  TGaxis::SetMaxDigits(2); // gStyle->SetTitleAlign(33);
  // gStyle->SetTitleX(.95);
  // gStyle->SetTitleY(.99);
  
  //histo make up settings
  Float_t Marker_Size = .8;
  Float_t Marker_Style[] = {    20,      24,      25,     21,     21,     25,       28,      28,     20};
  Color_t color[]        = {kBlack, kOrange, kCyan+2, kRed-1, kRed+1, kGray+2, kAzure-7, kBlue+1, kBlack};

  // open input file
  TFile * f = TFile::Open(projectionFile.Data());
  if (!f) return 1;

  TList * lUnlikePM = (TList *) f->Get("hUnlikePM");
  TList * lLikePP   = (TList *) f->Get("hLikePP");
  TList * lLikeMM   = (TList *) f->Get("hLikeMM");
  TList * lMixingPM = (TList *) f->Get("hMixingPM");
  
  // get axes
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  const Int_t nPtBins = ptbins->GetNbins();
  TAxis *centbins = (TAxis*)f->Get("centbins");
  const Int_t nCentBins = centbins->GetNbins();

  Bool_t nextBin = kFALSE;

  //---------------------
  // create output
  //---------------------
  TString folderName = Form("norm%3.2f-%3.2f", emNormInf,emNormSup);
  gSystem->Exec(Form("mkdir %s", folderName.Data()));
  TString fout_name = projectionFile;
  fout_name.ReplaceAll("proj_", "sub_");
  fout_name.Prepend(Form("%s/", folderName.Data()));
  
  TFile * fout = new TFile(fout_name.Data(),"recreate");
  ptbins->Write("ptbins");
  centbins->Write("centbins");
  
  Float_t rangeInf = emNormInf, rangeSup = emNormSup;
  Double_t normfactor_em[1],normfactor_ls[1];
  Float_t ptinf=0.0, ptsup=0.0, centinf=0.0, centsup=0.0;
  Int_t treept=-1, treecent=-1;
  
  TTree *ntree = new TTree("ntree","normFactors");
  ntree->Branch("factor_em", &normfactor_em[0], "factor/D");
  ntree->Branch("factor_ls", &normfactor_ls[0], "factor/D");
  ntree->Branch("inf", &rangeInf, "inf/F");
  ntree->Branch("sup", &rangeSup, "sup/F");
  ntree->Branch("pt_inf", &ptinf, "pt_inf/F");
  ntree->Branch("pt_sup", &ptsup, "pt_sup/F");
  ntree->Branch("pt_bin", &treept, "pt_bin/I");
  ntree->Branch("cent_inf", &centinf, "cent_inf/F");
  ntree->Branch("cent_sup", &centsup, "cent_sup/F");
  ntree->Branch("cent_bin", &treecent, "cent_bin/I");

  //canvas for display
  TCanvas * cdisplay[2][nCentBins];
  Char_t prefix[2][5]={"sub","dist"};

  for (Int_t cc=0; cc<nCentBins; cc++){ 
    for (Int_t j=0; j<2; j++){ 
      cdisplay[j][cc] = NULL; 
      if (!enableCanvas) continue;
    
      //display all cent and all pt
      if (icent<0 && ipt<0) {
	cdisplay[j][cc] = new TCanvas(Form("cent%i_%s",cc,prefix[j]), Form("%sB - cent bin %i", prefix[j], cc), 800, 1500);
	cdisplay[j][cc]->Divide(3,4);
      }
      
      //display selected cent and all pt
      if (icent>=0 && cc==icent) {
	if (ipt<0) {
	  cdisplay[j][cc] = new TCanvas(Form("cent%i_%s",icent, prefix[j]), Form("%sB - cent bin %i", prefix[j], icent), 1200, 900);
	  cdisplay[j][cc]->Divide(3,4);
	} else {
	  cdisplay[j][cc] = new TCanvas(Form("cent%i_%s",icent, prefix[j]), Form("%sB - cent bin %i", prefix[j], icent), 800, 600);
	}
      }
      
      //display selected pt for all cent
      if (ipt>=0) {
	if (icent<0) {
	  cdisplay[j][cc] = new TCanvas(Form("cent%i_%s",icent, prefix[j]), Form("%sB - cent bin %i", prefix[j], icent), 800, 600);
	} else {
	  if (cc==icent){
	    cdisplay[j][cc] = new TCanvas(Form("cent%i_%s",icent, prefix[j]), Form("%sB - cent bin %i", prefix[j], icent), 800, 600);
	  }
	}
      } //if pt
    }// loop on bg
  }// loop on cent


  
  for (Int_t icentbin=0;icentbin<nCentBins;icentbin++){
    //if only one bin selected skip the others
    if ((icent>=0) && (icentbin!=icent)) continue;

    TString pngSaveName = fout_name;
    pngSaveName.ReplaceAll("sub_", Form("sub_c%i_", icentbin));
    pngSaveName.ReplaceAll(".root","");
  
    for (Int_t iptbin=istartpt; iptbin<nPtBins; iptbin++){
      //if only one bin selected skip the others      
      if ((ipt>0) && (iptbin!=ipt)) continue;
      
      Printf("*************************************************");
      Printf("************** cent %i - pt %i ******************", icentbin, iptbin);
      Printf("*************************************************");

      Double_t lowPt = ptbins->GetBinLowEdge(iptbin+1);
      Double_t upPt = ptbins->GetBinLowEdge(iptbin+2);
      Double_t lowC = centbins->GetBinLowEdge(icentbin+1);
      Double_t upC = centbins->GetBinLowEdge(icentbin+2);
  
      //assign variables for norm factors tree
      ptinf=lowPt;
      ptsup=upPt;
      centinf=lowC;
      centsup=upC;
      treept=iptbin; 
      treecent=icentbin;

      //get input histos
      TString hs1_name   = Form("hUnlikePM_ptBin%02i_centBin%02i",iptbin,icentbin);
      TString hb1em_name = Form("hMixingPM_ptBin%02i_centBin%02i",iptbin,icentbin);
      TString hb1ls_name = Form("hLikePP_ptBin%02i_centBin%02i",iptbin,icentbin);
      TString hb2ls_name = Form("hLikeMM_ptBin%02i_centBin%02i",iptbin,icentbin);
    
      // signal
      TH1D * hs = (TH1D*) ((TH1D*) lUnlikePM->FindObject(hs1_name.Data()))->Clone(Form("Signal_ptBin%02i_centBin%02i",iptbin,icentbin));
      hs->SetTitle(Form("%4.2f < #it{p}_{T} < %5.2f GeV/#it{c} (%2.0f-%2.0f%%); #it{M}_{KK} (GeV/#it{c}); counts / (%4.3f GeV/#it{c})",lowPt,upPt,lowC,upC, desiredIMbinWidth));
      if (hs) Printf("Signal+background SE phi: %s",hs->GetName());
      else return 2;
      
      Beautify(hs, color[EHistStyle::kSig], 1, 2, Marker_Style[EHistStyle::kSig], Marker_Size);
      TH1D * hs_bis = (TH1D*) hs->Clone();
      hs_bis->SetTitle(Form("%4.2f < #it{p}_{T} < %5.2f GeV/#it{c} (%2.0f-%2.0f%%); #it{M}_{KK} (GeV/#it{c}); counts / (%4.3f GeV/#it{c})",lowPt,upPt,lowC,upC, desiredIMbinWidth));
      TH1D * hs_tris = (TH1D*) hs->Clone();
      hs_tris->SetTitle(Form("%4.2f < #it{p}_{T} < %5.2f GeV/#it{c} (%2.0f-%2.0f%%); #it{M}_{KK} (GeV/#it{c}); counts / (%4.3f GeV/#it{c})",lowPt,upPt,lowC,upC, desiredIMbinWidth));
      // mixed-event background
      TH1D * hbem = (TH1D*) ((TH1D*) lMixingPM->FindObject(hb1em_name.Data()))->Clone(Form("Mixing_ptBin%02i_centBin%02i",iptbin,icentbin));
      hbem->SetTitle(Form("EM bg: %4.2f<#it{p}_{T}<%5.2f GeV/#it{c} (%2.0f-%2.0f%%)",lowPt,upPt,lowC,upC));
      if (hbem) Printf("EM background for phi: %s",hbem->GetName());
      else return 3;
      Beautify(hbem, color[EHistStyle::kMEB], 1, 2, Marker_Style[EHistStyle::kMEB], Marker_Size);
      
      // like-sign background
      TH1D * hb1ls = (TH1D*) ((TH1D*) lLikePP->FindObject(hb1ls_name.Data()))->Clone(Form("LikePP_ptBin%02i_centBin%02i",iptbin,icentbin));
      TH1D * hb2ls = (TH1D*) ((TH1D*) lLikeMM->FindObject(hb2ls_name.Data()))->Clone(Form("LikeMM_ptBin%02i_centBin%02i",iptbin,icentbin));
      TH1D * hbls = (TH1D*) BgLike(hb1ls, hb2ls);
      hbls->SetNameTitle(Form("Like_ptBin%02i_centBin%02i",iptbin,icentbin), 
			 Form("LS bg: %4.2f<#it{p}_{T}<%5.2f GeV/#it{c} (%2.0f-%2.0f%%)",lowPt,upPt,lowC,upC));
      if (hbls) Printf("LS background for phi: %s",hbls->GetName());
      else return 4;
      Beautify(hb1ls, color[EHistStyle::kLSBPP], 1, 1, Marker_Style[EHistStyle::kLSBPP], Marker_Size);
      Beautify(hb2ls, color[EHistStyle::kLSBMM], 1, 1, Marker_Style[EHistStyle::kLSBMM], Marker_Size);
      Beautify(hbls, color[EHistStyle::kLSB], 1, 1, Marker_Style[EHistStyle::kLSB], Marker_Size);

      Int_t inorm1 =  hbem->GetXaxis()->FindBin(emNormInf);
      Int_t inorm2 =  hbem->GetXaxis()->FindBin(emNormSup);
      
      //MEB make up
      TH1D * hbem_norm;
      if ((emNormInf<=0) || (emNormSup<=0)) {
	hbem_norm = (TH1D*) hbem->Clone(Form("norm_%s", hbem->GetName()));
	Double_t norm_best_factor = BestNormalization(hs, hbem_norm, 0.005, 0.500, 0.001, inorm1, inorm2);
	hbem_norm->Scale(norm_best_factor);
	Printf(":::: EM Best normalization by factor = %5.3f", norm_best_factor);
	normfactor_em[0]=norm_best_factor;
      } else {
	hbem_norm = (TH1D*) normValuesInterval(hbem, hs, emNormInf, emNormSup, scaleEMNormFactor, normfactor_em); 
	Printf(":::: EM Normalization in [%3.2f,%3.2f] by factor = %5.3f", emNormInf, emNormSup, normfactor_em[0]);
      }
      hbem_norm->SetName(Form("norm_%s", hbem->GetName()));
      Beautify(hbem_norm, color[EHistStyle::kMEBnorm], 1, 2, Marker_Style[EHistStyle::kMEBnorm], Marker_Size);

      //LSB make up
      TH1D * hbls_norm = (TH1D*) hbls->Clone(Form("norm_%s", hbls->GetName()));
      if (doNormLS) {
	if ((emNormInf<=0) || (emNormSup<=0)) {
	  Double_t norm_best_factor = BestNormalization(hs, hbls_norm, 0.005, 0.500, 0.001, inorm1, inorm2);
	  hbls_norm->Scale(norm_best_factor);
	  Printf(":::: LS Best normalization by factor = %5.3f", norm_best_factor);
	  normfactor_ls[0] = norm_best_factor;      
	} else {
	  hbls_norm = (TH1D*) normValuesInterval(hbls, hs, emNormInf, emNormSup, scaleEMNormFactor, normfactor_ls); 
	  Printf(":::: LS Normalization in [%3.2f,%3.2f] by factor = %5.3f", emNormInf, emNormSup, normfactor_ls[0]);
	}
	hbls_norm->SetName(Form("norm_%s", hbls->GetName()));
      }
      Beautify(hbls_norm, color[EHistStyle::kLSBnorm], 1, 2, Marker_Style[EHistStyle::kLSBnorm], Marker_Size);

      //get invariant mass binning and rebin if requested
      Double_t invMassBinWidth = ((TAxis*) hs->GetXaxis())->GetBinWidth(1);
      Int_t rebinFactor = 2; //desiredIMbinWidth / invMassBinWidth;
      
      if (rebinFactor>0) {
	Printf(":::: REBINNING: Minv bin width: %4.3f --> %4.3f via rebin by factor %i", invMassBinWidth, desiredIMbinWidth, rebinFactor);
	hs->Rebin(rebinFactor);
	hs_bis->Rebin(rebinFactor);
	hs_tris->Rebin(rebinFactor);
	hbem->Rebin(rebinFactor);
	hbem_norm->Rebin(rebinFactor);
	hbls->Rebin(rebinFactor);
	hbls_norm->Rebin(rebinFactor);
      }

      //draw
      if (isDisplayDist && cdisplay[0][icentbin]){
	if (icent>0 && ipt>0) cdisplay[0][icentbin]->cd();
	else cdisplay[0][icentbin]->cd(iptbin+1-istartpt);
	hs_tris->Draw("same");
	hs_tris->SetTitleSize(18);
	hs_tris->GetYaxis()->SetTitleSize(0.07);
	hs_tris->GetYaxis()->SetTitleOffset(0.87);
	hs_tris->GetYaxis()->SetLabelSize(0.06);
	hs_tris->GetYaxis()->SetNdivisions(509);
	hs_tris->GetXaxis()->SetTitleSize(0.07);
	hs_tris->GetXaxis()->SetLabelSize(0.06);
	hbls_norm->Draw("hist same");
	hbem_norm->Draw("hist same");
      }


      /***************************************************************************
	 EVENT MIXING BACKGROUND
	 background subtraction formula: sub = s-b
	 s = ls_pm + ls_mp      
	 b = (em_pm + em_mp)*0.5      
      ***************************************************************************/
      TH1D * sub_em =(TH1D*) subtractBackgnd(hs, hbem_norm);
      if (!sub_em) return 5;
      sub_em->SetName(Form("sub_%s",hbem_norm->GetName()));
      sub_em->SetTitle(Form("MEB, %4.2f < #it{p}_{T} < %5.2f GeV/#it{c} (%2.0f-%2.0f%%); #it{M}_{KK} (GeV/#it{c}); counts / (%4.3f GeV/#it{c})",lowPt,upPt,lowC,upC, desiredIMbinWidth));
      Beautify(sub_em, color[EHistStyle::kMEBsub], 1, 1, Marker_Style[EHistStyle::kMEBsub], Marker_Size);
      SetLabels(sub_em, 0.05, 0.05, 1.1, 1.1);
      //fill normalization tree
      ntree->Fill();


      /***************************************************************************
	LIKE SIGN BACKGROUND
	background subtraction formula: sub = s-b
	s = ls_pm + ls_mp   
	b = Sqrt (ls_pp * ls_mm)
      ***************************************************************************/      
      //subtract
      TH1D * sub_ls = (TH1D*) subtractBackgnd(hs_bis, hbls_norm);
      if (!sub_ls) return 6;
      sub_ls->SetName(Form("sub_%s",hbls_norm->GetName()));
      sub_ls->SetTitle(Form("%4.2f < #it{p}_{T} < %5.2f GeV/#it{c} (%2.0f-%2.0f%%); #it{M}_{KK} (GeV/#it{c}); counts / (%4.3f GeV/#it{c})",lowPt,upPt,lowC,upC, desiredIMbinWidth));
      Beautify(sub_ls, color[EHistStyle::kLSBsub], 1, 1, Marker_Style[EHistStyle::kLSBsub], Marker_Size);
      SetLabels(sub_ls, 0.05, 0.05, 1.1, 1.1);

      sub_ls->GetXaxis()->SetRangeUser(sub_ls->GetXaxis()->GetBinLowEdge(1), 1.1);
      sub_ls->GetYaxis()->SetRangeUser(sub_ls->GetMinimum()*1.1, sub_ls->GetMaximum()*1.5);
      sub_em->GetXaxis()->SetRangeUser(sub_em->GetXaxis()->GetBinLowEdge(1), 1.1);

      sub_ls->GetYaxis()->SetTitleSize(0.07);
      sub_ls->GetYaxis()->SetTitleOffset(0.87);
      sub_ls->GetYaxis()->SetLabelSize(0.06);
      sub_ls->GetYaxis()->SetNdivisions(509);
      sub_ls->GetXaxis()->SetTitleSize(0.07);
      sub_ls->GetXaxis()->SetLabelSize(0.06);
	
      //draw
      if (cdisplay[1][icentbin]){
	if (icent>0 && ipt>0) cdisplay[1][icentbin]->cd();
	else cdisplay[1][icentbin]->cd(iptbin+1-istartpt);
	if (!displayEMonly) sub_ls->Draw("same");
	sub_em->Draw("same");
      }
      
      //save to file
      fout->cd();
      hs->Write();
      hbem->Write();
      hbls->Write();
      hbem_norm->Write();
      hbls_norm->Write();
      sub_em->Write();
      sub_ls->Write();
          
      /*
      cout << "Continue? (1/0)" << endl;
      cin >> nextBin;
      switch (nextBin) 
	{
	case 1:
	  continue;
	default:
	  return 0;
	}
      */
    } //loop over pt
    fout->cd();
    cdisplay[1][icentbin]->Write();
    cdisplay[1][icentbin]->Print(Form("%s.png", pngSaveName.Data()));
    cdisplay[1][icentbin]->Print(Form("%s.eps", pngSaveName.Data()));
    cdisplay[0][icentbin]->Write();
    pngSaveName.ReplaceAll("sub_", "dist_");
    cdisplay[0][icentbin]->Print(Form("%s.png",pngSaveName.Data()));
    cdisplay[0][icentbin]->Print(Form("%s.eps", pngSaveName.Data()));
  }//loop over centrality

  fout->cd();
  ntree->Write();
  Printf(":::: Succesfully saved output file %s",fout->GetName());
  return 0;
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
	bg =  2.0*TMath::Sqrt(y1*y2);//multiply by two if signal is not divided by two bg = 2.0 * TMath::Sqrt(y1*y2);
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
	sig = (y1+y2); //divide by 2 if kHalfSum: done in the main function
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
         bg = (y1+y2);
         err = e2*e2 + e1*e1;
         out->SetBinContent(i, bg);
         out->SetBinError(i, TMath::Sqrt(err));
      }
   }   
   return out;
}

/*****************************************/
void HistoMakeUp(TH1D*histo, Int_t markerStyle)
{
  if (!histo) return;
  TString hname = Form("%s",histo->GetName());

  // histo->Rebin(2);

  TAxis * xaxis = (TAxis*)histo->GetXaxis();
  xaxis->SetTitle("#it{M}_{KK} (GeV/#it{c}^{2})");

  TAxis * yaxis = (TAxis*)histo->GetYaxis();
  yaxis->SetTitle("Counts / (0.001 GeV/#it{c}^{2})");

  Color_t color;
  if (hname.Contains("Signal") || hname.Contains("UnlikePM")){
    color = kBlack;
    markerStyle = 20;
  }
  if (hname.Contains("Mixing")){
    color = kOrange-3;
    markerStyle = 24;
  }
  if (hname.Contains("Like")) {
    color = kAzure+1;
    markerStyle = 25;
  }
  if (hname.Contains("Mixing") && hname.Contains("norm")) {
    color = kRed;
    markerStyle = 24;
  }
  if (hname.Contains("Like") && hname.Contains("norm")) {
    color = kBlue;
    markerStyle = 25;
  }
  if (hname.Contains("Like") && hname.Contains("sub")) {
    color = kBlack;
    markerStyle = 20;
  }
  if (hname.Contains("Mixing") && hname.Contains("sub")) {
    color = kGray;
    markerStyle = 2;
  }
  histo->SetLineColor(color);
  histo->SetLineWidth(1);
  histo->SetMarkerColor(color);
  histo->SetFillColor(kWhite);
  histo->SetFillStyle(0);
  if (markerStyle>=0) histo->SetMarkerStyle(markerStyle);
 
  return;
}

//-----------------------------------------------------------------------
Double_t BestNormalization(TH1D *hSig, TH1D *hBg, Double_t normMin, Double_t normMax, Double_t normStep, Int_t firstBin, Int_t lastBin)
{
//
// search best factor which leaves all subtractions positive within error
//
   Double_t norm, bestNorm = normMin;
   Double_t min, bestMin  = 1E20;
   TH1D *htmp = NULL;
   for (norm = normMax; norm >= normMin; norm -= normStep) {
      htmp = (TH1D*)hSig->Clone("tmp");
      htmp->Add(hBg, -norm);
      min = 1E20;
      for (Int_t k = firstBin; k <= lastBin; k++) {
         Double_t y  = htmp->GetBinContent(k);
         if (htmp->GetBinError(k) < 1E-6) continue;
         y -= htmp->GetBinError(k);
         if (y < min) min = y;
      }
      if (min < 0.0)  
         continue;
      else if (min < bestMin) {
         bestMin  = min;
         bestNorm = norm;
      }
   }
   Double_t IntS = htmp->Integral(firstBin, lastBin);
   Double_t IntB = hBg->Integral(firstBin, lastBin);
   Printf(":::: Normalization: best -> intS = %f, intB = %f, norm = %8.4f", IntS, IntB, bestNorm);
   return bestNorm;
}


//---------------------------------------------------------------------------------
TH1D* normalize2Integral(TH1D* hist, TH1D*histRef, Double_t factor)
{
  
  if ((!hist)||(!histRef)){
    printf("normalize2Integral - ERROR: invalid histogram address passed to normalization function\n");
    return 0;
  }
  
  Int_t nbins=hist->GetXaxis()->GetNbins();
  Int_t nbinsRef=histRef->GetXaxis()->GetNbins();
  
  if (nbins!=nbinsRef){
    printf("normalize2Integral - ERROR: histgrams have different binning. Doing nothing.\n");
    return 0;    
  }
  
  TH1D * cloneH = (TH1D*) hist->Clone();
  Double_t integralS= histRef->Integral();
  Double_t integralB= hist->Integral();
  
  if (integralS==0){
    printf("normalize2Integral - INFO: signal + backgnd histogram integral = 0. Skipping normalization.\n");
    return 0;
  }
  
  if (integralB==0){
    printf("normalize2Integral - INFO: backgnd histogram integral = 0. Skipping normalization.\n");
    return 0;
  }
  
  if (factor<=0) factor=1.;
  Double_t normFactor= factor*integralS/integralB;
  cloneH->Scale(normFactor);
  return cloneH;
}

//---------------------------------------------------------------------------------
Float_t GetRangeValuesNormalizationFactor(TH1D* hist, TH1D* histRef, Double_t valueMin, Double_t valueMax){
  
  /*
    estimate normalization factor to be used on hist
    by using histRef as reference 
    in *values* range given by valueMin - valueMax
    
  */

  if ((!hist)||(!histRef)){
    printf("GetRangeValuesNormalizationFactor - ERROR: invalid histgram address passed to normalization function\n");
    return 0;
  }
  
  Int_t nbins=hist->GetXaxis()->GetNbins();
  Int_t nbinsRef=histRef->GetXaxis()->GetNbins();
  
  if (nbins!=nbinsRef){
    printf("GetRangeValuesNormalizationFactor - ERROR: histgrams have different binning. Doing nothing.\n");
    return 0;    
  }
  
  Int_t ibinmin,ibinmax;
  Double_t invMassBinWidth = ((TAxis*) hist->GetXaxis())->GetBinWidth(1);
  
  if ((valueMin<=kSmallNumber) || (valueMin>=valueMax)) {
    valueMin = hist->GetXaxis()->GetBinLowEdge(1);
    ibinmin=1;
    printf("GetRangeValuesNormalizationFactor - INFO:  min value used for normalization range as default.\n");
  } else {
    ibinmin= 1 + ( valueMin-(hist->GetXaxis()->GetBinLowEdge(1)) )/invMassBinWidth;
  } 
  
  if ((valueMax>=kBigNumber) || (valueMax<valueMin)) {
    valueMax = hist->GetXaxis()->GetBinUpEdge(nbins);
    ibinmax=nbins;
    printf("GetRangeValuesNormalizationFactor - INFO:  max value used for normalization range as default.\n");  
  } else {    
    ibinmax = nbins - ( (hist->GetXaxis()->GetBinLowEdge(nbins)) - valueMax )/invMassBinWidth;  
  }
  
  // printf("--------------- Normalizing histogram in interval [ %5.2f - %5.2f ]\n",valueMin,valueMax);
  //  printf("                corresponding to (x-axis) bin interval [ %i - %i ]\n",ibinmin,ibinmax);
  
  Double_t integralS= histRef->Integral(ibinmin,ibinmax);
  Double_t integralB= hist->Integral(ibinmin,ibinmax);
  
  if (integralS==0){
    printf("GetRangeValuesNormalizationFactor - INFO: reference histogram integral = 0. Skipping normalization.\n");
    return 0;
  }
  if (integralB==0){
    printf("GetRangeValuesNormalizationFactor - INFO: histogram to be normalized has integral = 0. Skipping normalization.\n");
    return 0;
  }
  
  Double_t normFactor= integralS/integralB;
  //  printf("GetRangeValuesNormalizationFactor - INFO: Normalization factor estimate = %6.3f\n",normFactor);
  return normFactor;
  
}

//---------------------------------------------------------------------------------
TH1D* normValuesInterval(TH1D* hist, TH1D* histRef, Double_t valueMin, Double_t valueMax, Double_t factor, Double_t * outputFactor)
{
  
  if ((!hist)||(!histRef)){
    printf("normValuesInterval - ERROR: invalid histgram address passed to normalization function\n");
    return 0;
  }
    printf("normValuesInterval - INFO: Normalizing background in IM interval [ %5.2f - %5.2f ]\n",valueMin,valueMax);
  
  TH1D * cloneH = (TH1D*) hist->Clone();  
  Double_t normFactor= GetRangeValuesNormalizationFactor(hist,histRef,valueMin,valueMax);
  printf("normValuesInterval: normalized histo %s with norm factor = (norm)%6.3f * (custom)%6.3f = (tot) %6.3f\n",cloneH->GetName(),normFactor,factor,normFactor*factor);  
  if (factor<=0) factor=1.;
  normFactor=normFactor*factor;
  cloneH->Scale(normFactor);
  outputFactor[0] = normFactor; 
  return cloneH;
  
}

//---------------------------------------------------------------------------------
TH1D* subtractBackgnd(TH1D* hist, TH1D* hist2)
{
  //hist - hist2
  if ((!hist)||(!hist2)){
    printf("subtractBackgnd - ERROR: invalid histogram address passed to base subtraction function\n");
    return 0x0;
  }
  
  Int_t nbins=hist->GetXaxis()->GetNbins();
  Int_t nbinsRef=hist2->GetXaxis()->GetNbins();
  
  if (nbins!=nbinsRef){
    printf("subtractBackgnd - ERROR: histograms have different binning. Doing nothing.\n");
    return 0x0;    
  }
  
  hist->Add(hist2,-1);
  return hist;
}
