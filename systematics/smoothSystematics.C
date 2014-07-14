void smoothSystematics(Bool_t plotSys=0,
		       TString suffix = "Fcn",
		       TString path = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/systUncert/", 
		       TString fileprefix = "sysFcn_LS_",
		       TString MBpath = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/0to100_binB/systUncert/",
		       TString MBfile = "0to100_smooth_Fcn_systematics_.root",
		       TString filebins = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/proj_2424_tpc2s_tof3sveto.root") 
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  Float_t maxYsys = 15.0;
  Color_t color[8]={ kRed+1, kPink+6, kGreen+1, kAzure+1, kBlue+3, kBlack, kMagenta, kOrange};  
  Color_t colorsys[9]={ kRed, kBlue, kGreen+1, kMagenta+1, kOrange, kCyan+2, kPink+2, kGray+2, kViolet+1};  
  //combined
  // kTeal+3, kSpring+5, kBlue-3, kCyan-3, kAzure-4, kBlack, //tpc
  // kRed+2, kOrange+6, kViolet-6, kMagenta, kBlue+2, kBlack}; //tof
  
  //get bins
  Int_t npt_axis = 0, ncent_axis=0; 
  TFile *f=TFile::Open(filebins.Data());
  if (!ptbins || !centbins) return;
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  npt_axis = ptbins->GetNbins();  
  TAxis *centbins = (TAxis*)f->Get("centbins");
  ncent_axis = centbins->GetNbins();
  f->Close();
  const Int_t npt = npt_axis;
  const Int_t ncent = ncent_axis;
  Double_t pt[npt+1];
  for (Int_t k=0; k<npt+1;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
  }
  Double_t cent[ncent+1]; 
  for (Int_t k=0; k<ncent+1;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
  }

  //access minimum bias file
  TFile * fin100 = TFile::Open(Form("%s/%s", MBpath.Data(), MBfile.Data()));
  if (!fin100) Printf("cannot open min bias file.");
  TString histSystName = "hSystVsPtPercentageOfCentral";
  TString MBhistSystName;
  if (suffix.Contains("Norm")||suffix.Contains("PID")) MBhistSystName = "hSystVsPtPercentageOfCentral_0_smooth2_iter1";
  else MBhistSystName = "hSystVsPtPercentageOfCentral_100_smooth2_iter1"; ;
  TH1D * hsys100 = (TH1D*) fin100->Get(Form("%s", MBhistSystName.Data()));//hSystVsPtPercentageOfCentral_100_smooth2_iter1
  if (!hsys100) Printf("cannot find min bias syst.");
  hsys100->SetLineColor(color[5]);
  hsys100->SetMarkerColor(color[5]);
  hsys100->GetYaxis()->SetRangeUser(0.0, maxYsys);
  hsys100->SetTitle(Form("%s sys., min. bias", suffix.Data()));

  //access multiplicity dependent files
  //Printf("Pt bins %i", nbins-1);
  TFile * fin[5];
  TH1D * hsys[5];
  TH1D * hratiosys[5];
  TH1D * hsys_smoothed[5];
  TH1D * hsys_smoothed2[5];
  TH1D * hsys_smoothed3[5];
  TH1D * hsys_smoothed4[5];
  TH1D * hsys_smoothed5[5];
  TH1D * hsys_smoothed6[5];
  TH1D * hsys_smoothed7[5];
  TH1D * hsys_smoothed8[5];
  TH1D * hsys_smoothed9[5];

  //create canvases
  TCanvas * csys = new TCanvas("csys","systematics", 800, 800); 
  csys->cd();
  hsys100->Draw();
  TCanvas * csyspt = new TCanvas("csyspt","systematics per pt bin", 900, 900); 
  csyspt->Divide(7, 3);
  TCanvas * csysratio = new TCanvas("csysratio","systematics ratio to min bias", 900, 500); 
  
  TCanvas * csyssmooth1 = new TCanvas("csyssmooth1","systematics - smoothing #1", 800, 800); 
  TCanvas * csyssmooth2 = new TCanvas("csyssmooth2","systematics - smoothing #2", 800, 800); 
  TCanvas * csyssmooth3 = new TCanvas("csyssmooth3","systematics - smoothing #3", 800, 800); 
  TCanvas * csyssmooth4 = new TCanvas("csyssmooth4","systematics - smoothing #4", 800, 800); 
  TCanvas * csyssmooth5 = new TCanvas("csyssmooth5","systematics - smoothing #5", 800, 800); 
  TCanvas * csyssmooth6 = new TCanvas("csyssmooth6","systematics - smoothing #6", 800, 800); 
  TCanvas * csyssmooth7 = new TCanvas("csyssmooth7","systematics - smoothing #7", 800, 800); 
  TCanvas * csyssmooth8 = new TCanvas("csyssmooth8","systematics - smoothing #8", 800, 800); 
  TCanvas * csyssmooth9 = new TCanvas("csyssmooth9","systematics - smoothing #9", 800, 800); 
  
  TCanvas * csyssmooth[6];
  for (Int_t ic=0;ic<ncent;ic++){
    csyssmooth[ic] = new TCanvas(Form("csyssmooth_%i",ic),Form("systematics multi bin %i",ic), 800, 800); 
    csyssmooth[ic]->Divide(4,3);  
  }

  //get input
  for (Int_t ic=0;ic<ncent;ic++){
    fin[ic] = TFile::Open(Form("%s/%s%i.root", path.Data(), fileprefix.Data(), ic));
    if (!fin[ic]) Printf("cannot open input file for bin %i.",ic);
    
    hsys[ic] = (TH1D*) fin[ic]->Get(Form("%s_%i", histSystName.Data(), ic));
    if (!hsys[ic]) Printf("cannot find syst. for bin %i.",ic);  
    hsys[ic]->SetLineColor(color[ic]);
    hsys[ic]->SetMarkerColor(color[ic]);
    hsys[ic]->GetYaxis()->SetRangeUser(0.0, maxYsys);
    hsys[ic]->SetTitle(Form("%s sys. %i-%i", suffix.Data(), ic*20, (ic+1)*20));
    //display original plots
    csys->cd();
    hsys[ic]->Draw("same");

    //ratio to min bias systematics
    hratiosys[ic] = (TH1D*) hsys[ic]->Clone(Form("ratio2mb_%i",ic));
    hratiosys[ic]->SetTitle(Form("ratio2mb_%i",ic));
    hratiosys[ic]->Divide(hratiosys[ic], hsys100, 1.0, 1.0,"B");
    hratiosys[ic]->GetYaxis()->SetTitle("ratio to min.bias");
    csysratio->cd();
    hratiosys[ic]->Draw((ic==0?"":"same"));
    
    //prepare copy of syst uncert. plot
    hsys_smoothed[ic] = (TH1D*) hsys[ic]->Clone(Form("%s_smooth1", hsys[ic]->GetName()));
    hsys_smoothed[ic]->SetTitle(Form("Smooth 1, use mean wrt multiplicity bins if #Delta>#sigma/2 (%i-%i%%)", ic*20, (ic+1)*20));
    hsys_smoothed[ic]->SetLineStyle(2);
    hsys_smoothed[ic]->SetLineWidth(2);
    if (plotSys)hsys_smoothed[ic]->SetLineColor(colorsys[0]);
    // hsys_smoothed[ic]->GetYaxis()->SetLabelSize(0.08);
    // hsys_smoothed[ic]->GetXaxis()->SetLabelSize(0.08);

    hsys_smoothed2[ic] = (TH1D*) hsys[ic]->Clone(Form("%s_smooth2", hsys[ic]->GetName()));
    hsys_smoothed2[ic]->SetTitle(Form("Smooth 2, use min. bias if #Delta>#sigma/2 (%i-%i%%)", ic*20, (ic+1)*20));
    hsys_smoothed2[ic]->SetLineStyle(2);
    hsys_smoothed2[ic]->SetLineWidth(2);
    if (plotSys)hsys_smoothed2[ic]->SetLineColor(colorsys[1]);
    // hsys_smoothed2[ic]->GetYaxis()->SetLabelSize(0.08);
    // hsys_smoothed2[ic]->GetXaxis()->SetLabelSize(0.08)

    hsys_smoothed3[ic] = (TH1D*) hsys[ic]->Clone(Form("%s_smooth3", hsys[ic]->GetName()));
    hsys_smoothed3[ic]->SetTitle(Form("Smooth 3, shift by RMS of residuals wrt fit vs multi bin (%i-%i%%)", ic*20, (ic+1)*20));
    hsys_smoothed3[ic]->SetLineStyle(3);
    hsys_smoothed3[ic]->SetLineWidth(2);
    if (plotSys)hsys_smoothed3[ic]->SetLineColor(colorsys[2]);

    // hsys_smoothed3[ic]->GetYaxis()->SetLabelSize(0.08);
    // hsys_smoothed3[ic]->GetXaxis()->SetLabelSize(0.08);

    hsys_smoothed4[ic] = (TH1D*) hsys[ic]->Clone(Form("%s_smooth4", hsys[ic]->GetName()));
    hsys_smoothed4[ic]->SetTitle(Form("Smooth 4, shift by RMS of residuals wrt fit vs p_{T} (%i-%i%%)", ic*20, (ic+1)*20));
    hsys_smoothed4[ic]->SetLineStyle(4);
    hsys_smoothed4[ic]->SetLineWidth(1);
    if (plotSys)hsys_smoothed4[ic]->SetLineColor(colorsys[3]);
    // hsys_smoothed4[ic]->GetYaxis()->SetLabelSize(0.08);
    // hsys_smoothed4[ic]->GetXaxis()->SetLabelSize(0.08);

    hsys_smoothed5[ic] = (TH1D*) hsys[ic]->Clone(Form("%s_smooth5", hsys[ic]->GetName()));
    hsys_smoothed5[ic]->SetTitle(Form("Smooth 5, as 4 but only if #Delta>#sigma/2 (%i-%i%%)", ic*20, (ic+1)*20));
    hsys_smoothed5[ic]->SetLineStyle(5);
    hsys_smoothed5[ic]->SetLineWidth(2);
   if (plotSys) hsys_smoothed5[ic]->SetLineColor(colorsys[4]);
    // hsys_smoothed5[ic]->GetYaxis()->SetLabelSize(0.08);
    // hsys_smoothed5[ic]->GetXaxis()->SetLabelSize(0.08);

    hsys_smoothed6[ic] = (TH1D*) hsys[ic]->Clone(Form("%s_smooth6", hsys[ic]->GetName()));
    hsys_smoothed6[ic]->SetTitle(Form("Smooth 6, as 1 but use median (%i-%i%%)", ic*20, (ic+1)*20));
    hsys_smoothed6[ic]->SetLineStyle(6);
    hsys_smoothed6[ic]->SetLineWidth(2);
    if (plotSys)hsys_smoothed6[ic]->SetLineColor(colorsys[5]);
    // hsys_smoothed6[ic]->GetYaxis()->SetLabelSize(0.08);
    // hsys_smoothed6[ic]->GetXaxis()->SetLabelSize(0.08);

    hsys_smoothed7[ic] = (TH1D*) hsys[ic]->Clone(Form("%s_smooth7", hsys[ic]->GetName()));
    hsys_smoothed7[ic]->SetTitle(Form("Smooth 7, average adjacent multiplicity bins (%i-%i%%)", ic*20, (ic+1)*20));
    hsys_smoothed7[ic]->SetLineStyle(7);
    hsys_smoothed7[ic]->SetLineWidth(2);
    if (plotSys)hsys_smoothed7[ic]->SetLineColor(colorsys[6]);
    // hsys_smoothed7[ic]->GetYaxis()->SetLabelSize(0.08);
    // hsys_smoothed7[ic]->GetXaxis()->SetLabelSize(0.08);

    hsys_smoothed8[ic] = (TH1D*) hsys[ic]->Clone(Form("%s_smooth8", hsys[ic]->GetName()));
    hsys_smoothed8[ic]->SetTitle(Form("Smooth 8, average with adjacent multiplicity bins (%i-%i%%)", ic*20, (ic+1)*20));
    hsys_smoothed8[ic]->SetLineStyle(8);
    hsys_smoothed8[ic]->SetLineWidth(2);
   if (plotSys) hsys_smoothed8[ic]->SetLineColor(colorsys[7]);

    // hsys_smoothed8[ic]->GetYaxis()->SetLabelSize(0.08);
    // hsys_smoothed8[ic]->GetXaxis()->SetLabelSize(0.08);

    hsys_smoothed9[ic] = (TH1D*) hsys[ic]->Clone(Form("%s_smooth9", hsys[ic]->GetName()));
    hsys_smoothed9[ic]->SetTitle(Form("Smooth 9, average with adjacent p_{T} bins (%i-%i%%)", ic*20, (ic+1)*20));
    hsys_smoothed9[ic]->SetLineStyle(9);
    hsys_smoothed9[ic]->SetLineWidth(2);
    if (plotSys) hsys_smoothed9[ic]->SetLineColor(colorsys[8]);
    // hsys_smoothed9[ic]->GetYaxis()->SetLabelSize(0.08);
    // hsys_smoothed9[ic]->GetXaxis()->SetLabelSize(0.08);
  }
  



  //TH1D * hSmooth1 = new TH1D("smoothed1","smoothed syst. uncert. method 1", npt+1, 0., 15.);
  TH1D * hsys_pt[22]; 
  TH1D * hsys_pt_cent[22];  
  
  for (Int_t ipt=0; ipt<npt; ipt++) {
    Double_t sysMB = hsys100->GetBinContent(ipt+1);
    hsys_pt[ipt] = new TH1D(Form("hsys_pt%i",ipt),Form("pt bin %i, smoothed syst. uncert. method 1",ipt), 400, 0., 40.);
    hsys_pt_cent[ipt] = new TH1D(Form("hsys_pt_cent%i",ipt),"pt bin %i- cent. dependence - smoothed syst. uncert. method 1", 6, 0., 6.);
    hsys_pt_cent[ipt]->GetXaxis()->SetBinLabel(1, "0-20%");
    hsys_pt_cent[ipt]->GetXaxis()->SetBinLabel(2, "20-40%");
    hsys_pt_cent[ipt]->GetXaxis()->SetBinLabel(3, "40-60%");
    hsys_pt_cent[ipt]->GetXaxis()->SetBinLabel(4, "60-80%");
    hsys_pt_cent[ipt]->GetXaxis()->SetBinLabel(5, "80-100%");
    hsys_pt_cent[ipt]->GetXaxis()->SetBinLabel(6, "MB");
    
    Printf("********** PT BIN %i ********** SMOOTHING STRATEGY 1", ipt);
    Printf("minimum bias:    sys = %4.2f", sysMB);
    hsys_pt_cent[ipt]->SetBinContent(6, sysMB);
    hsys_pt_cent[ipt]->GetYaxis()->SetRangeUser(0., maxYsys);
    
    //fill histo with values for all centralities except min bias
    for (Int_t ic=0;ic<ncent;ic++){
      hsys_pt[ipt]->Fill(hsys[ic]->GetBinContent(ipt+1));
      hsys_pt_cent[ipt]->SetBinContent(ic+1, hsys[ic]->GetBinContent(ipt+1));
      Printf("multi bin %i:     sys = %4.2f", ic, hsys[ic]->GetBinContent(ipt+1));
    }
    //get mean across all centralities
    Double_t mean_sys_cent = hsys_pt[ipt]->GetMean();
    Double_t rms_sys_cent  = hsys_pt[ipt]->GetRMS();
    Double_t median_sys_cent = Median(hsys_pt[ipt]);
    Printf("mean sys.:     sys = %4.2f", mean_sys_cent);
    Printf("rms sys.:      rms = %4.2f", rms_sys_cent);
    Printf("median sys.:   sys = %4.2f", median_sys_cent);
    TLine * lmean_vert = new TLine( mean_sys_cent, 0., mean_sys_cent, hsys_pt[ipt]->GetMaximum());
    lmean_vert->SetLineColor(kRed);
    TLine * lmean_horiz = new TLine( 0.0, mean_sys_cent, 5., mean_sys_cent);
    lmean_horiz->SetLineColor(kRed);
    
    //draw plot and line corresponding to mean uncert
    csyspt->cd(ipt+1);
    hsys_pt_cent[ipt]->GetYaxis()->SetRangeUser(0,maxYsys);
    hsys_pt_cent[ipt]->Draw();
    lmean_horiz->Draw("same");
    // if (ipt>0) hsys_pt[ipt]->Draw("same");
    // else hsys_pt[ipt]->Draw();
    // lmean_vert->Draw("same");
    
    //Smoothing 1- substitute with mean the point that deviates more than 1 rms from the mean wrt multi bin
    // minimum bias is not considered for the mean
    //if uncertainty is 0.0, set it to the mean
    for (Int_t ic=0;ic<ncent;ic++){ 
      if ((hsys_pt_cent[ipt]->GetBinContent(ic+1) <= 0.5) ||
	  ((rms_sys_cent>0) && 
	   (TMath::Abs(hsys_pt_cent[ipt]->GetBinContent(ic+1)-mean_sys_cent)/rms_sys_cent > 0.33))) {
	hsys_smoothed[ic]->SetBinContent(ipt+1, median_sys_cent);
      }
    }//end smooth 1
    
    //smoothing 6 - as n.1 but use median instead of mean
    for (Int_t ic=0;ic<ncent;ic++){ 
      if ((hsys_pt_cent[ipt]->GetBinContent(ic+1) <= 0.5) ||
	  ((rms_sys_cent>0) && 
	   (TMath::Abs(hsys_pt_cent[ipt]->GetBinContent(ic+1)-median_sys_cent)/rms_sys_cent > 1.0))) {
	hsys_smoothed6[ic]->SetBinContent(ipt+1, median_sys_cent);
      }
    }//end smooth 6
   
    
    //Smoothing 2 - substitute with min bias value the point that deviates more than 1 rms from the mean wrt multi bin
    // minimum bias is not considered for the mean
    for (Int_t ic=0;ic<ncent;ic++){ 
      if ((hsys_pt_cent[ipt]->GetBinContent(ic+1) <= 0.5) ||
	  ((rms_sys_cent>0) && 
	   (TMath::Abs(hsys_pt_cent[ipt]->GetBinContent(ic+1)-mean_sys_cent)/rms_sys_cent > 0.5))) {
	hsys_smoothed2[ic]->SetBinContent(ipt+1, sysMB);
      }
    }//end smooth 2
    
    
    //Smoothing 3 - for each pt bin: fit vs centrality --> p0, fitresidual = |sys-p0|
    //get rms of fit residuals, shift sys by this rms
    // minimum bias is not considered for the fit
    TF1 * f1 = new TF1("f1","pol0", 0., 6.);
    TFitResultPtr fitresult = hsys_pt_cent[ipt]->Fit(f1,"OQWL");
    Double_t fitp0 = f1->GetParameter(0);
    Double_t fitp0err = f1->GetParError(0);
    TH1D * hfitresiduals = new TH1D("hfitresiduals","hfitresiduals", 200, -10.0, 10.0);
    //fill histogram with residuals wrt the fit for each multi
    for (Int_t ic=0;ic<ncent;ic++){ 
      hfitresiduals->Fill(hsys_pt_cent[ipt]->GetBinContent(ic+1)-fitp0);
    }
    Double_t fitresiduals_rms = hfitresiduals->GetRMS();
    Printf("fit sys.:       p0 = %4.2f", fitp0);
    Printf("rms residuals: rms = %4.2f", fitresiduals_rms);   
    for (Int_t ic=0;ic<ncent;ic++){ 
      if (hsys_pt_cent[ipt]->GetBinContent(ic+1) >= fitp0)
	hsys_smoothed3[ic]->SetBinContent(ipt+1, hsys_pt_cent[ipt]->GetBinContent(ic+1)-fitresiduals_rms);
      else
	hsys_smoothed3[ic]->SetBinContent(ipt+1, hsys_pt_cent[ipt]->GetBinContent(ic+1)+fitresiduals_rms);
    }//end smooth 3
    
  } // end loop on pt bins
  
  
  //Smoothing 4 - smooth vs pt by shifting each uncert by the rms of the residuals wrt the fit 
  //Smoothing 5 - as the smoothing n. 4 but only if distance from fit is > than 1 rms
  
  for (Int_t ic=0;ic<ncent;ic++){ 
    TF1 * f2 = new TF1("f1","pol0", 0., 15.);
    TFitResultPtr fitresult = hsys[ic]->Fit(f2,"OQWL");
    Double_t fit2p0 = f2->GetParameter(0);
    Double_t fit2p0err = f2->GetParError(0);
    TH1D * hfitresiduals2 = new TH1D("hfitresiduals2","hfitresiduals", 200, -10.0, 10.0);
    //fill histogram with residuals wrt the fit for each multi
    for (Int_t ipt=0;ipt<npt;ipt++){ 
      hfitresiduals2->Fill(hsys[ic]->GetBinContent(ipt+1)-fit2p0);
    }
    Double_t fit2residuals_rms = hfitresiduals2->GetRMS();
    for (Int_t ipt=0;ipt<npt;ipt++){ 
      
      if (hsys[ic]->GetBinContent(ipt+1) >= fit2p0) {
	hsys_smoothed4[ic]->SetBinContent(ipt+1, hsys[ic]->GetBinContent(ipt+1)-fit2residuals_rms);
	if (TMath::Abs(hsys[ic]->GetBinContent(ipt+1)-fit2p0)>fit2residuals_rms) 
	  hsys_smoothed5[ic]->SetBinContent(ipt+1, hsys[ic]->GetBinContent(ipt+1)-fit2residuals_rms);
      } else {
	hsys_smoothed4[ic]->SetBinContent(ipt+1, hsys[ic]->GetBinContent(ipt+1)+fit2residuals_rms);
	if (TMath::Abs(hsys[ic]->GetBinContent(ipt+1)-fit2p0)>fit2residuals_rms)
	  hsys_smoothed5[ic]->SetBinContent(ipt+1, hsys[ic]->GetBinContent(ipt+1)+fit2residuals_rms);
      }
    }
  }//end smooth 4


  //Smoothing 7 - substitute for each bin the average between 2 neigbouring bins in centrality
  for (Int_t ipt=0;ipt<npt;ipt++){ 
    for (Int_t ic=0;ic<ncent;ic++){ 
      Double_t prev = 0.0;
      Double_t seq = 0.0;
      if (ic==0) prev = hsys_pt_cent[ipt]->GetBinContent(6);
      else prev = hsys_pt_cent[ipt]->GetBinContent(ic);
      seq = hsys_pt_cent[ipt]->GetBinContent(ic+2);
      Double_t average = (prev+seq)*0.5;
      hsys_smoothed7[ic]->SetBinContent(ipt+1, average);
    }
  }
  
  //Smoothing 8 - substitute for each bin the average between itself + 2 neighbouring bins in multi
  for (Int_t ipt=0;ipt<npt;ipt++){ 
    for (Int_t ic=0;ic<ncent;ic++){ 
      Double_t prev = 0.0;
      Double_t seq = 0.0;
      Double_t self = hsys_pt_cent[ipt]->GetBinContent(ic+1);
      if (ic==0) prev = hsys_pt_cent[ipt]->GetBinContent(6);
      else prev = hsys_pt_cent[ipt]->GetBinContent(ic);
      seq = hsys_pt_cent[ipt]->GetBinContent(ic+2);
      Double_t average = (prev+seq+self)/3.0;
      hsys_smoothed8[ic]->SetBinContent(ipt+1, average);
    }
  }
  
  //Smoothing 9 - substitute for each bin the average among itself and the 2 neighbouring bins in pt
  for (Int_t ic=0;ic<ncent;ic++){ 
    for (Int_t ipt=0;ipt<npt;ipt++){ 
      Double_t prev = hsys[ic]->GetBinContent(ipt);
      Double_t self = hsys[ic]->GetBinContent(ipt+1);
      Double_t seq = hsys[ic]->GetBinContent(ipt+2);
      Double_t average = 0.0;
      if (ipt == 0)
	average = (self+seq)/2.0;
      else 
	if (ipt==npt-1) 
	  average = (prev+self)/2.0;
	else 
	  average = (prev+seq+self)/3.0;
      hsys_smoothed9[ic]->SetBinContent(ipt+1, average);
    }
  }
  

  
  csyssmooth1->cd();
  hsys_smoothed[0]->GetYaxis()->SetRangeUser(0, maxYsys);
  hsys_smoothed[0]->Draw();
  for (Int_t ic=1;ic<ncent;ic++){ 
    hsys_smoothed[ic]->Draw("same");
  }
  gPad->BuildLegend(0.2,0.5,0.88,0.89)->SetFillColor(kWhite);
  
  csyssmooth2->cd();
  hsys_smoothed2[0]->GetYaxis()->SetRangeUser(0, maxYsys);
  hsys_smoothed2[0]->Draw();
  for (Int_t ic=1;ic<ncent;ic++){ 
    hsys_smoothed2[ic]->Draw("same");
  }
  gPad->BuildLegend(0.2,0.5,0.88,0.89)->SetFillColor(kWhite);

  csyssmooth3->cd();
  hsys_smoothed3[0]->GetYaxis()->SetRangeUser(0, maxYsys);
  hsys_smoothed3[0]->Draw();
  for (Int_t ic=1;ic<ncent;ic++){ 
    hsys_smoothed3[ic]->Draw("same");
  }
  gPad->BuildLegend(0.2,0.5,0.88,0.89)->SetFillColor(kWhite);

  csyssmooth4->cd();
  hsys_smoothed4[0]->GetYaxis()->SetRangeUser(0, maxYsys);
  hsys_smoothed4[0]->Draw();
  for (Int_t ic=1;ic<ncent;ic++){ 
    hsys_smoothed4[ic]->Draw("same");
  }
  gPad->BuildLegend(0.2,0.5,0.88,0.89)->SetFillColor(kWhite);

  csyssmooth5->cd();
  hsys_smoothed5[0]->GetYaxis()->SetRangeUser(0, maxYsys);
  hsys_smoothed5[0]->Draw();
  for (Int_t ic=1;ic<ncent;ic++){ 
    hsys_smoothed5[ic]->Draw("same");
  }
  gPad->BuildLegend(0.2,0.5,0.88,0.89)->SetFillColor(kWhite);
  
  csyssmooth6->cd();
  hsys_smoothed6[0]->GetYaxis()->SetRangeUser(0, maxYsys);
  hsys_smoothed6[0]->Draw();
  for (Int_t ic=1;ic<ncent;ic++){ 
    hsys_smoothed6[ic]->Draw("same");
  }
  gPad->BuildLegend(0.2,0.5,0.88,0.89)->SetFillColor(kWhite);

  csyssmooth7->cd();
  hsys_smoothed7[0]->GetYaxis()->SetRangeUser(0, maxYsys);
  hsys_smoothed7[0]->Draw();
  for (Int_t ic=1;ic<ncent;ic++){ 
    hsys_smoothed7[ic]->Draw("same");
  }
  gPad->BuildLegend(0.2,0.5,0.88,0.89)->SetFillColor(kWhite);

  csyssmooth8->cd();
  hsys_smoothed8[0]->GetYaxis()->SetRangeUser(0, maxYsys);
  hsys_smoothed8[0]->Draw();
  for (Int_t ic=1;ic<ncent;ic++){ 
    hsys_smoothed8[ic]->Draw("same");
  }
  gPad->BuildLegend(0.2,0.5,0.88,0.89)->SetFillColor(kWhite);
    
  csyssmooth9->cd();
  hsys_smoothed9[0]->GetYaxis()->SetRangeUser(0, maxYsys);
  hsys_smoothed9[0]->Draw();
  for (Int_t ic=1;ic<ncent;ic++){ 
    hsys_smoothed9[ic]->Draw("same");
  }
  gPad->BuildLegend(0.2,0.5,0.88,0.89)->SetFillColor(kWhite);
    

  TCanvas * csyssmooth[5];
  for (Int_t ic=0;ic<ncent;ic++){
    csyssmooth[ic] = new TCanvas(Form("csyssmooth_%i",ic),Form("systematics multi bin %i",ic), 1200, 800); 
    csyssmooth[ic]->Divide(2,5);  
    csyssmooth[ic]->cd(1);  hsys[ic]->Draw();   
    csyssmooth[ic]->cd(2);  hsys_smoothed[ic]->Draw(); 
    csyssmooth[ic]->cd(3);  hsys_smoothed2[ic]->Draw(); 
    csyssmooth[ic]->cd(4);  hsys_smoothed3[ic]->Draw(); 
    csyssmooth[ic]->cd(5);  hsys_smoothed4[ic]->Draw(); 
    csyssmooth[ic]->cd(6);  hsys_smoothed5[ic]->Draw(); 
    csyssmooth[ic]->cd(7);  hsys_smoothed6[ic]->Draw(); 
    csyssmooth[ic]->cd(8);  hsys_smoothed7[ic]->Draw(); 
    csyssmooth[ic]->cd(9);  hsys_smoothed8[ic]->Draw(); 
    csyssmooth[ic]->cd(10); hsys_smoothed9[ic]->Draw();    
  }
  for (Int_t ic=0;ic<ncent;ic++){
    for (Int_t is=1;is<11;is++){
      csyssmooth[ic]->cd(is); 
      gPad->BuildLegend(0.2,0.75,0.88,0.89)->SetFillColor(kWhite);
    }
  }
  
  csysratio->cd();
  gPad->BuildLegend(0.2,0.5,0.88,0.89)->SetFillColor(kWhite);

  TFile * fout = new TFile(Form("smooth_%s%s.root",fileprefix.Data(),suffix.Data()), "recreate");
  fout->cd();
  for (Int_t ic=0;ic<ncent;ic++){ 
    hsys[ic]->Write();
    hsys_smoothed[ic]->Write();
    hsys_smoothed2[ic]->Write();
    hsys_smoothed3[ic]->Write();  
    hsys_smoothed4[ic]->Write();
    hsys_smoothed5[ic]->Write();
    hsys_smoothed6[ic]->Write();
    hsys_smoothed7[ic]->Write();
    hsys_smoothed8[ic]->Write();
    hsys_smoothed9[ic]->Write();
  }
  csys->Write();
  csyspt->Write();
  csysratio->Write();
  return;
}


double Median(const TH1D * h1) { 

   int n = h1->GetXaxis()->GetNbins();  
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray(); 
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]); 
}



void smoothSystematics100(Bool_t plotSys=0,
		       TString suffix = "Fcn",
		       TString MBpath = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/0to100_binB/systFcn/",
		       TString MBfile = "systematics_EM.root",
		       TString filebins = "/Users/bellini/alice/resonances/kstar_pA5.02TeV/LF_pPb_8-9/multi_binB/proj_2424_tpc2s_tof3sveto.root") 
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  Float_t maxYsys = 15.0;
  Color_t color[8]={ kRed+1, kPink+6, kGreen+1, kAzure+1, kBlue+3, kBlack, kMagenta, kOrange};  
  Color_t colorsys[9]={ kRed, kBlue+1, kGreen+1, kMagenta+1, kOrange, kCyan+2, kPink+2, kGray+2, kViolet+1};  
  
  //get bins
  Int_t npt_axis = 0, ncent_axis=0; 
  TFile *f=TFile::Open(filebins.Data());
  if (!ptbins || !centbins) return;
  TAxis *ptbins = (TAxis*)f->Get("ptbins");
  npt_axis = ptbins->GetNbins();  
  TAxis *centbins = (TAxis*)f->Get("centbins");
  ncent_axis = centbins->GetNbins();
  f->Close();
  const Int_t npt = npt_axis;
  const Int_t ncent = ncent_axis;
  Double_t pt[npt+1];
  for (Int_t k=0; k<npt+1;k++){
    pt[k]=ptbins->GetBinLowEdge(k+1);
  }
  Double_t cent[ncent+1]; 
  for (Int_t k=0; k<ncent+1;k++){
    cent[k]=centbins->GetBinLowEdge(k+1);
  }

  //access minimum bias file
  TFile * fin100 = TFile::Open(Form("%s/%s", MBpath.Data(), MBfile.Data()));
  if (!fin100) Printf("cannot open min bias file.");
  TString histSystName = "hSystVsPtPercentageOfCentral";
  TH1D * hsys = (TH1D*) fin100->Get(Form("%s_0", histSystName.Data()));
  if (!hsys) Printf("cannot find min bias syst.");
  hsys->SetLineColor(color[5]);
  hsys->SetMarkerColor(color[5]);
  hsys->GetYaxis()->SetRangeUser(0.0, maxYsys);
  hsys->SetTitle(Form("%s sys., min. bias", suffix.Data()));
  
  TH1D * hsys_smoothed3[2];
  TH1D * hsys_smoothed2;
  TH1D * hsys_smoothed4;
  TH1D * hsys_smoothed5;
  TH1D * hsys_smoothed1;

  //create canvases
 TCanvas * csys = new TCanvas("csys","systematics", 800, 800); 
 TCanvas * csysFinal = new TCanvas("csysFinal","systematics", 800, 800); 
 
  Int_t ic = 0;
  for (Int_t j=0;j<2;j++){
    hsys_smoothed3[j] = (TH1D*) hsys->Clone(Form("%s_smooth2_iter%i", hsys->GetName(),j));
    hsys_smoothed3[j]->SetTitle(Form("Smooth 3, 2nd iteration of %i (%i-%i%%)", j+1, 0, 100));
    hsys_smoothed3[j]->SetLineStyle(1);
    hsys_smoothed3[j]->SetLineWidth(2);
    if (plotSys) {
      hsys_smoothed3[j]->SetLineColor(colorsys[j]);
      hsys_smoothed3[j]->SetMarkerColor(colorsys[j]);
    }
  }

  hsys_smoothed2 = (TH1D*) hsys->Clone(Form("%s_smooth2", hsys->GetName()));
  hsys_smoothed2->SetTitle(Form("Smooth 2, average with adjacent (%i-%i%%)", 0, 100));
  hsys_smoothed2->SetLineStyle(1);
  hsys_smoothed2->SetLineWidth(2);
  if (plotSys) {
    hsys_smoothed2->SetLineColor(colorsys[2]);
    hsys_smoothed2->SetMarkerColor(colorsys[2]);
  }

  hsys_smoothed4 = (TH1D*) hsys->Clone(Form("%s_smooth4", hsys->GetName()));
  hsys_smoothed4->SetTitle(Form("Smooth 4, shift by residuals to fit (%i-%i%%)", 0, 100));
  hsys_smoothed4->SetLineStyle(4);
  hsys_smoothed4->SetLineWidth(2);
  if (plotSys) {
    hsys_smoothed4->SetLineColor(colorsys[3]);
    hsys_smoothed4->SetMarkerColor(colorsys[3]);
  }
  hsys_smoothed5 = (TH1D*) hsys->Clone(Form("%s_smooth5", hsys->GetName()));
  hsys_smoothed5->SetTitle(Form("Smooth 5, shift as 4 if #Delta>1#sigma (%i-%i%%)", 0, 100));
  hsys_smoothed5->SetLineStyle(5);
  hsys_smoothed5->SetLineWidth(2);
  if (plotSys) {
    hsys_smoothed5->SetLineColor(colorsys[4]);
    hsys_smoothed5->SetMarkerColor(colorsys[4]);
  }
  hsys_smoothed1 = (TH1D*) hsys->Clone(Form("%s_smooth1", hsys->GetName()));
  hsys_smoothed1->SetTitle(Form("Smooth 1, average of adjacent (%i-%i%%)", 0, 100));
  hsys_smoothed1->SetLineStyle(9);
  hsys_smoothed1->SetLineWidth(2);
  if (plotSys) {
    hsys_smoothed1->SetLineColor(colorsys[8]);
    hsys_smoothed1->SetMarkerColor(colorsys[3]);
  }
  
  //Smoothing 1 - substitute for each bin the average among itself and the 2 neighbouring bins in pt
  for (Int_t ipt=0;ipt<npt;ipt++){ 
    Double_t prev = hsys->GetBinContent(ipt);
    Double_t self = hsys->GetBinContent(ipt+1);
    Double_t seq = hsys->GetBinContent(ipt+2);
    Double_t average2 = 0.0;
    Double_t average3 = 0.0;
    if (ipt == 0){
      average2 = (self+seq)/2.0;
      average3 = (self+seq)/2.0;
    }    else {
      if (ipt==npt-1) { 
	average2 = (prev+self)/2.0;
	average3 = (prev+self)/2.0;
      }      else { 
	average2 = (prev+seq)/2.0;
	average3 = (prev+self+seq)/3.0;
      }
      hsys_smoothed1->SetBinContent(ipt+1, average2);
      hsys_smoothed2->SetBinContent(ipt+1, average3);
    }
  }
  
  //Smoothing 2 - second iteration for each bin the average among itself and the 2 neighbouring bins in pt  
  for (Int_t ipt=0;ipt<npt;ipt++){ 
    Double_t prev = hsys_smoothed2->GetBinContent(ipt);
    Double_t self = hsys_smoothed2->GetBinContent(ipt+1);
    Double_t seq = hsys_smoothed2->GetBinContent(ipt+2);
    Double_t average2 = 0.0;
    Double_t average3 = 0.0;
    if (ipt == 0){
      average2 = (self+seq)/2.0;
      average3 = (self+seq)/2.0;
    }    else {
      if (ipt==npt-1) { 
	average2 = (prev+self)/2.0;
	average3 = (prev+self)/2.0;
      }      else { 
	average2 = (prev+seq)/2.0;
	average3 = (prev+self+seq)/3.0;
      }
      hsys_smoothed3[0]->SetBinContent(ipt+1, average2);
      hsys_smoothed3[1]->SetBinContent(ipt+1, average3);
    }
  }

 //Smoothing 4 - smooth vs pt by shifting each uncert by the rms of the residuals wrt the fit 
  //Smoothing 5 - as the smoothing n. 4 but only if distance from fit is > than 1 rms  
  TF1 * f2 = new TF1("f1","pol0", 0.4, 15.);
  TFitResultPtr fitresult = hsys->Fit(f2,"0QWLR");
  Double_t fit2p0 = f2->GetParameter(0);
  Double_t fit2p0err = f2->GetParError(0);
  TH1D * hfitresiduals2 = new TH1D("hfitresiduals2","hfitresiduals", 200, -10.0, 10.0);
  //fill histogram with residuals wrt the fit for each multi
  for (Int_t ipt=0;ipt<npt;ipt++){ 
    hfitresiduals2->Fill(hsys->GetBinContent(ipt+1)-fit2p0);
  }
  Double_t fit2residuals_rms = hfitresiduals2->GetRMS();
  for (Int_t ipt=0;ipt<npt;ipt++){ 
    if (hsys->GetBinContent(ipt+1) >= fit2p0) {
      hsys_smoothed4->SetBinContent(ipt+1, hsys->GetBinContent(ipt+1)-fit2residuals_rms);
      if (TMath::Abs(hsys->GetBinContent(ipt+1)-fit2p0)>fit2residuals_rms) 
	hsys_smoothed5->SetBinContent(ipt+1, hsys->GetBinContent(ipt+1)-fit2residuals_rms);
    } else {
      hsys_smoothed4->SetBinContent(ipt+1, hsys->GetBinContent(ipt+1)+fit2residuals_rms);
      if (TMath::Abs(hsys->GetBinContent(ipt+1)-fit2p0)>fit2residuals_rms)
	hsys_smoothed5->SetBinContent(ipt+1, hsys->GetBinContent(ipt+1)+fit2residuals_rms);
    }
  }//end smooth 4
  

  csys->cd();  
  hsys->Draw();
  hsys->GetYaxis()->SetRangeUser(0, maxYsys);
  hsys_smoothed1->Draw("same");
  hsys_smoothed2->Draw("same");
  hsys_smoothed3[0]->Draw("same");
  hsys_smoothed3[1]->Draw("same");
  hsys_smoothed4->Draw("same");
  hsys_smoothed5->Draw("same");
  gPad->BuildLegend(0.2,0.64,0.89,0.89)->SetFillColor(kWhite);
     
  csysFinal->cd();
  hsys->Draw();
  hsys_smoothed3[1]->Draw("same");
  hsys_smoothed3[1]->SetLineWidth(3);
  TLegend * leg = new TLegend(0.5,0.7,0.88,0.89);
  leg->SetHeader(Form("%s systematics",suffix.Data()));
  leg->AddEntry(hsys, "before smoothing","l");
  leg->AddEntry(hsys_smoothed3[1], "after smoothing","l");
  leg->Draw();
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  TString file_name = Form("smooth_%s_%s", suffix.Data(), MBfile.Data());
  TFile * fout = new TFile(Form("0to100_%s",file_name.Data()), "recreate");
  fout->cd();
  hsys->Write();
  hsys_smoothed1->Write();
  hsys_smoothed2->Write();
  hsys_smoothed3[0]->Write();
  hsys_smoothed3[1]->Write();
  hsys_smoothed4->Write();
  hsys_smoothed5->Write();
  csysFinal->Print(file_name.ReplaceAll(".root",".png"));  
  csys->Print(file_name.Prepend("summary_"));  
  
  return;
}

