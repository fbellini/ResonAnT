//<dN/dETA> pPb 5.02 TeV - 0-20, 20-40, 40-60, 60-80, 80-100
double x_pPb_small[5]={ 35.6, 23.2, 16.1, 9.8, 4.4};
double ux_pPb_small[5]={ 0.8,  0.5,  0.4, 0.2, 0.1};

//<dN/dETA> pPb 5.02 TeV - 0-5, 5-10, 10-20, 20-40, 40-60, 60-80, 80-100
double x_pPb[7]  = { 45.15, 36.22, 30.46, 23.24, 16.08, 9.82, 4.41};
double ux_pPb[7] = {  1,     0.8,   0.67,  0.51,  0.35, 0.21, 0.1};

//Phi in p-Pb, mean pt
double phi_MeanPt_pPb502[7] = {1.436918, 1.442072, 1.421458, 1.357387, 1.309693, 1.242375, 1.054900};
double phi_MeanPt_pPb502_stat[7] = {0.008972, 0.009296, 0.007526, 0.005608, 0.006364, 0.007829, 0.009562};
double phi_MeanPt_pPb502_sys[7] = {0.028202, 0.025423, 0.023951, 0.025181, 0.030687, 0.024283, 0.030287};

//K* in p-Pb, mean pt, FINAL 01 july 2015
double Ks_MeanPt_pPb502[5] = {1.3789, 1.2995, 1.2105, 1.1075, 0.9423};
double Ks_MeanPt_pPb502_stat[5] = {0.0110, 0.0099, 0.0091, 0.0094, 0.0089};//stat
double Ks_MeanPt_pPb502_sys[5] = { 0.0181, 0.0174, 0.0158, 0.0203, 0.0148};
 
//K* in p-Pb, dN/dy, FINAL 01 july 2015
double Ks_dNdy_pPb502[5] = { 0.616, 0.426, 0.302, 0.187*0.995, 0.089*0.945 };
double Ks_dNdy_pPb502_stat[5] = { 0.008, 0.005, 0.004, 0.003*0.995, 0.001*0.945 };
double Ks_dNdy_pPb502_sys[5] = { 0.05, 0.04, 0.03, 0.014*0.995, 0.008*0.945 };

//phi in p-Pb, dN/dy, FINAL 01 july 2015
double phi_dNdy_pPb502[7] = { 0.378, 0.288, 0.244, 0.185, 0.123, 0.070*0.995, 0.032*0.945 };
double phi_dNdy_pPb502_stat[7] = { 0.004, 0.003, 0.002, 0.001, 0.001, 0.001*0.995, 0.0005*0.945 };
double phi_dNdy_pPb502_sys[7] = { 0.034, 0.021, 0.02, 0.02, 0.01, 0.006*0.995, 0.004*0.945 };

Double_t XimYield_Final[7] =     {0.1175, 0.0923, 0.0746, 0.0545, 0.0356, 0.0196, 0.0068};
Double_t XimYieldStat_Final[7] = {0.0014, 0.0011, 0.0007, 0.0004, 0.0005, 0.0003, 0.0002};
Double_t XimYieldSyst_Final[7] = {0.0078, 0.0061, 0.0050, 0.0036, 0.0024, 0.0014, 0.0005};
Double_t XipYield_Final[7] =     {0.1173, 0.0939, 0.0754, 0.0555, 0.037, 0.0200, 0.0072};
Double_t XipYieldStat_Final[7] = {0.0014, 0.0011, 0.0007, 0.0005, 0.0004, 0.0003, 0.0002};
Double_t XipYieldSyst_Final[7] = {0.0076, 0.0075, 0.0057, 0.0045, 0.0038, 0.0016, 0.0009};

Double_t OmMYield_Final[7] =     {0.0126, 0.011,  0.0082, 0.0062, 0.0036, 0.0018, 0.00055};
Double_t OmMYieldStat_Final[7] = {0.0007, 0.0006, 0.0004, 0.0002, 0.0002, 0.0001, 0.00007};
Double_t OmMYieldSyst_Final[7] = {0.0016, 0.0014, 0.0011, 0.0008, 0.0005, 0.0003, 0.00011};
Double_t OmPYield_Final[7] =     {0.0134, 0.0104, 0.0085, 0.0058, 0.0036, 0.0024, 0.00071};
Double_t OmPYieldStat_Final[7] = {0.0008, 0.0006, 0.0005, 0.0002, 0.0002, 0.0002, 0.00008};
Double_t OmPYieldSyst_Final[7] = {0.0017, 0.0014, 0.0011, 0.0007, 0.0005, 0.0004, 0.00014};

Double_t XimMeanPt_Final[7] = {1.574, 1.529, 1.515, 1.468, 1.394, 1.280, 1.1546};
Double_t XimMeanPtStat_Final[7] = {0.010, 0.009, 0.007, 0.006, 0.012, 0.009, 0.0175};
Double_t XimMeanPtSyst_Final[7] = {0.037, 0.036, 0.034, 0.032, 0.030, 0.032, 0.0368};
// Double_t XimMeanPtUncorrSyst[7] = {0.012, 0.014, 0.016, 0.011, 0.016, 0.016, 0.017};
Double_t XipMeanPt_Final[7] = {1.549, 1.516, 1.502, 1.453, 1.358, 1.270, 1.1182};
Double_t XipMeanPtStat_Final[7] = {0.010, 0.009, 0.007, 0.006, 0.007, 0.009, 0.0172};
Double_t XipMeanPtSyst_Final[7] = {0.037, 0.036, 0.034, 0.032, 0.029, 0.041, 0.0398};
// Double_t XipMeanPtUncorrSyst[7] = {0.012, 0.014, 0.016, 0.011, 0.015, 0.019, 0.018};

Double_t OmMMeanPt_Final[7] = {1.818, 1.829, 1.771, 1.668, 1.563, 1.544, 1.350};
Double_t OmMMeanPtStat_Final[7] = {0.056, 0.062, 0.058, 0.040, 0.033, 0.055, 0.110};
Double_t OmMMeanPtSyst_Final[7] = {0.076, 0.079, 0.082, 0.082, 0.067, 0.084, 0.106};
// Double_t OmMMeanPtUncorrSyst[7] = {0.025, 0.030, 0.038, 0.029, 0.035, 0.040, 0.048};
Double_t OmPMeanPt_Final[7] = {1.761, 1.847, 1.739, 1.707, 1.593, 1.382, 1.187};
Double_t OmPMeanPtStat_Final[7] = {0.060, 0.062, 0.064, 0.036, 0.040, 0.078, 0.097};
Double_t OmPMeanPtSyst_Final[7] = {0.075, 0.086, 0.073, 0.061, 0.058, 0.076, 0.096};
// Double_t OmPMeanPtUncorrSyst[7] = {0.025, 0.032, 0.034, 0.022, 0.031, 0.036, 0.043};

// Double_t XimMeanPt_Final[7] = {1.574, 1.529, 1.515, 1.468, 1.394, 1.287, 1.1544};
// Double_t XimMeanPtStat_Final[7] = {0.010, 0.009, 0.007, 0.006, 0.012, 0.010, 0.0177};
// Double_t XimMeanPtSyst_Final[7] = {0.052, 0.052, 0.054, 0.055, 0.057, 0.063, 0.0660};
// Double_t XipMeanPt_Final[7] = {1.549, 1.516, 1.502, 1.453, 1.358, 1.271, 1.1185};
// Double_t XipMeanPtStat_Final[7] = {0.010, 0.009, 0.007, 0.006, 0.007, 0.009, 0.0175};
// Double_t XipMeanPtSyst_Final[7] = {0.068, 0.059, 0.061, 0.062, 0.078, 0.065, 0.0860};

// Double_t OmMMeanPt_Final[7] = {1.818, 1.830, 1.771, 1.669, 1.564, 1.543, 1.350};
// Double_t OmMMeanPtStat_Final[7] = {0.056, 0.060, 0.057, 0.040, 0.033, 0.054, 0.116};
// Double_t OmMMeanPtSyst_Final[7] = {0.156, 0.159, 0.157, 0.147, 0.146, 0.171, 0.209};
// Double_t OmPMeanPt_Final[7] = {1.761, 1.846, 1.739, 1.707, 1.592, 1.382, 1.244};
// Double_t OmPMeanPtStat_Final[7] = {0.060, 0.062, 0.064, 0.036, 0.040, 0.075, 0.077};
// Double_t OmPMeanPtSyst_Final[7] = {0.152, 0.163, 0.150, 0.139, 0.145, 0.154, 0.186};


void SaveGraphFromArrays(Double_t * arrayX, 
			 Double_t * arrayXerr,
			 Double_t * arrayY, 
			 Double_t * arrayYstat,
			 Double_t * arrayYsys,
			 Int_t npoints = 7,
			 TString quantity = "Yield",
			 TString particle = "phi",
			 Bool_t save2file = 1)
{
  if (!arrayX || !arrayXerr || !arrayY || !arrayYstat || !arrayYsys) return;

  TGraphErrors *g_stat = new TGraphErrors(npoints, arrayX, arrayY, arrayXerr, arrayYstat);
  g_stat->SetMarkerStyle(20);
  g_stat->SetMarkerColor(kBlue);
  g_stat->SetLineColor(kBlue);
  g_stat->SetMarkerSize(1.);
  g_stat->SetLineWidth(1);
  g_stat->SetName(Form("g%s_stat_%s_sum",quantity.Data(), particle.Data()));

  TGraphErrors *g_sys = new TGraphErrors(npoints, arrayX, arrayY, arrayXerr, arrayYsys);
  g_sys->SetMarkerStyle(20);
  g_sys->SetMarkerColor(kBlue);
  g_sys->SetLineColor(kBlue);
  g_sys->SetMarkerSize(1.);
  g_sys->SetLineWidth(1);
  g_sys->SetName(Form("g%s_sys_%s_sum",quantity.Data(), particle.Data()));
  

  if (save2file) {
    TString filename = "pPb_rsn_dNdy_graph.root";
    if (quantity.Contains("Mean")) filename.ReplaceAll("dNdy","mean");
    filename.ReplaceAll("rsn",particle.Data());
    TFile * fout = new TFile(filename.Data(), "update");
    fout->cd();
    g_stat->Write();
    g_sys->Write();
    fout->Close();
    delete fout;
  }
  return;
}
