void GetPhiPhaseSpace
(
   const char *nameData = "allMerged_RsnOut.root",
   const char *outName  = "phi_PhaseSpace",
   const char *listName = "RsnOut_tpc2sPtDep_tof2sveto5smism",
   Bool_t isMC = kTRUE
)
{

  Bool_t kComputeEfficiency=kFALSE;
  const Int_t kNhistosData=1;
  TString macroDir ="/Users/bellini/alice/macros/ResonAnT/phiXeXe";

   // initial setup
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat("1111");
  gStyle->SetTextFont(42);
   
  //output file
  TFile *fout = TFile::Open(Form("%s.root", outName), "RECREATE");
  // output list
  TList *out = new TList();
  
  // output lists - define once and then clear them
  TList lUnlikePM;
  lUnlikePM.SetName(Form("%s_phi", (!isMC ? "Data" : "MC")));
  
   // open input file
   TFile *fileData;// = 0x0;//, *fileMC = 0x0;
   TList *listData;// = 0x0;//, *listMC = 0x0;
   fileData = TFile::Open(nameData);
   if (fileData && fileData->IsOpen()) listData = (TList*)fileData->Get(listName);
   
   //get input 
   TH3F* hInput[kNhistosData]= {0};
   if (listData) {
     hInput[ 0] = (TH3F*)listData->FindObject("PhiXeXeMC_PhaseSpacetpc2sPtDep_tof2sveto5smism");
     hInput[ 0]->SetName("phi_phaseSpace");
   }
   if (! hInput[ 0]) return;
   // define binning in pT    
   Double_t   pt[] = {0.0, 0.3, 0.5, 0.7, 0.9, 1.10, 1.30, 1.50, 2.00, 3.00, 4.00, 5.0, 7.0, 10.0};
   Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;
   TAxis *ptbins = new TAxis(npt, pt);
   Double_t cent[]={ 0.0, 100.};
   //Double_t cent[]={ 0.0, 20.0, 40.0, 60.0, 80.0, 90.0};   
   Int_t   ncent  = sizeof(cent) / sizeof(cent[0]) - 1;
   TAxis *centbins = new TAxis(ncent, cent);
   
   Color_t color[2][21];
   for (Int_t ibin = 0; ibin < npt; ibin++) {
     if (ibin<7) {
       color[0][ibin] = kRed+2-ibin;
       color[1][ibin] = kBlue+2-ibin;
     } else {
       if (ibin<14) {
	 color[0][ibin] = kOrange-8+ibin;
	 color[1][ibin] = kTeal-8+ibin;
       } else {
	 color[0][ibin] = kMagenta-19+ibin;
	 color[1][ibin] = kAzure-19+ibin;
       }
     }
   }

   const Int_t knpt = npt;
   TH2D * p_yx[kNhistosData][knpt];
   TH1D * p_1st[kNhistosData][knpt];
   TH1D * p_2nd[kNhistosData][knpt];
   TH1D * px[kNhistosData];
   TH1D * py[kNhistosData];
   TH1D * pz[kNhistosData];
   TH1D * hPeakPt_1st[kNhistosData] ;
   TH1D * hPeakPt_2nd[kNhistosData];
   TH1D * hMeanPt_1st[kNhistosData] ;
   TH1D * hMeanPt_2nd[kNhistosData];

   TCanvas * c1 = new TCanvas("Ks","Ks", 1200, 900); c1->Divide(5,4);   
   TCanvas * c2 = new TCanvas("aKs","aKs", 1200, 900); c2->Divide(5,4);   
   TCanvas * c3 = new TCanvas("peak_pt","peak_pt", 800, 800);    
   TCanvas * c4 = new TCanvas("mean_pt","mean_pt", 800, 800);    
   
   //TProfile2D * prof_yx[kNhistosData][knpt];

   for (Int_t i = 0; i < kNhistosData; i++) {
     if (!hInput[i]) continue;
     hInput[i]->GetZaxis()->SetTitle("resonance p_{T} (GeV/c)");
     hInput[i]->GetYaxis()->SetTitle("2nd daughter (#pi) p_{T} (GeV/c)");
     hInput[i]->GetXaxis()->SetTitle("1st daughter (K) p_{T} (GeV/c)");
     px[i] = (TH1D*) hInput[i]->Project3D("x"); px[i]->SetName(Form("pt_%s_1daughter",(i==0?"phi":"phi")));
     py[i] = (TH1D*) hInput[i]->Project3D("y"); py[i]->SetName(Form("pt_%s_2daughter", (i==0?"phi":"phi")));
     pz[i] = (TH1D*) hInput[i]->Project3D("z"); pz[i]->SetName(Form("pt_%s_mother", (i==0?"phi":"phi")));

     hPeakPt_1st [i] = new TH1D(Form("hPeakPt_%s",(i==0?"Kplus":"Kminus")),Form("Peak p_{T} for 1st daughter of phi"), npt, pt );
     hPeakPt_2nd [i] = new TH1D(Form("hPeakPt_%s",(i==0?"Kminus":"Kplus")),Form("Peak p_{T} for 2nd daughter phi"), npt, pt );
     hMeanPt_1st [i] = new TH1D(Form("hMeanPt_%s",(i==0?"Kplus":"Kminus")),Form("Mean p_{T} for 1st daughter of phi"), npt, pt );
     hMeanPt_2nd [i] = new TH1D(Form("hMeanPt_%s",(i==0?"Kminus":"Kplus")),Form("Mean p_{T} for 2nd daughter phi"), npt, pt );
      
     for (Int_t ibin = 0; ibin < npt; ibin++) 
       {
	 hInput[i]->GetZaxis()->SetRangeUser(pt[ibin],pt[ibin+1]);
	 Printf("Projecting along z in range %4.2f-%4.2f", pt[ibin],pt[ibin+1]);
	 p_yx[i][ibin] = (TH2D*) hInput[i]->Project3D("yx");
	 p_1st[i][ibin] = (TH1D*) hInput[i]->Project3D("x");
	 p_2nd[i][ibin] = (TH1D*) hInput[i]->Project3D("y");
	 
	 //	 prof_yx[i][ibin] = (TProfile2D*) hInput[i]->Project3DProfile("yx");
	 p_yx[i][ibin]->SetNameTitle(Form("%s_ps_pt%i", (i==0?"phi":"phi"), ibin),Form("phase space %4.2f<p_{T}(K*)<%4.2f", pt[ibin],pt[ibin+1]));
	 p_1st[i][ibin]->SetNameTitle(Form("%s_1std_pt%i", (i==0?"phi":"phi"), ibin),Form("1st daughter, %4.2f<p_{T}(K*)<%4.2f", pt[ibin],pt[ibin+1]));
	 p_2nd[i][ibin]->SetNameTitle(Form("%s_2nsd_pt%i", (i==0?"phi":"phi"), ibin),Form("2nd daughter, %4.2f<p_{T}(K*)<%4.2f", pt[ibin],pt[ibin+1]));
	 //prof_yx[i][ibin]->SetNameTitle(Form("%s_ps_prof_pt%i", (i==0?"Ks":"aKs"), ibin),Form("phase space %4.2f<p_{T}(K*)<%4.2f", pt[ibin],pt[ibin+1]));
	 Int_t max = p_1st[i][ibin]->GetMaximumBin();
	 Double_t maxpt = p_1st[i][ibin]->GetBinCenter(max);
	 hPeakPt_1st[i]->SetBinContent(ibin+1, maxpt);
	 hMeanPt_1st[i]->SetBinContent(ibin+1, p_1st[i][ibin]->GetMean());
	 max = p_2nd[i][ibin]->GetMaximumBin();
	 maxpt = p_2nd[i][ibin]->GetBinCenter(max);
	 hPeakPt_2nd[i]->SetBinContent(ibin+1, maxpt);
	 hMeanPt_2nd[i]->SetBinContent(ibin+1, p_2nd[i][ibin]->GetMean());
	 //makeup
	 p_1st[i][ibin]->GetYaxis()->SetRangeUser(1, 1e6);
	 p_1st[i][ibin]->SetLineColor(color[0][ibin]); 	 p_1st[i][ibin]->SetMarkerColor(color[0][ibin]); p_1st[i][ibin]->SetFillColor(color[0][ibin]); p_1st[i][ibin]->SetFillStyle(3001);

	 p_2nd[i][ibin]->GetYaxis()->SetRangeUser(1, 1e6);
	 p_2nd[i][ibin]->SetLineColor(color[1][ibin]); 	 p_2nd[i][ibin]->SetMarkerColor(color[1][ibin]); p_2nd[i][ibin]->SetFillColor(color[1][ibin]); p_2nd[i][ibin]->SetFillStyle(1);

	 // save into a file
	 fout->cd();   
	 p_yx[i][ibin]->Write();
	 p_1st[i][ibin]->Write();
	 p_2nd[i][ibin]->Write();

	 //prof_yx[i][ibin]->Write();
	 if (i==0) c1->cd(ibin+1); 
	 else c2->cd(ibin+1);
	 gPad->SetLogy();
	 if (ibin==0)  {
	   p_1st[i][ibin]->Draw("hist");
	   p_2nd[i][ibin]->Draw("hist same");
	 } else {
	   p_1st[i][ibin]->Draw("hist same");
	   p_2nd[i][ibin]->Draw("hist same");
	 }
       }
     hPeakPt_2nd[i]->SetLineColor(color[i][0]);
     hPeakPt_2nd[i]->SetMarkerColor(color[i][0]);
     hPeakPt_2nd[i]->SetLineStyle(2);
     hPeakPt_1st[i]->SetLineColor(color[i][0]);
     hPeakPt_1st[i]->SetMarkerColor(color[i][0]);
     hPeakPt_1st[i]->SetLineStyle(1);

     hMeanPt_2nd[i]->SetLineColor(color[i][0]);
     hMeanPt_2nd[i]->SetMarkerColor(color[i][0]);
     hMeanPt_2nd[i]->SetLineStyle(2);
     hMeanPt_1st[i]->SetLineColor(color[i][0]);
     hMeanPt_1st[i]->SetMarkerColor(color[i][0]);
     hMeanPt_1st[i]->SetLineStyle(1);
   
     fout->cd();   
     hPeakPt_2nd[i]->Write();
     hPeakPt_1st[i]->Write();
     hMeanPt_2nd[i]->Write();
     hMeanPt_1st[i]->Write();
     px[i]->Write(); 
     py[i]->Write(); 
     pz[i]->Write();
     c3->cd();   
     if (i==0) hPeakPt_1st[i]->Draw("");     
     else hPeakPt_1st[i]->Draw("same");     
     hPeakPt_2nd[i]->Draw("same");
     c4->cd();   
     if (i==0) hMeanPt_1st[i]->Draw("");     
     else hMeanPt_1st[i]->Draw("same");     
     hMeanPt_2nd[i]->Draw("same");
   }
   hPeakPt_1st[0]->GetXaxis()->SetTitle("p_{T} (#phi) (GeV/c)");
   hPeakPt_1st[0]->GetYaxis()->SetTitle("daughter p_{T} (GeV/c)");
   hPeakPt_1st[0]->GetYaxis()->SetRangeUser(0.0, 6.0);
   hMeanPt_1st[0]->GetXaxis()->SetTitle("p_{T} (#phi) (GeV/c)");
   hMeanPt_1st[0]->GetYaxis()->SetTitle("daughter p_{T} (GeV/c)");
   hMeanPt_1st[0]->GetYaxis()->SetRangeUser(0.0, 6.0);

   c3->cd();   
   TLegend* leg = (TLegend*) gPad->BuildLegend(0.11,0.70, 0.55,0.89);
   leg->SetFillColor(kWhite);
   c3->Print("phi_phaseSpace.pdf");

   return;
}

 
//---------------------------------------------
void FillListWithProj(TList * listIn, TList * listOut, Int_t iBin)
{
  //selects all type of histos per centrality bin and saves them in separate lists  
  TString listName = listOut->GetName();
  TString type;
  
  if (listName.Contains("UnlikePM")) type = "UnlikePM";
  if (listName.Contains("UnlikeMP")) type = "UnlikeMP";
  
  Int_t nhisto = listIn->GetEntries();
  for (Int_t ih=0;ih<nhisto;ih++){
    TH1F * dummy = (TH1F *)listIn->At(ih);
    TString hName = dummy->GetName();
    if (hName.Contains(type.Data())){
      if (hName.Contains(Form("bin%02d",iBin))){
	listOut->AddLast(dummy);
	Printf("Histo found with name: %s",hName.Data());
      }
    }
  } 
  
  Printf("============== Output");
  listOut->ls();
  return;
}






 
