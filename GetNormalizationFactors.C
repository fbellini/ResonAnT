/* fbellini, created on 20 aug 2012 */

Int_t  GetNormalizationFactors(TString fileName, TString listNameSuffix = "_tof2s")
{
  gSystem->Execute("SetGraphicStyle.C","0");
  TString hEventsName = "hEventStat";
  TString hAccEventVsCentName = "hAEventsVsMulti";
  TString listName = "RsnOut";
  listName.Append(listNameSuffix.Data());
  
  TFile * fin = TFile::Open(fileName.Data(),"READ");
  if (!fin) return 1;
  TList * lin = (TList*) fin->Get(listName.Data());
  if (!lin) return 2;
  TH1D * hEvents = (TH1D*) lin->FindObject(hEventsName.Data());
  if (!hEvents) return 3;
  TH1D * hAccEventVsCent = (TH1D*) lin->FindObject(hAccEventVsCentName.Data());
  if (!hAccEventVsCent) return 4;

  Int_t counters[4] = {0,0,0,0};
  for (Int_t j=0; j<4;j++){
    counters[j] = hEvents->GetBinContent(j+1);
  }
  Printf("Event counters: \n total CINT1B = %i \n total V0AND = %i \n candle events = %i \n total ACCEPTED = %i \n ACCEPTED/CINT1B = %8.5f", 
	 counters[0],counters[1],counters[2],counters[3], (Float_t)counters[3]/counters[0]);
  
  Float_t centLowE[6] = {0.0, 20.0, 40.0, 60.0, 80.0, 101.0};
  Double_t centCounters[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  TH1D * hCentDraw[5];
  
  for (Int_t i=0; i<5;i++){
    Int_t lowBin = hAccEventVsCent->GetXaxis()->FindBin(centLowE[i]);
    Int_t upBin = hAccEventVsCent->GetXaxis()->FindBin(centLowE[i+1]-1);
    // Printf("lower value = %4.2f --> lower bin = %i \n upper value = %4.2f --> upper bin = %i", centLowE[i], lowBin, centLowE[i+1], upBin );
    centCounters[i] = (Double_t) hAccEventVsCent->Integral(lowBin, upBin);
    Printf("Accepted events in centrality bin (%2.0f -%3.0f)%% = %e", centLowE[i], centLowE[i+1], centCounters[i]);
    hCentDraw[i] = new TH1D(Form("Cent%i",i),Form("Cent%i",i), 
			    hAccEventVsCent->GetXaxis()->GetNbins(),
			    hAccEventVsCent->GetXaxis()->GetBinLowEdge(1),
			    hAccEventVsCent->GetXaxis()->GetBinUpEdge(hAccEventVsCent->GetXaxis()->GetNbins()));
    for (Int_t ibin = lowBin; ibin<upBin+1; ibin++){
      hCentDraw[i]->SetBinContent(ibin, hAccEventVsCent->GetBinContent(ibin));
    }

  }
  
  //  Color_t color[5] = {kOrange+6, kOrange, kGreen-5, kCyan-6, kAzure-9};
  Color_t color[5] = {kOrange-9, kYellow-8, kGreen-8, kCyan-6, kAzure-9};
  //hAccEventVsCent->Rebin();
  hAccEventVsCent->SetXTitle("centrality (%)");
  hAccEventVsCent->GetYaxis()->SetRangeUser(1,1e6);
  hAccEventVsCent->SetYTitle("accepted events");
  hAccEventVsCent->SetLineWidth(2);
  hAccEventVsCent->Draw();
  for (Int_t i=0; i<5;i++){
    hCentDraw[i]->SetLineStyle(2);
    hCentDraw[i]->SetFillStyle(3001);
    hCentDraw[i]->SetFillColor(color[i]);
    hCentDraw[i]->Draw("Hsame");
  }

  return 0;
}
