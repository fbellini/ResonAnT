/* fbellini, created on 20 aug 2012 */

Int_t  GetNormalizationFactors(TString fileName, TString listNameSuffix = "_tpc2s_tof3sveto")
{
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
  
  Float_t centLowE[5] = {0.0, 10., 30.0, 60.0, 90.0};
  Double_t centCounters[4] = {0.0, 0.0, 0.0, 0.0};
  TH1D * hCentDraw[4];
  
  for (Int_t i=0; i<4;i++){
    Int_t lowBin = hAccEventVsCent->GetXaxis()->FindBin(centLowE[i]);
    Int_t upBin = hAccEventVsCent->GetXaxis()->FindBin(centLowE[i+1]-1);
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
  
  Printf(Form("Normalization array: norm = { %e, %e, %e %e}", centCounters[0], centCounters[1], centCounters[2], centCounters[3]));

  Color_t color[4] = {kRed+2-9, kSpring-2, kBlue+2, kGray+1};
  hAccEventVsCent->SetXTitle("centrality (%)");
  hAccEventVsCent->GetYaxis()->SetRangeUser(1,1e6);
  hAccEventVsCent->SetYTitle("accepted events");
  hAccEventVsCent->SetLineWidth(2);
  hAccEventVsCent->Draw();
  for (Int_t i=0; i<4;i++){
    hCentDraw[i]->SetLineStyle(2);
    hCentDraw[i]->SetFillStyle(1001);
    hCentDraw[i]->SetFillColorAlpha(color[i],0.2);
    hCentDraw[i]->Draw("Hsame");
  }

  return 0;
}
