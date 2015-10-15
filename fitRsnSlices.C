Int_t fitRsnSlices(
	       TString filename = "train3_aod139.root",
	       TString listname = "RsnOut_tpc2s",
	       TString hname = "cutQ_bit5.TOF_nsigmaPi_vsP_p_pi")
{

  TFile * f = TFile::Open(filename.Data());
  if (!f) return 1;

  TList * list = (TList*) f->Get(listname.Data());
  if (!list) return 2;

  TH2D * h = (TH2D*) list->FindObject(hname.Data());
  if (!h) return 3;
  h->Draw("colz");     

  TF1 fg("fg","gaus",-2.,2.); // optimal fit range: [-2,2] sigma for TPC, [-3.,1.] for TOF
  TLine l;
  TObjArray arr;
  fg.SetParameters(1,0,1);
  h->FitSlicesY(&fg,0,-1,0,"NQR",&arr);

  TH1 *hM=(TH1*)arr.At(1);
  hM->SetMarkerStyle(20);
  hM->SetMarkerSize(.5);
  hM->DrawClone("sames");

  TH1 *hS=(TH1*)arr.At(2);
  hS->SetMarkerStyle(20);
  hS->SetMarkerSize(.5);
  hS->SetMarkerColor(kRed);
  hS->SetLineColor(kRed);
  hS->DrawClone("same");

  l.SetLineColor(kBlack);
  l.DrawLine(.2,0,20,0);
  l.SetLineColor(kRed);
  l.DrawLine(.2,1,20,1);


  return 0; 

}
