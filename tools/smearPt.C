
void smearPt(Int_t N = 1E6.)
{
  Float_t min = -10.;
  Float_t max = 10.;
  UInt_t seed = 12345678;
  TRandom3 * rnd3 = new TRandom3(seed);

  TCanvas * c1 = new TCanvas();

  //1) Histogram some gaussian random numbers
  TH1D * hist = new TH1D("hist","hist", 1000, min, max);
  for (int j=0; j<N; j++){    
    Float_t x = rnd3->Gaus(0., 1.);
    hist->Fill(x);
  }
  c1->cd();
  hist->Draw();
}
