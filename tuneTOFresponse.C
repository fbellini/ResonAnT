void systTOFresponse(Float_t minSigma = -2.0, Float_t maxSigma=2.0)
{
  gStyle->SetOptTitle(0);
  tuneTOFresponse(0.9, minSigma,  maxSigma, kBlue, 0);
  tuneTOFresponse(0.8, minSigma,  maxSigma, kRed, 1 );
  tuneTOFresponse(0.85, minSigma,  maxSigma, kMagenta-4, 1);
  tuneTOFresponse(0.95, minSigma,  maxSigma, kGreen+1, 1);
  tuneTOFresponse(1.0, minSigma,  maxSigma, kBlack, 1);
  return;
} 
Float_t tuneTOFresponse(Float_t paramTail = 0.9, Float_t minSigma, Float_t maxSigma, Color_t color = kBlue, Bool_t drawsame = 0) 
{
  // TOF response
  gStyle->SetOptTitle(0);
  fTOFResponseF = new TF1("fTOFprob","[0]*TMath::Exp(-(x-[1])*(x-[1])/2/[2]/[2])* (x < [1]+[3]*[2]) + (x > [1]+[3]*[2])*[0]*TMath::Exp(-(x-[1]-[3]*[2]*0.5)*[3]/[2])",-7,7);
  fTOFResponseF->SetParameter(0,1);
  fTOFResponseF->SetParameter(1,-0.1);
  fTOFResponseF->SetParameter(2,1);
  fTOFResponseF->SetParameter(3, paramTail); //tail param
  fTOFResponseF->SetParameter(0,1./fTOFResponseF->Integral(-7,7));
  fTOFResponseF->SetLineColor(color);
  fTOFResponseF->GetXaxis()->SetTitle("n#sigma");
  Float_t int7sigma = fTOFResponseF->Integral(-7.,7.); //deve fare 1
  Float_t int3sigma = fTOFResponseF->Integral(-3.,3.); //-> ti da la percentuale per il 3sigma cut
  Float_t int2sigma = fTOFResponseF->Integral(-2.,2.); //-> ti da la percentuale per il 2sigma cut
  Float_t intNsigma = fTOFResponseF->Integral(minSigma, maxSigma);
  fTOFResponseF->Draw((drawsame?"same":"")); //disegna la risposta nel MC
  
  TString text7sigma = Form("I(-7.0, 7.0) = %6.3f",int7sigma);
  TString text3sigma = Form("I(-3.0, 3.0) = %6.3f",int3sigma);
  TString text2sigma = Form("I(-2.0, 2.0) = %6.3f",int2sigma);
  TString textNsigma = Form("#tau_{TOF}= %4.2f: I(%2.1f, %2.1f) = %6.3f",paramTail, minSigma,maxSigma,intNsigma);

  Printf("TOF response tune \n#################### \n exp tail starting point = %4.2f sigma", paramTail);
  Printf("%s \n %s \n %s \n %s\n",text7sigma.Data(), text3sigma.Data(), text2sigma.Data(), textNsigma.Data());
 
  TPaveText * pave = new TPaveText(0.12,0.18,0.45,0.25,"NDC");
  pave->SetBorderSize(0);
  pave->SetFillColor(kWhite);
  pave->SetTextColor(color);
  pave->InsertText(textNsigma.Data());
  pave->Draw();

  TLine * line = new TLine(paramTail, 0.01,paramTail, 0.55);
  line->SetLineStyle(7);
  line->SetLineWidth(1);
  line->SetLineColor(color);
  //for the most recent version:
  //$ALICE_ROOT/STEER/STEERBase/AliTOFPIDResponse.cxx
  line->Draw();

  // fTOFResponseF->GetXaxis()->SetRangeUser(-3.0, 7.0);
  // fTOFResponseF->GetYaxis()->SetRangeUser(1e-3, 1.0);
 
  return; 
}
