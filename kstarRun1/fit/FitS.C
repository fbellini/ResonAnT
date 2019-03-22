enum EFit { kBW = 0,
	    kFitRelBW,
	    kFitRelBWpoly2,
	    kFitRelBWpoly3,
	    kFitRelBCpoly2 };

void FitS(TH1D * histo, Float_t Pt = -1.0, Double_t xmin = 0.8, Double_t xmax = 1.0, Int_t fitfcnID = 1, TString opt = "MWHI")
{
  if (!histo) return;
  if (Pt<0.0) return;
  Double_t pdgMass = 0.89545;
  Double_t pdgWidth = 0.0478;
   
  TF1 *fitfcn; Int_t parN = -1;
  switch (fitfcnID) {
  case kFitRelBW :
    parN = 6;
    fitfcn = new TF1("fitfcn",FitFunRelBW,xmin,xmax,parN); 
    fitfcn->SetParNames("Const","Slope","Yield","Mass","Width","Pt");
    fitfcn->SetParameters(60000,-50000,10000,pdgMass,pdgWidth,Pt); 
    fitfcn->FixParameter(4, pdgWidth);//width fix to PDG
    //funFit->SetParLimits(4,0.03,0.08);//constraint with in range
    fitfcn->SetParLimits(3, pdgMass*0.5, pdgMass*1.5);
    fitfcn->FixParameter(5, Pt);
    break;

  case kFitRelBWpoly2 :
    parN = 7;
    fitfcn = new TF1("fitfcn",FitFunRelBWPoly2,xmin,xmax, parN); 
    fitfcn->SetParNames("Yield","Mass","Width","Pt", "Const","p1","p2");
    fitfcn->SetParameters(histo->Integral()*0.5, pdgMass, pdgWidth, Pt, -50000, 0.6, -35000); 
    fitfcn->SetParLimits(0, histo->Integral()*0.01, histo->Integral());
    fitfcn->SetParLimits(1, pdgMass*0.5, pdgMass*1.5);
    fitfcn->FixParameter(2, pdgWidth);//width fix to PDG
    fitfcn->FixParameter(3, Pt);  
    // fitfcn->SetParLimits(4, -5.e5, 1.e7);
    // fitfcn->SetParLimits(5, 2., 2.);
    fitfcn->SetParLimits(6, -1.e6, 1.e3);

    // fitfcn = new TF1("fitfcn",FitFunRelBWPoly2,xmin,xmax,7); 
    // fitfcn->SetParNames("Const","p1","p2","Yield","Mass","Width","Pt");
    // fitfcn->SetParameters(60000,-50000,10000,-500, pdgMass,pdgWidth,Pt); 
    // fitfcn->FixParameter(5, pdgWidth);//width fix to PDG
    // //funFit->SetParLimits(4,0.03,0.08);//constraint with in range
    // fitfcn->SetParLimits(4, pdgMass*0.5, pdgMass*1.5);
    // fitfcn->FixParameter(6, Pt);
    break;

  case kFitRelBWpoly3 :
    parN = 8;
    fitfcn = new TF1("fitfcn",FitFunRelBWPoly3,xmin,xmax,parN); 
    fitfcn->SetParNames("Yield","Mass","Width","Pt", "Const","p1","p2","p3");
    fitfcn->SetParameters(histo->Integral()*0.5, pdgMass, pdgWidth, Pt, -50000, 150000, -35000, 10000); 
    fitfcn->FixParameter(2, pdgWidth);//width fix to PDG
    fitfcn->SetParLimits(1, pdgMass*0.5, pdgMass*1.5);
    fitfcn->FixParameter(3, Pt);  
    // fitfcn = new TF1("fitfcn",FitFunRelBWPoly2,xmin,xmax,7); 
    // fitfcn->SetParNames("Const","p1","p2","Yield","Mass","Width","Pt");
    // fitfcn->SetParameters(60000,-50000,10000,-500, pdgMass,pdgWidth,Pt); 
    // fitfcn->FixParameter(5, pdgWidth);//width fix to PDG
    // //funFit->SetParLimits(4,0.03,0.08);//constraint with in range
    // fitfcn->SetParLimits(4, pdgMass*0.5, pdgMass*1.5);
    // fitfcn->FixParameter(6, Pt);
    break;
  
  case kFitRelBCpoly2 :
    
    break;

  default :
    fitfcn = new TF1("fitfcn", KstarFun, xmin, xmax,4); 
    fitfcn->SetParNames("Yield","Mass","Width","Pt");
    fitfcn->SetParameters(histo->Integral()*0.5, pdgMass, pdgWidth, Pt); 
    fitfcn->SetParLimits(0, 0.0, histo->Integral());//constraint with in range
    fitfcn->SetParLimits(1, pdgMass*0.5, pdgMass*1.5);//constraint with in range
    fitfcn->FixParameter(2, pdgWidth);//width fix to PDG
    fitfcn->FixParameter(3, Pt);
    parN = 4;
  }

  fitfcn->SetNpx(500);
  fitfcn->SetLineWidth(2);
  fitfcn->SetLineColor(kBlue+2);
  
  TFitResultPtr result = histo->Fit("fitfcn","SEMBR","ep"); //add V=verbose or Q=quiet 
  Int_t fitstatus = result;
  result->Print("V");

  //mass parameter & error
  Double_t fMass   [2]; 
  fMass[0]=result->Parameter(1); 
  fMass[1]=result->ParError(1);
  
  //width parameter & error
  Double_t fGamma  [2]; 
  fGamma[0]=result->Parameter(2); 
  fGamma[1]=result->ParError(2);
  
  //Chi2 and dof
  Double_t fChi2   [2];
  fChi2[0]=fitfcn->GetChisquare(); 
  fChi2[1]=fitfcn->GetNDF();

  //total fit function parameters & error
  Double_t fFcn[7]; Double_t fFcnErr[7];  
  for (Int_t j=0;j<parN;j++) {
    fFcn[j]=result->Parameter(j);
    fFcnErr[j]=result->ParError(j);
    Printf("Fit param %i = %s = %e +/- %e", j, fitfcn->GetParName(j) , fFcn[j], fFcnErr[j]);
  }

  TF1 *backFcn;
  if (fitfcnID == kFitRelBW)       backFcn = new TF1("backFcn", Poly1, xmin, xmax, 2);//,"myFitFcn", "Bg");
  if (fitfcnID == kFitRelBWpoly2 ) backFcn = new TF1("backFcn", Poly2, xmin, xmax, 3);//, "myFitFcn", "Bg");
  if (fitfcnID == kFitRelBWpoly3 ) backFcn = new TF1("backFcn", Poly3, xmin, xmax, 4);//, "myFitFcn", "Bg");
  backFcn->SetParameters(&fFcn[4]);
  backFcn->SetParErrors(&fFcnErr[4]);
  backFcn->SetLineColor(kRed);
  backFcn->SetLineStyle(2);
  backFcn->SetLineWidth(2);
  backFcn->Draw("same"); 

  Float_t min2G = pdgMass - 2*pdgWidth;
  Float_t max2G = pdgMass + 2*pdgWidth;
  Int_t ibinmin2G=histo->GetXaxis()->FindBin(min2G);
  Int_t ibinmax2G=histo->GetXaxis()->FindBin(max2G);
  Int_t ibinmin=histo->GetXaxis()->FindBin(xmin);
  Int_t ibinmax=histo->GetXaxis()->FindBin(xmax);
  
  //integral of total fit function in fit range & error
  Double_t fFcnInt [2] ; 
  fFcnInt[0] = fitfcn->Integral(xmin, xmax);
  fFcnInt[1] = fitfcn->IntegralError(xmin, xmax);
  
  //integral of histogram in fit range & error
  Double_t fHistInt[2]= {0.0,0.0}; 
  fHistInt[0] = histo->IntegralAndError(ibinmin, ibinmax, fHistInt[1], "width");
  


  if (opt.Contains("M")) Printf("Mass          = %7.5f +- %7.5f (rel. err. %6.3f%%)", fMass  [0], fMass  [1], fMass  [1]*100./fMass  [0] );
  if (opt.Contains("W")) Printf("Width         = %7.5f +- %7.5f (rel. err. %6.3f%%)", fGamma [0], fGamma [1], fGamma [1]*100./fGamma [0]);
  //if (opt.Contains("S")) printf("Sigma         = %7.5f +- %7.5f", fSigma [0], fSigma [1]);
  if (opt.Contains("H")) Printf("Hist integral = %e +- %e (rel. err. %6.3f%%)", fHistInt[0], fHistInt[1], fHistInt[1]*100./fHistInt[0] );
  if (opt.Contains("I")) Printf("Func integral = %e +- %e (rel. err. %6.3f%%)", fFcnInt [0], fFcnInt [1], fFcnInt[1]*100./fFcnInt[0] );

  return;
}





Double_t Poly1(Double_t *x, Double_t *par)
{
  return par[0]+par[1]*x[0];
}
Double_t Poly2(Double_t *x, Double_t *par)
{
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}

Double_t Poly3(Double_t *x, Double_t *par)
{
  return par[0]*x[0]*x[0]*x[0] + par[1]*x[0]*x[0] + par[2]*x[0] + par[3];
}

Double_t arg(Double_t x0, Double_t x1, Double_t x2)
{
  return pow(x0*x0-x1*x1-x2*x2,2.0)-4.*x1*x1*x2*x2;
}

Double_t PS(Double_t m, Double_t pT, Double_t T)
{
  Double_t mT = sqrt(m*m+pT*pT);
  return m/mT*exp(-mT/T);
}

Double_t bw(Double_t m, Double_t m0, Double_t Gamma)
{
  return m*m0*Gamma/(pow(m*m-m0*m0,2.0)+m0*m0*Gamma*Gamma);
}

Double_t bw1(Double_t *x, Double_t *par)
{
  const Double_t MassK = 0.49368;
  const Double_t MassPi = 0.13957;
  Double_t Gamma = par[2]*pow(par[1]/x[0],4.0);
  Gamma *= pow(arg(x[0],MassK,MassPi)/arg(par[1],MassK,MassPi),1.5);
  //Double_t Gamma = par[2];
  return bw(x[0],par[1],Gamma);
}

Double_t bw2(Double_t *x, Double_t *par)
{
  const Double_t MassK = 0.49368;
  const Double_t MassPi = 0.13957;
  Double_t Gamma = par[2]*pow(par[1]/x[0],2.0);
  Gamma *= pow(arg(x[0],MassK,MassPi)/arg(par[1],MassK,MassPi),1.5);
  Gamma *= pow((pow(par[1]*LambdaPi,2)+arg(par[1],MassK,MassPi))/(pow(x[0]*LambdaPi,2)+arg(x[0],MassK,MassPi)),2);
  //Double_t Gamma = par[2];
  return bw(x[0],par[1],Gamma);
}

Double_t KstarFun(Double_t *x, Double_t *par)
{
  const Double_t Temp = 0.154;
  return par[0]*bw1(x, par)*PS(x[0], par[3], Temp)*1.e6;
  
}
Double_t FitFunRelBW(Double_t *x, Double_t *par)
{
  return KstarFun(x,&par[2])+Poly1(x,par);
}

Double_t FitFunRelBWPoly2(Double_t *x, Double_t *par)
{
  return KstarFun(x, par)+Poly2(x, &par[4]);
}

Double_t FitFunRelBWPoly3(Double_t *x, Double_t *par)
{
  return KstarFun(x, par)+Poly3(x, &par[4]);
}
