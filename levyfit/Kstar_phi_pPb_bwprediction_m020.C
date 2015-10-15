void Kstar_phi_pPb_bwprediction(Int_t cb = 0 //centrality bin
				)
{
  myOptions();
  gROOT->ForceStyle();
  
  double x,unorm;
  int j;

  TCanvas* c=new TCanvas("c","",10,10,400,500);
  c->SetFillColor(0);
  c->Draw();
  c->cd();
  c->cd();

  TPad* p1=new TPad("p1","",0.,0.4,1,1);
  myPadSetUp(p1,0.19,0.01,0.01,0.);
  p1->SetFillColor(0);
  p1->SetLogy();
  p1->Draw();
  p1->SetTicky(1);

  TPad* p2=new TPad("p2","",0.,0.,1,0.4);
  myPadSetUp(p2,0.19,0.,0.01,0.27);
  p2->SetFillColor(0);
  p2->Draw();
  p2->SetTicky(1);

  double xmax=3.05;
  double tr=1.5;

  TH1F* h1=new TH1F("h1","",1,0.,xmax);
  h1->SetXTitle("");
  h1->GetXaxis()->SetTitleSize(0.14);
  h1->GetXaxis()->SetTitleOffset(0.35);
  h1->GetXaxis()->SetLabelSize(0.08);
  h1->SetNdivisions(509,"x");
  h1->SetYTitle("d^{2}#it{N}/(d#it{p}_{T}d#it{y}) (GeV/#it{c})^{-1}");
  h1->GetYaxis()->SetTitleSize(0.08);
  h1->GetYaxis()->SetTitleOffset(1.15);
  h1->GetYaxis()->SetLabelSize(0.08);
  h1->GetYaxis()->SetLabelOffset(0.009);
  h1->SetNdivisions(505,"y");

  TH1F* h2=new TH1F("h2","",1,0.,xmax);
  h2->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  h2->GetXaxis()->SetTitleSize(tr*0.08);
  h2->GetXaxis()->SetTitleOffset(1.);
  h2->GetXaxis()->SetLabelSize(tr*0.08);
  h2->SetNdivisions(509,"x");
  h2->SetYTitle("Data/Prediction  ");
  h2->GetYaxis()->SetTitleSize(tr*0.08);
  h2->GetYaxis()->SetTitleOffset(0.81);
  h2->GetYaxis()->SetLabelSize(tr*0.08);
  h2->GetYaxis()->SetLabelOffset(tr*0.01);
  h2->SetNdivisions(505,"y");
  h2->SetMinimum(0.);
  h2->SetMaximum(2.15);

  TLine* one=new TLine(0.,1.,xmax,1.); one->SetLineStyle(2);

  TF1* dummy=new TF1("dummy","pol0",0.,1.);
  dummy->SetLineColor(1);
  dummy->SetLineWidth(2);

  //-----
  //Get Kstar spectrum
  //-----
  int nk=12;

  TGraphErrors* ks_data = new TGraphErrors(nk);
  TBox* ks_data_sys[20];
  get_Kstar_spectrum(ks_data,ks_data_sys,nk);
  ks_data->SetLineColor(2); ks_data->SetMarkerColor(2); ks_data->SetMarkerStyle(21);
  for(j=0;j<nk;j++){ks_data_sys[j]->SetLineColor(2); ks_data_sys[j]->SetFillStyle(0);}

  //-----
  //Get Kstar prediction
  //-----
  TF1* ks_model=new TF1("ks_model",x_blast,0.01,3.05,5);
  get_Kstar_prediction(ks_model,unorm);
  ks_model->SetLineColor(2);
  ks_model->SetLineStyle(2);

  TF1* ks_model2=(TF1*) ks_model->Clone("ks_model2");
  ks_model2->SetLineStyle(1);
  ks_model2->SetRange(0.3,xmax);

  TGraphErrors* ks_norm=new TGraphErrors(306);
  for(j=0;j<ks_norm->GetN();j++){
    x=0.005+0.01*j;
    ks_norm->SetPoint(j,x,ks_model->Eval(x));
    ks_norm->SetPointError(j,0.005,unorm*ks_model->Eval(x));
  }
  ks_norm->SetFillColor(TColor::GetColor("#ddaaaa"));

  TBox* ks_norm2=new TBox(0.,1.-unorm,0.2,1.+unorm);
  ks_norm2->SetFillColor(TColor::GetColor("#ddaaaa"));

  //-----
  //Get spectrum/prediction ratio
  //-----
  TGraphErrors* kt2=new TGraphErrors(nk);
  TBox* ky2[20];
  get_Kstar_ratio(kt2,ky2, nk);
  kt2->SetLineColor(2); kt2->SetMarkerColor(2); kt2->SetMarkerStyle(21);
  for(j=0;j<nk;j++){ky2[j]->SetLineColor(2); ky2[j]->SetFillStyle(0);}


  //-----

  int np=10;

  TGraphErrors* pt1=new TGraphErrors(np);
  TBox* py1[20];
  get_phi_spectrum(pt1,py1);
  pt1->SetLineColor(4); pt1->SetMarkerColor(4); pt1->SetMarkerStyle(4);
  for(j=0;j<np;j++){py1[j]->SetLineColor(4); py1[j]->SetFillStyle(0);}
   
  TF1* pf1=new TF1("pf1",x_blast,0.01,3.05,5);
  get_phi_prediction(pf1,unorm);
  pf1->SetLineColor(4);
  pf1->SetLineStyle(2);

  TF1* pf2=(TF1*) pf1->Clone("pf2");
  pf2->SetLineStyle(1);
  pf2->SetRange(0.3,xmax);

  TGraphErrors* pn1=new TGraphErrors(306);
  for(j=0;j<pn1->GetN();j++){
    x=0.005+0.01*j;
    pn1->SetPoint(j,x,pf1->Eval(x));
    pn1->SetPointError(j,0.005,unorm*pf1->Eval(x));
  }
  pn1->SetFillColor(33);

  TBox* pn2=new TBox(0.2,1.-unorm,0.4,1.+unorm);
  pn2->SetFillColor(33);

  TGraphErrors* pt2=new TGraphErrors(np);
  TBox* py2[20];
  get_phi_ratio(pt2,py2);
  pt2->SetLineColor(4); pt2->SetMarkerColor(4); pt2->SetMarkerStyle(4);
  for(j=0;j<np;j++){py2[j]->SetLineColor(4); py2[j]->SetFillStyle(0);}

  //-----

  p1->cd();

  h1->SetMinimum(0.007); h1->SetMaximum(0.7);
  h1->Draw();

  ks_norm->Draw("3same");
  pn1->Draw("3same");
  ks_model->Draw("same");
  ks_model2->Draw("same");
  pf1->Draw("same");
  pf2->Draw("same");
  for(j=0;j<nk;j++) ks_data_sys[j]->Draw();
  for(j=0;j<np;j++) py1[j]->Draw();
  ks_data->Draw("pzsame");
  pt1->Draw("pzsame");
  p1->RedrawAxis();

  TLatex* t1=new TLatex(0.96,0.94,"V0A Multiplicity Event Class 0-20% (Pb Side)");
  t1->SetTextAlign(31);
  t1->SetNDC();
  t1->SetTextSize(0.05);
  t1->Draw();

  TLatex* t2=new TLatex(0.96,0.88,"ALICE Preliminary");
  t2->SetTextAlign(31);
  t2->SetNDC();
  t2->SetTextSize(0.05);
  // t2->Draw();

  TLegend* l1=new TLegend(0.22,0.15,0.58,0.25);
  l1->SetNColumns(2);
  myLegendSetUp(l1,0.08);
  l1->AddEntry(ks_data," K*^{0}","p");
  l1->AddEntry(pt1," #phi","p");
  l1->Draw();

  TLegend* l2=new TLegend(0.22,0.07,0.48,0.1);
  myLegendSetUp(l2,0.08);
  l2->AddEntry(dummy,"blast-wave predictions","l");
  l2->Draw();

  //-----

  p2->cd();

  h2->Draw();
  ks_norm2->Draw();
  pn2->Draw();
  one->Draw();
  for(j=0;j<nk;j++) ky2[j]->Draw();
  for(j=0;j<np;j++) py2[j]->Draw();
  kt2->Draw("pzsame");
  pt2->Draw("pzsame");
  p2->RedrawAxis();

  TLatex* t3=new TLatex(0.205,0.88,"ALICE, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, -0.5 < #it{y} < 0");
  t3->SetNDC();
  t3->SetTextSize(tr*0.05);
  t3->Draw();

  TLatex* t4=new TLatex(0.205,0.33,"uncertainties: stat.(bars), sys.(boxes), norm.(shaded)");
  t4->SetNDC();
  t4->SetTextSize(tr*0.05);
  t4->Draw();

  //-----
  c->Print("Kstar_phi_pPb_bw_m020.pdf");
  c->Print("Kstar_phi_pPb_bw_m020.png");
  //  c->Print("Kstar_phi_pPb_bw_m020.gif");
  c->Print("Kstar_phi_pPb_bw_m020.eps");

  return;
}


void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  return;
}

void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}

void myGraphSetUp(TGraphErrors *currentGraph=0, Float_t currentMarkerSize = 1.0,
		  int currentMarkerStyle=20, int currentMarkerColor=0,
		  int currentLineStyle=1, int currentLineColor=0){
  currentGraph->SetMarkerSize(currentMarkerSize);
  currentGraph->SetMarkerStyle(currentMarkerStyle);
  currentGraph->SetMarkerColor(currentMarkerColor);
  currentGraph->SetLineStyle(currentLineStyle);
  currentGraph->SetLineColor(currentLineColor);
  return;
}

void myOptions(Int_t lStat=0){
  // Set gStyle
  int font = 42;
  // From plain
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.0,"xyz");  
  gStyle->SetTitleSize(0.06,"xyz");  
  gStyle->SetMarkerSize(1); 
  gStyle->SetPalette(1,0); 
  if (lStat){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    }
  else {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }
}

void get_Kstar_spectrum(TGraphErrors* g,TBox** b, Int_t useFirstNbins = 12)
{
  double dx=0.1;

  TFile * fin = TFile::Open("$HOME/alice/resonances/kstar_pA5.02TeV/output_LF5455/Levy/finalCand_kstar_pPb_smoothSys.root");
  if (!fin) {Printf("Cannot open input file with spectrum"); return;}
  TH1D * hstat = (TH1D*) fin->Get("hKstar_0");
  if (!hstat) {Printf("Error: invalid histogram: hKstar_0"); return;}
  TH1D * hsys = (TH1D*) fin->Get("hKstar_0_sys");
  if (!hsys) {Printf("Error: invalid histogram: hKstar_0_sys"); return;}

  for (Int_t i = 0; i<useFirstNbins; i++) {
    Double_t bincenter = hsys->GetXaxis()->GetBinCenter(i+1);
    Double_t value = hstat->GetBinContent(i+1);
    Double_t statunc = hstat->GetBinError(i+1);
    Double_t sysunc = hsys->GetBinError(i+1);
    g->SetPoint(i, bincenter, value); g->SetPointError(i, dx, statunc);
    b[i] = new TBox(bincenter-dx, value-sysunc, bincenter+dx, value+sysunc);
    b[i]->SetLineColor(hstat->GetLineColor());
    b[i]->SetFillStyle(0);
  }

  return;
}


void get_Kstar_prediction(TF1* g,double& a){
  g->SetParameters(8.958100000e-01,5.769672186e+03,1.474500703e-01,1.160745429e+00,8.330513291e-01);
  a=8.044801887e-02;
  return;
}


void get_Kstar_ratio(TGraphErrors* g,TBox** b, Int_t useFirstNbins = 12){
  double dx=0.1;

  TFile * fin = TFile::Open("$HOME/alice/resonances/kstar_pA5.02TeV/output_LF5455/Levy/finalCand_kstar_pPb_smoothSys.root");
  if (!fin) {Printf("Cannot open input file with spectrum"); return;}
  TH1D * hstat = (TH1D*) fin->Get("hKstar_0");
  if (!hstat) {Printf("Error: invalid histogram: hKstar_0"); return;}
  TH1D * hsys = (TH1D*) fin->Get("hKstar_0_sys");
  if (!hsys) {Printf("Error: invalid histogram: hKstar_0_sys"); return;}
  TF1* ks_model=new TF1("ks_model",x_blast,0.01,3.05,5);
  double unorm;
  get_Kstar_prediction(ks_model,unorm);

  hstat->Sumw2();
  hsys->Sumw2();

  hstat->Divide(ks_model);
  hsys->Divide(ks_model);

  for (Int_t i = 0; i<useFirstNbins; i++) {
    Double_t bincenter = hsys->GetXaxis()->GetBinCenter(i+1);
    Double_t value = hstat->GetBinContent(i+1);
    Double_t statunc = hstat->GetBinError(i+1);
    Double_t sysunc = hsys->GetBinError(i+1);
    g->SetPoint(i, bincenter, value); g->SetPointError(i, dx, statunc);
    b[i] = new TBox(bincenter-dx, value-sysunc, bincenter+dx, value+sysunc);
    b[i]->SetLineColor(hstat->GetLineColor());
    b[i]->SetFillStyle(0);
  }
  return;
}

//   g->SetPoint(0,0.10,7.630541921e-01); g->SetPointError(0,0.,1.401433628e-01);
//   b[0]=new TBox(0.10-dx,6.734424233e-01,0.10+dx,8.526659577e-01);
//   g->SetPoint(1,0.30,9.282789826e-01); g->SetPointError(1,0.,9.244309036e-02);
//   b[1]=new TBox(0.30-dx,8.537163287e-01,0.30+dx,1.002841638e+00);
//   g->SetPoint(2,0.50,8.237383366e-01); g->SetPointError(2,0.,5.294784265e-02);
//   b[2]=new TBox(0.50-dx,7.391020283e-01,0.50+dx,9.083746424e-01);
//   g->SetPoint(3,0.70,7.456531525e-01); g->SetPointError(3,0.,4.411097403e-02);
//   b[3]=new TBox(0.70-dx,6.711934432e-01,0.70+dx,8.201128604e-01);
//   g->SetPoint(4,0.90,9.045421481e-01); g->SetPointError(4,0.,4.652333243e-02);
//   b[4]=new TBox(0.90-dx,8.322859034e-01,0.90+dx,9.767983970e-01);
//   g->SetPoint(5,1.10,9.727353454e-01); g->SetPointError(5,0.,4.245869354e-02);
//   b[5]=new TBox(1.10-dx,8.961306736e-01,1.10+dx,1.049340019e+00);
//   g->SetPoint(6,1.30,9.296782017e-01); g->SetPointError(6,0.,4.021879334e-02);
//   b[6]=new TBox(1.30-dx,8.563731462e-01,1.30+dx,1.002983255e+00);
//   g->SetPoint(7,1.50,9.891349077e-01); g->SetPointError(7,0.,3.830089979e-02);
//   b[7]=new TBox(1.50-dx,9.084523618e-01,1.50+dx,1.069817451e+00);
//   g->SetPoint(8,1.70,9.670242667e-01); g->SetPointError(8,0.,3.791159758e-02);
//   b[8]=new TBox(1.70-dx,8.870233968e-01,1.70+dx,1.047025132e+00);
//   g->SetPoint(9,1.90,1.015060544e+00); g->SetPointError(9,0.,3.921354599e-02);
//   b[9]=new TBox(1.90-dx,9.132531583e-01,1.90+dx,1.116867934e+00);
//   g->SetPoint(10,2.25,1.168171525e+00); g->SetPointError(10,0.,2.688888595e-02);
//   b[10]=new TBox(2.25-dx,1.077770874e+00,2.25+dx,1.258572170e+00);
//   g->SetPoint(11,2.75,1.454425693e+00); g->SetPointError(11,0.,3.781229034e-02);
//   b[11]=new TBox(2.75-dx,1.348034367e+00,2.75+dx,1.560817013e+00);



void get_phi_spectrum(TGraphErrors* g,TBox** b){
  double dx=0.08,x,y;
  g->SetPoint(0,0.50,1.271879971e-01); g->SetPointError(0,0.,5.701878953e-03);
  b[0]=new TBox(0.50-dx,1.184087237e-01,0.50+dx,1.359672705e-01);
  g->SetPoint(1,0.70,1.464068890e-01); g->SetPointError(1,0.,2.726064742e-03);
  b[1]=new TBox(0.70-dx,1.334941499e-01,0.70+dx,1.593196280e-01);
  g->SetPoint(2,0.90,1.529591680e-01); g->SetPointError(2,0.,1.966534210e-03);
  b[2]=new TBox(0.90-dx,1.395936767e-01,0.90+dx,1.663246593e-01);
  g->SetPoint(3,1.10,1.498859972e-01); g->SetPointError(3,0.,1.816893637e-03);
  b[3]=new TBox(1.10-dx,1.377456058e-01,1.10+dx,1.620263886e-01);
  g->SetPoint(4,1.30,1.324973851e-01); g->SetPointError(4,0.,1.838837191e-03);
  b[4]=new TBox(1.30-dx,1.229800154e-01,1.30+dx,1.420147549e-01);
  g->SetPoint(5,1.50,1.148620397e-01); g->SetPointError(5,0.,1.882160245e-03);
  b[5]=new TBox(1.50-dx,1.066616969e-01,1.50+dx,1.230623825e-01);
  g->SetPoint(6,1.70,9.847348928e-02); g->SetPointError(6,0.,1.794616264e-03);
  b[6]=new TBox(1.70-dx,9.146592580e-02,1.70+dx,1.054810528e-01);
  g->SetPoint(7,1.90,7.767432928e-02); g->SetPointError(7,0.,1.533340896e-03);
  b[7]=new TBox(1.90-dx,7.218880672e-02,1.90+dx,8.315985184e-02);
  g->SetPoint(8,2.25,5.361141264e-02); g->SetPointError(8,0.,6.813266920e-04);
  b[8]=new TBox(2.25-dx,4.978517769e-02,2.25+dx,5.743764760e-02);
  g->SetPoint(9,2.75,2.912878618e-02); g->SetPointError(9,0.,4.014970092e-04);
  b[9]=new TBox(2.75-dx,2.710958081e-02,2.75+dx,3.114799154e-02);

  return;
}


void get_phi_prediction(TF1* g,double& a){
  g->SetParameters(1.019455000e+00,4.669971137e+03,1.474500703e-01,1.160745429e+00,8.330513291e-01);
  a=8.044801887e-02;
  return;
}


void get_phi_ratio(TGraphErrors* g,TBox** b){
  double dx=0.08;
  g->SetPoint(0,0.50,8.687492609e-01); g->SetPointError(0,0.,3.894631094e-02);
  b[0]=new TBox(0.50-dx,8.087830096e-01,0.50+dx,9.287155129e-01);
  g->SetPoint(1,0.70,8.891825080e-01); g->SetPointError(1,0.,1.655638661e-02);
  b[1]=new TBox(0.70-dx,8.107587248e-01,0.70+dx,9.676062875e-01);
  g->SetPoint(2,0.90,9.340893626e-01); g->SetPointError(2,0.,1.200920961e-02);
  b[2]=new TBox(0.90-dx,8.524691239e-01,0.90+dx,1.015709599e+00);
  g->SetPoint(3,1.10,9.941294789e-01); g->SetPointError(3,0.,1.205067567e-02);
  b[3]=new TBox(1.10-dx,9.136074707e-01,1.10+dx,1.074651484e+00);
  g->SetPoint(4,1.30,1.007074594e+00); g->SetPointError(4,0.,1.397647373e-02);
  b[4]=new TBox(1.30-dx,9.347357973e-01,1.30+dx,1.079413391e+00);
  g->SetPoint(5,1.50,1.041674376e+00); g->SetPointError(5,0.,1.706915628e-02);
  b[5]=new TBox(1.50-dx,9.673061445e-01,1.50+dx,1.116042609e+00);
  g->SetPoint(6,1.70,1.100156307e+00); g->SetPointError(6,0.,2.004964488e-02);
  b[6]=new TBox(1.70-dx,1.021867059e+00,1.70+dx,1.178445556e+00);
  g->SetPoint(7,1.90,1.097273946e+00); g->SetPointError(7,0.,2.166088967e-02);
  b[7]=new TBox(1.90-dx,1.019782178e+00,1.90+dx,1.174765715e+00);
  g->SetPoint(8,2.25,1.184418321e+00); g->SetPointError(8,0.,1.505231439e-02);
  b[8]=new TBox(2.25-dx,1.099886641e+00,2.25+dx,1.268950003e+00);
  g->SetPoint(9,2.75,1.341770172e+00); g->SetPointError(9,0.,1.849430611e-02);
  b[9]=new TBox(2.75-dx,1.248758763e+00,2.75+dx,1.434781580e+00);

  return;
}


Double_t blast(Double_t *x,Double_t *par){
  static TF1* fint=0;
  if(!fint) fint=new TF1("fint",blast_integrand,0.,1.,5);

  fint->SetParameters(par[0],par[2],par[3],par[4],x[0]);
  return fint->Integral(0.,1.)*par[1];
}

Double_t  x_blast(Double_t *x,Double_t *par){
  return x[0]*blast(x,par);
}


Double_t blast_integrand(const Double_t *x,const Double_t *par){
  Double_t x0=x[0];
  Double_t m=par[0];
  Double_t t=fabs(par[1]);
  Double_t n=par[2];
  Double_t beta_max=par[3];
  Double_t pt=par[4];

  //Keep beta within reasonable limits.
  Double_t beta=beta_max*TMath::Power(x0,n);
  if(beta>0.9999999999999999) beta=0.9999999999999999;

  Double_t mt=TMath::Sqrt(m*m+pt*pt);
  Double_t rho0=TMath::ATanH(beta);
  Double_t a0=pt*TMath::SinH(rho0)/t;
  if(a0>700.) a0=700.;//avoid floating point exception
  Double_t a1=mt*TMath::CosH(rho0)/t;
  return x0*mt*TMath::BesselI0(a0)*TMath::BesselK1(a1);
}
