/*--------------------------------------------
   //GetIntegral in a +/- 3sigma range from peak
   ----------------------------------------------*/
    // Define signal range in which events counts are to be defined
    x.setRange("fitrange", xrangem, xrangeM);
    x.setRange("extrapolate", 1.2, 1.30);
    x.setRange("fullRange", 0.63, 1.30);
    // Associated nsig as expected number of events
    RooRealVar sevts("sevts","number of interesting events",histo_integral*0.5, 0., histo_integral);
    RooExtendPdf esig("esig","extended signal p.d.f",model,sevts,"fullRange");

    // Perform unbinned extended ML fit to data in fitrange
    RooFitResult* r2 =  esig.chi2FitTo(data, Extended(kTRUE), SumW2Error(kFALSE), Save(), Range("fitrange")); 

    //Get fitted function
    RooArgList pars(* esig.getParameters(RooArgSet(x)));   
    RooArgSet prodSet(esig); //prodSet.add(nsig);
    RooProduct unNormPdf("fitted Function", "fitted Function", prodSet);
    TF1 * f2 = unNormPdf.asTF(RooArgList(x), pars);
    
    float nsig1 = ((RooRealVar*) pars.find("sevts"))->getVal();
    float dnsig1 = ((RooRealVar*) pars.find("sevts"))->getError();
    // float mean1 = ((RooRealVar*) pars.find("mean"))->getVal();
    // float dmean1 = ((RooRealVar*) pars.find("mean"))->getError();
    // float sigma1 = ((RooRealVar*) pars.find("sigma"))->getVal();
    // float dsigma1 = ((RooRealVar*) pars.find("sigma"))->getError();
    // cout << " nsig = " << nsig1 << " +- " << dnsig1 << endl;
    // cout << " mean = " << mean1 << " +- " << dmean1 << endl;
    // cout << " sigma = " << sigma1 << " +- " << dsigma1 << endl;
    Double_t integ2_full = f2->Integral(0.63,1.30);
    Double_t integ2 = nsig1*f2->Integral(xrangem, xrangeM, 0)/integ2_full;
    new TCanvas(); f2->DrawClone();
    Double_t dinteg2 = nsig1*f2->IntegralError(xrangem, xrangeM, 0, r2->covarianceMatrix().GetMatrixArray())/integ2_full;
