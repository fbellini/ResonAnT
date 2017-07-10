void checkRsnQA(TString filename="/Users/fbellini/alice/resonances/kstar_pA5.02TeV/LF_pPb17-21/train1920.root"){}

TH1 * projectCutPlot(TString filename="/Users/fbellini/alice/resonances/kstar_pA5.02TeV/LF_pPb17-21/train1920.root",
                    TString foldername="RsnOut_cut5",
                    TString plot="dcaxy",
                    TString projopt="y",
                    TString newname = "cut5")
{
    TFile * fin = TFile::Open(filename.Data());
    if (!fin) return 0x0;
    TList * list = (TList*) fin->Get(foldername.Data());
    if (!list) return 0x0;
    TH2F * hist = (TH2F*) list->FindObject(getPlotName(plot.Data(),5));
    if (!hist) return;
    
    TH1F * proj = 0x0;
    if (projopt.Contains("x")) proj = (TH1F*) hist->ProjectionX(Form("%s_px",newname.Data()));
    else proj = (TH1F*) hist->ProjectionY(Form("%s_py",newname.Data()));

    return proj;
    
}

TString getPlotName(TString what="dcaxy", Int_t bit=5)
{
    
    TString speciename = "cutQ";
    TString suffix;
/*
    if (what.Contains("mom")) suffix = "P_p";
    if (what.Contains("pt")) suffix = "Pt_pt"; // at: 0x7fdbd2cfaa00
    if (what.Contains("eta")) suffix = "Eta_eta"; // at: 0x7fdbd2cfadd0
    if (what.Contains("phi")) suffix = "Phi_phi"; // at: 0x7fdbd2cfb1a0
    if (what.Contains("tpcout")) suffix = "PhiOuterTPC_phiOuterTPC"; // at: 0x7fdbd2cfb570
*/
    if (what.Contains("phi")) suffix = "PhiVsPt_pt_phi"; // at: 0x7fdbd2cfb940
    if (what.Contains("dcaxy")) suffix = "DCAxyVsPt_pt_DCAxy"; // at: 0x7fdbd4800000
    if (what.Contains("dcaz")) suffix = "DCAzVsP_p_DCAz"; // at: 0x7fdbd48003f0
    if (what.Contains("itscls")) suffix = "ITSclsVsPt_pt_ITScls"; // at: 0x7fdbd48007e0
//    if (what.Contains("tpccls")) suffix = "TPCclsVsPt_pt_TPCcls"; // at: 0x7fdbd4800bd0
    if (what.Contains("tpccls")) suffix = "TPCclsVsPtpc_pTPC_TPCcls"; // at: 0x7fdbd4800fc0
    if (what.Contains("crossedrows")) suffix = "TPCcrossedRowsVsPtpc_pTPC_TPCcrossedRows"; // at: 0x7fdbd48013b0
//    if (what.Contains("")) suffix = "TPCcrossedRows2FclsVsPtpc_pTPC_TPCcrossedRows2Fcls"; // at: 0x7fdbd48017a0
    if (what.Contains("itschi2")) suffix = "ITSchi2VsPt_pt_ITSchi2"; // at: 0x7fdbd4801b90
//    if (what.Contains("tpcchi2")) suffix = "TPCchi2VsPt_pt_TPCchi2"; // at: 0x7fdbd4802350
    if (what.Contains("tpcchi2")) suffix = "TPCchi2VsPtpc_pTPC_TPCchi2"; // at: 0x7fdbd4802740
    if (what.Contains("dedx")) suffix = "dEdx_VsPtpc_pTPC_sTPC"; // at: 0x7fdbd4802e60
    if (what.Contains("nstpcpi")) suffix = "TPC_nsigmaPi_VsPtpc_pTPC_pi"; // at: 0x7fdbd4803250
    if (what.Contains("nstpck")) suffix = "TPC_nsigmaK_VsPtpc_pTPC_K"; // at: 0x7fdbd4803640
    if (what.Contains("nstpcpro")) suffix = "TPC_nsigmaPro_VsPtpc_pTPC_p"; // at: 0x7fdbd4803a30
    if (what.Contains("allsigmapi")) suffix = "TOFnsigmaPi_TPCnsigmaPi_pi_pi"; // at: 0x7fdbd4803e20
    if (what.Contains("allsigmak")) suffix = "TOFnsigmaK_TPCnsigmaK_K_K"; // at: 0x7fdbd4804240
    if (what.Contains("allsigmapro")) suffix = "TOFnsigmaP_TPCnsigmaP_p_p"; // at: 0x7fdbd48046b0
    if (what.Contains("nstofpi")) suffix = "TOF_nsigmaPi_vsP_p_pi"; // at: 0x7fdbd4804b20
    if (what.Contains("nstofk")) suffix = "TOF_nsigmaK_vsP_p_K"; // at: 0x7fdbd4804f60
    if (what.Contains("nstofpro")) suffix = "TOF_nsigmaPro_vsP_p_p"; // at: 0x7fdbd48053a0
    if (what.Contains("deltapi")) suffix = "TOF_deltaPi_vsP_p_Dpi"; // at: 0x7fdbd48057e0
    if (what.Contains("deltak")) suffix = "TOF_deltaK_vsP_p_DK"; // at: 0x7fdbd4805c20
    if (what.Contains("deltapro")) suffix = "TOF_deltaPro_vsP_p_Dp"; // at: 0x7fdbd4806060

    //cutK24_2.0sigma.TOF_deltaPro_vsP_p_Dp
    TString name = Form("%s_bit%i.%s", speciename.Data(), bit, suffix.Data());
    return name;
    
}
