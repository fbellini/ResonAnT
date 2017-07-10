void checkSimulatedProd()
{

  //Printf("check Kinematics.root: %s", (CheckGenerated()?"Not OK":"OK"));
  // Printf("check decays: %s", (checkDecay()?"Not OK":"OK"));
  Printf("check decays: %s", (checkAliDecay()?"Not OK":"OK"));
  return;
}


Int_t CheckGenerated()
{
  TFile * fin = TFile::Open("Kinematics.root","read");
  if (!fin) return 1;

  TH1I * hallpdg = new TH1I("hallpdg","hallpdg", 1000, 0, 1000);
  TH1I * hpdg = new TH1I("hpdg","hpdg", 5, 1, 6);
  hpdg->GetXaxis()->SetBinLabel(1, "#phi(1020)");
  hpdg->GetXaxis()->SetBinLabel(2, "#Lambda(1520)");
  hpdg->GetXaxis()->SetBinLabel(3, "#bar{#Lambda}(1520)");
  hpdg->GetXaxis()->SetBinLabel(4, "f0");
  hpdg->GetXaxis()->SetBinLabel(5, "f2");

  for (Int_t evCounter = 0; evCounter<20; evCounter++){
    TDirectoryFile * dir = (TDirectoryFile*) fin->Get(Form("Event%i",evCounter));
    if (!dir) return 2;
    TTree * tree = (TTree*) dir->Get("TreeK");
    if (!tree) return 3;

    Int_t pdgcode = -1;
    TParticle * part = new TParticle();
    tree->SetBranchAddress("Particles",&part);

    for (Int_t ip = 0; ip<tree->GetEntries(); ip++) {
      tree->GetEntry(ip);

      Int_t pdgcode = part->GetPdgCode();
      switch (pdgcode) {
      case 333 :
	hpdg->Fill(1.5);
	break;
      case 3124 :
	hpdg->Fill(2.5);
	break;
      case -3124 :
	hpdg->Fill(3.5);
	break;
      case 9010221:
	hpdg->Fill(4.5);
	break;
      case 225 :
	hpdg->Fill(5.5);
	break;
      default;
      }
      hallpdg->Fill(pdgcode);
    }
  }

  TCanvas * c = new TCanvas("cpdg","cpdg", 700, 500);
  c->cd();
  hpdg->Draw();
  TCanvas * call = new TCanvas("callpdg","callpdg", 700, 500);
  call->cd();
  hallpdg->Draw();

  return 0;
}


Int_t checkDecay(Int_t checkpdg=9010221)
{
  TFile * fin = TFile::Open("Kinematics_002.root","read");
  if (!fin) return 1;

  TH1I * hdecayed = new TH1I("hdecayed","hdecayed", 5, 1, 6);
  hdecayed->GetXaxis()->SetBinLabel(1, "#phi(1020)");
  hdecayed->GetXaxis()->SetBinLabel(2, "#Lambda(1520)");
  hdecayed->GetXaxis()->SetBinLabel(3, "#bar{#Lambda}(1520)");
  hdecayed->GetXaxis()->SetBinLabel(4, "f0");
  hdecayed->GetXaxis()->SetBinLabel(5, "f2");

  for (Int_t evCounter = 0; evCounter<20; evCounter++){
    TDirectoryFile * dir = (TDirectoryFile*) fin->Get(Form("Event%i",evCounter));
    if (!dir) return 2;
    TTree * tree = (TTree*) dir->Get("TreeK");
    if (!tree) return 3;

    Int_t pdgcode = -1;
    TParticle * part = new TParticle();
    tree->SetBranchAddress("Particles",&part);

    for (Int_t ip = 0; ip<tree->GetEntries(); ip++) {
      tree->GetEntry(ip);
      Int_t pdgcode = part->GetPdgCode();
      Int_t daughters[2] = { part->GetDaughter(0), part->GetDaughter(1)};
      //do they decay?
      if (pdgcode==333 && daughters[0]>0 && daughters[1]>0) hdecayed->Fill(1.5);
      if (pdgcode==3124 && daughters[0]>0 && daughters[1]>0) hdecayed->Fill(2.5);
      if (pdgcode==-3124 && daughters[0]>0 && daughters[1]>0) hdecayed->Fill(3.5);
      if (pdgcode==9010221 && daughters[0]>0 && daughters[1]>0) hdecayed->Fill(4.5);
      if (pdgcode==225 && daughters[0]>0 && daughters[1]>0) hdecayed->Fill(5.5);
    }
  }
  TCanvas * cd2k = new TCanvas("cd2k","cd2k", 700, 500);
  cd2k->cd();
  hdecayed->Draw();
  return 0;
}


Int_t checkAliDecay(Int_t checkpdg=9010221)
{

  // TFile * fin = TFile::Open("Kinematics_002.root","read");
  // if (!fin) return 1;

  TChain * chain = new TChain("esdTree");
  chain->Add("AliESDs.root");

  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  AliAnalysisTaskSE * task = new AliAnalysisTaskSE("mytask");
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("input0",
                      TTree::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("output0",
                      TList::Class(), AliAnalysisManager::kOutputContainer);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,0,coutput1);
  mgr->AddTask(task);
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
  }

  AliInputEventHandler *mceh = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mceh);

  AliMCEvent* mcEvent = 0x0;

  for (Int_t evCounter = 0; evCounter<20; evCounter++){
    if (!mceh->LoadEvent(evCounter)) {Printf("Failed loading event %i", evCounter); continue;}
    mcEvent = ( AliMCEvent* ) mceh->MCEvent();
    if (!mcEvent) {Printf("Failed retrieve MC event"); continue;}
    AliStack* stack = 0x0;
    stack = mcEvent->Stack();
    if (!stack) {Printf("Failed retrieve stack"); continue;}
    cout << "Ntrk = " << stack->GetNtrack() << endl;

    continue;

    Int_t pdgcode = -1;
    TParticle * part = new TParticle();
    tree->SetBranchAddress("Particles",&part);

    for (Int_t ip = 0; ip<tree->GetEntries(); ip++) {
      tree->GetEntry(ip);
      Int_t pdgcode = part->GetPdgCode();
      Int_t daughters[2] = { part->GetDaughter(0), part->GetDaughter(1)};

      //check specific decayed particle
      if (checkpdg>0 && pdgcode==checkpdg) {
	Printf(":::: Particle %i found with PDG = %i", ip, pdgcode);
	Printf("daughter[0] = %i        daughter[1] = %i", daughters[0], daughters[1]);
	TParticle * d1 = (TParticle *) tree->Particle(daughters[0]);
	TParticle * d2 = (TParticle *) tree->Particle(daughters[1]);
	if (!d1) {Printf("First daughter  %i not found",daughters[0]); continue;}
	if (!d2) {Printf("Second daughter  %i not found",daughters[1]); continue;}
	Printf("Particle = %i PDG = %i      PDG_d1 = %i      PDG_d2 = %i", ip, pdgcode, d1->GetPdgCode(), d2->GetPdgCode());
	Printf("d1->GetMother() %i", d1->GetFirstMother());
	Printf("d2->GetMother() %i", d2->GetFirstMother());
      }
    }
  }

  return 0;
}
