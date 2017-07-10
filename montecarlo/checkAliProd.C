void checkAliProd()
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
      case 9010221ls
 :
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

Int_t checkAliDecay(Int_t checkpdg=9010221)
{

  // open run loader and load gAlice, kinematics and header
  TString gAliceFileName("galice.root");
  AliRunLoader* runLoader = AliRunLoader::Open(gAliceFileName);
  if (!runLoader) {
    Error("CheckESD", "getting run loader from file %s failed", 
	    gAliceFileName);
    return kFALSE;
  }
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  if (!gAlice) {
    Error("CheckESD", "no galice object found");
    return kFALSE;
  }
  runLoader->LoadKinematics();
  runLoader->LoadHeader();
  
  // open the ESD file
  TString esdFileName = "AliESDs.root";
  TFile* esdFile = TFile::Open(esdFileName);
  if (!esdFile || !esdFile->IsOpen()) {
    Error("CheckESD", "opening ESD file %s failed", esdFileName);
    return kFALSE;
  }
  
  AliESDEvent * esd = new AliESDEvent;
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("CheckESD", "no ESD tree found");
    return kFALSE;
  }
  esd->ReadFromTree(tree);


  // loop over events
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    runLoader->GetEvent(iEvent);

    // select simulated primary particles, V0s and cascades
    AliStack* stack = runLoader->Stack();
    Int_t nParticles = stack->GetNtrack();
    TArrayF vertex(3);
    runLoader->GetHeader()->GenEventHeader()->PrimaryVertex(vertex);
    TObjArray selParticles;
    TObjArray selV0s;
    TObjArray selCascades;
    for (Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
      TParticle* particle = stack->Particle(iParticle);
      if (!particle) continue;
      if (particle->Pt() < 0.001) continue;
      if (TMath::Abs(particle->Eta()) > 0.8) continue;
      TVector3 dVertex(particle->Vx() - vertex[0], 
		       particle->Vy() - vertex[1],
		       particle->Vz() - vertex[2]);
      if (dVertex.Mag() > 0.0001) continue;
      Int_t pdgcode = particle->GetPdgCode();

      if ( (checkpdg>0) && (pdgcode!=checkpdg) ) continue;
	TParticle* particle = stack->Particle(iParticle);
	if (!particle) continue;
	if (particle->Pt() < 0.001) continue;
	if (TMath::Abs(particle->Eta()) > 0.8) continue;
	TVector3 dVertex(particle->Vx() - vertex[0], 
			 particle->Vy() - vertex[1],
			 particle->Vz() - vertex[2]);
	if (dVertex.Mag() > 0.0001) continue;
	Int_t pdgcode = particle->GetPdgCode();

	if ( (checkpdg>0) && (pdgcode!=checkpdg) ) continue;
	// switch (TMath::Abs(pdgcode)) {
	// case kElectron:
	// case kMuonMinus:
	// case kPiPlus:
	// case kKPlus:
	// case kProton: 
	// case kPhi:
	// case kLambdaStar :
	// case f0 :
	// default: break;
	// }
	Int_t daughters[2] = { particle->GetDaughter(0), particle->GetDaughter(1) }; 
      
	//check specific decayed particle
	Printf(":::: Particle %i found with PDG = %i", ip, pdgcode);
	Printf("daughter[0] = %i        daughter[1] = %i", daughters[0], daughters[1]); 
	TParticle * d1 = (TParticle *) tree->Particle(daughters[0]);
	TParticle * d2 = (TParticle *) tree->Particle(daughters[1]);

	// switch (TMath::Abs(pdgcode)) {
	// case kElectron:
	// case kMuonMinus:
	// case kPiPlus:
	// case kKPlus:
	// case kProton: 
	// case kPhi:
	// case kLambdaStar :
	// case f0 :
	// default: break;
	// }
	Int_t daughters[2] = { particle->GetDaughter(0), particle->GetDaughter(1) }; 
      
	//check specific decayed particle
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
  return 0;
}
