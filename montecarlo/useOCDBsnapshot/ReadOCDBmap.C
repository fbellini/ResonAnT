Int_t run = 122374;
ReadOCDBmap_REC(Char_t *filename, Char_t *snapname = NULL)
{
  TGrid::Connect("alien://");
  // get CDB map and list from AliESDs UserInfo
  TFile *fin = TFile::Open(filename);
  TTree *tin = (TTree *)fin->Get("esdTree");
  TList *ui = (TList *)tin->GetUserInfo();
  TMap *cdbm = (TMap *)ui->FindObject("cdbMap");
  TList *cdbl = (TList *)ui->FindObject("cdbList");

  ReadOCDBmap(cdbm, cdbl, snapname);
}
//"ocdb.sim.root"
ReadOCDBmap_SIM(Char_t *filename, Char_t *snapname = NULL )
{
  gSystem->Load("libPythia6.so");
  TGrid::Connect("alien://");
  // get CDB map and list from galice
  TFile *fin = TFile::Open(filename);
  TMap *cdbm = (TMap *)fin->Get("cdbMap");
  TList *cdbl = (TList *)fin->Get("cdbList");

  ReadOCDBmap(cdbm, cdbl, snapname);
}
  
ReadOCDBmap(TMap *cdbm, TList *cdbl, Char_t *snapname = NULL)
{
  
  cdbm->Print();
  cdbl->Print();

  if (!snapname) return;
  
  // instance CDB manager
  AliCDBManager *cdb = AliCDBManager::Instance();
  // set default storage
  TPair *pair = (TPair *)cdbm->FindObject("default");
  TObjString *ostr = pair->Value();
  cdb->SetDefaultStorage(ostr->String().Data());
  // set run
  cdb->SetRun(run);
  // get objects
  for (Int_t iobj = 0; iobj < cdbl->GetEntries(); iobj++) {
    
    TObjString *ostr = (TObjString *)cdbl->At(iobj);
    TObjArray *oa = ostr->String().Tokenize(";");

    // get path
    TObjString *ostr_path = (TObjString *)oa->At(0);
    TObjArray *oa_path = ostr_path->String().Tokenize(" ");
    ostr_path = (TObjString *)oa_path->At(1);
    TString str = ostr_path->String();
    str = str.Strip(TString::kLeading, '\"');
    str = str.Strip(TString::kTrailing, '\"');
    TString path = str.Data();  

    // get range
    TObjString *ostr_range = (TObjString *)oa->At(1);
    TObjArray *oa_range = ostr_range->String().Tokenize(" ");
    ostr_range = (TObjString *)oa_range->At(2);
    TString str = ostr_range->String();
    str = str.Strip(TString::kLeading, '[');
    str = str.Strip(TString::kTrailing, ']');
    oa_range = str.Tokenize(",");
    Int_t runMin = ((TObjString *)oa_range->At(0))->String().Atoi();
    Int_t runMax = ((TObjString *)oa_range->At(1))->String().Atoi();
    AliCDBRunRange runRange(runMin, runMax);
    
    // get version
    TObjString *ostr_version = (TObjString *)oa->At(2);
    TObjArray *oa_version = ostr_version->String().Tokenize(" ");
    ostr_version = (TObjString *)oa_version->At(1);
    TString str = ostr_version->String();
    oa_version = str.Tokenize("_");
    ostr_version = (TObjString *)oa_version->At(0);
    str = ostr_version->String();
    str = str.Strip(TString::kLeading, 'v');
    Int_t version = str.Atoi();
    
    // set specific storage
    printf("--->>> Setting specific storage for %s (version = %d)\n", path.Data(), version);
    TPair *pair = (TPair *)cdbm->FindObject(path.Data());
    if (pair) {// && !path.EqualTo("GRP/GRP/Data")) {
      TString str = ((TObjString *)pair->Value())->String();
      TObjArray *oa = str.Tokenize("?");
      str = ((TObjString *)oa->At(0))->String() + ((TObjString *)oa->At(2))->String();
      printf("--->>> %s\n", str.Data());
      cdb->SetSpecificStorage(path.Data(), str.Data(), version);
    }
    //    else if (path.EqualTo("GRP/GRP/Data")) { // R+HACK
    //      TString str = "alien://DBFolder=/alice/data/2010/OCDB";
    //      printf("--->>> %s\n", str.Data());
    //      cdb->SetSpecificStorage(path.Data(), str.Data(), version);
    //    }
    else {
      pair = (TPair *)cdbm->FindObject("default");
      TString str = ((TObjString *)pair->Value())->String();
      TObjArray *oa = str.Tokenize("?");
      str = ((TObjString *)oa->At(0))->String() + ((TObjString *)oa->At(2))->String();
      printf("--->>> %s\n", str.Data());
      cdb->SetSpecificStorage(path.Data(), str.Data(), version);
    }
    // get object
    cdb->Get(path.Data());
  }
  // dump to snapshot file
  cdb->DumpToSnapshotFile(snapname, kFALSE);
}
