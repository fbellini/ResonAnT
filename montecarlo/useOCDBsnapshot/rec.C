void rec()
{

  AliReconstruction reco;

  reco.SetRunReconstruction("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");

  //
  // switch off cleanESD, write ESDfriends and Alignment data
  
  reco.SetCleanESD(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  reco.SetFractionFriends(.1);

  //
  // ITS Efficiency and tracking errors

  reco.SetRunPlaneEff(kTRUE);
  reco.SetUseTrackingErrorsForAlignment("ITS");

  reco.SetCDBSnapshotMode("ocdb.rec.root");

  reco.SetRunQA(":");

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
