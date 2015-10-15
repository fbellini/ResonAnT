void sim(Int_t nev=4)
{

  AliSimulation simulator;

  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  simulator.SetMakeDigitsFromHits("ITS TPC");
  
  simulator.SetCDBSnapshotMode("ocdb.sim.root");
  
  //
  // Vertex and Mag.field from OCDB
  simulator.UseVertexFromCDB();
  simulator.UseMagFieldFromGRP();

  simulator.SetRunQA(":");

  //
  // The rest
  
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
