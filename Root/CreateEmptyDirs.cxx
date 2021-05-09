void CreateEmptyDirs()
{
  gSystem->Exec("mkdir -p sim");
  gSystem->Exec("mkdir -p hist");
  gSystem->Exec("mkdir -p img/PassingThroughDetector");
  gSystem->Exec("mkdir -p img/Limits");
  gSystem->Exec("mkdir -p img/SupressionFactor");
  
  
}
