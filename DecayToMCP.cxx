
void DecayToMCP(TString meson ="pi0",TString horn = "fhc")
{
  // decay parameters
  Double_t mCPmass   = 0.01;   // GeV
  Double_t mCPcharge = 0.01;   // electron charges

  TDatabasePDG pdg;
  double kMesonMass;
  double kLifetime;
  double c = 3e8;
  if ( meson == "pi0" ) {
    kMesonMass = pdg.GetParticle(111)->Mass();
    kLifetime = 8.4e-17;
  }
  if ( meson == "eta" ) {
    kMesonMass = pdg.GetParticle(221)->Mass();
    kLifetime = 5.0e-19;
  }

  /////////
  // input
  TString mode = horn+"_"+meson+"s";
  TFile *f = new TFile(mode+"_tree.root");
  TTreeReader reader(mode,f);
  TTreeReaderValue<TLorentzVector> Mom4v(reader,"Mom");
  TTreeReaderValue<TLorentzVector> Pos4v(reader,"Pos");
  TTreeReaderValue<Float_t> weight_pi0(reader,"Weight");
  TTreeReaderValue<Int_t> evnt(reader,"Event");


  //////////
  // decay
  Double_t fmM = sqrt(1.-(mCPmass*mCPmass)/(kMesonMass*kMesonMass));
  Double_t branching_ratio = mCPcharge*mCPcharge * 0.01174 * fmM;
  cout << "Branching ratio " << branching_ratio << endl;
  if ( 2*mCPmass > kMesonMass ) {
    cout << "Forbidden decay, breaking" << endl;
    return;
  }

  //////////
  // output
  //TString out_filename = "mCP_q_"+mCPcharge+"_m_"+mCPmass+"_"+mode+".root";
  TString out_filename
    = Form("mCP_q_%0.3f_m_%0.3f_%s.root",mCPcharge,mCPmass,mode.Data());
  auto g = TFile::Open(out_filename,"recreate");
  //tree->SetDirectory(g); // remember this line!

  TLorentzVector mCPMom;
  TLorentzVector mCPPos;
  TTree * tree = new TTree("mCP",Form("Millicharged particles with mass %0.3f, charge %0.3f",mCPmass, mCPcharge));
  tree->Branch("Mom",&mCPMom); // particle tlorentz4v
  tree->Branch("Pos",&mCPPos); // particle tlorentz4v
  Float_t weight_final;    // weight from tgenphasespace
  tree->Branch("Weight",&weight_final);
  Int_t event;
  tree->Branch("Event",&event);


  
  //////////////
  // events loop
  int events = 0;
  while ( reader.Next() ) {
    events++;
    // for testing
    /*
    if ( events < 300 ) {
      cout << "Meson Px: " << Mom4v->Px() << endl;
    } else break;
    */
    Double_t masses[3] = { 0, mCPmass, mCPmass};
    TGenPhaseSpace decay;
    decay.SetDecay(*Mom4v, 3, masses);
    Double_t weight_decay = decay.Generate();

    TLorentzVector *gamma = decay.GetDecay(0);
    TLorentzVector *mCP1 = decay.GetDecay(1);
    TLorentzVector *mCP2 = decay.GetDecay(2);
    
    event = *evnt;
    weight_final = (*weight_pi0)*(weight_decay)*(branching_ratio);

    // fills two times, so two events, but variable event has event
    // number
    mCPMom = *mCP1;
    mCPPos = *Pos4v;
    tree->Fill();
    mCPMom = *mCP2;
    tree->Fill();

  }

  cout << "done " << events << " events" << endl;
  // write output
  g->Write();
}
