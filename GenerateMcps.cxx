#include "Root/CrossSection.cxx"
#include "Root/DiffCrossSection.cxx"
#include "Root/CreateEmptyDirs.cxx"

const bool kTest = false;

void GenerateMcps(TString meson ="eta",TString horn = "rhc", Double_t mCPmass = 0.01, Double_t mCPcharge = 0.01)
{
  CreateEmptyDirs();
  bool kWriteFile = true;
  if ( kTest ) {
  kWriteFile = false;
 }
  // decay parameters
  //Double_t mCPmass   = 0.01;   // GeV
  //Double_t mCPcharge = 0.01;   // electron charges

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
  
  std::cout << "To generate mCP decaying from "
	    << Form("%ss in %s horn mode",meson.Data(),horn.Data())
	    <<std::endl;

  /////////
  // input
  TString mode = horn+"_"+meson+"s";
  TFile *f = new TFile("sim/"+mode+"_tree.root");
  TTreeReader reader(mode,f);
  TTreeReaderValue<TLorentzVector> Mom4v(reader,"Mom");
  TTreeReaderValue<TLorentzVector> Pos4v(reader,"Pos");
  TTreeReaderValue<Float_t> weight_pi0(reader,"Weight");
  TTreeReaderValue<Int_t> evnt(reader,"Event");
  Long64_t TotalEntries = reader.GetEntries();

  //////////
  // decay
  // cross section dedicated macro
  Double_t cross_section = CrossSection(mCPmass,meson);
  if ( 2*mCPmass > kMesonMass ) {
    std::cout << "ERROR: Kinematically forbidden decay. Aborting." << std::endl;
    return;
  }

  //////////
  // output

  // output file
  TString out_filename
    = Form("sim/mCP_q_%0.3f_m_%0.3f_%s.root",mCPcharge,mCPmass,mode.Data());
  TFile g(out_filename,"recreate");
  if ( kWriteFile == false ) {
    std::cout << "NO OUTPUT FILE IS BEING GENERATED" << std::endl;
    g.Close();
  }
  //tree->SetDirectory(g); // remember this line in the future!

  // metadata tree
  TTree meta("Metadata","variables used");
  meta.Branch("mCPmass",&mCPmass);
  meta.Branch("mCPcharge",&mCPcharge);
  meta.Branch("CrossSection",&cross_section);
  meta.Branch("Mother",&meson);
  meta.Branch("HornMode",&horn);
  meta.Fill();
  if ( kWriteFile ) {
    meta.Write();
  } else {
    std::cout << "NO OUTPUT FILE IS BEING GENERATED" << std::endl;
  }

  // main mcp particle tree
  TLorentzVector mCPMom;
  TLorentzVector mCPPos;
  TTree * tree = new TTree("mCP",Form("Millicharged particles with mass %0.3f, charge %0.3f",mCPmass, mCPcharge));
  tree->Branch("Mom",&mCPMom); // particle tlorentz4v
  tree->Branch("Pos",&mCPPos); // particle tlorentz4v
  Float_t weight_meson;      // weight from the files
  Float_t weight_decay;    // weight from tgenphasespace
  tree->Branch("WeightMeson",&weight_meson);
  tree->Branch("WeightDecay",&weight_decay);
  Int_t event;
  tree->Branch("Event",&event);
  Double_t cross_section_differential;
  tree->Branch("DiffCrossSection",&cross_section_differential);
  

  
  //////////////
  // events loop

  std::cout << "Generating mCPs for "
	    << Form("mass %0.3f GeV, charge %0.3f epsilon",mCPmass,mCPcharge)
	    << std::endl;
  
  int events = 0;
  while ( reader.Next() ) {
    events++;

    // print progress at regular intervals
    if ( events % 100000 == 0 ) {
      std::cout << Form("Completed %04.1f%% of %lli events",((Double_t)events/TotalEntries)*100,TotalEntries) << "\r" << std::flush;
    }

    // do only 5 events if testing
    if ( kTest && events < 5 ) {
      std::cout << "Meson Px: " << Mom4v->Px() << std::endl;
    } else if ( kTest && events >= 5 ) {
      break;
    }    


    // generate decay
    Double_t masses[3] = { 0, mCPmass, mCPmass};
    TGenPhaseSpace decay;
    decay.SetDecay(*Mom4v, 3, masses);
    weight_decay = decay.Generate();
    Double_t wtmax = decay.GetWtMax();
    
    //cout << Form("weight decay: %.3f, weightmax: %.3f",weight_decay,wtmax) << endl;
    
    // collect decay products
    TLorentzVector *gamma = decay.GetDecay(0);
    TLorentzVector *mCP1 = decay.GetDecay(1);
    TLorentzVector *mCP2 = decay.GetDecay(2);


    cross_section_differential =
      DiffCrossSection(*mCP1,*mCP2,mCPcharge,meson,kMesonMass);

    //cout << "diferential cross section " << cross_section_differential << endl;

    // fill metadata variables
    event = *evnt;
    weight_meson = *weight_pi0;

    // fills two times, so two events for each decay
    // but branch "event" keeps track of the event
    mCPMom = *mCP1;
    mCPPos = *Pos4v; // both mCPs are generated in the same position,
		     // so fill position once
    tree->Fill();
    mCPMom = *mCP2;
    tree->Fill(); // fill again for the second mCP particle
  }

  // finish message
  std::cout << "Finished looping. Generated: " << events << " pairs" << std::endl;
  // write output
  if ( kWriteFile ) {
    g.Write();
    std::cout << "Wrote mCPs to file: " << out_filename <<  std::endl;
  } else {
    std::cout << "NO OUTPUT FILE IS BEING GENERATED" << std::endl;
  }
}
