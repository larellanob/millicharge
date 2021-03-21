// Reads mCP file in sim/ folder and returns "geometric acceptance"
// and Number of mCPs through ArgoNeuT

// uses cross section from eq (2) in arXiv:1812.03998v2
// or the PhysRevD associated with it
#include "CrossSection.cxx"

Double_t AcceptedArgoneut(TString fstr, int WEIGHT = 1, bool DUNE = false)
{
  // opens files containing mCP particles
  TFile *f  = new TFile(fstr);
  TTreeReader reader("mCP",f);
  TTreeReaderValue<TLorentzVector> Mom(reader,"Mom");
  TTreeReaderValue<Float_t> weight_decay(reader,"WeightDecay");
  TTreeReaderValue<Float_t> weight_meson(reader,"WeightMeson");
  TTreeReaderValue<Int_t> event(reader,"Event");
  TTreeReaderValue<Double_t> diff_xsec(reader,"DiffCrossSection");
  // loops over events
  int events = 0;
  Double_t value1 = 0; // result stored here
  Double_t value2 = 0; // result stored here
  Double_t value3 = 0; // result stored here
  Double_t value4 = 0; // result stored here
  Double_t sum_weight_decay = 0;
  Double_t sum_weight_meson = 0;
  Double_t sum_diff_xsec = 0;
  while ( reader.Next() ) {
    events++;

    // transverse momentum and angle
    Double_t Pt = sqrt(Mom->Px()*Mom->Px() + Mom->Py()*Mom->Py());
    Double_t Th = atan2(Pt,Mom->Pz());
    
    /*
    if ( events < 5 ) {
      cout  << Pt << endl;
    } else break;
    */

    // if they pass, compute integral
    sum_weight_decay += *weight_decay;
    sum_weight_meson += *weight_meson;
    sum_diff_xsec += *diff_xsec;
    // "angular acceptance"
    // this is a 'circular detector' approximation used on Friday 12 March 2021
    /*
    if ( Mom->Theta() <  0.00019344230 ) {
      //value2+= (*weight_meson);
      //value3+= (*weight_meson)*(*weight_decay);
      value1++;
    }
    */

    // "cartesian acceptance"
    // closer 'rectangular detector' approximation
    // they use 975m as the distance to the target in the paper
    Double_t xdev = atan2(0.235,975.); // x deviation
    Double_t ydev = atan2(0.2  ,975.); // y deviation

    if ( DUNE ) {
      // DUNE geometry instead of ArgoNeuT
      xdev = atan2(0.5,1040.);
      ydev = atan2(0.5,1040.);
    }
    
    /*
    Double_t xdev = atan2(0.235,1033.); // x deviation
    Double_t ydev = atan2(0.2  ,1033.); // y deviation
    */
    if ( abs(atan2(Mom->Px(),Mom->Pz())) < xdev && abs(atan2(Mom->Py(),Mom->Pz())) < ydev ) {
      value4+= (*weight_meson)*(*weight_decay)*(*diff_xsec);
      value3+= (*weight_meson)*(*weight_decay);
      value2+= (*weight_meson);
      value1++;
    }

    /*
    cout << "pt over pz " << atan2(Pt,Mom->Pz()) << endl;
    cout << "geo " << Mom->Theta() << endl;
    cout << "px over pz " << atan2(Mom->Px(),Mom->Pz()) << endl;
    cout << "py over pz " << atan2(Mom->Py(),Mom->Pz()) << endl;
    */
  }
  cout << "----------------------------------" << endl;
  cout << "== AcceptedArgoneut.cxx output: ==" << endl;
  if ( DUNE ) {
    cout << "**********************************" << endl;
    cout << "***Using DUNE geometry and POTs***" << endl;
    cout << "**********************************" << endl;
  }
  cout << Form("finished looping %i events",events) << endl;
  
  TTreeReader meta("Metadata",f);
  TTreeReaderValue<Double_t> mass(meta,"mCPmass");
  TTreeReaderValue<Double_t> charge(meta,"mCPcharge");
  // the generated mCP files have one cross section but that's not
  // necesairly what needs to be used
  TTreeReaderValue<Double_t> xsec(meta,"CrossSection"); 
  TTreeReaderValue<TString> meson(meta,"Mother");
  meta.Next();

  // cross section correction
  TDatabasePDG pdg;
  double mesonmass = 0;
  if ( *meson == "pi0" ) {
    mesonmass = pdg.GetParticle(111)->Mass();
  } else if ( *meson == "eta" ) {
    mesonmass = pdg.GetParticle(221)->Mass();
  }
  //Double_t truexsec = CrossSection(*mass,mesonmass); // recalculate xsec using CrossSection
  Double_t truexsec = *xsec; // no need to recalculate

  cout << "TTree cross section: " << *xsec << endl;
  cout << "Recalculated cross section: " << truexsec << endl;
  cout << "Unweighted mCPs passing: " << value1 << endl;
  
  /*
  cout << Form("value: %.3f\nPOT normalized value: %.3f\n(old) xsec: %.7f\nsum weight decays: %.3f\nsum weight meson: %.3f",
	       value1,(1e20/500000.)*value1/(sum_weight_decay*0.5),*xsec,sum_weight_decay/2.,sum_weight_meson/2.) << endl;
  */
  Double_t result;
  Double_t POT_norm = 1e20/500000.; // POT normalization
  if ( DUNE ) {
    POT_norm = 1e21/500000.;
  }
  if ( WEIGHT == 1 ) {
    result = POT_norm*(truexsec)*value1; // result uses truexsec
  } else if ( WEIGHT == 2 ) {
    result = POT_norm*(truexsec)*value2; // result uses truexsec
  } else if ( WEIGHT == 3 ) {
    result = POT_norm*(truexsec)*value3*(events/sum_weight_decay); // result uses truexsec
  } else if ( WEIGHT == 4 ) {
    result = POT_norm*value4*(events/sum_weight_decay);
  }
  cout << "geometrical acceptance " << value1/(events) << endl;
  cout << "result " << result << endl;
  cout << "----------------------------------" << endl;
  return result;
  
}
