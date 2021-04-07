// Reads mCP file in sim/ folder and returns "geometric acceptance"
// and Number of mCPs through ArgoNeuT

// uses cross section from eq (2) in arXiv:1812.03998v2
// or the PhysRevD associated with it
#include "CrossSection.cxx"

#include "ubooneGeo.cxx"


Double_t AcceptedArgoneut(TString fstr, int WEIGHT = 1, TString detector = "uboone")
{
  // opens files containing mCP particles
  TFile *f  = new TFile(fstr);
  TTreeReader reader("mCP",f);
  TTreeReaderValue<TLorentzVector> Mom(reader,"Mom");
  TTreeReaderValue<TLorentzVector> Pos(reader,"Pos");
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


  Double_t xdev;
  Double_t ydev;

  if ( detector == "dune" ) {
    // DUNE geometry
    // fermini group: Phys. Rev. D 100, 015043
    xdev = 0.5/574.;
    ydev = 0.5/574.;
  }

  if ( detector == "argoneut" ) {
    // argoneut group: arXiv:1902.03246
    xdev = 0.235/975.;
    ydev = 0.2/975.;
  }

  if ( detector == "duneOrnella" ) {
    // DUNE geometry
    // argoneut group: arXiv:1902.03246
    xdev = 1.5/574.;
    ydev = 2.0/574.;
  }

  // uboone position wrt numi target (in cm)
  const TVector3 uboonepos(31387.58422,
			   3316.402543,
			   60100.2414); // proton on target in uboone coords * cm
  Double_t theta;
  if ( detector == "naiveuboone" ) {
    // naive approximation to uboone geometry
    // circle of radius 5m a distance 679m from target
    theta = atan2(5,679);
  }

  
  // event loop
  while ( reader.Next() ) {
    events++;

    // transverse momentum and angle
    Double_t Pt = sqrt(Mom->Px()*Mom->Px() + Mom->Py()*Mom->Py());
    Double_t Th = atan2(Pt,Mom->Pz());
    
    /*
    if ( events < 5 ) {
      //cout  << Pt << endl;
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


    if ( ( detector ==  "dune" || detector == "argoneut" || detector == "duneOrnella" ) &&
	 //abs(atan2(Mom->Px(),Mom->Pz())) < xdev &&
	 //abs(atan2(Mom->Py(),Mom->Pz())) < ydev ) {
	 abs(Mom->Px()/Mom->Pz()) < xdev &&
	 abs(Mom->Py()/Mom->Pz()) < ydev ) {
      value4+= (*weight_meson)*(*weight_decay)*(*diff_xsec);
      value3+= (*weight_meson)*(*weight_decay);
      value2+= (*weight_meson);
      value1++;
    }

    if ( detector ==  "naiveuboone" &&
	 Mom->Angle(uboonepos) <  theta ) {
      value4+= (*weight_meson)*(*weight_decay)*(*diff_xsec);
      value3+= (*weight_meson)*(*weight_decay);
      value2+= (*weight_meson);
      value1++;
    }

    // uses ubooneGeo.cxx to check uboone mCP hit
    if ( detector == "uboone" && ubooneGeo(Pos->Vect(),Mom->Vect()) == true) {
      value4+= (*weight_meson)*(*weight_decay)*(*diff_xsec);
      value3+= (*weight_meson)*(*weight_decay);
      value2+= (*weight_meson);
      value1++;
    }
    
  }
  
  TTreeReader meta("Metadata",f);
  TTreeReaderValue<Double_t> mass(meta,"mCPmass");
  TTreeReaderValue<Double_t> charge(meta,"mCPcharge");
  // the generated mCP files have one cross section but that's not
  // necesairly what needs to be used
  TTreeReaderValue<Double_t> xsec(meta,"CrossSection"); 
  TTreeReaderValue<TString> meson(meta,"Mother");
  meta.Next();


  cout << "============ AcceptedArgoneut.cxx output: ============" << endl;
  cout << "******************************************************" << endl;
  cout << Form("Decays of %s into mCPs of mass %.3f and charge %.3f",
	       meson->Data(),*mass,*charge) << endl;
  cout << "         Using "+detector+" geometry and POTs" << endl;
  cout << "******************************************************" << endl;
  cout << Form("finished looping %i events",events) << endl;



  // cross section correction
  TDatabasePDG pdg;
  double mesonmass = 0;
  if ( *meson == "pi0" ) {
    mesonmass = pdg.GetParticle(111)->Mass();
  } else if ( *meson == "eta" ) {
    mesonmass = pdg.GetParticle(221)->Mass();
  }
  Double_t decay_factor_zhenliu = CrossSection(*mass,mesonmass,"zhenliu");
  Double_t decay_factor_physrevd = CrossSection(*mass,mesonmass,"physrevd");
  Double_t decay_factor_naive = CrossSection(*mass,mesonmass,"naive");
  Double_t decay_factor_t2k = CrossSection(*mass,mesonmass,"t2k");
  Double_t truexsec = *xsec; // no need to recalculate

  cout << "Simulated particles entering detector geometry: " << value1 << endl;
  cout << "TTree cross section: " << *xsec << endl;
  cout << "Recalculated cross section: " << truexsec << endl;
  cout << "Unweighted mCPs passing: " << value1 << endl;
  
  /*
  cout << Form("value: %.3f\nPOT normalized value: %.3f\n(old) xsec: %.7f\nsum weight decays: %.3f\nsum weight meson: %.3f",
	       value1,(1e20/500000.)*value1/(sum_weight_decay*0.5),*xsec,sum_weight_decay/2.,sum_weight_meson/2.) << endl;
  */
  Double_t result;

  // POT normalization
  Double_t POT_norm;
  if ( detector == "argoneut" ) {
    POT_norm = 1e20/500000.;
  } else if ( detector == "dune" ) {
    POT_norm = 1e21/500000.;
  } else if ( detector == "duneOrnella" ) {
    POT_norm = 3*1e22/500000.;
  } else if ( detector == "uboone" ) {
    POT_norm = 1e21/500000.;
  }

  
  if ( WEIGHT == 1 ) {
    result = POT_norm*(truexsec)*value1; // result uses truexsec
  } else if ( WEIGHT == 2 ) {
    result = POT_norm*(truexsec)*value2; // result uses truexsec
  } else if ( WEIGHT == 3 ) {
    result = POT_norm*(truexsec)*value3*(events/sum_weight_decay); // result uses truexsec
  } else if ( WEIGHT == 4 ) {
    result = POT_norm*value4*(events/sum_weight_decay); // differential cross section
  } else if ( WEIGHT == 5 ) {
    //result = POT_norm*value3*(events/sum_weight_decay)*decay_factor_zhenliu;
    result = POT_norm*value3*(events/sum_weight_decay)*decay_factor_physrevd;
    cout << "USING DECAY FACTOR decay_factor_physrevd" << endl;
  } else if ( WEIGHT == 6 ) {
    result = decay_factor_t2k;
  }
  //cout << "geometrical acceptance " << value1/(events) << endl;
  cout << "geometrical acceptance " << value3/(sum_weight_decay) << endl;
  cout << "result " << result << endl;
  cout << "----------------------------------" << endl;
  return result;
  
}
