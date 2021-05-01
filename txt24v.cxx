void txt24v(TString txtfile) {
  // files originally located in 
  // /afs/hep.man.ac.uk/d/lartpc-RanD/millicharge/meson_flux/;
  
  TString meson, horn;
   if ( txtfile.Contains("pi0") ) {
    meson = "pi0";
  } else if ( txtfile.Contains("eta") ) {
    meson = "eta";
  }
  if ( txtfile.Contains("fhc") ) {
    horn = "fhc";
  } else if ( txtfile.Contains("rhc") ) {
    horn = "rhc";
  }

  if ( meson == "" || horn == "" ) {
    std::cout << "Can't understand file input" << std::endl;
    return;
  }
  
  std::cout
    << "Generating tree from text files for "+meson+"s in horn mode "+horn
    << std::endl;
     
  TString mode = horn+"_"+meson+"s";
  ifstream in;
  in.open(Form("%s",txtfile.Data()));
  
  auto f = TFile::Open("sim/"mode+"_tree.root","RECREATE");
  Float_t Px,Py,Pz,E,x,y,z,t,weight;
  Int_t event;

  TLorentzVector mom(0.0,0.0,0.0,0.0);
  TLorentzVector pos(0.0,0.0,0.0,0.0);
  //TLorentzVector pos(x,y,z,t);
  TTree *tree = new TTree(mode,meson+"s from NuMI in "+mode+" mode");
  tree->Branch("Mom","TLorentzVector",&mom);
  tree->Branch("Pos","TLorentzVector",&pos);
  tree->Branch("Weight",&weight);
  tree->Branch("Event",&event);
  
  int nlines = -1;
  while (1) {
    nlines++;

    /*
    if  (nlines > 9 ) {
      break;
    }
    */
    in >> Px >> Py >> Pz >> E >> x >> y >> z >> t >> weight >> event;
    if (!in.good()) break; // stop reading when file is over
    if ( nlines < 1 ) {
      cout << "first five events" << endl;
    }
    if (nlines < 5) {
      printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
    }
    mom.SetPxPyPzE(Px,Py,Pz,E);
    pos.SetXYZT(x,y,z,t);
    tree->Fill();
    
  }
  cout << "events written: " << nlines << endl;
  in.close();
  f->Write();
  
}
