#include "Root/UbooneAcceptanceChecker.cxx"
#include "Root/CreateEmptyDirs.cxx"

void FilterAccepted(TString fstr = "sim/mCP_q_0.010_m_0.010_fhc_pi0s.root",
		    int WEIGHT = 5,
		    TString detector = "uboone",
		    Bool_t kTest = false)
{
  CreateEmptyDirs();
  TFile *f  = new TFile(fstr);
  TTree *mcp = f->Get<TTree>("mCP");
  ROOT::RDataFrame df("mCP",fstr.Data());
  ROOT::RDataFrame dfmetadata("Metadata",fstr.Data());

  
  // Detector cuts
  Double_t xdev;
  Double_t ydev;
  if ( detector == "dune" ) {
    // DUNE from fermini group: Phys. Rev. D 100, 015043
    xdev = 0.5/574.;
    ydev = 0.5/574.;
  }

  if ( detector == "argoneut" ) {
    // argoneut group: arXiv:1902.03246
    xdev = 0.235/975.;
    ydev = 0.2/975.;
  }

  if ( detector == "duneArgo" ) {
    // DUNE from argoneut group: arXiv:1902.03246
    xdev = 1.5/574.;
    ydev = 2.0/574.;
  }

  Double_t theta;
  if ( detector == "naiveuboone" ) {
    // naive approximation to uboone geometry
    // circle of radius 5m a distance 679m from target
    theta = atan2(5,679);
  }
  // event loop

  auto cut_geometry = [detector,xdev, ydev] (TLorentzVector pos, TLorentzVector mom)
		      {
			if ( detector == "uboone" ) {
			  return UbooneAcceptanceChecker(pos.Vect(),mom.Vect()) > 0;
			} else if ( detector ==  "argoneut" ||
				    detector == "dune" ||
				    detector == "duneArgo" ||
				    detector == "naiveboone" ) {
			  return abs(mom.Px()/mom.Pz()) < xdev && abs(mom.Py()/mom.Pz()) < ydev;
			} else {
			  std::cout << "WARNING: don't know what to do with this geometry" << std::endl;
			  return false;
			}
		      };

  ROOT::RDF::RNode rangedDF = df;
  
  if ( kTest == true ) {
    rangedDF = df.Range(0,1000);
    std::cout << "INFO: Running in test mode, stopping at 1000 events" << std::endl;
  }


  // metadata values
  // (one entry in each branch of tree contained in dfmetadata)
  auto meson  = dfmetadata.Take<TString>("Mother").GetValue()[0];
  auto horn   = dfmetadata.Take<TString>("HornMode").GetValue()[0];
  auto charge = dfmetadata.Take<Double_t>("mCPcharge").GetValue()[0];
  auto mass   = dfmetadata.Take<Double_t>("mCPmass").GetValue()[0];


  TString filtered_output_filename = Form("sim/mCP_%s_q_%.3f_m_%.3f_%s_%ss.root",
					  detector.Data(),charge,mass,horn.Data(),meson.Data());
  auto df2 = rangedDF.Filter(cut_geometry,{"Pos","Mom"}).Snapshot("mCP",filtered_output_filename.Data());
  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "UPDATE";
  dfmetadata.Snapshot("Metadata",filtered_output_filename.Data(),dfmetadata.GetColumnNames(),opts);

  // number of entries collected
  auto df2_count = df2->Count();
  auto df_count = rangedDF.Count();
  
  // cout accepted event numbers, if any
  std::cout << "List of events passing through detector space: " << std::endl;
  if ( *df2_count > 0 ) {
    auto b1Vec = df2->Take<int>("Event");
    if (auto b1List = df2->Take<int, std::list<int>>("Event") ) {
      for (auto b1_entry : *b1List)
	std::cout << b1_entry << " ";
      std::cout << std::endl;
    }
  }

  //std::cout << "The type of mesonCol is " << mesonColClass->GetName() << std::endl;
  //std::cout << mesonCol[0] << endl;
  std::cout << "============= FilterAccepted.cxx output: =============" << std::endl;
  std::cout << "******************************************************" << std::endl;
  std::cout << Form("Filtered decays of %s into mCPs of mass GeV %.3f, charge %.3f",
		    meson.Data(),mass,charge) << std::endl;
  std::cout << "using "+detector+" geometry and POTs" << std::endl;
  std::cout << "******************************************************" << std::endl;
  //std::cout << Form("finished looping %i events",events) << std::endl;
  std::cout << "Passing events: " << *df2_count << ", Parsed events: " << *df_count << std::endl;
  
}
