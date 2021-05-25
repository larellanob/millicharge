#include "Root/DetectorInteraction.cxx"
#include "Root/CreateEmptyDirs.cxx"
#include "Root/UbooneAcceptanceChecker.cxx"

void GenerateElectrons
(
 TString mcp_file = "sim/mCP_uboone_q_0.010_m_0.010_fhc_pi0s.root",
 TString detector = "uboone"
)
{
  CreateEmptyDirs();

  // INPUT
  ROOT::RDataFrame df_mcp("mCP",mcp_file.Data());
  ROOT::RDataFrame df_metadata("Metadata",mcp_file.Data());
  auto meson  = df_metadata.Take<TString>("Mother").GetValue()[0];
  auto horn   = df_metadata.Take<TString>("HornMode").GetValue()[0];
  auto charge_mcp = df_metadata.Take<Double_t>("mCPcharge").GetValue()[0];
  auto mass_mcp   = df_metadata.Take<Double_t>("mCPmass").GetValue()[0];

  // PARAMETERS
  int n_hits = 5;
  Double_t recoil_threshold = 1000; // detection threshold in keV
  TString electron_output_filename
    = Form("sim/e_%s_q_%.3f_m_%.3f_%s_%ss_%ihits",
	   detector.Data(),charge_mcp,mass_mcp,horn.Data(),meson.Data(),n_hits);

  // OUTPUT
  auto outFileName = electron_output_filename;
  ofstream out_hepevt;
  out_hepevt.open(Form("%s.txt",outFileName.Data()));
  
  // detector half-dimensions in cm
  const TVector3 det_half_dims(0.5*(246.35-10.),0.5*(107.47+105.53),.5*(1026.8-10.1));
  
  // ~= {0.,0.,46536.3525}; // detector centre in det coords
  const TVector3 det_centre = TVector3(10.,-105.53,10.1) + det_half_dims;
  
  // numi to uboone frame
   /* beam rotation, from /uboone/app/users/zarko/windowRotation/correctWindow.c */
  const TVector3 beampos(-31387.58422,
			 -3316.402543,
			 -60100.2414); // proton on target in uboone coords * cm
  
  
  // defines TRotation c++ 'lambda' named rot
  const TRotation rot
    = []()
      {
	TRotation R;
	// Rotation matrix using the 0,0,0 position for MicroBooNE (beam to det input)
	const TVector3 Rot_row_x = {  0.92103853804025681562,
				      0.022713504803924120662,
				      0.38880857519374290021  };
	const TVector3 Rot_row_y = {  4.6254001262154668408e-05,
				      0.99829162468141474651,
				      -0.058427989452906302359 };
	const TVector3 Rot_row_z = { -0.38947144863934973769,
				     0.053832413938664107345,
				     0.91946400794392302291  };
	
	R.RotateAxes(Rot_row_x, Rot_row_y, Rot_row_z); // Also inverts so now to det to beam
	R.Invert(); // go back to beam to det
	return R;
      }();
  
  

  ///////////////////////////
  // electron recoil
  int event_counter = 0;
    
  std::vector<ROOT::Math::PxPyPzEVector> electron_mom_vector;
  std::vector<ROOT::Math::PxPyPzEVector> electron_pos_vector;
  std::vector<Double_t> weight;
  std::vector<TString> hepevt_line;


  /// electron recoil energy-momemtum
  auto electron_recoil = [ mass_mcp,
			   charge_mcp,
			   n_hits,
			   rot,
			   beampos,
			   &event_counter,
			   &out_hepevt,
			   &electron_mom_vector,
			   &weight,
			   &hepevt_line
			   ]( TLorentzVector mom )
			 {
			   event_counter++;
			   out_hepevt << Form("%i",event_counter) << " " << n_hits << "\n";
			   TVector3 mom3v = mom.Vect();
			   mom3v = rot*mom3v;
			   mom.SetVect(mom3v);
			   electron_mom_vector.clear();
			   weight.clear();
			   hepevt_line.clear();
			   Double_t e_mass = 0.000511;
			   TLorentzVector e_rest(0.0,0.0,0.0,e_mass);
			   TLorentzVector W = mom+e_rest;
			   TGenPhaseSpace recoil;
			   Double_t masses[2] = {mass_mcp,e_mass};
			   recoil.SetDecay(W,2,masses);
			   
			   // NOTE: weight will be 1 in a two-body decay
			   Double_t recoil_weight;
			   // n electrons per each mcp
			   for ( int i = 0; i < n_hits; i++ ) {
			     recoil_weight = recoil.Generate();
			     weight.push_back(recoil_weight);
			     TLorentzVector *mcp_recoil = recoil.GetDecay(0);
			     TLorentzVector *aux = recoil.GetDecay(1);
			     // sample electron recoil energy from cross section
			     Double_t Er,pxr,pyr,pzr;
			     Er = DetectorInteraction(mcp_recoil->T(),mass_mcp,charge_mcp)/1000.;
			     Double_t lambda = aux->T()/Er;
			     Double_t rho = sqrt((lambda*lambda*aux->P()*aux->P())/(aux->E()*aux->E()-lambda*lambda*e_mass*e_mass));
			     
			     Er = aux->T()/lambda;
			     pxr = aux->Px()/rho;
			     pyr = aux->Py()/rho;
			     pzr = aux->Pz()/rho;
			     
			     // adjusted electron
			     ROOT::Math::PxPyPzEVector adj_e_recoil(pxr,pyr,pzr,Er);
			     ROOT::Math::PxPyPzEVector e_recoil(aux->X(),aux->Y(),aux->Z(),aux->T());
			     hepevt_line.push_back(Form("1 11 0 0 0 0 %f %f %f %f %f ",pxr,pyr,pzr,Er,e_mass));
			     // the position needs time, use beta from the momentum 4vector
			     electron_mom_vector.push_back(adj_e_recoil);
			   }
			   return electron_mom_vector;
			 };
  
  /// electron recoil time-position
  auto electron_spawn = [ n_hits,
			  rot,
			  beampos,
			  det_half_dims,
			  det_centre,
			  &out_hepevt,
			  &event_counter,
			  &electron_pos_vector,
			  &weight,
			  &hepevt_line
			  ]( TLorentzVector pos, TLorentzVector mom )
			{
			  electron_pos_vector.clear();
			  TVector3 pos3v = pos.Vect();
			  TVector3 mom3v = mom.Vect();
			  pos3v = (rot*pos3v)+beampos;
			  mom3v = rot*mom3v;
			  auto entry_exit_points = intersects(pos3v,mom3v,det_centre,det_half_dims);
			  TVector3 travel;
			  Int_t entry_index;
			  // will use z coordinate to determine which are entry and exit points through uboone
			  Double_t entry_z = min(entry_exit_points[0].Z(),entry_exit_points[1].Z());
			  Double_t exit_z = max(entry_exit_points[0].Z(),entry_exit_points[1].Z());
			  TVector3 entry_pos;
			  TVector3 exit_pos;
			  if ( entry_z == entry_exit_points[0].Z() ) {
			    entry_index = 0;
			    entry_pos = entry_exit_points[0];
			    exit_pos = entry_exit_points[1];
			  } else if ( exit_z == entry_exit_points[0].Z() ) {
			    entry_index = 1;
			    entry_pos = entry_exit_points[1];
			    exit_pos = entry_exit_points[0];
			  }
			  travel = exit_pos - entry_pos;
			  // TIMING
			  // arrival time to detector in nanoseconds
			  Double_t t_in = 0;
			  t_in = pos.T() + ((entry_pos-pos3v).Mag()/(mom.Beta()*2.99792e10))*1e9;
			  // this is giving us the 'same' hit distribution every time (uses same random seed)
			  TRandom randy;
			  for ( int i = 0; i < n_hits; i++ ) {
			    Double_t random_z = randy.Uniform(0,exit_z-entry_z);
			    Double_t lambda = random_z/travel.Z();
			    // recoil electron position
			    TVector3 hit_pos = entry_pos + lambda*travel;
			    // recoil electron time
			    Double_t hit_time = t_in + ((lambda*travel).Mag()/(mom.Beta()*2.99792e10))*1e9;
			    hit_time += 5627.5; // time offset (in ns)
			    ROOT::Math::XYZTVector hit_4v(hit_pos.X(),hit_pos.Y(),hit_pos.Z(),hit_time);
			    electron_pos_vector.push_back(hit_4v);
			    // write to hepevt file
			    hepevt_line[i]+=Form("%f %f %f %f\n",hit_pos.X(),hit_pos.Y(),hit_pos.Z(),hit_time);
			  }
			  // write to hepevt file
			  for ( int i = 0; i < n_hits; i++ ) {
			    out_hepevt << hepevt_line[i];
			  }
			  return electron_pos_vector;
			};
  
  auto dnew = df_mcp.Define("Mom_e",electron_recoil,{"Mom"})
    .Define("Pos_e",electron_spawn,{"Pos","Mom"})
    .Define("Weight_e",[&weight]() {return weight;})
    ;


  TString out_root = (TString)outFileName+".root";
  dnew.Snapshot("ee",out_root.Data(),{"Mom_e","Pos_e","Weight_e"})
    ;
  
}
