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
  //TFile *f = new TFile
  //TFile*f = new TFile
  ROOT::RDataFrame df_mcp("mCP",mcp_file.Data());
  ROOT::RDataFrame df_metadata("Metadata",mcp_file.Data());
  auto meson  = df_metadata.Take<TString>("Mother").GetValue()[0];
  auto horn   = df_metadata.Take<TString>("HornMode").GetValue()[0];
  auto charge_mcp = df_metadata.Take<Double_t>("mCPcharge").GetValue()[0];
  auto mass_mcp   = df_metadata.Take<Double_t>("mCPmass").GetValue()[0];

  Double_t recoil_threshold = 1000; // detection threshold in keV
  TString electron_output_filename
    = Form("sim/e_%s_q_%.3f_m_%.3f_%s_%ss.root",
	   detector.Data(),charge_mcp,mass_mcp,horn.Data(),meson.Data());

  
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
  
  

  
  // electron recoil lambda
  int n_hits = 5;
  int event_counter = 0;
  
  //ROOT::VecOps::RVec<TLorentzVector> electron_mom;
  std::vector<ROOT::Math::PxPyPzEVector> electron_mom_vector;
  std::vector<ROOT::Math::PxPyPzEVector> electron_pos_vector;
  std::vector<Double_t> weight;
  auto electron_recoil = [ mass_mcp,
			   charge_mcp,
			   n_hits,
			   rot,
			   beampos,
			   &event_counter,
			   &electron_mom_vector,
			   &weight]( TLorentzVector mom )
			 {
			   event_counter++;
			   TVector3 mom3v = mom.Vect();
			   mom3v = rot*mom3v;
			   mom.SetVect(mom3v);
			   
			   electron_mom_vector.clear();
			   weight.clear();
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

			     // the position needs time, use beta from the momentum 4vector
			     electron_mom_vector.push_back(adj_e_recoil);
			   }
			   return electron_mom_vector;
			 };

  auto electron_spawn = [ n_hits,
			  rot,
			  beampos,
			  det_half_dims,
			  det_centre,
			  &event_counter,
			  &electron_pos_vector,
			  &weight]( TLorentzVector pos, TLorentzVector mom )
			{
			  electron_pos_vector.clear();
			  
			  TVector3 pos3v = pos.Vect();
			  TVector3 mom3v = mom.Vect();
			  pos3v = (rot*pos3v)+beampos;
			  mom3v = rot*mom3v;
			  auto entry_exit_points = intersects(pos3v,mom3v,det_centre,det_half_dims);
			  // not sure why it is [0]-[1], I thought it would be [1]-[0]
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
			    travel = entry_exit_points[1] - entry_exit_points[0];
			  } else if ( exit_z == entry_exit_points[0].Z() ) {
			    entry_index = 1;
			    entry_pos = entry_exit_points[1];
			    exit_pos = entry_exit_points[0];
			    travel = entry_exit_points[0] - entry_exit_points[1];
			  }
			  
			  // TIMING
			  // arrival time to detector in seconds
			  Double_t t_in = 0;
			  t_in = pos.T() + (entry_pos-pos3v).Mag()/(mom.Beta()*2.99792e10);
			  
			  TRandom randy;
			  //std::cout << "greatest hits: " << std::endl;
			  for ( int i = 0; i < n_hits; i++ ) {
			    Double_t random_z = randy.Uniform(entry_z,exit_z);
			    Double_t lambda = random_z/travel.Z();
			    //std::cout << lambda << std::endl;
			    TVector3 hit_pos = lambda*travel;
			    Double_t hit_time = t_in + (hit_pos-entry_pos).Mag()/(mom.Beta()*2.99792e10);
			    ROOT::Math::XYZTVector hit_4v(hit_pos.X(),hit_pos.Y(),hit_pos.Z(),hit_time);
			    //TLorentzVector hit4v (hit_pos,hit_time);
			    //hit4v.Print();
			    //hit_pos.Print();
			    electron_pos_vector.push_back(hit_4v);
			  }
			  //pos.SetVect(pos3v);
			  //pos = (rot*pos) + beampos;
			  // the position needs time, use beta from the momentum 4vector
			  

			  return electron_pos_vector;
			};
  
  auto dnew = df_mcp.Define("Mom_e",electron_recoil,{"Mom"})
    //.Define("Pos_e",[&electron_pos_vector](){return electron_pos_vector; })
    .Define("Pos_e",electron_spawn,{"Pos","Mom"})
    .Define("Weight_e",[&weight]() {return weight;})
    .Define("E",[&electron_mom_vector](){
		  std::vector<Double_t> E; 
		  for ( int i = 0; i < electron_mom_vector.size(); i++ ) {
		    E.push_back(electron_mom_vector[i].E());
		  }
		  return E;
		})
    ;

  //auto h1 = dnew.Histo1D("E[0]");

  auto outFileName = "test.root";
  //dnew.Snapshot<std::vector< ROOT::Math::PxPyPzEVector>>("ee",outFileName,{"Mom_e","Pos_e"});
  dnew.Snapshot("ee",outFileName,{"Mom_e","Pos_e","Weight_e","E"});
  //auto fileName = "df002_dataModel.root";
  //auto treeName = "myTree";
  //fill_tree(fileName, treeName);
  
}
