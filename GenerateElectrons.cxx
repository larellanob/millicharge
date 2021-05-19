#include "Root/DetectorInteraction.cxx"
#include "Root/CreateEmptyDirs.cxx"

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

  
  TString electron_output_filename
    = Form("sim/e_%s_q_%.3f_m_%.3f_%s_%ss.root",
	   detector.Data(),charge_mcp,mass_mcp,horn.Data(),meson.Data());

  // electron recoil lambda

  //ROOT::VecOps::RVec<TLorentzVector> electron_mom;
  std::vector<ROOT::Math::PxPyPzEVector> electron_mom;
  std::vector<ROOT::Math::PxPyPzEVector> electron_pos;
  std::vector<Double_t> weight;
  auto electron_recoil = [mass_mcp,
			  charge_mcp,
			  &electron_mom,
			  &electron_pos,
			  &weight](TLorentzVector pos,
				   TLorentzVector mom)
			 {
			   electron_mom.clear();
			   electron_pos.clear();
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
			   int n = 3;
			   for ( int i = 0; i < n; i++ ) {
			     recoil_weight = recoil.Generate();
			     weight.push_back(recoil_weight);
			     TLorentzVector *mcp_recoil = recoil.GetDecay(0);
			     TLorentzVector *aux = recoil.GetDecay(1);
			     // sample electron recoil energy from cross section
			     Double_t Er,pxr,pyr,pzr;
			     Er = DetectorInteraction(aux->T(),mass_mcp,charge_mcp)/1000.;
			     Double_t lambda = aux->T()/Er;
			     Double_t rho = sqrt((lambda*lambda*aux->P()*aux->P())/(aux->E()*aux->E()-lambda*lambda*e_mass*e_mass));
			    
			     Er = aux->T()/lambda;
			     pxr = aux->Px()/rho;
			     pyr = aux->Py()/rho;
			     pzr = aux->Pz()/rho;
			     
			     ROOT::Math::PxPyPzEVector adje_recoil(pxr,pyr,pzr,Er);
			     ROOT::Math::PxPyPzEVector e_recoil(aux->X(),aux->Y(),aux->Z(),aux->T());
			     
			     ROOT::Math::PxPyPzEVector e_pos(aux->X()+1.0,aux->Y(),aux->Z(),aux->T());
			     electron_mom.push_back(e_recoil);
			     electron_pos.push_back(e_pos);
			     //return *aux;
			   }
			   return electron_mom;
			 };
  

  auto dnew = df_mcp.Define("Mom_e",electron_recoil,{"Pos","Mom"})
    .Define("Pos_e",[&electron_pos](){return electron_pos; })
    .Define("Weight_e",[&weight]() {return weight;})
    .Define("E",[&electron_mom](){
		  std::vector<Double_t> E; 
		  for ( int i = 0; i < electron_mom.size(); i++ ) {
		    E.push_back(electron_mom[i].E());
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
