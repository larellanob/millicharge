Double_t FormFactor(Double_t x, Double_t mass_mes ) {
  Double_t api0 = 0.11; // dimensionless
  Double_t leta = 0.72; // GeV
  if ( mass_mes < 0.3 ) { // pi0
    return pow((1.0+api0*x/(mass_mes*mass_mes)),2);
  } else { // eta
    return pow((1.0-x/(leta*leta)),-2);
  }
}

Double_t I3ZhenLiu(Double_t *t, Double_t *par ) {
  double q2 = t[0]; // integration variable
  double Ml = par[0]; // mCP mass (lepton)
  double MP = par[1]; // meson mass (Pseudoscalar meson)
  double Ml2 = Ml*Ml;
  double MP2 = MP*MP;
  return (1.0/q2) *  sqrt(1.0-((4.0*Ml2)/q2)) * (1.0+2.0*(Ml2/q2))
    * pow((1.0-(q2/MP2)),3);
}

Double_t I3ZhenLiuFormFactors(Double_t *t, Double_t *par ) {
  double q2 = t[0]; // integration variable
  double Ml = par[0]; // mCP mass (lepton)
  double MP = par[1]; // meson mass (Pseudoscalar meson)
  double Ml2 = Ml*Ml;
  double MP2 = MP*MP;
  return (1.0/q2) *  sqrt(1.0-((4.0*Ml2)/q2)) * (1.0+2.0*(Ml2/q2))
    * pow((1.0-(q2/MP2)),3) * FormFactor(q2,MP);
}


Double_t DecayFactor(Double_t mass_mcp, TString meson = "", TString mode = "zhenliu" )
{
  // this might generate a warning of 'object already instantiated
  if ( meson == "" ) {
    std::cout << "ERROR: no meson for CrossSection.cxx" << std::endl;
    return 0;
  }
  TDatabasePDG pdgCrossSection;
  Double_t mass_mes = pdgCrossSection.GetParticle(meson.Data())->Mass();
  Double_t x = (mass_mcp*mass_mcp)/(mass_mes*mass_mes);
  if ( x == 0 ) {
    cout << "ERROR: mCP mass 0?" << endl;
    //return 0;
  }
  Double_t alpha = 1./137.;

  // zhen liu mode (ornella paper)
  Double_t mass_ele = pdgCrossSection.GetParticle(11)->Mass();
  TF1 *fZhenMCP = new TF1("I3ZhenMCP",I3ZhenLiu,4.0*mass_mcp*mass_mcp,mass_mes*mass_mes,2);
  TF1 *fZhenEle = new TF1("I3ZhenEle",I3ZhenLiu,4.0*mass_ele*mass_ele,mass_mes*mass_mes,2);
  fZhenMCP->SetParameters(mass_mcp,mass_mes);
  fZhenEle->SetParameters(mass_ele,mass_mes);

  //Double_t iZhenMCP = fZhenMCP->Integral(2.0*mass_mcp,mass_mes);
  //Double_t iZhenEle = fZhenEle->Integral(2.0*mass_ele,mass_mes);
  Double_t iZhenMCP = fZhenMCP->Integral(4.0*mass_mcp*mass_mcp,mass_mes*mass_mes);
  Double_t iZhenEle = fZhenEle->Integral(4.0*mass_ele*mass_ele,mass_mes*mass_mes);
  

  // t2k calculation (includes form factors)
  Double_t BRMesonPhotonPhoton;
  Double_t BRMesonDalitz;
  //Double_t FormFactor;
  if ( mass_mes == pdgCrossSection.GetParticle(111)->Mass() ) {
    // pi0 -> gamma gamma
    BRMesonPhotonPhoton = 0.98823;
    BRMesonDalitz = 0.01174; // pdg
  } else if ( mass_mes == pdgCrossSection.GetParticle(221)->Mass() ) {
    // eta -> gamma gamma
    BRMesonPhotonPhoton = 0.3941;
    BRMesonDalitz = 0.0069; // pdg
  }
  else {
    cout << "ERROR CrossSection.cxx: ";
    cout << "COULD NOT COMPUTE CROSS SECTION" << endl;
    return 0;
  }
  
  //return i; // to check the integrals

  if ( mode == "zhenliu" ) {
    return BRMesonDalitz*(iZhenMCP/iZhenEle);
  } else {
    cout << "ERROR CrossSection.cxx: no mode selected" << endl;
    return 0;
  }
  //return BRMesonPhotonPhoton*i2; // for testing
  // factors 2 c_m and N_{POT} don't enter BR calculation
  // eq. (2) is for total number of mCP produced
}
