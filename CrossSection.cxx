Double_t I3(Double_t *t, Double_t *p ) {
  double z = t[0]; // integration variable
  double x = p[0]; // masses as parameter
  return (2.0)/(3.0*TMath::Pi()) * sqrt(1.0-((4.0*x)/z)) * (1-z)/(z*z)
    * (12.0*x*x*x + 6.0*x*x*(3.0*z-2.0) + x*(5.0*z-2)*(z-1.0)+z*(z-1.0)*(z-1.0));
  
}

Double_t I3PhysRevD(Double_t *t, Double_t *p ) {
  double z = t[0]; // integration variable
  double x = p[0]; // masses as parameter
  return (2.0)/(3.0*TMath::Pi()) * sqrt(1.0-((4.0*x)/z))
    *(((1.0-z)*(1.0-z)*(1.0-z)*(2.0*x + z))/(z*z));
  
}

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
    * pow((1.0-(q2/MP2)),3) * FormFactor(q2,MP);
}


Double_t CrossSection(Double_t mass_mcp= 0.0001, Double_t mass_mes = 9.0, TString mode = "zhenliu" )
{
  // this might generate a warning of 'object already instantiated
  TDatabasePDG pdgCrossSection;
  //mass_mes = pdgCrossSection.GetParticle(221)->Mass();
  //mass_mcp= pdgCrossSection.GetParticle(11)->Mass();
  Double_t x = (mass_mcp*mass_mcp)/(mass_mes*mass_mes);
  //Double_t x = (mass_mcp*mass_mcp*mass_mcp*mass_mcp)/(mass_mes*mass_mes*mass_mes*mass_mes);
  if ( x == 0 ) {
    cout << "ERROR: mCP mass 0?" << endl;
    //return 0;
  }
  TF1 *f = new TF1("I3",I3,4.0*x,1.0,1);
  TF1 *f2 = new TF1("I3PhysRevD",I3PhysRevD,4.0*x,1.0,1);
  f->SetParameter(0,x);
  f2->SetParameter(0,x);
  //f2->Draw();
  Double_t i = f->Integral(4*x,1.0);
  Double_t i2 = f2->Integral(4*x,1.0);
  Double_t epsilon = 0.01;
  Double_t alpha = 1./137.;
  //Double_t alpha = 1.;


  Double_t mass_ele = pdgCrossSection.GetParticle(11)->Mass();
  TF1 *fZhenMCP = new TF1("I3ZhenMCP",I3ZhenLiu,4.0*mass_mcp*mass_mcp,mass_mes*mass_mes,2);
  TF1 *fZhenEle = new TF1("I3ZhenEle",I3ZhenLiu,4.0*mass_ele*mass_ele,mass_mes*mass_mes,2);
  fZhenMCP->SetParameters(mass_mcp,mass_mes);
  fZhenEle->SetParameters(mass_ele,mass_mes);

  //Double_t iZhenMCP = fZhenMCP->Integral(2.0*mass_mcp,mass_mes);
  //Double_t iZhenEle = fZhenEle->Integral(2.0*mass_ele,mass_mes);
  Double_t iZhenMCP = fZhenMCP->Integral(4.0*mass_mcp*mass_mcp,mass_mes*mass_mes);
  Double_t iZhenEle = fZhenEle->Integral(4.0*mass_ele*mass_ele,mass_mes*mass_mes);
  
  
  Double_t BRMesonPhotonPhoton;
  Double_t BRMesonDalitz;
  //Double_t FormFactor;
  if ( mass_mes == pdgCrossSection.GetParticle(111)->Mass() ) {
    // pi0 -> gamma gamma
    BRMesonPhotonPhoton = 0.98823;
    BRMesonDalitz = 0.0174; // pdg
    //FormFactor = pow(1.0+ 0.11*((mass_mcp*mass_mcp)/(mass_mes*mass_mes)),2);
    //return 0.023*epsilon*epsilon;
  } else if ( mass_mes == pdgCrossSection.GetParticle(221)->Mass() ) {
    // eta -> gamma gamma
    BRMesonPhotonPhoton = 0.3941;
    BRMesonDalitz = 0.0069; // pdg
    //FormFactor = pow(1.0- ((mass_mcp*mass_mcp)/(0.72*0.72)),-2);
    //return 0.014*epsilon*epsilon;
  }
  else {
    cout << "ERROR CrossSection.cxx: ";
    cout << "COULD NOT COMPUTE CROSS SECTION" << endl;
    return 0;
  }
  
  //return i; // to check the integrals

  if ( mode == "naive" ) {
    // eq. (2) in arXiv 1812.03998v2
    return epsilon*epsilon*BRMesonPhotonPhoton*alpha*i;
  }  else if ( mode == "physrevd" ) {
    // eq. (2) in PhysRevD is different from arXiv!!
    return epsilon*epsilon*BRMesonPhotonPhoton*alpha*i2;
  } else if ( mode == "zhenliu" ) {
    //return 2.0*BRMesonDalitz*(iZhenMCP/iZhenEle)*epsilon*epsilon;
    /*
    cout << "Integral check: Dalitz from gamma gamma:" << endl;
    cout << (2.0*alpha)/(3.0*TMath::Pi()) * iZhenEle * BRMesonPhotonPhoton << endl;
    cout << "CrossSection.cxx Result: "
	 << 2.0*BRMesonDalitz*(iZhenMCP/iZhenEle)*epsilon*epsilon << endl;
    */
    //cout << "FORMFACTOR " << FormFactor << endl;
    
    //return 2.0*BRMesonDalitz*(iZhenMCP/iZhenEle)*epsilon*epsilon*FormFactor;
    return 2.0*BRMesonPhotonPhoton*((2.0*alpha)/(3.0*TMath::Pi()))*(iZhenMCP)*epsilon*epsilon;//*FormFactor;
    // the 2.0 factor is not due to particle/antiparticle
  }
  else {
    cout << "ERROR CrossSection.cxx: no mode selected" << endl;
    return 0;
  }
  //return BRMesonPhotonPhoton*i2; // for testing
  // factors 2 c_m and N_{POT} don't enter BR calculation
  // eq. (2) is for total number of mCP produced
}
