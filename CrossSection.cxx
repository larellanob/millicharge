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

Double_t CrossSection(Double_t massmcp= 0.0001, Double_t massmeson = 9.0 )
{
  // this might generate a warning of 'object already instantiated
  TDatabasePDG pdg;
  //massmeson = pdg.GetParticle(221)->Mass();
  Double_t x = (massmcp*massmcp)/(massmeson*massmeson);
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
  Double_t BRMesonPhotonPhoton;
  if ( massmeson == pdg.GetParticle(111)->Mass() ) {
    // pi0 -> gamma gamma
    BRMesonPhotonPhoton = 0.98823;
  } else if ( massmeson == pdg.GetParticle(221)->Mass() ) {
    // eta -> gamma gamma
    BRMesonPhotonPhoton = 0.3941;
  }
  Double_t alpha = 1./137.;
  //Double_t alpha = 1.;

  // eq. (2) in arXiv 1812.03998v2
  //return epsilon*epsilon*BRMesonPhotonPhoton*alpha*i;
  //return i; // to check the integral

  // eq. (2) in PhysRevD is different from arXiv!!
  return epsilon*epsilon*BRMesonPhotonPhoton*alpha*i2;
  //return BRMesonPhotonPhoton*i2; // for testing
  // factors 2 c_m and N_{POT} don't enter BR calculation
  // eq. (2) is for total number of mCP produced
}
