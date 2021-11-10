//#include "../Root/DecayFactor.cxx"
#include "../Root/UbooneAcceptanceChecker.cxx"
#include "../Root/DetectorInteraction.cxx" // only used when generating electrons

class TmCP
{
private:
  //double fProbability;    // final result, probability of n hit interaction with detector
  //TVector3 fUboone_entry; // entry point to the microboone detector (det. coordinates in cm)
  //TVector3 fUboone_exit;  // exit point of microboone detector (det. coordinates in cm)
  
  //int fNhits;                  // number of (true) hits for the signal event

  // kinematic parameters
  double fEchi;
  double fMchi;
  double fCharge;
  double fEmax;
  double fCross_section;       // this will be calculated from another function
  double fMean_path;           // mean free path of mcp

  // detector parameters  
  double fEmin;
  double fDetector_resolution = 0.005; // detector resolution being considered (in m)
  double fDistance_travelled; // distance travelled in detector 
  bool fVerbose = false;
  
public:

  TmCP(TLorentzVector mcp_mom, TLorentzVector mcp_pos, double Emin, double mass, double charge);

  double Emax();
  double CrossSection();
  double MeanPath();
  double GetProbability(int nhits);

  
};


// constructor
TmCP::TmCP(
	   TLorentzVector mcp_mom,
	   TLorentzVector mcp_pos,
	   double Emin,
	   double mass,
	   double charge
	   )
{
  //fNhits = nhits; // number of hits being considered for signal
  fEmin  = Emin;  // low-energy recoil threshold
  fCharge = charge;
  fEchi = mcp_mom.E()*1000;
  fMchi = mass*1000;
  fEmax = Emax();

  
  // fMchi need to get it somehow, as well as fCharge
  // make 100% sure to get fMchi in MeV!!!!
  fDistance_travelled = UbooneAcceptanceChecker(mcp_pos.Vect(),mcp_mom.Vect())*0.01; // cm to m
  fCross_section = CrossSection();
  fMean_path = MeanPath();

  if ( fVerbose ) {
    std::cout << Form("mass %.3f, charge %.6f, Echi (MeV) %.3f ",fMchi, fCharge, fEchi) << std::endl;
    std::cout << Form("meanpath %.3f, travelled %.3f, det reso %.3f", fMean_path, fDistance_travelled, fDetector_resolution) << std::endl;
  }
}


double TmCP::Emax() {
  double m_e         = 0.511;
  double result = ((fEchi*fEchi-fMchi*fMchi)*m_e)/(fMchi*fMchi+2*fEchi*m_e+m_e*m_e);
  if ( fVerbose ) {
    std::cout << "Emax " << result << std::endl;
    std::cout << "fMchi " << fMchi << std::endl;
  }
  return result;
}

double TmCP::CrossSection() {
  double alpha       = 1./137.;
  double m_e         = 0.511;   // mass of electron in MeV
  double constant    = TMath::Pi()*alpha*alpha*fCharge*fCharge;
  // please make sure all of these are in MeV
  double denominator = (fEmax*fEmin*(fEchi*fEchi-fMchi*fMchi)*m_e*m_e);
  double first_term  = m_e*(fEmax-fEmin)*(2*fEchi*fEchi + fEmax*fEmin);
  double second_term = fEmax*fEmin*(fMchi*fMchi+m_e*(2*fEchi+m_e))*log(fEmax/fEmin);
  double result      = constant*(first_term-second_term)/denominator;
  if ( fVerbose ) {
    std::cout << "-------- " << std::endl;
    std::cout << fCharge << std::endl;
    std::cout << "const " << constant << std::endl;
    std::cout << "denom " << denominator << std::endl;
    std::cout << "first " << first_term << std::endl;
    std::cout << "second " << second_term << std::endl;
    std::cout << "CrossSection " << result << std::endl;
  }
  return result;
  
}
  
double TmCP::MeanPath() {
  double Z = 18.;
  double NA = 6.022e23;
  double rho = 1.3954;
  double m_a = 39.948;
  double n_det = NA*rho/m_a; // number density in 1/cm^3
  // now need to convert cross section from 1/MeV^2 to cm^2
  // 1 MeV = 8.06554e9 cm^-1
  double xsec_cm = fCross_section/(8.06554e9*8.06554e9);
  double result = 1./(Z*n_det*xsec_cm); // mean free path in cm
  result *= 0.01;                       // mean free path in m
  if ( fVerbose ) {
    std::cout << "xsec_cm " <<  fCross_section/(8.06554e9*8.06554e9) << std::endl;
    std::cout << "MeanPath " << result << std::endl;
  }
  return result;
}

double TmCP::GetProbability(int nhits) {
  if ( fMean_path < 0 ) {
    if ( fVerbose ) {
      std::cout << "GetProbability:: Mean path < 0! " << std::endl;
    }
    return 0;
  }
  double subdivisions = fDistance_travelled/fDetector_resolution;
  //double result = TMath::Binomial(subdivisions,nhits) * pow(fDetector_resolution/fMean_path,nhits) * pow(1.-fDetector_resolution/fMean_path,subdivisions-nhits);
  double result = TMath::Binomial(subdivisions,nhits) * pow(fDetector_resolution/fMean_path,nhits);
  

  if ( fVerbose ) {
    std::cout << "binomial " << TMath::Binomial(subdivisions,nhits) << std::endl;
    std::cout << "pow " << pow(fDetector_resolution/fMean_path,nhits) << std::endl;
    std::cout << "Probability " << result << std::endl;
  }
  return result;
}
