// gives differential cross section as shown in
// arXiv 2010.07941v1

Double_t DiffCrossSection(TLorentzVector mcp,
			  TLorentzVector mcpbar,
			  Double_t epsilon,
			  TString meson,
			  Double_t mesonmass)
{
  Double_t alpha = 1./137.;
  TLorentzVector g = (mcp)+(mcpbar); // virtual photon resulting in mcp pair
  Double_t s = g.M2();
  //g.Print();
  // boost vector to virtual photon rest frame
  /* BoostVector() already implements this
  TVector3 boo(
	       g.Beta()*cos(g.Phi())*sin(g.Theta()),
	       g.Beta()*sin(g.Phi())*sin(g.Theta()),
	       g.Beta()*cos(g.Theta())
	       );
  */
  TVector3 boo = g.BoostVector();
  // check: if you perform this boost and print it, g should be in rest frame
  //g.Boost(-boo);
  //g.Print();
  TLorentzVector mcp_in_grest = mcp; // copy of mcp tlorentzvector
  mcp_in_grest.Boost(-boo); // mcp in photon rest frame
  //mcp_in_grest.Print();
  Double_t theta = mcp_in_grest.Angle(boo);

  Double_t factor1 = (epsilon*epsilon*alpha)/(4.0*TMath::Pi()*s);
  Double_t factor2 = pow(1.0 - (s/(mesonmass*mesonmass)),3);
  Double_t factor3 = sqrt(1.0 - (4.0*mcp.M()*mcp.M())/s);
  Double_t factor4 = 2.0 - (1.0 - (4.0*mcp.M()*mcp.M())/s)*pow(sin(theta),2);
  Double_t factor5;
  if ( meson == "pi0" ) {
    factor5 = 0.98823;
  } else if ( meson == "eta" ) {
    factor5 = 0.3941;
  }


  return factor1*factor2*factor3*factor4*factor5;
  
}
  
