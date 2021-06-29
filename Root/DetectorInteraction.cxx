Double_t DetDiffCrossSection(Double_t *x, Double_t *par ) {
  Double_t Er = x[0];
  double Echi = par[0];
  double mchi = par[1];
  double epsilon = par[2];
  Double_t alpha = 1/137.;
  Double_t mele = 0.00051099891;
  
  // constants in front
  Double_t con = TMath::Pi()*alpha*alpha*epsilon*epsilon;
  // numerator
  Double_t num = 2*Echi*Echi*mele + Er*Er*mele - Er*(mchi*mchi + mele*(2*Echi+mele));
  // denominator
  Double_t den = Er*Er*(Echi*Echi-mchi*mchi)*mele*mele;
  return (con*num)/den;
  
}

Double_t DetTotalCrossSection(Double_t *Er, Double_t *par) {
  double Emax = Er[0];
  double Echi = par[0];
  double mchi = par[1];
  double epsilon = par[2];
  double Emin = par[3];
  Double_t alpha = 1/137.;
  Double_t mele = 0.00051099891;

  // constants in front
  Double_t con = TMath::Pi()*alpha*alpha*epsilon*epsilon;
  // numerator
  Double_t num = mele*(Emax-Emin)*(2*Echi*Echi+Emax*Emin) - Emax*Emin*(mchi*mchi + mele*(2*Echi+mele))*log(Emax/Emin);
  // denominator
  Double_t den = Emax*Emin*(Echi*Echi - mchi*mchi) * mele*mele;
  return (con*num)/den;
}

// returns probability distribution of recoil energy
Double_t DetectorInteraction(Double_t Echi, Double_t mchi, Double_t epsilon, Double_t threshold = 1000, Bool_t kPrint = false)
{
  // minimum recoil energy in GeV, given by experiment detection
  // threshold is given in keV
  Double_t Emin = threshold*1e-6;
  
  // electron mass
  Double_t mele = 0.00051099891;
  
  // maximal kinematically allowed recoil energyu
  Double_t Emax = ((Echi*Echi-mchi*mchi)*mele)/(mchi*mchi + 2*Echi*mele + mele*mele);
  if ( Emax < Emin ) {
    return 0;
  }
  Double_t *Emaxptr = &Emax;

  TF1 * dif = new TF1("DetDiffCrossSection",DetDiffCrossSection,Emin,Emax,3);
  //TF1 * tot = new TF1("DetTotalCrossSection",DetTotalCrossSection,Emin,Emax,4);
  dif->SetParameters(Echi,mchi,epsilon);
  //tot->SetParameters(Echi,mchi,epsilon,Emin);

  if ( kPrint ) {
    auto *c1 = new TCanvas();
    TTimeStamp ts;
    dif->Draw("L");
    dif->GetXaxis()->SetTitle("Energy (GeV)");
    TString savefile = Form("./img/cross_sect/DiffCrossSect_%i_%i.png",ts.GetTime(),ts.GetNanoSec());
    c1->SaveAs(savefile);
  }
  
  // cross section
  // generates random seed
  //gRandom= new TRandom3(0);
  Double_t recoil = dif->GetRandom();
  ///////////////////////////
  ////////// MeV ////////////
  ///////////////////////////
  return recoil*1000.;
  // returns recoil energy in MeV
  ///////////////////////////
  ////////// MeV ////////////
  ///////////////////////////

}


