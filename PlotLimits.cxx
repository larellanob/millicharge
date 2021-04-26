// Reads mCP file in sim/ folder and returns "geometric acceptance"
// and Number of mCPs through ArgoNeuT

// uses cross section from eq (2) in arXiv:1812.03998v2
// or the PhysRevD associated with it
#include "Root/CrossSection.cxx"

#include "Root/UbooneAcceptanceChecker.cxx"

#include "Root/DetectorInteraction.cxx"


Double_t mean_path(Double_t epsilon, Double_t Emin )
{
  // returns mean path in km
  return pow((0.01)/epsilon,2)*(Emin/0.001)*1.;
}

Double_t prob_hit(Double_t mean_path, Double_t Length, Int_t hits)
{
  return 1./TMath::Factorial(hits) * pow(Length/mean_path,hits);
}

Double_t nbkghits(Double_t empty, Int_t hits, double frames)
{
  // usage example
  //cout << nbkghits(0.88,1,3.26e6) << endl;
  //return 0;
  
  // empty = percentage of empty frames
  if ( empty > 1 ) {
    cout << "ERROR: can't have that many empty frames " << endl;
    return 0;
  }
  double phit = -log(empty);
  double sum = 0;
  double oldsum = -1;
  int n = hits;

  // summing, convergence for 0.00001
  while ( abs(oldsum-sum) > 0.00001 ) {
    cout << "sum " << sum << endl;
    oldsum = sum;
    sum += TMath::Binomial(n,hits)*pow(phit,n)/TMath::Factorial(n);
    n++;
  }
  return frames * exp(-phit) * sum;
}

void PlotLimits(TString fstr = "sim/mCP_uboone_q_0.010_m_0.010_fhc_etas.root",
		int WEIGHT = 5,
		TString detector = "uboone")
{
  cout << "initial " << endl;
  // opens files containing mCP particles
  TFile *f  = new TFile(fstr);
  TTreeReader reader("mCP",f);
  TTreeReaderValue<TLorentzVector> Mom(reader,"Mom");
  TTreeReaderValue<TLorentzVector> Pos(reader,"Pos");
  TTreeReaderValue<Float_t> weight_decay(reader,"WeightDecay");
  TTreeReaderValue<Float_t> weight_meson(reader,"WeightMeson");
  TTreeReaderValue<Int_t> event(reader,"Event");
  TTreeReaderValue<Double_t> diff_xsec(reader,"DiffCrossSection");
  // loops over events
  int events = 0;
  Double_t value1 = 0; // result stored here
  Double_t value2 = 0; // result stored here
  Double_t value3 = 0; // result stored here
  Double_t value4 = 0; // result stored here
  Double_t sum_weight_decay = 0;
  Double_t sum_weight_meson = 0;
  Double_t sum_diff_xsec = 0;

  // uboone position wrt numi target (in cm)
  const TVector3 uboonepos(31387.58422,
			   3316.402543,
			   60100.2414); // proton on target in uboone coords * cm

  TTreeReader reader_meta("Metadata",f);
  TTreeReaderValue<Double_t> mass(reader_meta,"mCPmass");
  TTreeReaderValue<Double_t> charge(reader_meta,"mCPcharge");
  // the generated mCP files have one cross section but that's not
  // necesairly what needs to be used
  TTreeReaderValue<Double_t> xsec(reader_meta,"CrossSection"); 
  TTreeReaderValue<TString> meson(reader_meta,"Mother");
  TTreeReaderValue<TString> mode(reader_meta,"HornMode");
  reader_meta.Next();


  TH1 * AccE = new TH1F("AccE", "mCPs through detector;Energy (GeV);Entries (unweighted)", 30,0,12);
  TH1 * AccL = new TH1F("AccL", "Detector distance crossed;Distance (cm);Entries (unweighted)", 30,0,600);
  TH1 * AccEw = new TH1F("AccEw", "mCPs through detector;Energy (GeV);Entries (weighted)", 30,0,12);
  TH1 * AccLw = new TH1F("AccLw", "Detector distance crossed;Distance (cm);Entries (weighted)", 30,0,600);


  // vertical axis of limits plot
  Float_t charges[29] = {
			 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4,
			 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 
			 1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2, 
			 1e-1, 2e-1
  };

  // different limits for 1-hit, 2-hit, 3-hit, 4-hit
  Float_t probs1hit[29];
  Float_t probs2hit[29];
  Float_t probs3hit[29];
  Float_t probs4hit[29];

  for ( int i = 0; i < 29; i++ ) {
    probs1hit[i] = 0;
    probs2hit[i] = 0;
    probs3hit[i] = 0;
    probs4hit[i] = 0;
  }
  
  // event loop
  while ( reader.Next() ) {
    events++;

    // event weights
    sum_weight_decay += *weight_decay;
    sum_weight_meson += *weight_meson;
    sum_diff_xsec += *diff_xsec;
    
    
    Double_t travelled = UbooneAcceptanceChecker(Pos->Vect(),Mom->Vect());
    
    value4+= (*weight_meson)*(*weight_decay)*(*diff_xsec);
    value3+= (*weight_meson)*(*weight_decay);
    value2+= (*weight_meson);
    value1++;
    AccEw->Fill(Mom->E(),(*weight_meson)*(*weight_decay));
    AccLw->Fill(travelled,(*weight_meson)*(*weight_decay));
    AccE->Fill(Mom->E());
    AccL->Fill(travelled);
    
    // DetectorInteraction returns the energy of the recoil electron
    // will be useful later
    //DetectorInteraction(Mom->E(),Mom->M(),*charge);
    

    // here I have to loop ever charges, not just use *charge
    // travelled from cm to km
    travelled *= 0.00001;
    for ( int i = 0; i < 29; i++ ) {
      Double_t meanpath = mean_path(charges[i],0.001);
      /*
	if ( travelled/meanpath > 0.01 ) {
	cout << Form("WARNING: mean path (%.3f) close to travel length (%.3f) for charge %.3f",meanpath,travelled,charges[i]) << endl;
	}
      */
	
      probs1hit[i] += (*weight_meson)*(*weight_decay)*prob_hit(meanpath,travelled,1);
      probs2hit[i] += (*weight_meson)*(*weight_decay)*prob_hit(meanpath,travelled,2);
      probs3hit[i] += (*weight_meson)*(*weight_decay)*prob_hit(meanpath,travelled,3);
      probs4hit[i] += (*weight_meson)*(*weight_decay)*prob_hit(meanpath,travelled,4);
    }
  }

  
  cout << "============ PlotLimits.cxx output: ============" << endl;
  cout << "******************************************************" << endl;
  cout << Form("Decays of %s into mCPs of mass %.3f and charge %.3f",
	       meson->Data(),*mass,*charge) << endl;
  cout << "         Using "+detector+" geometry and POTs" << endl;
  cout << "******************************************************" << endl;
  cout << Form("finished looping %i events",events) << endl;
  

  // Accepted histograms
  TString AccHistoFilename =
    Form("hist/Acc_%s_q_%0.3f_m_%0.3f_%s_%ss.root",
	 detector.Data(),*charge,*mass,mode->Data(),meson->Data());
  TFile *AccFileOut = new TFile(AccHistoFilename,"recreate");
  AccEw->Write();
  AccLw->Write();
  AccE->Write();
  AccL->Write();

  
  // cross section correction
  TDatabasePDG pdg;
  double mesonmass = 0;
  if ( *meson == "pi0" ) {
    mesonmass = pdg.GetParticle(111)->Mass();
  } else if ( *meson == "eta" ) {
    mesonmass = pdg.GetParticle(221)->Mass();
  }
  Double_t decay_factor_zhenliu = CrossSection(*mass,mesonmass,"zhenliu");
  Double_t decay_factor_physrevd = CrossSection(*mass,mesonmass,"physrevd");
  Double_t decay_factor_naive = CrossSection(*mass,mesonmass,"naive");
  Double_t decay_factor_t2k = CrossSection(*mass,mesonmass,"t2k");
  Double_t truexsec = *xsec; // no need to recalculate

  cout << "Simulated particles entering detector geometry: " << value1 << endl;
  cout << "TTree cross section: " << *xsec << endl;
  cout << "Recalculated cross section: " << truexsec << endl;
  cout << "Unweighted mCPs passing: " << value1 << endl;
  
  /*
  cout << Form("value: %.3f\nPOT normalized value: %.3f\n(old) xsec: %.7f\nsum weight decays: %.3f\nsum weight meson: %.3f",
	       value1,(1e20/500000.)*value1/(sum_weight_decay*0.5),*xsec,sum_weight_decay/2.,sum_weight_meson/2.) << endl;
  */
  Double_t result;

  // POT normalization
  Double_t POT_norm;
  if ( detector == "argoneut" ) {
    POT_norm = 1e20/500000.;
  } else if ( detector == "dune" ) {
    POT_norm = 1e21/500000.;
  } else if ( detector == "duneOrnella" ) {
    POT_norm = 3*1e22/500000.;
  } else if ( detector == "uboone" ) {
    POT_norm = 1e21/500000.;
  }

  
  if ( WEIGHT == 1 ) {
    result = POT_norm*(truexsec)*value1; // result uses truexsec
  } else if ( WEIGHT == 2 ) {
    result = POT_norm*(truexsec)*value2; // result uses truexsec
  } else if ( WEIGHT == 3 ) {
    result = POT_norm*(truexsec)*value3*(events/sum_weight_decay); // result uses truexsec
  } else if ( WEIGHT == 4 ) {
    result = POT_norm*value4*(events/sum_weight_decay); // differential cross section
  } else if ( WEIGHT == 5 ) {
    if ( detector == "argoneut" || detector == "duneOrnella" ) {
      result = POT_norm*value3*(events/sum_weight_decay)*decay_factor_zhenliu;
      cout << "USING DECAY FACTOR decay_factor_zhenliu" << endl;
    } else if ( detector == "uboone" || detector == "dune" || detector == "naiveuboone" ) { 
      result = POT_norm*value3*(events/sum_weight_decay)*decay_factor_physrevd;
      cout << "USING DECAY FACTOR decay_factor_physrevd" << endl;
    } else {
      result = 0;
      cout << "USING DECAY FACTOR NULL"  << endl;
    }
  } else if ( WEIGHT == 6 ) {
    result = decay_factor_t2k;
  }

  // DO LIMITS
  // UNCOMMENT
  TFile *LimitsFileUpdate = new TFile("hist/Lim_uboone.root","update");
  TH2F * lim1hit;
  lim1hit = (TH2F*)gDirectory->Get("lim1hit");
  TH2F * lim2hit;
  lim2hit = (TH2F*)gDirectory->Get("lim2hit");
  TH2F * lim3hit;
  lim3hit = (TH2F*)gDirectory->Get("lim3hit");
  TH2F * lim4hit;
  lim4hit = (TH2F*)gDirectory->Get("lim4hit");

  // this are the axes of the limit plot as reference
  /*
  const Float_t xaxis[29] = {9.5,15,25,35,45,55,65,75,85,95,
		       150,250,350,450,550,650,750,850,950,
		       1500,2500,3500,4500,5500,6500,7500,8500,9500,
		       10500
  };
  
  const Float_t yaxis[30] = {
		       0.5e-4, 1.5e-4, 2.5e-4, 3.5e-4, 4.5e-4, 5.5e-4, 6.5e-4, 7.5e-4, 8.5e-4, 9.5e-4,
		       1.5e-3, 2.5e-3, 3.5e-3, 4.5e-3, 5.5e-3, 6.5e-3, 7.5e-3, 8.5e-3, 9.5e-3, 
		       1.5e-2, 2.5e-2, 3.5e-2, 4.5e-2, 5.5e-2, 6.5e-2, 7.5e-2, 8.5e-2, 9.5e-2, 
		       1.5e-1, 2.5e-1
  };


  TH2 * lim1hit = new TH2F("lim1hit",";m_{#chi} (GeV);#epsilon",28,xaxis,29,yaxis);
  TH2 * lim2hit = new TH2F("lim2hit",";m_{#chi} (GeV);#epsilon",28,xaxis,29,yaxis);
  TH2 * lim3hit = new TH2F("lim3hit",";m_{#chi} (GeV);#epsilon",28,xaxis,29,yaxis);
  TH2 * lim4hit = new TH2F("lim4hit",";m_{#chi} (GeV);#epsilon",28,xaxis,29,yaxis);
  */


  // add points for this mass to the limits plot

  for ( int i = 0; i < 29; i++ ) {
    lim1hit->Fill(*mass*1000.,charges[i],POT_norm*(events/sum_weight_decay)*decay_factor_physrevd*probs1hit[i]);
    lim2hit->Fill(*mass*1000.,charges[i],POT_norm*(events/sum_weight_decay)*decay_factor_physrevd*probs2hit[i]);
    lim3hit->Fill(*mass*1000.,charges[i],POT_norm*(events/sum_weight_decay)*decay_factor_physrevd*probs3hit[i]);
    lim4hit->Fill(*mass*1000.,charges[i],POT_norm*(events/sum_weight_decay)*decay_factor_physrevd*probs4hit[i]);
  }
  lim1hit->Write();
  lim2hit->Write();
  lim3hit->Write();
  lim4hit->Write();

  auto c1 = new TCanvas();
  c1->SetLogy();
  c1->SetLogx();
  c1->SetLogz();
  lim1hit->Draw("colz");
  c1->SaveAs("img/Lim_1hit_uboone.png");
  lim2hit->Draw("colz");
  c1->SaveAs("img/Lim_2hit_uboone.png");
  lim3hit->Draw("colz");
  c1->SaveAs("img/Lim_3hit_uboone.png");
  lim4hit->Draw("colz");
  c1->SaveAs("img/Lim_4hit_uboone.png");
  lim1hit->Draw("colz text");
  c1->SaveAs("img/Lim_1hit_uboone_text.png");
  lim2hit->Draw("colz text");
  c1->SaveAs("img/Lim_2hit_uboone_text.png");
  lim3hit->Draw("colz text");
  c1->SaveAs("img/Lim_3hit_uboone_text.png");
  lim4hit->Draw("colz text");
  c1->SaveAs("img/Lim_4hit_uboone_text.png");
  
  //cout << "geometrical acceptance " << value1/(events) << endl;
  cout << "geometrical acceptance " << value3/(sum_weight_decay) << endl;
  cout << "result " << result << endl;
  cout << "----------------------------------" << endl;

  return result;


}
