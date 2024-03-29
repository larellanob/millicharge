#include "Root/CrossSection.cxx"
#include "Root/UbooneAcceptanceChecker.cxx"
#include "Root/DetectorInteraction.cxx"
#include "Root/CreateEmptyDirs.cxx"

Double_t emax ( double m_chi, double E_chi ) {
  double m_e         = 0.511;
  return ((E_chi*E_chi-m_chi*m_chi)*m_e)/(m_chi*m_chi+2*E_chi*m_e+m_e*m_e);
}

Double_t det_interaction(double Ermin,
			 double Ermax,
			 double m_chi,
			 double charge){
  double alpha       = 1./137.;
  double m_e         = 0.511;
  double constant    = TMath::Pi()*alpha*alpha*charge*charge;
  double denominator = (Ermax*Ermin*(E_chi*E_chi-m_chi*m_chi)*m_e*m_e);
  double first_term  = m_e*(Ermax-Ermin)*(2*E_chi*E_chi + Ermax*Ermin);
  double second_term = Ermax*Ermin*(m_chi*m_chi+m_e*(2*E_chi+m_e))*log(Ermax/Ermin);
  double result      = constant*(first_term-second_term)/denominator;
  return result;
}

Double_t mean_path(Double_t epsilon, Double_t Emin ) {
  return pow((0.01)/epsilon,2)*(Emin/0.001)*1.;
}

Double_t mean_path2(Double_t epsilon, Double_t Emin, Double_t Emax, Double_t m_chi )
{
  // returns mean path in km
  double Z  = 18.;
  double NA = 6.022e23;
  double rho = 1.3954;
  double m_a = 39.948;
  double n_det = NA*rho/m_a;

  double xsec_cm = det_interaction(Emin,Emax,m_chi,epsilon)/(8.06554e9*8.06554e9);
  double result = 1./(Z*n_det*xsec_cm); // mean free path in cm
  result *= 0.01;                       // mean fre path in m
  return result;
}

Double_t prob_hit(Double_t mean_path, Double_t Length, Int_t hits)
{
  return 1./TMath::Factorial(hits) * pow(Length/mean_path,hits);
}

Double_t prob_hit_fix(Double_t mean_path, Double_t Length, Int_t hits)
{
  Int_t subdivisions = 1000;
  Double_t DeltaL = Length/subdivisions;
  return TMath::Binomial(subdivisions,hits) * pow(DeltaL/mean_path,hits)* pow(DeltaL/mean_path,subdivisions-hits);
}
Double_t prob_hit_fix2(Double_t mean_path, Double_t Length, Int_t hits, Int_t subdivisions = 1000)
{
  //Int_t subdivisions = 1000;
  Double_t DeltaL = Length/subdivisions;
  return TMath::Binomial(subdivisions,hits) * pow(DeltaL/mean_path,hits);
}


Double_t nbkghits(Double_t empty, Int_t hits, double frames)
{
  // usage example
  //std::cout << nbkghits(0.88,1,3.26e6) << std::endl;
  //return 0;
  
  // empty = percentage of empty frames
  if ( empty > 1 ) {
    std::cout << "ERROR: can't have that many empty frames " << std::endl;
    return 0;
  }
  double phit = -log(empty);
  double sum = 0;
  double oldsum = -1;
  int n = hits;

  // summing, convergence for 0.00001
  while ( abs(oldsum-sum) > 0.00001 ) {
    std::cout << "sum " << sum << std::endl;
    oldsum = sum;
    sum += TMath::Binomial(n,hits)*pow(phit,n)/TMath::Factorial(n);
    n++;
  }
  return frames * exp(-phit) * sum;
}

void PlotSignal(TString fstr = "sim/mCP_uboone_q_0.010_m_0.010_fhc_etas.root",
		TString detector = "uboone")
{
  CreateEmptyDirs();
  std::cout << "initial " << std::endl;
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


  // vertical axis of limits plot
  Float_t charges[29] = {
			 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4,
			 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 
			 1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2, 
			 1e-1, 2e-1
  };

  // horizontal axis bins
  // {9.5,15,25,35,45,55,65,75,85,95,
  //  150,250,350,450,550,650,750,850,950,
  //  1500,2500,3500,4500,5500,6500,7500,8500,9500,
  //  10500 }
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

    if ( detector == "argoneut_published" ) {
      break;
    }


    events++;
    
    // event weights
    sum_weight_decay += *weight_decay;
    sum_weight_meson += *weight_meson;
    sum_diff_xsec += *diff_xsec;
    
    
    Double_t travelled = UbooneAcceptanceChecker(Pos->Vect(),Mom->Vect());
    if ( detector == "argoneut" ){
      // travelled in cm
      travelled = 90.;
    }
    value4+= (*weight_meson)*(*weight_decay)*(*diff_xsec);
    value3+= (*weight_meson)*(*weight_decay);
    value2+= (*weight_meson);
    value1++;
    
    // DetectorInteraction returns the energy of the recoil electron
    // will be useful later
    //DetectorInteraction(Mom->E(),Mom->M(),*charge);
    

    // here I have to loop ever charges, not just use *charge
    // travelled from cm to km
    travelled *= 0.00001;
    for ( int i = 0; i < 29; i++ ) {
      //Double_t meanpath = mean_path(charges[i],0.0008);
      Double_t meanpath = mean_path2(charges[i],0.0008,emax(*mass,Mom.E()),*mass);
      /*
      if ( travelled/meanpath > 0.01 ) {
	std::cout << Form("WARNING: mean path (%.3f) close to travel length (%.3f) for charge %.3f",meanpath,travelled,charges[i]) << std::endl;
	std::cout << "Prob 1hit 1 2 3 " << prob_hit(meanpath,travelled,1) << " " << prob_hit_fix(meanpath,travelled,1) << " " << prob_hit_fix2(meanpath,travelled,1) << std::endl;
	std::cout << "Prob 2hit 1 2 3 " << prob_hit(meanpath,travelled,2) << " " << prob_hit_fix(meanpath,travelled,2) << " " << prob_hit_fix2(meanpath,travelled,2) << std::endl;
	std::cout << "Prob 3hit 1 2 3 " << prob_hit(meanpath,travelled,3) << " " << prob_hit_fix(meanpath,travelled,3) << " " << prob_hit_fix2(meanpath,travelled,3) << std::endl;
	std::cout << "Prob 4hit 1 2 3 " << prob_hit(meanpath,travelled,4) << " " << prob_hit_fix(meanpath,travelled,4) << " "
		  << prob_hit_fix2(meanpath,travelled,4, 10) << " " << prob_hit_fix2(meanpath,travelled,4, 1000) << " " << prob_hit_fix2(meanpath,travelled,4, 1000000) << std::endl;
      }
      */
      probs1hit[i] += (*weight_meson)*(*weight_decay)*prob_hit(meanpath,travelled,1);
      probs2hit[i] += (*weight_meson)*(*weight_decay)*prob_hit(meanpath,travelled,2);
      probs3hit[i] += (*weight_meson)*(*weight_decay)*prob_hit(meanpath,travelled,3);
      probs4hit[i] += (*weight_meson)*(*weight_decay)*prob_hit(meanpath,travelled,4);
    }
  }

  
  std::cout << "============ PlotLimits.cxx output: ============" << std::endl;
  std::cout << "******************************************************" << std::endl;
  std::cout << Form("Decays of %s into mCPs of mass %.3f and charge %.3f",
	       meson->Data(),*mass,*charge) << std::endl;
  std::cout << "         Using "+detector+" geometry and POTs" << std::endl;
  std::cout << "******************************************************" << std::endl;
  std::cout << Form("finished looping %i events",events) << std::endl;
  

  // cross section correction
  Double_t decay_factor_zhenliu = CrossSection(*mass,*meson,"zhenliu");
  Double_t decay_factor_physrevd = CrossSection(*mass,*meson,"physrevd");
  Double_t decay_factor_naive = CrossSection(*mass,*meson,"naive");
  Double_t decay_factor_t2k = CrossSection(*mass,*meson,"t2k");
  Double_t truexsec = *xsec; // no need to recalculate

  std::cout << "Simulated particles entering detector geometry: " << value1 << std::endl;
  std::cout << "TTree cross section: " << *xsec << std::endl;
  std::cout << "Recalculated cross section: " << truexsec << std::endl;
  std::cout << "Unweighted mCPs passing: " << value1 << std::endl;
  
  /*
  std::cout << Form("value: %.3f\nPOT normalized value: %.3f\n(old) xsec: %.7f\nsum weight decays: %.3f\nsum weight meson: %.3f",
	       value1,(1e20/500000.)*value1/(sum_weight_decay*0.5),*xsec,sum_weight_decay/2.,sum_weight_meson/2.) << std::endl;
  */
  Double_t result;

  // POT normalization
  Double_t POT_norm;
  if ( detector == "argoneut" ) {
    POT_norm = 1e20/500000.;
  } else if ( detector == "dune" ) {
    POT_norm = 1e21/500000.;
  } else if ( detector == "duneArgo" ) {
    POT_norm = 3*1e22/500000.;
  } else if ( detector == "uboone" ) {
    POT_norm = 1e21/500000.;
  }

  
  if ( detector == "argoneut" || detector == "duneArgo" ) {
    result = POT_norm*value3*(events/sum_weight_decay)*decay_factor_zhenliu;
    std::cout << "USING DECAY FACTOR decay_factor_zhenliu" << std::endl;
  } else if ( detector == "uboone" || detector == "dune" || detector == "naiveuboone" ) { 
    result = POT_norm*value3*(events/sum_weight_decay)*decay_factor_physrevd;
    std::cout << "USING DECAY FACTOR decay_factor_physrevd" << std::endl;
  } else {
    result = 0;
    std::cout << "USING DECAY FACTOR NULL"  << std::endl;
  }


  // DO LIMITS
  // UNCOMMENT
  std::cout << "all good" << std::endl;
  TFile *LimitsFileUpdate = new TFile("hist/Lim_"+detector+".root","update");
  TH2F * lim1hit;
  lim1hit = (TH2F*)gDirectory->Get("lim1hit");
  TH2F * lim2hit;
  lim2hit = (TH2F*)gDirectory->Get("lim2hit");
  TH2F * lim3hit;
  lim3hit = (TH2F*)gDirectory->Get("lim3hit");
  TH2F * lim4hit;
  lim4hit = (TH2F*)gDirectory->Get("lim4hit");
  //LimitsFileUpdate->Write();

  if ( detector == "argoneut_published" ) {
    std::cout << "det = argo published" << std::endl;
    Double_t meanpath = mean_path(0.01,0.001);
    Double_t travelled = 90.; // 90cm = z axis detector length
    travelled*=0.0001; // in km
    Double_t argo_flux;
    if ( *mass == 0.01 ) {
      if ( *meson == "pi0" )
	argo_flux = 29053078518.935688;
      if ( *meson == "eta")
	argo_flux = 1559272026.0810509;
    }
    if ( *mass == 0.02 ) {
      if ( *meson == "pi0" )
	argo_flux = 29053078518.935688;
      if ( *meson == "eta")
	argo_flux = 1559272026.0810509;
    }
    if ( *mass == 0.03 ) {
      if ( *meson == "pi0" )
	argo_flux = 15592720260.810507;
      if ( *meson == "eta")
	argo_flux = 1559272026.0810509;
    }
    if ( *mass == 0.05 ) {
      if ( *meson == "pi0" )
	argo_flux = 7389202410.398621;
      if ( *meson == "eta")
	argo_flux = 1659391704.1669834;
    }
    if ( *mass == 0.06 ) {
      if ( *meson == "pi0" )
	argo_flux = 1765939991.098824;
      if ( *meson == "eta")
	argo_flux = 1659391704.1669834;
    }
    if ( *mass == 0.1 ) {
      if ( *meson == "eta")
	argo_flux = 1559272026.0810509;
    }
    if ( *mass == 0.2 ) {
      if ( *meson == "eta")
	argo_flux = 836857701.5803154;
    }
    if ( *mass == 0.25 ) {
      if ( *meson == "eta")
	argo_flux = 129372153.23092696;
    }
    std::cout << "meanpath  " << meanpath << " ";
    std::cout << "travelled  " << travelled << " ";
    std::cout << "prob hit " << prob_hit(meanpath,travelled,1) << std::endl;
    lim1hit->Fill( *mass*1000., 0.01, argo_flux*prob_hit(meanpath,travelled,1) );
    lim2hit->Fill( *mass*1000., 0.01, argo_flux*prob_hit(meanpath,travelled,2) );
    lim3hit->Fill( *mass*1000., 0.01, argo_flux*prob_hit(meanpath,travelled,3) );
    lim4hit->Fill( *mass*1000., 0.01, argo_flux*prob_hit(meanpath,travelled,4) );
    lim1hit->Write();
    lim2hit->Write();
    lim3hit->Write();
    lim4hit->Write();

  }
  std::cout << "not good" << std::endl;
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
  if ( detector != "argoneut_published" ) {
    for ( int i = 0; i < 29; i++ ) {
      if ( i == 7 ) {
	std::cout << "charge bin " << i << " ";
	std::cout << "charge " << charges[i] << " ";
	std::cout << "potnorm " << POT_norm << " ";
	std::cout << "events/sum_weigtdecay " << (events/sum_weight_decay) << " ";
	std::cout << "decayfact " << decay_factor_physrevd << " ";
	std::cout << "prob 1-hit " << probs1hit[i] << std::endl;
	std::cout << "result fill " << POT_norm*(events/sum_weight_decay)*decay_factor_physrevd*probs1hit[i] << std::endl;
	//lim1hit->Fill(*mass*1000.,charges[i],1000000000000000);
	//break;
      }
      /*
      lim1hit->Fill(*mass*1000.,charges[i],POT_norm*(events/sum_weight_decay)*decay_factor_physrevd*probs1hit[i]*charges[i]*charges[i]);
      lim2hit->Fill(*mass*1000.,charges[i],POT_norm*(events/sum_weight_decay)*decay_factor_physrevd*probs2hit[i]*charges[i]*charges[i]);
      lim3hit->Fill(*mass*1000.,charges[i],POT_norm*(events/sum_weight_decay)*decay_factor_physrevd*probs3hit[i]*charges[i]*charges[i]);
      lim4hit->Fill(*mass*1000.,charges[i],POT_norm*(events/sum_weight_decay)*decay_factor_physrevd*probs4hit[i]*charges[i]*charges[i]);
      */
      std::cout << "about to fill " << std::endl;
      std::cout << "mass" << std::endl;
      std::cout << *mass << std::endl;
      std::cout << "charges[i]" << std::endl;
      std::cout << charges[i] << std::endl;
      std::cout << "POTNORM " << std::endl;
      std::cout << POT_norm << std::endl;
      std::cout << "decay factor " << std::endl;
      std::cout << decay_factor_physrevd << std::endl;
      std::cout << "probs 1hit " << std::endl;
      std::cout << probs1hit[i] << std::endl;
      std::cout << "probs 2hit " << std::endl;
      std::cout << probs2hit[i] << std::endl;
      std::cout << "probs 3hit " << std::endl;
      std::cout << probs3hit[i] << std::endl;
      std::cout << "probs 4hit " << std::endl;
      std::cout << probs4hit[i] << std::endl;
      
      lim1hit->Fill(*mass*1000.,charges[i],POT_norm*decay_factor_physrevd*probs1hit[i]*charges[i]*charges[i]);
      std::cout << "filling failing" << std::endl;
      lim2hit->Fill(*mass*1000.,charges[i],POT_norm*decay_factor_physrevd*probs2hit[i]*charges[i]*charges[i]);
      lim3hit->Fill(*mass*1000.,charges[i],POT_norm*decay_factor_physrevd*probs3hit[i]*charges[i]*charges[i]);
      lim4hit->Fill(*mass*1000.,charges[i],POT_norm*decay_factor_physrevd*probs4hit[i]*charges[i]*charges[i]);

    }
    lim1hit->Write();
    lim2hit->Write();
    lim3hit->Write();
    lim4hit->Write();
  }

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".2e");
  auto DrawLimits = [detector](TH2F* lim, int nhits)
		    {
		      auto c1 = new TCanvas();
		      c1->SetLogy();
		      c1->SetLogx();
		      c1->SetLogz();
		      lim->Draw("colz");
		      lim->GetXaxis()->CenterTitle();
		      lim->GetYaxis()->CenterTitle();
		      lim->SetMarkerSize(0.7);
		      TString outname
			= Form("img/Limits/Sig_%ihit_%s.png",
			       nhits,detector.Data());
		      c1->SaveAs(outname);
		      c1->SetCanvasSize(3000,1800);
		      outname
			  = Form("img/Limits/Sig_%ihit_%s_text.png",
				 nhits,detector.Data());
		      lim->Draw("colz text");
		      c1->SaveAs(outname);
		      delete c1;
		    };

  DrawLimits(lim1hit,1);
  DrawLimits(lim2hit,2);
  DrawLimits(lim3hit,3);
  DrawLimits(lim4hit,4);

  
  std::cout << "geometrical acceptance " << value3/(sum_weight_decay) << std::endl;
  std::cout << "result " << result << std::endl;
  std::cout << "----------------------------------" << std::endl;

  return result;


}
