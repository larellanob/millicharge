#include "Root/DecayFactor.cxx"
//#include "Root/UbooneAcceptanceChecker.cxx"
//#include "Root/DetectorInteraction.cxx"
#include "Root/CreateEmptyDirs.cxx"
#include "mCP/mCP.cxx"

void GenerateSensitivityHistograms(
		     TString fstr = "sim/mCP_uboone_q_0.010_m_0.060_fhc_pi0s.root",
		     double Ethres = 0.6  // minimum energy threshold in MeV
		     )
{
  CreateEmptyDirs();
  std::cout << "initial " << std::endl;
  // opens files containing mCP particles

  // INITIAL VALUES
  //Double_t Ethres = 0.6; // in MeV
  Double_t POT_norm = 1e21/500000.;  // POT normalization

  
  // reading input file which gives the flux
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

  TTreeReader reader_meta("Metadata",f);
  TTreeReaderValue<Double_t> mass_mcp(reader_meta,"mCPmass");
  TTreeReaderValue<Double_t> charge(reader_meta,"mCPcharge");
  // the generated mCP files have one cross section but that's not
  // necesairly what needs to be used
  TTreeReaderValue<Double_t> xsec(reader_meta,"CrossSection"); 
  TTreeReaderValue<TString> meson_name(reader_meta,"Mother");
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
  
  // event loop for MICROBOONE
  while ( reader.Next() ) {

    events++;

    /////////////////////
    /////////////////////
    /////////////////////
    /////////////////////
    // FLUX
    // event weights
    
    value3+= (*weight_meson)*(*weight_decay);
    value2+= (*weight_meson);
    value1++;
    
    // DetectorInteraction returns the energy of the recoil electron
    //DetectorInteraction(Mom->E(),Mom->M(),*charge);
    


    /////////////////////
    /////////////////////
    // PROBABILITY n-HITS
    // here I have to loop ever charges, not just use *charge

    // event characteristics
    
    // travelled from cm to m
    Double_t travelled = UbooneAcceptanceChecker(Pos->Vect(),Mom->Vect());
    travelled *= 0.01;

    if ( travelled < 0.025 ) { // need travelled > (det_resolution x nhits)
      continue;                // 0.025 > 0.005 * 5 (5 hits at 5mm resolution)  
    }
    for ( int i = 0; i < 29; i++ ) {

      // prob of n hits and sum
      TmCP eve(*Mom,*Pos,Ethres, *mass_mcp, charges[i]);
      probs1hit[i] += (*weight_meson)*(*weight_decay)*eve.GetProbability(1);
      probs2hit[i] += (*weight_meson)*(*weight_decay)*eve.GetProbability(2);
      probs3hit[i] += (*weight_meson)*(*weight_decay)*eve.GetProbability(3);
      probs4hit[i] += (*weight_meson)*(*weight_decay)*eve.GetProbability(4);

      /*
      probs1hit[i] += (*weight_meson)*(*weight_decay)*prob_hit(meanpath,travelled,1,det_resolution);
      probs2hit[i] += (*weight_meson)*(*weight_decay)*prob_hit(meanpath,travelled,2,det_resolution);
      probs3hit[i] += (*weight_meson)*(*weight_decay)*prob_hit(meanpath,travelled,3,det_resolution);
      probs4hit[i] += (*weight_meson)*(*weight_decay)*prob_hit(meanpath,travelled,4,det_resolution);
      */

    }
  }
  //////////////
  // finished looping over file
  
  std::cout << "============ PlotLimits.cxx output: ============" << std::endl;
  std::cout << "******************************************************" << std::endl;
  std::cout << Form("Decays of %s into mCPs of mass %.3f and charge %.3f",
	       meson_name->Data(),*mass_mcp,*charge) << std::endl;
  std::cout << "         Using geometry and POTs" << std::endl;
  std::cout << "******************************************************" << std::endl;
  std::cout << Form("finished looping %i events",events) << std::endl;
  

  // cross section correction
  Double_t decay_factor = DecayFactor(*mass_mcp,*meson_name,"zhenliu");

  std::cout << "Simulated particles entering detector geometry: " << value1 << std::endl;
  std::cout << "TTree cross section: " << *xsec << std::endl;
  std::cout << "Unweighted mCPs passing: " << value1 << std::endl;
  

  // renormalization using POTs and decay factor
  /*
  for ( int i = 0; i < 29; i++ ) {
    probs1hit[i] = probs1hit[i]*decay_factor*POT_norm;
    probs2hit[i] = probs2hit[i]*decay_factor*POT_norm;
    probs3hit[i] = probs3hit[i]*decay_factor*POT_norm;
    probs4hit[i] = probs4hit[i]*decay_factor*POT_norm;
  }
  */

  // DO LIMITS
  // UNCOMMENT
  std::cout << "all good" << std::endl;

  TString histfilename = Form("hist/Lim_uboone_%.1f.root",Ethres);
  TFile *LimitsFileUpdate = new TFile(histfilename,"update");
  TH2F * lim1hit;
  lim1hit = (TH2F*)gDirectory->Get("lim1hit");
  TH2F * lim2hit;
  lim2hit = (TH2F*)gDirectory->Get("lim2hit");
  TH2F * lim3hit;
  lim3hit = (TH2F*)gDirectory->Get("lim3hit");
  TH2F * lim4hit;
  lim4hit = (TH2F*)gDirectory->Get("lim4hit");
  //LimitsFileUpdate->Write();


  // add points for this mass to the limits plot
  for ( int i = 0; i < 29; i++ ) {
    if ( i == 10 ) {
      	std::cout << Form("Charge %.3f", charges[i]) << std::endl;
      	std::cout << Form("decay factor %.10f", decay_factor) << std::endl;
	std::cout << "Prob 1 hit: " << probs1hit[i] << std::endl;
	std::cout << "Prob 2 hit: " << probs2hit[i] << std::endl;
	std::cout << "Prob 3 hit: " << probs3hit[i] << std::endl;
	std::cout << "Prob 4 hit: " << probs4hit[i] << std::endl;

	std::cout << Form("filling mass %.3f, charge %.3f with %.3f",*mass_mcp*1000.,charges[i],POT_norm*decay_factor*probs2hit[i]*charges[i]*charges[i]) << std::endl;
    }
    // it has an epsilon^2 because i'm leaving that out of the decay factor
    lim1hit->Fill(*mass_mcp*1000.,charges[i],POT_norm*decay_factor*probs1hit[i]*charges[i]*charges[i]);
    lim2hit->Fill(*mass_mcp*1000.,charges[i],POT_norm*decay_factor*probs2hit[i]*charges[i]*charges[i]);
    lim3hit->Fill(*mass_mcp*1000.,charges[i],POT_norm*decay_factor*probs3hit[i]*charges[i]*charges[i]);
    lim4hit->Fill(*mass_mcp*1000.,charges[i],POT_norm*decay_factor*probs4hit[i]*charges[i]*charges[i]);
    
  }

  

  lim1hit->Write();
  lim2hit->Write();
  lim3hit->Write();
  lim4hit->Write();

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".2e");
  TString detector = "uboone";
  auto DrawLimits = [detector,Ethres](TH2F* lim, int nhits)
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
			= Form("img/Limits/Sig_class_%ihit_%s_%.1f.png",
			       nhits,detector.Data(),Ethres);
		      c1->SaveAs(outname);
		      c1->SetCanvasSize(3000,1800);
		      outname
			  = Form("img/Limits/Sig_class_%ihit_%s_%.1f_text.png",
				 nhits,detector.Data(),Ethres);
		      lim->Draw("colz text");
		      c1->SaveAs(outname);
		      delete c1;
		    };

  DrawLimits(lim1hit,1);
  DrawLimits(lim2hit,2);
  DrawLimits(lim3hit,3);
  DrawLimits(lim4hit,4);

  


}
