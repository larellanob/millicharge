#include "Root/UbooneAcceptanceChecker.cxx"
#include "Root/CrossSection.cxx"
#include "Root/DetectorInteraction.cxx"
#include "Root/CreateEmptyDirs.cxx"

void GenerateHistogramsDF(TString fstr,TString detector) {
  std::cout << "Generating histograms using RDataFrame for file " << fstr << std::endl;

  
  
}

void GenerateHistograms(TString fstr,TString detector) {
  CreateEmptyDirs();
  std::cout << "Generating histograms for file " << fstr << std::endl;
  TH1 * AccE = new TH1F("AccE",
			"mCPs through detector;Energy (GeV);Entries (unweighted)",
			30,0,12);
  TH1 * AccL = new TH1F("AccL",
			"Detector distance crossed;Distance (cm);Entries (unweighted)",
			30,0,600);
  TH1 * AccEw = new TH1F("AccEw",
			 "mCPs through detector;Energy (GeV);Entries (weighted)",
			 30,0,12);
  TH1 * AccLw = new TH1F("AccLw",
			 "Detector distance crossed;Distance (cm);Entries (weighted)",
			 30,0,600);

  TH1 * AcceEw = new TH1F("AcceEw",
			  "Electron recoil energies (one electron per mCP);Energy (MeV);Entries (weighted)",
			  30,0,300);
  TH1 * AcceE = new TH1F("AcceE",
			 "Electron recoil energies (one electron per mCP);Energy (MeV);Entries (unweighted)",
			 30,0,300);

  /*
    const int pi0_points = 5;
  const int pi0_points = 8;
  
  TGraph * sim_pi0_flux = new TGraph();
  TGraph * sim_eta_flux = new TGraph();
  */

  TH1 * AccFlux = new TH1F("AccFlux",
			   "mCPs passing though "+detector+";arbitrary;# mCP passing through",
			   1,0,1);
  
  // read accepted mcp file
  TFile *f  = new TFile(fstr);

  TTreeReader reader_meta("Metadata",f);
  TTreeReaderValue<Double_t> mass(reader_meta,"mCPmass");
  TTreeReaderValue<Double_t> charge(reader_meta,"mCPcharge");
  // the generated mCP files have one cross section but that's not
  // necesairly what needs to be used
  TTreeReaderValue<Double_t> xsec(reader_meta,"CrossSection"); 
  TTreeReaderValue<TString> meson(reader_meta,"Mother");
  TTreeReaderValue<TString> mode(reader_meta,"HornMode");
  reader_meta.Next();

  
  TTreeReader reader("mCP",f);
  TTreeReaderValue<TLorentzVector> Mom(reader,"Mom");
  TTreeReaderValue<TLorentzVector> Pos(reader,"Pos");
  TTreeReaderValue<Float_t> weight_decay(reader,"WeightDecay");
  TTreeReaderValue<Float_t> weight_meson(reader,"WeightMeson");
  TTreeReaderValue<Int_t> event(reader,"Event");
  TTreeReaderValue<Double_t> diff_xsec(reader,"DiffCrossSection");

  Int_t events = 0;
  // event loop
  while ( reader.Next() ) {
    events++;

    Double_t travelled = UbooneAcceptanceChecker(Pos->Vect(),Mom->Vect());
    AccEw->Fill(Mom->E(),(*weight_meson)*(*weight_decay));
    AccLw->Fill(travelled,(*weight_meson)*(*weight_decay));
    AccE->Fill(Mom->E());
    AccL->Fill(travelled);
    // returns in MeV!!
    Double_t recoil_electron_energy = DetectorInteraction(Mom->E(),Mom->M(),*charge);
    AcceEw->Fill(recoil_electron_energy,(*weight_meson)*(*weight_decay));
    AccFlux->Fill(0.5,(*weight_meson)*(*weight_decay));
  }

  // Flux normalization
  double POT;
  if ( detector == "argoneut" ) {
    POT = 1e20/500000.;
  } else if ( detector == "dune" ) {
    POT = 1e21/500000.;
  } else if ( detector == "duneOrnella" ) {
    POT = 3*1e22/500000.;
  } else if ( detector == "uboone" ) {
    POT = 1e21/500000.;
  }
  std::cout << "mCP Weights " << AccFlux->GetBinContent(1) << std::endl;
  double phase_space_supression = CrossSection(*mass,*meson,"physrevd");
  double epsilon = 0.01;
  AccFlux->Scale(POT*phase_space_supression*epsilon*epsilon);
  std::cout << "mass, meson " << *mass << " " << *meson << std::endl;
  std::cout << "PS supression " << phase_space_supression << std::endl;
  std::cout << "POT " << POT << std::endl;
  std::cout << "epsilon^2 " << epsilon*epsilon << std::endl;
  // then you get the value with
  // AccFlux->GetBinContent(1);
  
  // Accepted histograms
  TString AccHistoFilename =
    Form("hist/Acc_%s_q_%0.3f_m_%0.3f_%s_%ss.root",
	 detector.Data(),*charge,*mass,mode->Data(),meson->Data());
  TFile *AccFileOut = new TFile(AccHistoFilename,"recreate");
  AccEw->Write();
  AccLw->Write();
  AccE->Write();
  AccL->Write();
  AcceEw->Write();
  AcceE->Write();
  AccFlux->Write();
  delete AccEw;
  delete AccLw;
  delete AccE;
  delete AccL;
  delete AcceEw;
  delete AcceE;
  delete AccFlux;
  AccFileOut->Close();
}

void PlotAcceptedVariables(Bool_t Generate = true)
{
  CreateEmptyDirs();
  std::vector<TString> mass_points_pi0 =
    {
     "0.010",
     "0.020",
     "0.030",
     "0.050",
     "0.060"
    };

  std::vector<TString> mesons =
    {
     "pi0",
     "eta"
    };

  std::vector<TString> vars =
    {
     "E",
     "eEw"
    };

  std::vector<TString> detectors =
    {
     //"dune",
     //"uboone",
     "argoneut"
    };

  TString input_dir = "hist/";
  TFile *f = new TFile();
  TH1F *h1;
  TH1F *htitle;
  std::vector<TH1F *> h1v;
  
  if ( Generate == true ) {
    for ( TString det: detectors ) {
      for ( TString meson: mesons ) {
	for ( TString mass: mass_points_pi0 ) {
	  TString generate_file = "sim/mCP_"+det+"_q_0.010_m_"+mass+"_fhc_"+meson+"s.root";
	  GenerateHistograms(generate_file,det);
	}
      }
    }
  }
  
  gStyle->SetOptStat(0);
  for ( TString det: detectors ) {
    for ( TString meson: mesons ) {
      for ( TString var: vars ) {
	TLegend *leg;
	auto *c1 = new TCanvas();
	int plot_counter = 0;
	Double_t HistMaximum = 0;
	for ( TString mass: mass_points_pi0 ) {
	  f->Open(input_dir+"Acc_"+det+"_q_0.010_m_"+mass+"_fhc_"+meson+"s.root");
	  h1 = (TH1F*)gDirectory->Get("Acc"+var);
	  if ( plot_counter == 0 ) {
	    h1->Draw("hist L");
	    htitle = h1;
	  } else if ( plot_counter > 0 ) {
	    h1->Draw("hist same L");
	  }
	  if ( h1->GetMaximum() > HistMaximum ) {
	    HistMaximum = h1->GetMaximum();
	  }
	  h1->SetLineColor(plot_counter+1);
	  h1->SetLineWidth(2);
	  h1->SetTitle("mCP mass: "+mass+" GeV");
	  
	  plot_counter++;
	  
	}
	htitle->SetMaximum(HistMaximum+0.1*HistMaximum);
	leg = c1->BuildLegend(0.5,0.7,0.9,0.9,"","");
	if ( var == "E" || var == "Ew" ) {
	  htitle->SetTitle("mCPs from "+meson+" in "+det);
	} else if ( var == "eE" || var == "eEw" ) {
	  htitle->SetTitle("Electron recoil energy from "+meson+" in "+det);
	}
	c1->SaveAs("img/PassingThroughDetector/Acc_"+var+"_"+det+"_"+meson+"s.png");
      }
    }
  }
    
}
