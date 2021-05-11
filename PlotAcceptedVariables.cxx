#include "Root/UbooneAcceptanceChecker.cxx"
#include "Root/CrossSection.cxx"
#include "Root/DetectorInteraction.cxx"
#include "Root/CreateEmptyDirs.cxx"


void LatexText(Double_t x, Double_t y, int font, TString text)
{
  TLatex l2;
  l2.SetNDC();
  l2.SetTextFont(font);
  l2.DrawLatex(x,y,text);
}


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

void PlotAcceptedVariables(Bool_t Generate = false)
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

  std::vector<TString> mass_points_eta =
    {
     "0.010",
     "0.020",
     "0.030",
     "0.050",
     "0.060",
     "0.100",
     "0.200",
     "0.250"
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
      for ( TString mass: mass_points_pi0 ) {
	TString generate_file = "sim/mCP_"+det+"_q_0.010_m_"+mass+"_fhc_pi0s.root";
	GenerateHistograms(generate_file,det);
      }
      for ( TString mass: mass_points_eta ) {
	TString generate_file = "sim/mCP_"+det+"_q_0.010_m_"+mass+"_fhc_etas.root";
	GenerateHistograms(generate_file,det);
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

  /////////////////////////////////
  /////////////////////////////////
  // Flux compared to publications
  // i.e. number of particles passing through the detector
  /////////////////////////////////
  /////////////////////////////////

  // creating and filling the tgraphcs
  TGraph * sim_pi0_flux = new TGraph();
  TGraph * sim_eta_flux = new TGraph();

  TGraph * pub_pi0_flux = new TGraph();
  TGraph * pub_eta_flux = new TGraph();
  for ( auto det: detectors ) {
    TString pub_data_filename_pi0 = "data/passing_through/"+det+"/pi0.csv";
    TString pub_data_filename_eta = "data/passing_through/"+det+"/eta.csv";
    std::ifstream pub_data_file_pi0(pub_data_filename_pi0.Data());
    std::ifstream pub_data_file_eta(pub_data_filename_eta.Data());
    TString reader;
    // takes two column csv files
    while ( reader.ReadToDelim(pub_data_file_pi0,',') ) {
      Double_t x = reader.Atof();
      reader.ReadToDelim(pub_data_file_pi0,'\n');
      Double_t y = reader.Atof();
      pub_pi0_flux->AddPoint(x,y);
    }
    while ( reader.ReadToDelim(pub_data_file_eta,',') ) {
      Double_t x = reader.Atof();
      reader.ReadToDelim(pub_data_file_eta,'\n');
      Double_t y = reader.Atof();
      pub_eta_flux->AddPoint(x,y);
    }
  }


  // drawiing properties
  pub_pi0_flux->SetLineColor(kBlue);
  pub_eta_flux->SetLineColor(kOrange);
  pub_pi0_flux->SetLineWidth(5);
  pub_eta_flux->SetLineWidth(5);

  sim_pi0_flux->SetMarkerStyle(kFullTriangleUp);
  sim_eta_flux->SetMarkerStyle(kFullSquare);
  sim_pi0_flux->SetMarkerSize(1.5);
  sim_eta_flux->SetMarkerSize(1);
  sim_pi0_flux->SetMarkerColor(kBlue+2);
  sim_eta_flux->SetMarkerColor(kOrange-2);
  
  //gStyle->SetPadTickX(1);
  //gStyle->SetPadTickY(1);

  for ( TString det: detectors ) {
    TString det_formal;
    if ( det == "uboone" || det == "naiveuboone" ) {
      det_formal = "MicroBooNE";
    } else if ( det == "dune" || det == "duneOrnella" ) {
      det_formal = "DUNE";
    } else if ( det == "argoneut" ) {
      det_formal = "ArgoNeuT";
    } else if ( det == "t2k" ) {
      det_formal = "T2K";
    } else {
      det_formal = det;
    }
    
    for ( TString mass: mass_points_pi0 ) {
      f->Open(input_dir+"Acc_"+det+"_q_0.010_m_"+mass+"_fhc_pi0s.root");
      h1 = (TH1F*)gDirectory->Get("AccFlux");
      sim_pi0_flux->AddPoint(mass.Atof(),h1->GetBinContent(1));
    }
    
    for ( TString mass: mass_points_eta ) {
      f->Open(input_dir+"Acc_"+det+"_q_0.010_m_"+mass+"_fhc_etas.root");
      h1 = (TH1F*)gDirectory->Get("AccFlux");
      sim_eta_flux->AddPoint(mass.Atof(),h1->GetBinContent(1));
    }


    Double_t global_maximum;
    double max1 = TMath::MaxElement(sim_pi0_flux->GetN(),sim_pi0_flux->GetY());
    double max2 = TMath::MaxElement(sim_eta_flux->GetN(),sim_eta_flux->GetY());
    double max3 = TMath::MaxElement(pub_pi0_flux->GetN(),pub_pi0_flux->GetY());
    double max4 = TMath::MaxElement(pub_eta_flux->GetN(),pub_eta_flux->GetY());
    global_maximum = max(max1,max2);
    global_maximum = max(global_maximum,max3);
    global_maximum = max(global_maximum,max4);
    std::cout << global_maximum << std::endl;
    global_maximum *= 100;
    sim_eta_flux->SetMaximum(global_maximum);

    Double_t global_minimum;
    double min1 = TMath::MinElement(sim_pi0_flux->GetN(),sim_pi0_flux->GetY());
    double min2 = TMath::MinElement(sim_eta_flux->GetN(),sim_eta_flux->GetY());
    double min3 = TMath::MinElement(pub_pi0_flux->GetN(),pub_pi0_flux->GetY());
    double min4 = TMath::MinElement(pub_eta_flux->GetN(),pub_eta_flux->GetY());
    global_minimum = min(min1,min2);
    global_minimum = min(global_minimum,min3);
    global_minimum = min(global_minimum,min4);
    std::cout << global_minimum << std::endl;
    global_minimum /= 10;
    sim_eta_flux->SetMinimum(global_minimum);
    

    /*
    sim_eta_flux->Draw("A P");
    sim_pi0_flux->Draw("P same");

    pub_pi0_flux->Draw("same");
    pub_eta_flux->Draw("same");
    */


    //////////////
    // RATIO between simulation / published
    
    auto *c2 =  new TCanvas();  
    double ratiopi0[5];
    double ratioeta[8];
    
    
    cout << "Plotting and outputing ratio of simulation/published " << endl;
    
    TGraph *ratio_pi0 = new TGraph();
    TGraph *ratio_eta = new TGraph();
    for ( int i = 0; i < sim_pi0_flux->GetN(); i++ ) {
      ratio_pi0->AddPoint(sim_pi0_flux->GetPointX(i),sim_pi0_flux->GetPointY(i)/pub_pi0_flux->GetPointY(i));
    }
    for ( int i = 0; i < sim_eta_flux->GetN(); i++ ) {
      ratio_eta->AddPoint(sim_eta_flux->GetPointX(i),sim_eta_flux->GetPointY(i)/pub_eta_flux->GetPointY(i));
    }
    
    c2->SetLogy();
    c2->SetLogx();

    ratio_pi0->SetMarkerStyle(kFullTriangleUp);
    ratio_eta->SetMarkerStyle(kFullSquare);
    ratio_pi0->SetMarkerSize(1.5);
    ratio_eta->SetMarkerSize(1.2);
    ratio_pi0->SetMarkerColor(kBlue+2);
    ratio_eta->SetMarkerColor(kOrange-2);

    
    double ratiomin = min(TMath::MinElement(ratio_eta->GetN(),ratio_eta->GetY()), TMath::MinElement(ratio_pi0->GetN(),ratio_pi0->GetY()));
    double ratiomax = max(TMath::MaxElement(ratio_eta->GetN(),ratio_eta->GetY()), TMath::MaxElement(ratio_pi0->GetN(),ratio_pi0->GetY()));
    
    if ( det == "argoneut" ) {
      ratio_eta->SetMaximum(40);
      ratio_eta->SetMinimum(0.0001);
    } else if ( det == "dune" ) {
      ratio_eta->SetMaximum(40);
      ratio_eta->SetMinimum(0.1);
    }
    ratio_eta->SetTitle("mCP flux comparison "+det);
  
    ratio_eta->Draw("A P");
    ratio_pi0->Draw("same P");
    
    ratio_eta->GetYaxis()->SetTitle("(mCP our simul)/("+det+" published)");
    ratio_eta->GetXaxis()->SetTitle("m_{#chi} (GeV)");
    c2->SaveAs("img/DifferenceAccepted_"+det+".png");



    // ratio plot using TRatioPlot (duh)
    auto *c3 = new TCanvas();
    //sim_pi0_flux->GetY()[i]/ypi0[i];
    //auto rp = new TRatioPlot(sim_pi0_flux->GetHistogram(),pub_pi0_flux->GetHistogram());
    //rp->Draw();
    
    // top pad (plot)
    // xlow, ylow, xhigh, yhigh
    

    
    
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1.0,1.0);
    pad1->SetBottomMargin(0.05); // upper and lower plot are joined
    pad1->SetGridx();
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    pad1->SetLogx();
    sim_eta_flux->Draw("A P");
    pub_eta_flux->Draw("same");
    
    // y axis
    if ( det != "t2k" ) {
      sim_eta_flux->GetYaxis()->SetTitle("# mCP accepted by "+det_formal);
    } else if ( det == "t2k" ) {
      sim_eta_flux->GetYaxis()->SetTitle("mCP Br for #epsilon^{2} = 1");
    }
    
    // top pad range
    if ( det == "duneOrnella" ) {
      sim_eta_flux->SetMinimum(1e8);
      sim_eta_flux->SetMaximum(1e17);
    } else if ( det == "t2k" ) {
      sim_eta_flux->SetMinimum(1e-7);
      sim_eta_flux->SetMaximum(1e-2);
    }
      
    // drawing
    //sim_eta_flux->Draw("same P");
    pub_pi0_flux->Draw("same");
    sim_pi0_flux->Draw("same P");
  

    // legend has to be after drawing
    TLegend *ratioleg;
    sim_eta_flux->SetTitle("#eta Manc. simulation");
    sim_pi0_flux->SetTitle("#pi^{0} Manc. simulation");
    if ( det == "dune" ) {
      pub_pi0_flux->SetTitle("#pi^{0} Phys. Rev. D 100, 015043");
      pub_eta_flux->SetTitle("#eta Phys. Rev. D 100, 015043");
    } else if ( det == "argoneut" ||
		det == "duneOrnella" ||
		det == "uboone" ) {
      pub_pi0_flux->SetTitle("#pi^{0} arXiv:1902.03246");
      pub_eta_flux->SetTitle("#eta arXiv:1902.03246");
    } else if ( det == "t2k" ) {
      pub_pi0_flux->SetTitle("#pi^{0} arXiv 2103.11814");
      pub_eta_flux->SetTitle("#eta arXiv 2103.11814");
    } 
    ratioleg = pad1->BuildLegend(0.6,0.7,0.9,0.9,"","");  
    
    // title after legend is built
    if ( det == "dune" ) {
      sim_eta_flux->SetTitle("Validation with FerMINI group: "+det_formal+" detector");
    } else if ( det == "argoneut" || det == "duneOrnella" ) {
      sim_eta_flux->SetTitle("Validation with ArgoNeuT group: "+det_formal+" detector");
    } else if ( det == "t2k" ) {
      sim_eta_flux->SetTitle("Validation with T2K group: Branching ratio only");
    } else if ( det == "uboone" ) {
      sim_eta_flux->SetTitle("Comparison of MicroBooNE simulation (10^{21} POT fhc) with ArgoNeuT group");
    }
  
    // first ratio (pi0)
    c3->cd();   // Go back to the main canvas before defining pad2!!!!
    TPad *pad2 = new TPad("pad2","pad3",0,0.05,1.0,0.32);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.3);
    pad2->SetGridx();
    if ( det == "dune" ) {
      pad2->SetGridy();
    }
    pad2->Draw();
    pad2->cd();
    
    // per det adjustment
    pad2->SetLogx();
    
    if ( det == "argoneut" ||
	 det == "uboone" ||
	 det == "naiveuboone" ||
	 det == "duneOrnella"
	 ) {
      pad2->SetLogy();
    } if ( det == "dune" ) {
      ratio_eta->SetMaximum(5.5);
      ratio_eta->SetMinimum(-0.5);
    } if ( det == "t2k" ) {
      ratio_eta->SetMaximum(1.5);
      ratio_eta->SetMinimum(0.2);
    }
    


    // drawing
    ratio_eta->Draw("A P");
    ratio_pi0->Draw("same P");
    ratio_eta->GetYaxis()->SetTitle("Simulation/published");
    ratio_eta->SetTitle("");
    
    // y axis
    sim_eta_flux->GetYaxis()->SetTitleSize(0.05);
    sim_eta_flux->GetYaxis()->SetTitleOffset(0.7);
    ratio_eta->GetYaxis()->SetTitleSize(0.09);
    ratio_eta->GetYaxis()->SetTitleOffset(0.4);
    ratio_eta->GetYaxis()->SetLabelSize(0.1);
    sim_eta_flux->GetYaxis()->CenterTitle();
    //ratio_eta->GetYaxis()->CenterTitle();
    
    // x axis
    sim_eta_flux->GetXaxis()->SetLabelSize(0.0);
    sim_eta_flux->GetXaxis()->SetTitleSize(0.0);
    ratio_eta->GetXaxis()->SetLabelSize(0.12);
    ratio_eta->GetXaxis()->SetTitleSize(0.14);
    ratio_eta->GetXaxis()->CenterTitle();
    
    //ratio_eta->G
    Double_t xrange = ratio_eta->GetHistogram()->GetBinLowEdge(101);
    TLine *line1 = new TLine(0,1,xrange,1);
    TLine *line10 = new TLine(0,10,xrange,10);
    TLine *line01 = new TLine(0,0.1,xrange,0.1);
    line1->SetLineStyle(7);
    line10->SetLineStyle(7);
    line01->SetLineStyle(7);

    // line drawing
    if ( det != "dune" ) {
      line1->Draw("same");
      if ( ratio_eta->GetMaximum() > 10 ) {
	line10->Draw("same");
      }
      if ( ratio_eta->GetMinimum() < 0.1 ) {
	line01->Draw("same");
      }
    }
    ratio_eta->Draw("same P"); // redraw points on top of lines
    ratio_pi0->Draw("same P"); // redraw points on top of lines
    //c3->SaveAs("img/2padDifferenceAccepted_"+det+"_.pdf");
    
    c3->SaveAs("img/PassingThroughDetector/Acc_Flux_"+det+"_s.png");
    c3->SaveAs("img/PassingThroughDetector/Acc_Flux_"+det+"_s.pdf");
  }
}
