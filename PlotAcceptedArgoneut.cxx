#include "AcceptedArgoneut.cxx"

void LatexText(Double_t x, Double_t y, int font, TString text)
{
  TLatex l2;
  l2.SetNDC();
  l2.SetTextFont(font);
  l2.DrawLatex(x,y,text);
}

void PlotAcceptedArgoneut(int WEIGHT = 0, bool DUNE = false)
{


  // These are the datapoints obtained by Friday 12 March 2021. If
  // script run with WEIGHT != 0, it calculates the y points again
  // using AcceptedArgoneut()
  const int NSimulPi0 = 5;
  const int NSimulEta = 8;
  
  Double_t xpi0decay[NSimulPi0]
    = { 0.01,
	0.02,
	0.03,
	0.05,
	0.06
  };

  Double_t ypi0decay[NSimulPi0]
    = { 4.5881e+10,
	5.99157e+10,
	6.40965e+10,
	7.49126e+10,
	5.43682e+10
  };

  Double_t xetadecay[NSimulEta]
    = { 0.01,
	0.02,
	0.03,
	0.05,
	0.06,
	0.10,
	0.20,
	0.25
  };

  Double_t yetadecay[NSimulEta]
    = { 3.36867e+09,
	1.74123e+09,
	6.59142e+09,
	9.67536e+09,
	8.32365e+09,
	7.53147e+09,
	8.5745e+09,
	8.05567e+09
  };
    
  TString weight_tstring;
  weight_tstring.Form("%i",abs(WEIGHT));


  if ( WEIGHT > 0 ) {
    cout << "Calculating simulation plot points" << endl;
    cout << "Will be stored in ./hist folder, run with negative number next time" << endl;
    cout << "e.g. root -l -b -q PlotAcceptedArgoneut.cxx(-"+weight_tstring+")" << endl;
    for ( int i = 0; i < NSimulPi0; i++ ) {
      TString fstr = Form("sim/mCP_q_0.010_m_%0.3f_fhc_pi0s.root",xpi0decay[i]);
      if ( !DUNE ) {
	ypi0decay[i] = AcceptedArgoneut(fstr,WEIGHT);
      } else if ( DUNE ) {
	ypi0decay[i] = AcceptedArgoneut(fstr,WEIGHT,true);
      }
    }
    for ( int i = 0; i < NSimulEta; i++ ) {
      TString fstr = Form("sim/mCP_q_0.010_m_%0.3f_fhc_etas.root",xetadecay[i]);
      if ( !DUNE ) {
	yetadecay[i] = AcceptedArgoneut(fstr,WEIGHT);
      } else if ( DUNE ) {
	yetadecay[i] = AcceptedArgoneut(fstr,WEIGHT,true);
      }
    }
  }

  TGraph *from_pi0decay = new TGraph(NSimulPi0,xpi0decay,ypi0decay);
  TGraph *from_etadecay = new TGraph(NSimulEta,xetadecay,yetadecay);
  from_pi0decay->SetName("pi0_decay");
  from_etadecay->SetName("eta_decay");
  
  TFile g;
  TString outfile_graph_decay;

  TString detector;
  if ( !DUNE ) {
    detector = "argoneut";
  } else if ( DUNE ) {
    detector = "dune";
  }
  outfile_graph_decay =
    "hist/"+detector+"_acceptance_decay_weight"+weight_tstring+".root";


  if ( WEIGHT > 0 ) {
    g.Open(outfile_graph_decay,"recreate");
    from_pi0decay->Write();
    from_etadecay->Write();
  }
  g.Close();

  ////////////////////////////////////////////////////////
  // data points obtained using WebPlotDigitizer
  // from fig. 2 (left) in arXiv:1902.03246v2
  // ArgoNeuT
  Double_t xpi0[5] = { 0.010088954277555864,
		       0.020130110267451286,
		       0.029986313485755686,
		       //0.039810717055349755,
		       0.04967682673934071,
		       0.0598305612629714,
		         };

  Double_t ypi0[5] = { 29053078518.935688,
		       29053078518.935688,
		       15592720260.810507,
		       //10733953891.080954,
		       7389202410.398621,
		       1765939991.098824,
  };

  Double_t xeta[8] = { 0.01000000000000001,
		       0.020309176209047368,
		       0.029986313485755686,
		       //0.039810717055349755,
		       0.05011872336272725,
		       0.0598305612629714,
		       0.1,
		       0.1995262314968879,
		       0.24897391368871477
  };
  
  Double_t yeta[8] = { 1559272026.0810509,
		       1559272026.0810509,
		       1559272026.0810509,
		       //1659391704.1669834,
		       1659391704.1669834,
		       1659391704.1669834,
		       1559272026.0810509,
		       836857701.5803154,
		       129372153.23092696
  };

  ////////////// DUNE data points
  // from https://journals.aps.org/prd/pdf/10.1103/PhysRevD.100.015043
  if ( DUNE ) {
    xpi0[0] = 0.010192923275448693; 
    xpi0[1] = 0.02027934585187628; 
    xpi0[2] = 0.030291989724257083; 
    xpi0[3] = 0.05026258415779048; 
    xpi0[4] = 0.05969440255395195; 
    // xpi0[5] = 0.06382329961306857;

    ypi0[0] = 2848035868435805;
    ypi0[1] = 1123324032978031.1;
    ypi0[2] = 403701725859656.6;
    ypi0[3] = 15556761439304.787;
    ypi0[4] = 453487850812.85913;
    //ypi0[5] = 17475284000.0769;
    
    xeta[0] = 0.010192923275448693; 
    xeta[1] = 0.02027934585187628; 
    xeta[2] = 0.0300039493211703; 
    xeta[3] = 0.05026258415779048; 
    xeta[4] = 0.06026747377166441; 
    xeta[5] = 0.1000000000000001; 
    xeta[6] = 0.2008651365410111; 
    xeta[7] = 0.25023050302971983;
      
    yeta[0] = 335160265093884.75;
    yeta[1] = 253536449397011.66;
    yeta[2] = 191791026167249.28;
    yeta[3] = 120450354025878.86;
    yeta[4] = 109749876549305.9;
    yeta[5] = 43287612810830.62;
    yeta[6] = 1047615752789.6663;
    yeta[7] = 5722367659.35022;
  }
  //sizeof(xpi0decay)/sizeof(xpi0decay[0])
  //int pi0n = sizeof(xpi0)/sizeof(xpi0[0]);
  //int etan = sizeof(xeta)/sizeof(xeta[0]);
  int pi0n = 5;
  int etan = 8;
  TGraph *from_pi0 = new TGraph(pi0n,xpi0,ypi0);
  TGraph *from_eta = new TGraph(etan,xeta,yeta);
  from_pi0->SetName("pi0_paper");
  from_eta->SetName("eta_paper");
  
  TFile h("hist/"+detector+"_acceptance_decay_paper.root","recreate");
  
  from_pi0->Write();
  from_eta->Write();

  h.Close();


  TFile f2;
  if ( WEIGHT < 0 ) {
    f2.Open("hist/"+detector+"_acceptance_decay_weight"+weight_tstring+".root");
    from_pi0decay = (TGraph*)gDirectory->Get("pi0_decay");
    from_etadecay = (TGraph*)gDirectory->Get("eta_decay");
  }

  // DUNE plots are N_accepted/epsilon^2, so we have to multiply by 1/epsilon^2
  if ( DUNE ) {
    for ( int i = 0; i < NSimulPi0; i++ ) {
      from_pi0decay->GetY()[i] *= 1e4;
    }
    for ( int i = 0; i < NSimulEta; i++ ) {
      from_etadecay->GetY()[i] *=1e4;
    }
  }
  
  cout << "PlotAcceptedArgoneut.cxx: " << endl;
  cout << "Ploting the following " << endl;
  for ( int i = 0; i < 5; i++ ) {
    cout << ypi0[i] << " " << ypi0decay[i] << endl;
  }
  cout << endl;
  for ( int i = 0; i < 8; i++ ) {
    cout << yeta[i] << " " << yetadecay[i] << endl;
  }  cout << endl;

  // PROTIP these have to be before the TCanvas creation
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  auto *c1 = new TCanvas();
  c1->SetLogy();
  c1->SetLogx();

  from_pi0->SetLineColor(kBlue);
  from_eta->SetLineColor(kOrange);
  from_pi0->SetLineWidth(5);
  from_eta->SetLineWidth(5);
  from_eta->Draw();
  from_eta->SetMinimum(1e-2);
  from_eta->SetMaximum(1e12);

  Double_t GlobalMaximum;
  GlobalMaximum = max(from_eta->GetMaximum(),from_pi0->GetMaximum());
  GlobalMaximum = max(GlobalMaximum,from_pi0decay->GetMaximum());
  GlobalMaximum = max(GlobalMaximum,from_etadecay->GetMaximum());
  GlobalMaximum *= 10; 
  from_eta->SetMaximum(GlobalMaximum);
  cout << GlobalMaximum << endl;


  if ( DUNE ) {
    from_eta->SetMinimum(1e8);
    //from_eta->SetMaximum(1e16);
    GlobalMaximum *= 100000; // to make the plot look better
    from_eta->SetMaximum(GlobalMaximum);
  }
  from_pi0->Draw("same");

  from_pi0decay->SetMarkerStyle(kFullTriangleUp);
  from_etadecay->SetMarkerStyle(kFullSquare);
  from_pi0decay->SetMarkerSize(1.5);
  from_etadecay->SetMarkerSize(1);
  from_pi0decay->SetMarkerColor(kBlue+5);
  from_etadecay->SetMarkerColor(kOrange-5);
  from_pi0decay->Draw("same P");
  from_etadecay->Draw("same P");
  
  TLegend *myleg;

  if ( !DUNE ) {
    from_pi0decay->SetTitle("#pi^{0} our simulation");
    from_etadecay->SetTitle("#eta our simulation");
    from_pi0->SetTitle("#pi^{0} arXiv 1902.03246v2");
    from_eta->SetTitle("#eta arXiv 1902.03246v2");
    
    from_eta->GetYaxis()->SetTitle("Number of mCPs geometrically accepted by ArgoNeut");
    from_eta->GetXaxis()->SetTitle("m_{#chi} (GeV)");
    
    myleg = c1->BuildLegend(0.5,0.15,0.9,0.45,"","");
    myleg->SetFillStyle(0);
    
    LatexText(0.2,0.3,42,"#epsilon = 0.01");
    LatexText(0.2,0.2,42,"10^{20} POT");
  } else if ( DUNE ) {
    from_pi0decay->SetTitle("#pi^{0} our simulation");
    from_etadecay->SetTitle("#eta our simulation");
    from_pi0->SetTitle("#pi^{0} Phys. Rev. D 100, 015043");
    from_eta->SetTitle("#eta Phys. Rev. D 100, 015043");
    
    from_eta->GetYaxis()->SetTitle("N_{mCP} /#epsilon^{2} at 10^{21}POT");
    from_eta->GetXaxis()->SetTitle("m_{#chi} (GeV)");
    
    myleg = c1->BuildLegend(0.15,0.15,0.5,0.45,"","");
    myleg->SetFillStyle(0);
    
    LatexText(0.2,0.85,42,"\"DUNE ND\", 10^{21} POT, 574m, 1m x 1m detector");
  }
  //gStyle->SetOptTitle(0);

  if ( abs(WEIGHT) == 1 ) {
    from_eta->SetTitle("No weights");
  } else if ( abs(WEIGHT) == 2 ) {
    from_eta->SetTitle("Meson simulation weights");
  } else if ( abs(WEIGHT) == 3 ) {
    from_eta->SetTitle("Meson simulation and normalized TGenPhaseSpace weights");
  } else if ( abs(WEIGHT) == 4 ) {
    from_eta->SetTitle("Differential Branching fraction from arXiv 2010.07941v1");
  } else if ( abs(WEIGHT) == 5 ) {
    from_eta->SetTitle("Branching fraction from Zhen Liu email");
  }

  c1->SaveAs("img/Accepted_"+detector+"_"+weight_tstring+".png");

  //////////////
  // RATIO between simulation / published

  auto *c2 =  new TCanvas();  
  double ratiopi0[5];
  double ratioeta[8];

  
  cout << "Plotting and outputing ratio of simulation/published " << endl;
  for ( int i = 0; i < 5; i++ ) {
    ratiopi0[i] = from_pi0decay->GetY()[i]/ypi0[i];
    cout << ypi0[i] << " " << from_pi0decay->GetY()[i] << " " << ratiopi0[i];
      cout << endl;
  }
  cout << endl;
  for ( int i = 0; i < 8; i++ ) {
    ratioeta[i] = from_etadecay->GetY()[i]/yeta[i];
    cout << yeta[i] << " " << from_etadecay->GetY()[i] << " " << ratioeta[i];
      cout << endl;
  }  cout << endl;
  TGraph *ratio_pi0 = new TGraph(5,xpi0decay,ratiopi0);
  TGraph *ratio_eta = new TGraph(8,xetadecay,ratioeta);

  c2->SetLogy();
  c2->SetLogx();

  ratio_pi0->SetMarkerStyle(kFullTriangleUp);
  ratio_eta->SetMarkerStyle(kFullTriangleUp);
  ratio_pi0->SetMarkerSize(1.5);
  ratio_eta->SetMarkerSize(1.5);
  ratio_pi0->SetMarkerColor(kBlue+2);
  ratio_eta->SetMarkerColor(kOrange-2);
  
  double ratiomin = min(ratio_eta->GetMinimum(), ratio_pi0->GetMinimum());
  double ratiomax = max(ratio_eta->GetMaximum(), ratio_pi0->GetMaximum());
  cout << "etamin, etamax " <<  ratio_eta->GetMinimum()<< " " << ratio_eta->GetMaximum() << endl;
  cout << "pi0min, pi0max " <<  ratio_pi0->GetMinimum()<< " " << ratio_pi0->GetMaximum() << endl;
  cout << "ratiomin, ratiomax " << ratiomin << " " << ratiomax << endl;

  if ( detector == "argoneut" ) {
    ratio_eta->SetMaximum(40);
    ratio_eta->SetMinimum(0.0001);
  } else if ( detector == "dune" ) {
    ratio_eta->SetMaximum(40);
    ratio_eta->SetMinimum(0.1);
  }
  ratio_eta->SetTitle("mCP flux comparison "+detector);
  
  ratio_eta->Draw("A P");
  ratio_pi0->Draw("same P");
    
  ratio_eta->GetYaxis()->SetTitle("(mCP our simul)/("+detector+" published)");
  ratio_eta->GetXaxis()->SetTitle("m_{#chi} (GeV)");
  c2->SaveAs("img/DifferenceAccepted_"+detector+"_"+weight_tstring+".png");



  /*
** pi0
0.020130110267451286, 29053078518.935688
0.010088954277555864, 29053078518.935688
0.04967682673934071, 7389202410.398621
0.0598305612629714, 1765939991.098824
0.029986313485755686, 15592720260.810507
0.039810717055349755, 10733953891.080954
*** eta
0.01000000000000001, 1559272026.0810509
0.020309176209047368, 1559272026.0810509
0.029986313485755686, 1559272026.0810509
0.039810717055349755, 1659391704.1669834
0.0598305612629714, 1659391704.1669834
0.05011872336272725, 1659391704.1669834
0.1, 1559272026.0810509
0.1995262314968879, 836857701.5803154
0.24897391368871477, 129372153.23092696

   */
}
