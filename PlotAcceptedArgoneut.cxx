#include "AcceptedArgoneut.cxx"

void LatexText(Double_t x, Double_t y, int font, TString text)
{
  TLatex l2;
  l2.SetNDC();
  l2.SetTextFont(font);
  l2.DrawLatex(x,y,text);
}

void PlotAcceptedArgoneut(int WEIGHT = 5, TString detector = "uboone")
{

  TString detector_formal;
  if ( detector == "uboone" || detector == "naiveuboone" ) {
    detector_formal = "MicroBooNE";
  } else if ( detector == "dune" || detector == "duneOrnella" ) {
    detector_formal = "DUNE";
  } else if ( detector == "argoneut" ) {
    detector_formal = "ArgoNeuT";
  } else if ( detector == "t2k" ) {
    detector_formal = "T2K";
    WEIGHT = 6;
  } else {
    detector_formal = detector;
  }
    
  
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
      ypi0decay[i] = AcceptedArgoneut(fstr,WEIGHT,detector);
    }
    for ( int i = 0; i < NSimulEta; i++ ) {
      TString fstr = Form("sim/mCP_q_0.010_m_%0.3f_fhc_etas.root",xetadecay[i]);
      yetadecay[i] = AcceptedArgoneut(fstr,WEIGHT,detector);
    }
  }

  TGraph *from_pi0decay = new TGraph(NSimulPi0,xpi0decay,ypi0decay);
  TGraph *from_etadecay = new TGraph(NSimulEta,xetadecay,yetadecay);
  from_pi0decay->SetName("pi0_decay");
  from_etadecay->SetName("eta_decay");
  
  TFile g;
  TString outfile_graph_decay;

  /*
  TString detector;
  if ( !DUNE ) {
    detector = "argoneut";
  } else if ( DUNE ) {
    detector = "dune";
  }
  */
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
  /*
  Double_t xpi0[5] = { 0.010088954277555864,
		       0.020130110267451286,
		       0.029986313485755686,
		       //0.039810717055349755,
		       0.04967682673934071,
		       0.0598305612629714,
		         };
  */

  Double_t xpi0[5] = { 0.01,
		       0.02,
		       0.03,
		       //0.039810717055349755,
		       0.05,
		       0.06,
  };

  Double_t ypi0[5] = { 29053078518.935688,
		       29053078518.935688,
		       15592720260.810507,
		       //10733953891.080954,
		       7389202410.398621,
		       1765939991.098824,
  };
  /*
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
  */

  Double_t xeta[8] = { 0.01,
		       0.02,
		       0.03,
		       //0.039810717055349755,
		       0.05,
		       0.06,
		       0.1,
		       0.2,
		       0.25
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
  if ( detector == "dune" ) {
    /*
    xpi0[0] = 0.010192923275448693; 
    xpi0[1] = 0.02027934585187628; 
    xpi0[2] = 0.030291989724257083; 
    xpi0[3] = 0.05026258415779048; 
    xpi0[4] = 0.05969440255395195; 
    */
    // xpi0[5] = 0.06382329961306857;


    // re extraction 05-apr-21
    ypi0[0] = 4245286227393835;
    ypi0[1] = 1654263361624045.5;
    ypi0[2] = 543104899786251.3;
    ypi0[3] = 17640529950316.74;
    ypi0[4] = 288712186834.80786;
    
    /*
    ypi0[0] = 2848035868435805;
    ypi0[1] = 1123324032978031.1;
    ypi0[2] = 403701725859656.6;
    ypi0[3] = 15556761439304.787;
    ypi0[4] = 453487850812.85913;
    */
    //ypi0[5] = 17475284000.0769;
    /*
    xeta[0] = 0.010192923275448693; 
    xeta[1] = 0.02027934585187628; 
    xeta[2] = 0.0300039493211703; 
    xeta[3] = 0.05026258415779048; 
    xeta[4] = 0.06026747377166441; 
    xeta[5] = 0.1000000000000001; 
    xeta[6] = 0.2008651365410111; 
    xeta[7] = 0.25023050302971983;
    */

    // re extraction 05-apr-21
    yeta[0] = 543104899786251.3;
    yeta[1] = 420006308884957.9;
    yeta[2] = 298138872452931.3;
    yeta[3] = 211632028822338.97;
    yeta[4] = 163664123270640.62;
    yeta[5] = 69480104843934.05;
    yeta[6] = 1349685425291.2366;
    yeta[7] = 7900860950.738926;
    
    

    /*
    yeta[0] = 335160265093884.75;
    yeta[1] = 253536449397011.66;
    yeta[2] = 191791026167249.28;
    yeta[3] = 120450354025878.86;
    yeta[4] = 109749876549305.9;
    yeta[5] = 43287612810830.62;
    yeta[6] = 1047615752789.6663;
    yeta[7] = 5722367659.35022;
    */
  }

  if ( detector == "duneOrnella" ) {
    ypi0[0] = 2106344542324141.5;
    ypi0[1] = 1968419447286631.5;
    ypi0[2] = 1070068955693183.4;
    ypi0[3] = 508021804691307.44;
    ypi0[4] = 114504756993829.06;

    yeta[0] = 114504756993829.06;
    yeta[1] = 122527985738287.39;
    yeta[2] = 122527985738287.39;
    yeta[3] = 122527985738287.39;
    yeta[4] = 114504756993829.06;
    yeta[5] = 114504756993829.06;
    yeta[6] = 50802180469130.75;
    yeta[7] = 8161400793251.868;
  }

  if ( detector == "t2k" ) {
    ypi0[0] = 0.0031207061265801687;
    ypi0[1] = 0.0013029943385750237;
    ypi0[2] = 0.0005159928433650851;
    ypi0[3] = 0.000026624861440398474;
    ypi0[4] = 9.484435649458936e-7;

    yeta[0] = 0.002959813643842186;
    yeta[1] = 0.0020981623055342426;
    yeta[2] = 0.0016102620275609408;
    yeta[3] = 0.0010543589908346815;
    yeta[4] = 0.0008760496806274342;
    yeta[5] = 0.0003657786820891421;
    yeta[6] = 0.000014485034710121493;
    yeta[7] = 1.4873521072935117e-7;
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
  if ( detector == "dune" ) {
    for ( int i = 0; i < NSimulPi0; i++ ) {
      from_pi0decay->GetY()[i] *= 1e4;
    }
    for ( int i = 0; i < NSimulEta; i++ ) {
      from_etadecay->GetY()[i] *=1e4;
    }
  }
  //cout << scientific;
  cout << "PlotAcceptedArgoneut.cxx: \n" << endl;
  cout << "Ploting the following " << endl;
  cout << "pi0" << endl;
  cout << "mCP mass\tpublished\tour simulation" << endl;
  cout << "----------------------------------------------" << endl;
  for ( int i = 0; i < 5; i++ ) {
    cout << Form("%8.3f\t%e\t%.3f",xpi0[i], ypi0[i], ypi0decay[i]) << endl;
  }
  cout << endl;
  cout << "eta" << endl;
  cout << "mCP mass\tpublished\tour simulation" << endl;
  cout << "----------------------------------------------" << endl;
  for ( int i = 0; i < 8; i++ ) {
    cout << Form("%8.3f\t%e\t%.3f",xeta[i], yeta[i], yetadecay[i]) << endl;
    //cout << xeta[i] << "\t" << yeta[i] << "\t" << yetadecay[i] << endl;
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
  cout << "Global Maximum: " << GlobalMaximum << endl;


  if ( detector == "dune" ) {
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
  from_pi0decay->SetMarkerColor(kBlue+2);
  from_etadecay->SetMarkerColor(kOrange-2);
  from_pi0decay->Draw("same P");
  from_etadecay->Draw("same P");
  
  TLegend *myleg;

  if ( detector != "dune" ) {
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
  } else if ( detector == "dune" ) {
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
    if ( detector == "dune" || detector == "duneOrnella" ) {
      from_eta->SetTitle(detector_formal+" points, DUNE lines");
    } else if ( detector == "uboone" || detector == "argoneut" ) {
      from_eta->SetTitle(detector_formal+" points (simulation), ArgoNeuT lines (published)");
    }
  } else if ( abs(WEIGHT) == 5 && detector == "uboone" ) {
    cout << "what" << endl;
  }

  c1->SaveAs("img/Accepted_"+detector+"_"+weight_tstring+".png");
  cout << endl;
  //////////////
  // RATIO between simulation / published

  auto *c2 =  new TCanvas();  
  double ratiopi0[5];
  double ratioeta[8];

  
  cout << "Plotting and outputing ratio of simulation/published " << endl;

  cout << "pi0" << endl;
  cout << "mCP mass\tpublished\tour simulation\tratio" << endl;
  cout << "----------------------------------------------------------" << endl;
  for ( int i = 0; i < 5; i++ ) {
    ratiopi0[i] = from_pi0decay->GetY()[i]/ypi0[i];
    cout << Form("%8.3f\t%e\t%10.3f\t%.3f",xpi0[i], ypi0[i], from_pi0decay->GetY()[i], ratiopi0[i]) << endl;
    //cout << ypi0[i] << " " << from_pi0decay->GetY()[i] << " " << ratiopi0[i];
  }
  cout << "eta" << endl;
  cout << "mCP mass\tpublished\tour simulation\tratio" << endl;
  cout << "----------------------------------------------------------" << endl;
  for ( int i = 0; i < 8; i++ ) {
    ratioeta[i] = from_etadecay->GetY()[i]/yeta[i];
    cout << Form("%8.3f\t%e\t%10.3f\t%.3f",xeta[i], yeta[i], from_etadecay->GetY()[i], ratioeta[i]) << endl;
    //cout << yeta[i] << " " << from_etadecay->GetY()[i] << " " << ratioeta[i];
  }  cout << endl;
  TGraph *ratio_pi0 = new TGraph(5,xpi0decay,ratiopi0);
  TGraph *ratio_eta = new TGraph(8,xetadecay,ratioeta);

  c2->SetLogy();
  c2->SetLogx();

  ratio_pi0->SetMarkerStyle(kFullTriangleUp);
  ratio_eta->SetMarkerStyle(kFullSquare);
  ratio_pi0->SetMarkerSize(1.5);
  ratio_eta->SetMarkerSize(1.2);
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



  // ratio plot using TRatioPlot (duh)
  auto *c3 = new TCanvas();
    //from_pi0decay->GetY()[i]/ypi0[i];
  //auto rp = new TRatioPlot(from_pi0decay->GetHistogram(),from_pi0->GetHistogram());
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
  from_eta->Draw();
  //from_eta->SetMinimum(1e5);


  // y axis
  if ( detector != "t2k" ) {
    from_eta->GetYaxis()->SetTitle("# mCP accepted by "+detector_formal);
  } else if ( detector == "t2k" ) {
    from_eta->GetYaxis()->SetTitle("mCP Br for #epsilon^{2} = 1");
  }
  
  // top pad range
  if ( detector == "duneOrnella" ) {
    from_eta->SetMinimum(1e8);
    from_eta->SetMaximum(1e17);
  } else if ( detector == "t2k" ) {
    from_eta->SetMinimum(1e-7);
    from_eta->SetMaximum(1e-2);
  }

  // drawing
  from_etadecay->Draw("same P");
  from_pi0->Draw("same");
  from_pi0decay->Draw("same P");


  // legend has to be after drawing
  TLegend *ratioleg;
  from_etadecay->SetTitle("#eta Manc. simulation");
  from_pi0decay->SetTitle("#pi^{0} Manc. simulation");
  if ( detector == "dune" ) {
    from_pi0->SetTitle("#pi^{0} Phys. Rev. D 100, 015043");
    from_eta->SetTitle("#eta Phys. Rev. D 100, 015043");
  } else if ( detector == "argoneut" ||
	      detector == "duneOrnella" ||
	      detector == "uboone" ) {
    from_pi0->SetTitle("#pi^{0} arXiv:1902.03246");
    from_eta->SetTitle("#eta arXiv:1902.03246");
  } else if ( detector == "t2k" ) {
    from_pi0->SetTitle("#pi^{0} arXiv 2103.11814");
    from_eta->SetTitle("#eta arXiv 2103.11814");
  } 
  ratioleg = pad1->BuildLegend(0.6,0.7,0.9,0.9,"","");  

  // title after legend is built
  if ( detector == "dune" ) {
    from_eta->SetTitle("Validation with FerMINI group: "+detector_formal+" detector");
  } else if ( detector == "argoneut" || detector == "duneOrnella" ) {
    from_eta->SetTitle("Validation with ArgoNeuT group: "+detector_formal+" detector");
  } else if ( detector == "t2k" ) {
    from_eta->SetTitle("Validation with T2K group: Branching ratio only");
  } else if ( detector == "uboone" ) {
    from_eta->SetTitle("Comparison of MicroBooNE simulation (10^{21} POT fhc) with ArgoNeuT group");
  }
  
  // first ratio (pi0)
  c3->cd();   // Go back to the main canvas before defining pad2!!!!
  TPad *pad2 = new TPad("pad2","pad3",0,0.05,1.0,0.32);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.3);
  pad2->SetGridx();
  if ( detector == "dune" ) {
    pad2->SetGridy();
  }
  pad2->Draw();
  pad2->cd();

  // per detector adjustment
  pad2->SetLogx();

  if ( detector == "argoneut" ||
       detector == "uboone" ||
       detector == "naiveuboone" ||
       detector == "duneOrnella"
       ) {
    pad2->SetLogy();
  } if ( detector == "dune" ) {
    ratio_eta->SetMaximum(5.5);
    ratio_eta->SetMinimum(-0.5);
  } if ( detector == "t2k" ) {
    ratio_eta->SetMaximum(1.5);
    ratio_eta->SetMinimum(0.2);
  }
      


  // drawing
  ratio_eta->Draw("A P");
  ratio_pi0->Draw("same P");
  ratio_eta->GetYaxis()->SetTitle("Simulation/published");
  ratio_eta->SetTitle("");

  // y axis
  from_eta->GetYaxis()->SetTitleSize(0.05);
  from_eta->GetYaxis()->SetTitleOffset(0.7);
  ratio_eta->GetYaxis()->SetTitleSize(0.09);
  ratio_eta->GetYaxis()->SetTitleOffset(0.4);
  ratio_eta->GetYaxis()->SetLabelSize(0.1);
  from_eta->GetYaxis()->CenterTitle();
  //ratio_eta->GetYaxis()->CenterTitle();

  // x axis
  from_eta->GetXaxis()->SetLabelSize(0.0);
  from_eta->GetXaxis()->SetTitleSize(0.0);
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
  if ( detector != "dune" ) {
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

  // save
  c3->SaveAs("img/2padDifferenceAccepted_"+detector+"_"+weight_tstring+".pdf");


  
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
