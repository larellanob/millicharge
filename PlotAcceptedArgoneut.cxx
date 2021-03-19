#include "AcceptedArgoneut.cxx"

void LatexText(Double_t x, Double_t y, int font, TString text)
{
  TLatex l2;
  l2.SetNDC();
  l2.SetTextFont(font);
  l2.DrawLatex(x,y,text);
}

void PlotAcceptedArgoneut(int WEIGHT = 0)
{


  // These are the datapoints obtained by Friday 12 March 2021. If
  // script run with WEIGHT != 0, it calculates the y points again
  // using AcceptedArgoneut()
  Double_t xpi0decay[5] = { 0.01,
			    0.02,
			    0.03,
			    0.05,
			    0.06
  };

  Double_t xetadecay[8] = { 0.01,
			    0.02,
			    0.03,
			    0.05,
			    0.06,
			    0.10,
			    0.20,
			    0.25
  };

  Double_t ypi0decay[5] = { 4.5881e+10,
			    5.99157e+10,
			    6.40965e+10,
			    7.49126e+10,
			    5.43682e+10
  };

  Double_t yetadecay[8] = { 3.36867e+09,
			    1.74123e+09,
			    6.59142e+09,
			    9.67536e+09,
			    8.32365e+09,
			    7.53147e+09,
			    8.5745e+09,
			    8.05567e+09
  };
    

  if ( WEIGHT == 1 || WEIGHT == 2 || WEIGHT == 3 ) {
    for ( int i = 0; i < sizeof(xpi0decay)/sizeof(xpi0decay[0]); i++ ) {
      TString fstr = Form("sim/mCP_q_0.010_m_%0.3f_fhc_pi0s.root",xpi0decay[i]);
      ypi0decay[i] = AcceptedArgoneut(fstr,WEIGHT);
    }
    for ( int i = 0; i < sizeof(xetadecay)/sizeof(xetadecay[0]); i++ ) {
      TString fstr = Form("sim/mCP_q_0.010_m_%0.3f_fhc_etas.root",xetadecay[i]);
      yetadecay[i] = AcceptedArgoneut(fstr,WEIGHT);
    }
  }

  TGraph *from_pi0decay = new TGraph(5,xpi0decay,ypi0decay);
  TGraph *from_etadecay = new TGraph(8,xetadecay,yetadecay);
  from_pi0decay->SetName("pi0_decay");
  from_etadecay->SetName("eta_decay");
  
  TFile g;
  TString outfile_graph_decay;
  if ( WEIGHT == 1 ) {
    outfile_graph_decay =
      "hist/argoneut_acceptance_decay_weight1.root";
  } else if ( WEIGHT == 2 ) {
    outfile_graph_decay =
      "hist/argoneut_acceptance_decay_weight2.root";
  } else if ( WEIGHT == 3 ) {
    outfile_graph_decay =
      "hist/argoneut_acceptance_decay_weight3.root";
  }


  if ( WEIGHT == 1 || WEIGHT == 2 || WEIGHT == 3 ) {
    g.Open(outfile_graph_decay,"recreate");
    from_pi0decay->Write();
    from_etadecay->Write();
  }
  g.Close();

  ////////////////////////////////////////////////////////
  // data points obtained using WebPlotDigitizer
  // from fig. 2 (left) in arXiv:1902.03246v2
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

  
  TGraph *from_pi0 = new TGraph(5,xpi0,ypi0);
  TGraph *from_eta = new TGraph(8,xeta,yeta);
  from_pi0->SetName("pi0_paper");
  from_eta->SetName("eta_paper");
  
  TFile h("hist/argoneut_acceptance_decay_paper.root","recreate");

  from_pi0->Write();
  from_eta->Write();

  h.Close();


  TFile f2;
  if ( WEIGHT == -1 ) {
    f2.Open("hist/argoneut_acceptance_decay_weight1.root");
    from_pi0decay = (TGraph*)gDirectory->Get("pi0_decay");
    from_etadecay = (TGraph*)gDirectory->Get("eta_decay");
  } else if ( WEIGHT == -2 ) {
    f2.Open("hist/argoneut_acceptance_decay_weight2.root");
    from_pi0decay = (TGraph*)gDirectory->Get("pi0_decay");
    from_etadecay = (TGraph*)gDirectory->Get("eta_decay");
  } else if ( WEIGHT == -3 ) {
    f2.Open("hist/argoneut_acceptance_decay_weight3.root");
    from_pi0decay = (TGraph*)gDirectory->Get("pi0_decay");
    from_etadecay = (TGraph*)gDirectory->Get("eta_decay");
  }
      
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
  from_pi0->Draw("same");

  from_pi0decay->SetMarkerStyle(kFullTriangleUp);
  from_etadecay->SetMarkerStyle(kFullTriangleUp);
  from_pi0decay->SetMarkerColor(kBlue+2);
  from_etadecay->SetMarkerColor(kOrange-2);
  from_pi0decay->Draw("same P");
  from_etadecay->Draw("same P");
  
  TLegend *myleg;
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

  gStyle->SetOptTitle(0);


  TString weight_tstring;
  weight_tstring.Form("%i",abs(WEIGHT));
  c1->SaveAs("img/AcceptedArgoneut_"+weight_tstring+".png");


  //////////////
  
  auto *c2 =  new TCanvas();  
  double ratiopi0[5];
  double ratioeta[8];
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
  
  ratio_eta->Draw("A P");
  ratio_pi0->Draw("same P");
  ratio_eta->GetYaxis()->SetTitle("(mCP our simul)/(mCP published)");
  ratio_eta->GetXaxis()->SetTitle("m_{#chi} (GeV)");
  c2->SaveAs("img/DifferenceAcceptedArgoneut_"+weight_tstring+".png");
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
