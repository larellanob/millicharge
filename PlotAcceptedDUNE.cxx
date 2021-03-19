#include "AcceptedDUNE.cxx"

void LatexText(Double_t x, Double_t y, int font, TString text)
{
  TLatex l2;
  l2.SetNDC();
  l2.SetTextFont(font);
  l2.DrawLatex(x,y,text);
}


void PlotAcceptedDUNE(int WEIGHT = 0)
{

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
    


  if ( WEIGHT == 1 || WEIGHT == 2 || WEIGHT == 3 || WEIGHT == 4 ) {
    for ( int i = 0; i < sizeof(xpi0decay)/sizeof(xpi0decay[0]); i++ ) {
      TString fstr = Form("sim/mCP_q_0.010_m_%0.3f_fhc_pi0s.root",xpi0decay[i]);
      ypi0decay[i] = AcceptedDUNE(fstr,WEIGHT);
    }
    for ( int i = 0; i < sizeof(xetadecay)/sizeof(xetadecay[0]); i++ ) {
      TString fstr = Form("sim/mCP_q_0.010_m_%0.3f_fhc_etas.root",xetadecay[i]);
      yetadecay[i] = AcceptedDUNE(fstr,WEIGHT);
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
      "hist/dune_acceptance_decay_weight1.root";
  } else if ( WEIGHT == 2 ) {
    outfile_graph_decay =
      "hist/dune_acceptance_decay_weight2.root";
  } else if ( WEIGHT == 3 ) {
    outfile_graph_decay =
      "hist/dune_acceptance_decay_weight3.root";
  } else if ( WEIGHT == 4 ) {
    outfile_graph_decay =
      "hist/dune_acceptance_decay_weight4.root";
  }


  if ( WEIGHT == 1 || WEIGHT == 2 || WEIGHT == 3 || WEIGHT == 4 ) {
    g.Open(outfile_graph_decay,"recreate");
    from_pi0decay->Write();
    from_etadecay->Write();
  }
  g.Close();
  
  ////////////////////////////////////////////////////////
  // data points obtained using WebPlotDigitizer
  // https://journals.aps.org/prd/pdf/10.1103/PhysRevD.100.015043
  
  Double_t xpi0[6]
    = { 
       0.010192923275448693, 
       0.02027934585187628, 
       0.030291989724257083, 
       0.05026258415779048, 
       0.05969440255395195, 
       0.06382329961306857,
  };




  Double_t ypi0[6]
    = { 
       2848035868435805,
       1123324032978031.1,
       403701725859656.6,
       15556761439304.787,
       453487850812.85913,
       17475284000.0769,
  };
  
  Double_t xeta[8]
    = {
       0.010192923275448693, 
       0.02027934585187628, 
       0.0300039493211703, 
       0.05026258415779048, 
       0.06026747377166441, 
       0.1000000000000001, 
       0.2008651365410111, 
       0.25023050302971983
  };
  
  Double_t yeta[8]
    = {
       335160265093884.75,
       253536449397011.66,
       191791026167249.28,
       120450354025878.86,
       109749876549305.9,
       43287612810830.62,
       1047615752789.6663,
       5722367659.35022
  };

  /* eta
0.010192923275448693, 335160265093884.75
0.02027934585187628, 253536449397011.66
0.0300039493211703, 191791026167249.28
0.05026258415779048, 120450354025878.86
0.06026747377166441, 109749876549305.9
0.1000000000000001, 43287612810830.62
0.2008651365410111, 1047615752789.6663
0.25023050302971983, 5722367659.35022

  */
  
  /* pi0
0.010192923275448693, 2848035868435805
0.02027934585187628, 1123324032978031.1
0.030291989724257083, 403701725859656.6
0.05026258415779048, 15556761439304.787
0.05969440255395195, 453487850812.85913
0.06382329961306857, 17475284000.0769

   */
  
  TGraph *from_pi0 = new TGraph(6,xpi0,ypi0);
  TGraph *from_eta = new TGraph(8,xeta,yeta);
  from_pi0->SetName("pi0_paper");
  from_eta->SetName("eta_paper");
  
  TFile h("hist/dune_acceptance_decay_paper.root","recreate");

  from_pi0->Write();
  from_eta->Write();

  h.Close();


  TFile f2;
  if ( WEIGHT == -1 ) {
    f2.Open("hist/dune_acceptance_decay_weight1.root");
    from_pi0decay = (TGraph*)gDirectory->Get("pi0_decay");
    from_etadecay = (TGraph*)gDirectory->Get("eta_decay");
  } else if ( WEIGHT == -2 ) {
    f2.Open("hist/dune_acceptance_decay_weight2.root");
    from_pi0decay = (TGraph*)gDirectory->Get("pi0_decay");
    from_etadecay = (TGraph*)gDirectory->Get("eta_decay");
  } else if ( WEIGHT == -3 ) {
    f2.Open("hist/dune_acceptance_decay_weight3.root");
    from_pi0decay = (TGraph*)gDirectory->Get("pi0_decay");
    from_etadecay = (TGraph*)gDirectory->Get("eta_decay");
  } else if ( WEIGHT == -4 ) {
    f2.Open("hist/dune_acceptance_decay_weight4.root");
    from_pi0decay = (TGraph*)gDirectory->Get("pi0_decay");
    from_etadecay = (TGraph*)gDirectory->Get("eta_decay");
  }
  if ( WEIGHT == -1 || WEIGHT == -2 || WEIGHT == -3 || WEIGHT == -4 ) {
    for ( int i = 0; i < 5; i++ ) {
      from_pi0decay->GetY()[i] *= 1e4;
    }
    for ( int i = 0; i < 8; i++ ) {
      from_etadecay->GetY()[i] *=1e4;
    }
  }
  
  auto *c1 = new TCanvas();
  c1->SetLogy();
  c1->SetLogx();

  from_pi0->SetLineColor(kBlue);
  from_eta->SetLineColor(kOrange);
  from_pi0->SetLineWidth(5);
  from_eta->SetLineWidth(5);
  
  from_eta->Draw();
  
  from_eta->SetMinimum(1e8);
  from_eta->SetMaximum(1e16);
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
  from_pi0decay->SetTitle("#pi^{0} our simulation");
  from_etadecay->SetTitle("#eta our simulation");
  from_pi0->SetTitle("#pi^{0} Phys. Rev. D 100, 015043");
  from_eta->SetTitle("#eta Phys. Rev. D 100, 015043");

  from_eta->GetYaxis()->SetTitle("N_{mCP} /#epsilon^{2} at 10^{21}POT");
  from_eta->GetXaxis()->SetTitle("m_{#chi} (GeV)");
  
  myleg = c1->BuildLegend(0.15,0.15,0.5,0.45,"","");
  myleg->SetFillStyle(0);

  LatexText(0.2,0.85,42,"\"DUNE ND\", 10^{21} POT, 574m, 1m x 1m detector");

  if ( WEIGHT == -1 ) {
    //LatexText(0.2,0.3,42,"No weights");
    from_eta->SetTitle("No weights");
  } else if ( WEIGHT == -2 ) {
    //LatexText(0.2,0.3,42,"Meson simulation weights");
    from_eta->SetTitle("Meson simulation weights");
  } else if ( WEIGHT == -3 ) {
    //LatexText(0.2,0.3,42,"Meson simulation and TGenPhaseSpace weights");
    from_eta->SetTitle("Meson simulation and TGenPhaseSpace weights");
  } else if ( WEIGHT == -4 ) {
    //LatexText(0.2,0.3,42,"Meson simulation and TGenPhaseSpace weights");
    from_eta->SetTitle("Meson simulation and normalized TGenPhaseSpace weights");
  }


  
  //LatexText(0.2,0.2,42,"10^{21} POT");
  
  //gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TString weight_tstring;
  weight_tstring.Form("%i",abs(WEIGHT));
  c1->SaveAs("img/AcceptedDUNE_"+weight_tstring+".png");
    
}
