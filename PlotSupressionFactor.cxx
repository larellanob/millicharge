// Plots the phase space suppression factor

#include "Root/CrossSection.cxx"
#include "Root/CreateEmptyDirs.cxx"

void PlotSupressionFactor()
{
  CreateEmptyDirs();
  std::vector<double> masspoints_pi0 =
    {
     0.001,
     0.002,
     0.003,
     0.005,
     0.006
    };
  std::vector<double> masspoints_eta =
    {
     0.001,
     0.002,
     0.003,
     0.005,
     0.006,
     0.010,
     0.020,
     0.025
    };

  double charge = 0.01; // epsilon factor
  
  TGraph *zhenliu_pi0 = new TGraph(masspoints_pi0.size());
  TGraph *physrevd_pi0 = new TGraph(masspoints_pi0.size());
  TGraph *t2k_pi0 = new TGraph(masspoints_pi0.size());
  TGraph *naive_pi0 = new TGraph(masspoints_pi0.size());

  TGraph *zhenliu_eta = new TGraph(masspoints_eta.size());
  TGraph *physrevd_eta = new TGraph(masspoints_eta.size());
  TGraph *t2k_eta = new TGraph(masspoints_eta.size());
  TGraph *naive_eta = new TGraph(masspoints_eta.size());


  TDatabasePDG pdgPlotCrossSection;
  double mass_pi0 = pdgPlotCrossSection.GetParticle(111)->Mass();
  double mass_eta = pdgPlotCrossSection.GetParticle(221)->Mass();


  double global_minimum_pi0 = 100;
  double global_maximum_pi0 = -100;
  double global_minimum_eta = 100;
  double global_maximum_eta = -100;
  double xsec;
  
  for ( int i = 0; i < masspoints_pi0.size(); i++ ) {
    xsec = CrossSection(masspoints_pi0[i],mass_pi0,"zhenliu");
    zhenliu_pi0->SetPoint(i,masspoints_pi0[i],xsec);
    if ( global_minimum_pi0 > xsec ) {
      global_minimum_pi0 = xsec;
    }
    if ( global_maximum_pi0 < xsec ) {
      global_maximum_pi0 = xsec;
    }
    
    xsec = CrossSection(masspoints_pi0[i],mass_pi0,"physrevd");
    physrevd_pi0->SetPoint(i,masspoints_pi0[i],xsec);
    if ( global_minimum_pi0 > xsec ) {
      global_minimum_pi0 = xsec;
    }
    if ( global_maximum_pi0 < xsec ) {
      global_maximum_pi0 = xsec;
    }

    
    xsec = CrossSection(masspoints_pi0[i],mass_pi0,"naive");
    naive_pi0->SetPoint(i,masspoints_pi0[i],xsec);
    if ( global_minimum_pi0 > xsec ) {
      global_minimum_pi0 = xsec;
    }
    if ( global_maximum_pi0 < xsec ) {
      global_maximum_pi0 = xsec;
    }

    
    xsec = CrossSection(masspoints_pi0[i],mass_pi0,"t2k");
    t2k_pi0->SetPoint(i,masspoints_pi0[i],xsec);
    if ( global_minimum_pi0 > xsec ) {
      global_minimum_pi0 = xsec;
    }
    if ( global_maximum_pi0 < xsec ) {
      global_maximum_pi0 = xsec;
    }

  }

  for ( int i = 0; i < masspoints_eta.size(); i++ ) {
    xsec = CrossSection(masspoints_eta[i],mass_eta,"zhenliu");
    zhenliu_eta->SetPoint(i,masspoints_eta[i],xsec);
    if ( global_minimum_eta > xsec ) {
      global_minimum_eta = xsec;
    }
    if ( global_maximum_eta < xsec ) {
      global_maximum_eta = xsec;
    }
    
    xsec = CrossSection(masspoints_eta[i],mass_eta,"physrevd");
    physrevd_eta->SetPoint(i,masspoints_eta[i],xsec);
    if ( global_minimum_eta > xsec ) {
      global_minimum_eta = xsec;
    }
    if ( global_maximum_eta < xsec ) {
      global_maximum_eta = xsec;
    }

    
    xsec = CrossSection(masspoints_eta[i],mass_eta,"naive");
    naive_eta->SetPoint(i,masspoints_eta[i],xsec);
    if ( global_minimum_eta > xsec ) {
      global_minimum_eta = xsec;
    }
    if ( global_maximum_eta < xsec ) {
      global_maximum_eta = xsec;
    }

    
    xsec = CrossSection(masspoints_eta[i],mass_eta,"t2k");
    t2k_eta->SetPoint(i,masspoints_eta[i],xsec);
    if ( global_minimum_eta > xsec ) {
      global_minimum_eta = xsec;
    }
    if ( global_maximum_eta < xsec ) {
      global_maximum_eta = xsec;
    }

  }

  auto c1 = new TCanvas();
  //c1->SetLogy();
  //c1->SetLogx();
  
  zhenliu_pi0->SetLineColor(kRed);
  physrevd_pi0->SetLineColor(kBlue);
  naive_pi0->SetLineColor(kGreen);
  t2k_pi0->SetLineColor(kOrange);

  zhenliu_pi0->SetLineWidth(2);
  physrevd_pi0->SetLineWidth(2);
  naive_pi0->SetLineWidth(2);
  t2k_pi0->SetLineWidth(2);

  std::cout << global_minimum_pi0 << std::endl;
  std::cout << global_maximum_pi0 << std::endl;
  zhenliu_pi0->SetMinimum(global_minimum_pi0 - 0.1*global_minimum_pi0);
  zhenliu_pi0->SetMaximum(global_maximum_pi0 + 0.1*global_maximum_pi0);
  zhenliu_pi0->Draw("A L*");
  physrevd_pi0->Draw("same L*");
  naive_pi0->Draw("same L*");
  t2k_pi0->Draw("same L*");

  zhenliu_pi0->SetTitle("arXiv:1902.03246 (ArgoNeuT)");
  physrevd_pi0->SetTitle("Phys. Rev. D 100, 015043 (FerMINI)");
  naive_pi0->SetTitle("arXiv: 1812.03998 (FerMINI)");
  t2k_pi0->SetTitle("arXiv: 2103.11814 (T2K/T2HK)");
  
  TLegend *myleg;
  myleg = c1->BuildLegend(0.55,0.65,0.88,0.85,"","");
  
  zhenliu_pi0->SetTitle("#pi^{0} decays, #epsilon = 1");
  zhenliu_pi0->GetXaxis()->SetTitle("m_{#chi} (GeV)");
  zhenliu_pi0->GetYaxis()->SetTitle("Supression factor");
  zhenliu_pi0->GetXaxis()->CenterTitle();
  zhenliu_pi0->GetYaxis()->CenterTitle();

  c1->SaveAs("img/SupressionFactor/Factor_pi0_decays.png");
  c1->SaveAs("img/SupressionFactor/Factor_pi0_decays.pdf");

  auto c2 = new TCanvas();
  //c2->SetLogy();
  //c2->SetLogx();
  
  zhenliu_eta->SetLineColor(kRed);
  physrevd_eta->SetLineColor(kBlue);
  naive_eta->SetLineColor(kGreen);
  t2k_eta->SetLineColor(kOrange);

  zhenliu_eta->SetLineWidth(2);
  physrevd_eta->SetLineWidth(2);
  naive_eta->SetLineWidth(2);
  t2k_eta->SetLineWidth(2);

  std::cout << global_minimum_eta << std::endl;
  std::cout << global_maximum_eta << std::endl;
  zhenliu_eta->SetMinimum(global_minimum_eta - 0.1*global_minimum_eta);
  zhenliu_eta->SetMaximum(global_maximum_eta + 0.1*global_maximum_eta);
  zhenliu_eta->Draw("A L*");
  physrevd_eta->Draw("same L*");
  naive_eta->Draw("same L*");
  t2k_eta->Draw("same L*");

  zhenliu_eta->SetTitle("arXiv:1902.03246 (ArgoNeuT)");
  physrevd_eta->SetTitle("Phys. Rev. D 100, 015043 (FerMINI)");
  naive_eta->SetTitle("arXiv: 1812.03998 (FerMINI)");
  t2k_eta->SetTitle("arXiv: 2103.11814 (T2K/T2HK)");

  myleg = c2->BuildLegend(0.55,0.65,0.88,0.85,"","");
  
  zhenliu_eta->SetTitle("#eta decays, #epsilon = 1");
  zhenliu_eta->GetXaxis()->SetTitle("m_{#chi} (GeV)");
  zhenliu_eta->GetYaxis()->SetTitle("Supression factor");
  zhenliu_eta->GetXaxis()->CenterTitle();
  zhenliu_eta->GetYaxis()->CenterTitle();

  c2->SaveAs("img/SupressionFactor/Factor_eta_decays.png");
  c2->SaveAs("img/SupressionFactor/Factor_eta_decays.pdf");
}
