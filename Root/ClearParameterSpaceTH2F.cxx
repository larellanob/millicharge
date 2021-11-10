void ClearParameterSpaceTH2F(TString detector = "uboone", bool multithreshold = false)
{
  


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


  TH2 * lim1hit = new TH2F("lim1hit",";m_{#chi} (MeV);#epsilon",28,xaxis,29,yaxis);
  TH2 * lim2hit = new TH2F("lim2hit",";m_{#chi} (MeV);#epsilon",28,xaxis,29,yaxis);
  TH2 * lim3hit = new TH2F("lim3hit",";m_{#chi} (MeV);#epsilon",28,xaxis,29,yaxis);
  TH2 * lim4hit = new TH2F("lim4hit",";m_{#chi} (MeV);#epsilon",28,xaxis,29,yaxis);

  TFile *f;


  
  if ( multithreshold == false ) {
    f = new TFile("hist/Lim_"+detector+".root","recreate");
    lim1hit->Write();
    lim2hit->Write();
    lim3hit->Write();
    lim4hit->Write();
  } else if ( multithreshold == true ) {
    std::vector<TString> thresholds = {
      "0.1",
      "0.6",
      "0.8",
      "1.0"
    };
    for ( auto th: thresholds ) {
      f = new TFile("hist/Lim_"+detector+"_"+th+".root","recreate");
      lim1hit->Write();
      lim2hit->Write();
      lim3hit->Write();
      lim4hit->Write();
    }

    /*
    for ( auto th: thresholds ) {
      f = new TFile("hist/Lim_"+detector+"_"+th+".root","recreate");
      std::cout << "doing it" << std::endl;
      TH2F lim1hit_th("lim_1hit_"+th+"thres",";m_{#chi} (MeV);#epsilon",28,xaxis,29,yaxis);
      TH2F lim2hit_th("lim_2hit_"+th+"thres",";m_{#chi} (MeV);#epsilon",28,xaxis,29,yaxis);
      TH2F lim3hit_th("lim_3hit_"+th+"thres",";m_{#chi} (MeV);#epsilon",28,xaxis,29,yaxis);
      TH2F lim4hit_th("lim_4hit_"+th+"thres",";m_{#chi} (MeV);#epsilon",28,xaxis,29,yaxis);
      lim1hit_th.Write();
      lim2hit_th.Write();
      lim3hit_th.Write();
      lim4hit_th.Write();
    }
    */
  }
  
  
  
}
