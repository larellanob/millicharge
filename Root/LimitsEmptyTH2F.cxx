void LimitsEmptyTH2F()
{
  TFile *f = new TFile("hist/Lim_uboone.root","recreate");


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

  
  
		       //0.00005,0.00015,0.00025,0.00035,0.00045,0.00055,0.00065,0.00075,0.00085,0.00095,
  
  //TH2 * lim = new TH2F("lim",";m_{#chi} (GeV);#epsilon",28,xaxis,29,yaxis);
  TH2 * lim1hit = new TH2F("lim1hit",";m_{#chi} (GeV);#epsilon",28,xaxis,29,yaxis);
  TH2 * lim2hit = new TH2F("lim2hit",";m_{#chi} (GeV);#epsilon",28,xaxis,29,yaxis);
  TH2 * lim3hit = new TH2F("lim3hit",";m_{#chi} (GeV);#epsilon",28,xaxis,29,yaxis);
  TH2 * lim4hit = new TH2F("lim4hit",";m_{#chi} (GeV);#epsilon",28,xaxis,29,yaxis);
  lim1hit->Write();
  lim2hit->Write();
  lim3hit->Write();
  lim4hit->Write();

  /*
  TF2 *xyg = new TF2("xyg","xygaus",0,10,0,10);
  //amplitude, meanx,sigmax,meany,sigmay 
  xyg->SetParameters(10,200,50,0.01,0.02);
  //lim->FillRandom("xyg",10000000);
  
  lim->Write();
  auto c1 = new TCanvas();
  c1->SetLogy();
  c1->SetLogx();
  lim->Draw("colz text");
  c1->SaveAs("img/Lim_test.png");
  */
}
