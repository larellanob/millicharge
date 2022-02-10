#include "Root/CreateEmptyDirs.cxx"
void MultiDraw(std::vector<TGraph * >, std::vector<std::vector<std::vector<TGraph *>>>,TH2F*);

// this macro iterates over these three vectors

// comment values you are not interested in

// when trying to add new values you also have to add it to the
// previous run_Limits.sh list

std::vector<TString> Ethres_vector = {
  "0.1",
  "0.6",
  "0.8",
  "1.0"
};

std::vector<int> events = 
  {
    1,
    10,
    100
  };

std::vector<int> nhits_vector = {
  1,
  2,
  3,
  4
};
  


void PlotPublishedLimits()
{
  CreateEmptyDirs();

  // this method of getting the data might create some problems. the
  // order of file iteration chances (it's not alphabetical), so in
  // the final plot the legends might be off
  const char* dir = "/home/luciano/Physics/manchester/millicharge/data/limits_argoneut_experiment_stfcposter/";
  
  //char *dir  = gSystem->ExpandPathName(dir_tstring);
  
  void *dirpointer = gSystem->OpenDirectory(dir);

  std::vector<TString> full_files;
  std::vector<TString> files;
  std::vector<TString> experiments;
  const char *entry;
  TString dir_iteration;
  while ( ( entry = (char*)gSystem->GetDirEntry(dirpointer))) {
    dir_iteration = entry;
    if ( dir_iteration.EndsWith(".") || dir_iteration == "README.txt" ) {
      continue;
    }
    full_files.push_back(dir+dir_iteration);
    files.push_back(dir_iteration);
    experiments.push_back(dir_iteration.ReplaceAll(".dat",""));
    std::cout << dir+dir_iteration << std::endl;
  }

  std::cout << "Found limits for the following experiments:" << std::endl;
  for ( int i = 0; i < experiments.size(); i++ ) {
    std::cout << experiments[i] << std::endl;
  }

  std::vector<std::vector<Double_t>> experimentsx;
  std::vector<std::vector<Double_t>> experimentsy;
  std::vector<Double_t> datax;
  std::vector<Double_t> datay;

  for ( int i = 0; i < files.size(); i++ ) {
    Double_t x,y;
    std::cout << "opening file " << full_files[i] << std::endl;
    ifstream in;
    in.open(full_files[i]);
    datax.clear();
    datay.clear();
    while (1) {
      in >> x >> y;
      //std::cout << x << " " << y << std::endl;
      datax.push_back(x);
      datay.push_back(y);
      if (!in.good()) break;
    }
    experimentsx.push_back(datax);
    experimentsy.push_back(datay);
  }
  std::cout << experimentsx[0].size() << std::endl;


  std::vector<TGraph *> tgraphs;
  for ( int exp = 0; exp < experiments.size(); exp++ ) {
    TGraph *gr = new TGraph(experimentsx[exp].size());
    tgraphs.push_back(gr);
  }
  
  for ( int exp = 0; exp < experiments.size(); exp++ ) {
    std::cout << "Graph " << experiments[exp] << std::endl;
    for ( int point = 0; point < experimentsx[exp].size(); point++ ){
      //std::cout << "point x " << point << " " <<experimentsx[exp][point] << std::endl;
      tgraphs[exp]->SetPoint(point, experimentsx[exp][point], experimentsy[exp][point]);
    }
  }

  
  TH2F * limnhit;

  // 3d vector
  // events_graphs[low energy threshold][n of events][n of scatters]
  std::vector<std::vector<std::vector<TGraph *>>> events_graphs; 


  // x-bin numbers where calculated mass points are stored
  // 0.01 0.02 0.03 0.05 0.06 0.1 0.2 0.25 GeV
  std::vector<int> massbins =
    {
     1,
     2,
     3,
     5,
     6,
     10,
     11,
     12
    };

  int counter = -1;
  std::cout << "ok " << std::endl;

  std::vector<std::vector<TGraph*>> vec2;
  std::vector<TGraph*> vec1;
  
  for ( int eth = 0; eth < Ethres_vector.size(); eth++ ) {
    TFile *f = new TFile("hist/Lim_uboone_"+Ethres_vector[eth]+".root");
    vec2.clear();
    for ( int e = 0; e < events.size(); e++ ) {
      vec1.clear();
      for ( auto nhit: nhits_vector ) {
	if ( nhit == 1 ) {
	  limnhit = (TH2F*)gDirectory->Get("lim1hit");
	} else if ( nhit == 2 ) {
	  limnhit = (TH2F*)gDirectory->Get("lim2hit");
	} else if ( nhit == 3 ) {
	  limnhit = (TH2F*)gDirectory->Get("lim3hit");
	} else if ( nhit == 4 ) {
	  limnhit = (TH2F*)gDirectory->Get("lim4hit");
	}
	TGraph *uboonenhit = new TGraph(8);
	for ( int x: massbins ) {

	  counter++;
	  double centerx = limnhit->GetXaxis()->GetBinCenter(x);
	  bool findlimnhit = true;
	  
	  // searches the first point in th2f which is greater than
	  // threshold and sets the limit there (stops searching)
	  for ( int y = 0; y < 29; y++ ) {
	    // find nhit threshold

	    if ( limnhit->GetBinContent(limnhit->GetBin(x,y+1)) > events[e]  && findlimnhit ) {
	      
	      double meany
		= 0.5*(limnhit->GetYaxis()->GetBinCenter(y+1)
		       + limnhit->GetYaxis()->GetBinCenter(y));
	      double lowedge = limnhit->GetYaxis()->GetBinLowEdge(y+1);
	      uboonenhit->SetPoint(counter,centerx,lowedge);
	      //uboonenhit->SetPoint(counter,centerx,meany);
	      findlimnhit = false; // stop searching
	    }
	    
	    // if it's the last point and didn't find limit, set it at the
	    // highest charge (top of y axis)
	    if ( y == 28 ) {
	      double meany
		= 0.5*(limnhit->GetYaxis()->GetBinCenter(y+1)
		       + limnhit->GetYaxis()->GetBinCenter(y));
	      if ( findlimnhit == true ) {
		// reached the end
		std::cout << "Reached the end and didn't find charge limit" << std::endl;
		uboonenhit->SetPoint(counter,centerx,meany);
	      }
	    }
	  }
	}
	vec1.push_back(uboonenhit);
	
      }
      vec2.push_back(vec1);
    }
    events_graphs.push_back(vec2);
  }

  std::cout << "finished reading, about to Draw" << std::endl;
  MultiDraw(tgraphs, events_graphs, limnhit);
}

void MultiDraw(std::vector<TGraph * > published, std::vector<std::vector<std::vector<TGraph *>>> limits,TH2F* dummy)
{
  // draws all combinations of plotsn
  
  // e.g. keeping threshold and number of events fixed, make plot for
  // different number of scatters

  // three loops in total

  std::cout << "published size: " << published.size() << std::endl;
  std::cout << "limits size: " << limits.size() << std::endl;
  std::cout << "limits[0] size: " << limits[0].size() << std::endl;
  std::cout << "limits[0][0] size: " << limits[0][0].size() << std::endl;
  
  gStyle->SetOptStat(0);
  TString detector = "uboone";
  TString outname;

  // base histogram (axes and title)

  dummy->GetXaxis()->CenterTitle();
  dummy->GetYaxis()->SetTitle("Millicharge #epsilon = Q/e");
  dummy->GetXaxis()->SetTitle("Mass m_{#chi} (MeV)");
  dummy->GetYaxis()->CenterTitle();
  dummy->Reset();

  //dummy->SetTitle(Ethres+" MeV threshold");

  // experiment lines
  published[1]->SetLineColor(kBlue);
  published[1]->SetLineStyle(2);
  published[1]->SetLineWidth(5);

  published[0]->SetLineColor(kRed);
  published[0]->SetLineWidth(3);

  published[2]->SetLineColor(kTeal);
  published[2]->SetLineWidth(5);

  published[3]->SetLineColor(kGray);
  published[3]->SetLineWidth(3);

  published[4]->SetLineColor(kViolet);
  published[4]->SetLineWidth(3);

  published[5]->SetLineColor(kBlack);
  published[5]->SetLineWidth(3);

  published[0]->SetTitle("MilliQ");
  published[1]->SetTitle("MicroBooNE (estimate using BNB)");
  published[2]->SetTitle("ArgoNeuT");
  published[3]->SetTitle("LHC");
  published[4]->SetTitle("MiniBooNE");
  published[5]->SetTitle("LSND");

  gStyle->SetOptStat(0);


  std::vector<TGraph*>::iterator it3;
  std::vector<std::vector<TGraph*>>::iterator it2;
  std::vector<std::vector<std::vector<TGraph*>>>::iterator it1;
  
  std::cout << "loop over thresholds " << std::endl;
  //////////////////////
  // loop over thesholds
  bool draw_clean = true;
  for ( int i2 = 0; i2 < limits[0].size(); i2++ ) {
    for ( int i3 = 0; i3 < limits[0][0].size(); i3++ ) {
      auto c1 = new TCanvas();
      c1->SetLogy();
      c1->SetLogx();
      c1->SetLogz();

      dummy->Draw("text");
      published[3]->Draw("L");
      published[0]->Draw("L");
      published[1]->Draw("L"); // pheno estimate
      published[2]->Draw("L");
      published[4]->Draw("L");
      published[5]->Draw("L");
  
      TLegend *myleg = new TLegend(0.5,0.15,0.9,0.45);
      myleg->AddEntry(published[0]);
      myleg->AddEntry(published[1]); // pheno estimate
      myleg->AddEntry(published[2]);
      myleg->AddEntry(published[3]);
      myleg->AddEntry(published[4]);
      myleg->AddEntry(published[5]);

      // clean image
      if ( draw_clean ) {
	outname = Form("img/Limits/Limits.pdf");
	myleg->Draw();
	c1->SaveAs(outname);
	draw_clean = false;
      }
      
      dummy->SetTitle(Form("#muBooNE - %i-hit signal events, %i scatters",events[i2],nhits_vector[i3]));
      for ( int i1 = 0; i1 < limits.size(); i1++ ) {
	limits[i1][i2][i3]->SetLineColor(kGreen+i1);
	limits[i1][i2][i3]->SetLineWidth(5);
	limits[i1][i2][i3]->SetTitle(Form("%s MeV threshold",Ethres_vector[i1].Data()));
	limits[i1][i2][i3]->SetLineStyle(10-i1);
	limits[i1][i2][i3]->Draw("L");
	myleg->AddEntry(limits[i1][i2][i3]);
      }
      outname = Form("img/Limits/Limits_loop-thresholds_%ievents_%ihits.pdf",events[i2],nhits_vector[i3]);
      myleg->Draw();
      c1->SaveAs(outname);
    }
    
  }

  std::cout << "loop over nevents" << std::endl;
  //////////////////////
  // loop over events
  for ( int i1 = 0; i1 < limits.size(); i1++ ) {
    for ( int i3 = 0; i3 < limits[0][0].size(); i3++ ) {
      auto c1 = new TCanvas();
      c1->SetLogy();
      c1->SetLogx();
      c1->SetLogz();
      dummy->SetTitle(Form("#muBooNE - %s MeV threshold, %i scatters",Ethres_vector[i1].Data(),nhits_vector[i3]));
      dummy->Draw("text");
      published[3]->Draw("L");
      published[0]->Draw("L");
      //published[1]->Draw("L"); // pheno estimate
      published[2]->Draw("L");
      published[4]->Draw("L");
      published[5]->Draw("L");
  
      TLegend *myleg = new TLegend(0.5,0.15,0.9,0.45);
      myleg->AddEntry(published[0]);
      // myleg->AddEntry(published[1]); pheno estimate
      myleg->AddEntry(published[2]);
      myleg->AddEntry(published[3]);
      myleg->AddEntry(published[4]);
      myleg->AddEntry(published[5]);
      

      for ( int i2 = 0; i2 < limits[0].size(); i2++ ) {
	limits[i1][i2][i3]->SetLineColor(kGreen+i2);
	limits[i1][i2][i3]->SetLineWidth(5);
	limits[i1][i2][i3]->SetTitle(Form("%i signal events",events[i2]));
	limits[i1][i2][i3]->SetLineStyle(10-i2);
	limits[i1][i2][i3]->Draw("L");
	myleg->AddEntry(limits[i1][i2][i3]);
      }

      outname = Form("img/Limits/Limits_loop-events_%sthres_%ihits.pdf",Ethres_vector[i1].Data(),nhits_vector[i3]);
      myleg->Draw();
      c1->SaveAs(outname);

    }
    
  }

  std::cout << "loop over nhits" << std::endl;
  //////////////////////
  // loop over nhits
  for ( int i1 = 0; i1 < limits.size(); i1++ ) {
    for ( int i2 = 0; i2 < limits[0].size(); i2++ ) {    
      auto c1 = new TCanvas();
      c1->SetLogy();
      c1->SetLogx();
      c1->SetLogz();
      dummy->SetTitle(Form("#muBooNE %s MeV threshold,  %i-hit signal events",Ethres_vector[i1].Data(),events[i2]));      
      dummy->Draw("text");
      published[3]->Draw("L");
      published[0]->Draw("L");
      //published[1]->Draw("L"); // pheno estimate
      published[2]->Draw("L");
      published[4]->Draw("L");
      published[5]->Draw("L");
  
      TLegend *myleg = new TLegend(0.5,0.15,0.9,0.45);
      myleg->AddEntry(published[0]);
      // myleg->AddEntry(published[1]); pheno estimate
      myleg->AddEntry(published[2]);
      myleg->AddEntry(published[3]);
      myleg->AddEntry(published[4]);
      myleg->AddEntry(published[5]);
      
      for ( int i3 = 0; i3 < limits[0][0].size(); i3++ ) {
	limits[i1][i2][i3]->SetLineColor(kGreen+i3);
	limits[i1][i2][i3]->SetLineWidth(5);
	limits[i1][i2][i3]->SetTitle(Form("%i scatters",nhits_vector[i3]));
	limits[i1][i2][i3]->SetLineStyle(10-i3);
	limits[i1][i2][i3]->Draw("L");
	myleg->AddEntry(limits[i1][i2][i3]);
      }
      
      outname = Form("img/Limits/Limits_loop-scatters_%sthres_%ievents.pdf",Ethres_vector[i1].Data(),events[i2]);
      myleg->Draw();
      c1->SaveAs(outname);
    }
    
  }

  std::cout << "finished. " << std::endl;
  
}
  
