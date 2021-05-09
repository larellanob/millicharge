#include "Root/CreateEmptyDirs.cxx"

Double_t nbkghits(Double_t empty, Int_t hits, double frames)
{
  // usage example
  //std::cout << nbkghits(0.88,1,3.26e6) << std::endl;
  //return 0;
  
  // empty = percentage of empty frames
  if ( empty > 1 ) {
    std::cout << "ERROR: can't have that many empty frames " << std::endl;
    return 0;
  }
  double phit = -log(empty);
  double sum = 0;
  double oldsum = -1;
  int n = hits;

  // summing, convergence for 0.00001
  while ( abs(oldsum-sum) > 0.00001 ) {
    std::cout << "sum " << sum << std::endl;
    oldsum = sum;
    sum += TMath::Binomial(n,hits)*pow(phit,n)/TMath::Factorial(n);
    n++;
  }
  return frames * exp(-phit) * sum;
}


void PlotPublishedLimits()
{
  CreateEmptyDirs();
  const char* dir = "/home/luciano/Physics/neutrino/millicharge/data/limits_argoneut_experiment_12hits/";
  
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
  
  //TGraph *gr = new TGraph(e0_size,x,y);
  TFile *f = new TFile("hist/Lim_argoneut.root");
  auto c1 = new TCanvas();
  c1->SetLogy();
  c1->SetLogx();
  c1->SetLogz();
  //tgraphs[1]->Draw("A L");

  TH2F * lim1hit;
  lim1hit = (TH2F*)gDirectory->Get("lim1hit");

  TH2F * lim2hit;
  lim2hit = (TH2F*)gDirectory->Get("lim2hit");

  //nbkghits(Double_t empty, Int_t hits, double frames)
  Double_t percentage_empty = 0.88;
  Double_t frames = 3.26e6;
  Double_t threshold = nbkghits(percentage_empty,1,frames);
  Double_t threshold2hit = nbkghits(percentage_empty,2,frames);

  std::cout << "threshold " << threshold << std::endl;
  TGraph *uboone1hit = new TGraph(8);
  TGraph *uboone2hit = new TGraph(8);

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

  // argoneut numbers
  Double_t resolution_factor = (0.3/470)*(5.6/400); // put microboone resolution
  
  std::cout << "Resolution factor " << resolution_factor << std::endl;
  std::cout << "bkg 2 hits " << threshold2hit << std::endl;
  std::cout << "reduced 2 hit bkg " << resolution_factor*threshold2hit << std::endl;
  int counter = -1;
  std::cout << "ok " << std::endl;
  for ( int x: massbins ) {
    counter++;
    double centerx = lim1hit->GetXaxis()->GetBinCenter(x);
    bool findlim1hit = true;
    bool findlim2hit = true;

    // searches the first point in th2f which is greater than
    // threshold and sets the limit there (stops searching)
    for ( int y = 0; y < 29; y++ ) {
      // 1 hit
      if ( lim1hit->GetBinContent(lim1hit->GetBin(x,y+1)) > threshold  && findlim1hit ) {
	
	double meany
	  = 0.5*(lim1hit->GetYaxis()->GetBinCenter(y+1)
		 + lim1hit->GetYaxis()->GetBinCenter(y));
	uboone1hit->SetPoint(counter,centerx,meany);
	findlim1hit = false; // stop searching
      }
      // 2 hits
      if ( lim2hit->GetBinContent(lim2hit->GetBin(x,y+1))*resolution_factor > 3.5 && findlim2hit) {
	double meany
	  = 0.5*(lim2hit->GetYaxis()->GetBinCenter(y+1)
		 + lim2hit->GetYaxis()->GetBinCenter(y));
	uboone2hit->SetPoint(counter,centerx,meany);
	findlim2hit = false; // stop searching
      }
      // if it's the last point and didn't find limit, set it at the
      // highest charge (top of y axis)
      if ( y == 28 ) {
	double meany
	  = 0.5*(lim1hit->GetYaxis()->GetBinCenter(y+1)
		 + lim1hit->GetYaxis()->GetBinCenter(y));
	if ( findlim1hit == true ) {
	  // reached the end
	  std::cout << "Reached the end and didn't find charge limit" << std::endl;
	  uboone1hit->SetPoint(counter,centerx,meany);
	}
	if ( findlim2hit == true ) {
	  // reached the end
	  std::cout << "Reached the end and didn't find charge limit" << std::endl;
	  uboone2hit->SetPoint(counter,centerx,meany);
	}
      }
    }
    
  }
  gStyle->SetOptStat(0);
  lim1hit->SetTitle("1 hit");
  lim1hit->Draw("colz text");
  
  uboone1hit->SetLineColor(kYellow);
  uboone1hit->SetLineWidth(5);

  uboone2hit->SetLineColor(kGreen);
  uboone2hit->SetLineWidth(5);
  
  tgraphs[1]->SetLineColor(kBlue);
  tgraphs[1]->SetLineWidth(5);

  tgraphs[0]->SetLineColor(kRed);
  tgraphs[0]->SetLineWidth(5);

  tgraphs[1]->Draw("same L");
  uboone1hit->Draw("same L");

  TLegend *myleg;
  uboone1hit->SetTitle("ArgoNeuT 1 hit our simulation");
  tgraphs[1]->SetTitle("ArgoNeuT 1 hit published");
  myleg = c1->BuildLegend(0.6,0.15,0.9,0.35,"","");
  auto c2 = new TCanvas();
  c2->SetLogy();
  c2->SetLogx();
  c2->SetLogz();
  gStyle->SetOptStat(0);
  lim2hit->SetTitle("2 hits");
  lim2hit->Draw("colz");
  tgraphs[0]->Draw("same L");
  uboone2hit->Draw("same L");
  

  //TF1 * fun = new TF1("fun",[&](double*x, double *p){ return p[0]*tgraphs[1]->Eval(x[0]); }, 2.5, 2.7, 1);
  //fun->Draw();
  
}
