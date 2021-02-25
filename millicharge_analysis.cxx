bool gPlot = false;

void millicharge_analysis() {
  Double_t kMCPmass   = 0.01;    // GeV
  Double_t kMCPcharge = 0.01;   // electron charges

  TString kNM = "pi0";   // eta or pi0
  TString kMode = "fhc"; // fhc or rhc
  TFile *f = new TFile(kMode+"_"+kNM+"s_ntuple.root");

  TTreeReader reader(kMode+"_"+kNM+"s",f);
  TTreeReaderValue<Float_t> Px(reader,"Px");
  TTreeReaderValue<Float_t> Py(reader,"Py");
  TTreeReaderValue<Float_t> Pz(reader,"Pz");
  TTreeReaderValue<Float_t> E(reader,"E");
  TTreeReaderValue<Float_t> x(reader,"x");
  TTreeReaderValue<Float_t> y(reader,"y");
  TTreeReaderValue<Float_t> z(reader,"z");
  TTreeReaderValue<Float_t> t(reader,"t");
  TTreeReaderValue<Float_t> weight(reader,"weight");
  TTreeReaderValue<Float_t> event(reader,"event");

  TDatabasePDG pdg;
  double kNMmass;
  double kLifetime;
  double c = 3e8;
  if ( kNM == "pi0" ) {
    kNMmass = pdg.GetParticle(111)->Mass();
    kLifetime = 8.4e-17;
  }
  if ( kNM == "eta" ) {
    kNMmass = pdg.GetParticle(221)->Mass();
    kLifetime = 5.0e-19;
  }

  
  
  Double_t fmM = sqrt(1.-(kMCPmass*kMCPmass)/(kNMmass*kNMmass));
  Double_t branching_ratio = kMCPcharge*kMCPcharge * 0.01174 * fmM;

  cout << Form("mCP charge %.3f mass %.3f, branching ratio %e",kMCPcharge,kMCPmass,branching_ratio) << endl;
  
  int events = 0;
  TLorentzVector _NM;
  TLorentzVector _NM_pos;

  TH1 * NMpx = new TH1F("NMpx",";P_{x} (GeV)",100,-20,20);
  TH1 * NMpy = new TH1F("NMpy",";P_{y} (GeV)",100,-20,20);
  TH1 * NMpz = new TH1F("NMpz",";P_{z} (GeV)",100,0,20);
  TH1 * NMe = new TH1F("NMe",";E (GeV)",100,0,120);
  TH1 * mcp1px = new TH1F("mcp1px",";P_{x} (GeV)",100,-20,20);
  TH1 * mcp1py = new TH1F("mcp1py",";P_{y} (GeV)",100,-20,20);
  TH1 * mcp1pz = new TH1F("mcp1pz",";P_{z} (GeV)",100,0,20);
  TH1 * mcp1e = new TH1F("mcp1e","mCP_{1};E (GeV)",100,0,40);
  TH1 * mcp2px = new TH1F("mcp2px",";P_{x} (GeV)",100,-20,20);
  TH1 * mcp2py = new TH1F("mcp2py",";P_{y} (GeV)",100,-20,20);
  TH1 * mcp2pz = new TH1F("mcp2pz",";P_{z} (GeV)",100,0,20);
  TH1 * mcp2e = new TH1F("mcp2e","mCP_{2};E (GeV)",100,0,40);
  TH1 * NMtheta = new TH1F("NMtheta",";#theta (Rad)",100,0,1.57);
  TH1 * mcp1theta = new TH1F("mcp1theta",";#theta (Rad)",100,0,1.57);
  TH1 * mcp2theta = new TH1F("mcp2theta",";#theta (Rad)",100,0,1.57);
  TH1 * NMtheta2 = new TH1F("NMtheta2",";#theta (Rad)",100,0,1.57);
  TH1 * mcp1thetadiff = new TH1F("mcp1thetadiff",";#theta_{#pi^{0}}-#theta_{mCP_{1}} (Rad)",100,0,3.15);
  TH1 * mcp2thetadiff = new TH1F("mcp2thetadiff",";#theta_{#pi^{0}}-#theta_{mCP_{2}} (Rad)",100,0,3.15);
  TH1 * gammathetadiff = new TH1F("gammathetadiff",";#theta_{#pi^{0}}-#theta_{#gamma} (Rad)",100,0,3.15);
  TH1 * traveldistance = new TH1F("traveldistance",";Travel distance (m)",100,0,0.00005);

  TH1 * mcpe_argo = new TH1F("mcpe_argo","mCP_{1},mCP_{2};E (GeV)",100,0,40);

  
  Double_t NM_per_POT = 0;
  while ( reader.Next() ) {
    events++;
    /*
    if ( events > 1000 ) {
      break;
    }
    */
    _NM.SetPxPyPzE(*Px,*Py,*Pz,*E);
    _NM_pos.SetXYZT(*x,*y,*z,*t);
    // need to make NMs travel a bit before decaying
    double travel = _NM.Beta()*c*_NM.Gamma()*kLifetime;
    traveldistance->Fill(travel,*weight);
    
    
    Double_t masses[3] = { 0, kMCPmass, kMCPmass};
    TGenPhaseSpace dalitz;
    dalitz.SetDecay(_NM, 3, masses);
    
    Double_t decay_weight = dalitz.Generate();
    
    TLorentzVector *gamma = dalitz.GetDecay(0);
    TLorentzVector *MCP1 = dalitz.GetDecay(1);
    TLorentzVector *MCP2 = dalitz.GetDecay(2);
    
    if ( events < 5 ) {
      cout << Form("NM mass %.3f, mcp1 Px %.3f, decay weight %.4f",_NM.M(), MCP1->Px(),decay_weight) << endl;
    }
    // final weight
    Double_t mcp_weight = *weight * decay_weight * branching_ratio;
    // NMs that decay to mcp (all of them)
    //Double_t NM_weight = *weight * branching_ratio;
    Double_t NM_weight = *weight;
    NMpx->Fill(_NM.Px(),NM_weight);
    NMpy->Fill(_NM.Py(),NM_weight);
    NMpz->Fill(_NM.Pz(),NM_weight);
    mcp1px->Fill(MCP1->Px(),mcp_weight);
    mcp1py->Fill(MCP1->Py(),mcp_weight);
    mcp1pz->Fill(MCP1->Pz(),mcp_weight);
    mcp2px->Fill(MCP2->Px(),mcp_weight);
    mcp2py->Fill(MCP2->Py(),mcp_weight);
    mcp2pz->Fill(MCP2->Pz(),mcp_weight);
    NMe->Fill(_NM.E(),NM_weight);
    NM_per_POT += *weight;
    NMtheta->Fill(_NM.Theta(),NM_weight);
    /*
    if ( _NM.Theta() > 1.45 ) {
      cout << _NM.Px() << " " << _NM.Py() << " " << _NM.Pz() << endl;
    }
    */
    mcp1e->Fill(MCP1->E(),mcp_weight);
    mcp2e->Fill(MCP2->E(),mcp_weight);
    mcp1theta->Fill(MCP1->Theta(),mcp_weight);
    mcp2theta->Fill(MCP2->Theta(),mcp_weight);

    if ( MCP1->Theta() < 0.00019344230 ) {
      mcpe_argo->Fill(MCP1->E(),mcp_weight);
    }
    if ( MCP2->Theta() < 0.00019344230 ) {
      mcpe_argo->Fill(MCP2->E(),mcp_weight);
    }
    
    mcp1thetadiff->Fill(_NM.Angle(MCP1->Vect()),mcp_weight);
    mcp2thetadiff->Fill(_NM.Angle(MCP2->Vect()),mcp_weight);
    gammathetadiff->Fill(_NM.Angle(gamma->Vect()),mcp_weight);

    Double_t Pt = sqrt((*Px)*(*Px) + (*Py)*(*Py));
    Double_t th = atan2(Pt,*Pz);
    NMtheta2->Fill(th,NM_weight);
  }

  auto * c1 = new TCanvas();
  TString hfilename;

  if ( gPlot == true ) {
    NMpx->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)NMpx->GetName()+".png");
    NMpy->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)NMpy->GetName()+".png");
    NMpz->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)NMpz->GetName()+".png");
    NMe->Draw("hist");
    c1->SetLogy();
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)NMe->GetName()+".png");
    NMtheta->Draw("hist");
    c1->SetLogy(0);
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)NMtheta->GetName()+".png");
    NMtheta2->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)NMtheta2->GetName()+".png");
    
    mcp1px->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcp1px->GetName()+".png");
    mcp1py->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcp1py->GetName()+".png");
    mcp1pz->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcp1pz->GetName()+".png");
    mcp1e->Draw("hist");
    c1->SetLogy();
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcp1e->GetName()+".png");
    mcp1theta->Draw("hist");
    c1->SetLogy(0);
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcp1theta->GetName()+".png");
    mcp1thetadiff->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcp1thetadiff->GetName()+".png");
    
    mcp2px->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcp2px->GetName()+".png");
    mcp2py->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcp2py->GetName()+".png");
    mcp2pz->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcp2pz->GetName()+".png");
    mcp2e->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcp2e->GetName()+".png");
    mcp2theta->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcp2theta->GetName()+".png");
    mcp2thetadiff->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcp2thetadiff->GetName()+".png");
    
    gammathetadiff->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)gammathetadiff->GetName()+".png");
    
    traveldistance->Draw("hist");
    c1->SetLogy();
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)traveldistance->GetName()+".png");
    
    mcpe_argo->Draw("hist");
    c1->SaveAs("./img/"+kMode+"_"+kNM+"s_"+(TString)mcpe_argo->GetName()+".png");
  }
  
  THStack *mcpe = new THStack("mcpe",";E (GeV);Entries #times BR #times TGen");
  mcp1e->SetFillColor(kBlue-6);
  mcp2e->SetFillColor(kBlue-4);
  //mcp1e->Scale(ftot->GetParameter(0));
  //mcp2e->Scale(ftot->GetParameter(1));
  mcpe->Add(mcp1e);
  mcpe->Add(mcp2e);
  //mcp->Draw("EP");
  c1->SetLogy();
  mcpe->Draw("hist");
  /*
  mcp1e->SetMarkerStyle(8);
  mcp1e->Draw("EP same");    
  */
  TLegend *mylegsb;
  mylegsb = c1->BuildLegend(0.62,0.8,0.9,0.9,"","");
  mylegsb->SetFillStyle(0);
  c1->SaveAs("./img/stack_"+kMode+"_"+kNM+"s_mcpe.png");

  // # of mCPs entering argoneut as shown in fig. 2 of 1902.03246
  double accepted_argo;
  accepted_argo = (10e20/500000.)*mcpe_argo->Integral();
  
  cout << "NM_per_POT: " << NM_per_POT << endl;
  cout << "accepted in argoneut: " << accepted_argo << endl;
}
