// generates uboone geometry wrt numi and stores as a root file in geom/ dir

void geom_uboone(TString fstr ="sim/mCP_q_0.010_m_0.020_fhc_pi0s.root" )
{

  TGeoManager *geom = new TGeoManager("geom","My first 3D geometry");
  //TEveManager::Create();
  // name, a, z, rho
  TGeoMaterial *vacuum=new TGeoMaterial("vacuum",0,0,0);
  TGeoMaterial *Fe=new TGeoMaterial("Fe",55.845,26,7.87);

  // name, numed?, material
  TGeoMedium *Air=new TGeoMedium("Vacuum",0,vacuum);
  TGeoMedium *Iron=new TGeoMedium("Iron",1,Fe);


  // transformation vectors/matrices
  const TVector3 beampos(-31387.58422,
			 -3316.402543,
			 -60100.2414
			 ); // proton on target in uboone coords * cm
  
  //rotation from beam to det coordinates
  const TRotation rot
    = []()
      {
	TRotation R;
	// Rotation matrix using the 0,0,0 position for MicroBooNE (beam to det input)
	const TVector3 Rot_row_x = {  0.92103853804025681562,
				      0.022713504803924120662,
				      0.38880857519374290021  };
	const TVector3 Rot_row_y = {  4.6254001262154668408e-05,
				      0.99829162468141474651,
				      -0.058427989452906302359 };
	const TVector3 Rot_row_z = { -0.38947144863934973769,
				     0.053832413938664107345,
				     0.91946400794392302291  };
	
	R.RotateAxes(Rot_row_x, Rot_row_y, Rot_row_z); // Also inverts so now to det to beam
	R.Invert(); // go back to beam to det
	return R;
      }();
  TRotation R;
  const TVector3 Rot_row_x = {  0.92103853804025681562,
				0.022713504803924120662,
				0.38880857519374290021  };
  const TVector3 Rot_row_y = {  4.6254001262154668408e-05,
				0.99829162468141474651,
				-0.058427989452906302359 };
  const TVector3 Rot_row_z = { -0.38947144863934973769,
			       0.053832413938664107345,
			       0.91946400794392302291  };
  R.RotateAxes(Rot_row_x, Rot_row_y, Rot_row_z);
  //R.Invert();

  /*
  const TVector3 pos_argo(-31387.58422,
			 -3316.402543,
			 -60100.2414
			 ); // proton on target in uboone coords * cm
  */
  
  // georotation angles are input in degrees for some reason
  TGeoRotation * geo_rot = new TGeoRotation("geo_rot",R.GetXPhi()*180./TMath::Pi(),R.GetXTheta()*180./TMath::Pi(),R.GetXPsi()*180./TMath::Pi());
  TGeoRotation * geo_rot2 = new TGeoRotation("geo_rot",-R.GetXPhi()*180./TMath::Pi(),-R.GetXTheta()*180./TMath::Pi(),-R.GetXPsi()*180./TMath::Pi());
  
  cout << Form("Radians: Phi %.3f, Theta %.3f, Psi %.3f\n",R.GetXPhi(),R.GetXTheta(),R.GetXPsi());
  cout << Form("Degrees: Phi %.3f, Theta %.3f, Psi %.3f\n",R.GetXPhi()*180./TMath::Pi(),R.GetXTheta()*180./TMath::Pi(),R.GetXPsi()*180./TMath::Pi());
  
  // top geometry, the Universe
  TGeoVolume *top=geom->MakeBox("top",
				Air,
				1000000,
				1000000,
				1000000
				);
  geom->SetTopVolume(top);
  geom->SetTopVisible(0);
 
  // uboone detector
  //TGeoVolume *det_uboone=geom->MakeBox("det_uboone",Air,236.35,213.,1016.7);
  TGeoVolume *det_uboone=geom->MakeBox("det_uboone",
				       Air,
				       236.35/2.,
				       213./2.,
				       1016.7/2.
				       );
  top->AddNodeOverlap(det_uboone,0);
  // det_uboone->SetLineColor(kBlue);
  //det_uboone->SetFillStyle(0);
  det_uboone->SetLineColorAlpha(kGreen,0.0);
  det_uboone->SetLineStyle(10);
  // argoneut detector

  TVector3 pos_argo(0,0,97500); // pos in numi frame
  pos_argo.Print();
  pos_argo = beampos+rot*pos_argo; // pos in uboone frame
  pos_argo.Print();

  TGeoVolume *det_argoneut=geom->MakeBox("det_argoneut",
					 Air,
					 470./2.,
					 400./2.,
					 900./2.
					 );
  Double_t dist_argoneut = 97500; // distance in z to target (cm)
  TGeoCombiTrans * trans_argoneut = new  TGeoCombiTrans(pos_argo.X(),pos_argo.Y(),pos_argo.Z(),geo_rot);
  top->AddNodeOverlap(det_argoneut,0, trans_argoneut);
  det_argoneut->SetLineColor(kGreen);
  

   // name, medium, rmin, rmax, themin, themax, phimin, phimax
   TGeoVolume *numitarg=geom->MakeSphere("numitarg",Iron,0,350,0,180,0,360);
   //numitarg->SetLineColor(41);
   //det_uboone->AddNodeOverlap(numitarg,1,new TGeoCombiTrans(0,0,4000,new TGeoRotation("numitarg",0,0,0)));
   top->AddNodeOverlap(
		       numitarg,
		       0,
		       new TGeoCombiTrans(beampos.X(),
					  beampos.Y(),
					  beampos.Z(),
					  geo_rot)
		       );
   numitarg->SetLineColor(kGreen);

   TGeoVolume *bnbtarg=geom->MakeSphere("bnbtarg",Iron,0,350,0,180,0,360);
   //bnbtarg->SetLineColor(41);
   //det_uboone->AddNodeOverlap(bnbtarg,1,new TGeoCombiTrans(0,0,4000,new TGeoRotation("bnbtarg",0,0,0)));
   top->AddNodeOverlap(
		       bnbtarg,
		       0,
		       new TGeoCombiTrans(0,
					  0,
					  -47000,
					  0)
		       );
   bnbtarg->SetLineColor(kBlue);

   // numi beam
   TGeoVolume * numibeam = geom->MakeTube("numibeam",
					  Iron,
					  0,
					  50,
					  dist_argoneut/2.
					  );
   numibeam->SetLineColor(kGreen);
   TVector3 pos_numibeam(0,0,97500/2.); // pos in numi frame
   pos_numibeam = beampos+rot*pos_numibeam; // pos in uboone frame

   top->AddNodeOverlap(
		       numibeam,
		       1,
		       new TGeoCombiTrans(pos_numibeam.X(),
					  pos_numibeam.Y(),
					  pos_numibeam.Z(),
					  geo_rot2)
		       ); 

   // bnb beam
   TGeoVolume * bnbbeam = geom->MakeTube("bnbbeam",
					  Iron,
					  0,
					  50,
					  47000/2.
					  );
   bnbbeam->SetLineColor(kBlue);
   TVector3 pos_bnbbeam(0,0,-47000/2.); // pos in numi frame
   //pos_bnbbeam = beampos+rot*pos_bnbbeam; // pos in uboone frame

   top->AddNodeOverlap(
		       bnbbeam,
		       1,
		       new TGeoCombiTrans(pos_bnbbeam.X(),
					  pos_bnbbeam.Y(),
					  pos_bnbbeam.Z(),
					  0)
		       ); 


   
   /// event script (separate these two)
   /*
   auto tracksXYZ = new TEveStraightLineSet("StraightLinesXYZ");
   tracksXYZ->SetLineColor(kRed;
   tracksXYZ->SetLineWidth(2);


   TFile *f  = new TFile(fstr);
   TTreeReader reader("mCP",f);
   TTreeReaderValue<TLorentzVector> Mom(reader,"Mom");
   TTreeReaderValue<TLorentzVector> Pos(reader,"Pos");
   TTreeReaderValue<Float_t> weight_decay(reader,"WeightDecay");
   TTreeReaderValue<Float_t> weight_meson(reader,"WeightMeson");
   int events = 0;
    while ( reader.Next() ) {
    events++;
    if ( events < 6 ) {
      //cout  << Pt << endl;
    } else break;

      TEveVector trackDir(Mom->Px(), Mom->Py() ,Mom->Pz());
    Float_t length = 4000;
    TEveVector trackEnd = trackDir * length;
    tracksXYZ->AddLine(0., 0., 0., trackEnd.fX, trackEnd.fY, trackEnd.fZ );
    //tracksXYZ->AddLine(0., 0., 0., trackDir.fX, trackDir.fY, trackDir.fZ );
    }
    gEve->AddElement(tracksXYZ);

   */
  
   top->SetVisibility(0);

   geom->Export("geom/uboone_numi.root","uboone_numi");
   //top->Draw("ogl");
   geom->CloseGeometry();
}
