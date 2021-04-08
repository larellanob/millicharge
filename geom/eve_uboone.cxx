#include "../UbooneAcceptanceChecker.cxx"
const char* esd_geom_file_name =
   "/home/luciano/Physics/neutrino/millicharge/geom/uboone_numi.root";

void geomGentleTPC();

void eve_uboone(TString fstr ="../sim/mCP_q_0.010_m_0.020_fhc_pi0s.root" )
{

  TEveManager::Create();

  gGeoManager =
    gEve->GetGeometry("/home/luciano/Physics/neutrino/millicharge/geom/uboone_numi.root");

  TGeoNode *_wnode = (TGeoNode *)gGeoManager->GetListOfNodes()->At(0);
  TEveGeoTopNode *_node = new TEveGeoTopNode(gGeoManager, _wnode);
  gEve->AddGlobalElement(_node);
  
  //auto det_argoneut = new TEveGeoTopNode(gGeoManager,"det_argoneut");
  //auto det_argo = (TEveElement*) gGeoManager->GetVolume("det_argoneut");
  //det_argo->SetMainTransparency(0.5);
  //auto det_uboone = (TEveElement*) gGeoManager->GetVolume("det_uboone");
  //gEve->AddGlobalElement(det_argo);
  //gEve->AddGlobalElement(det_uboone);
  gEve->Redraw3D(kTRUE);

  auto v = gEve->GetDefaultGLViewer();
  
  //v->GetClipSet()->SetClipType(TGLClip::EType(1));
  //v->ColorSet().Background().SetColor(kMagenta+4);
  //  v->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kFALSE, 0);
  v->RefreshPadEditor(v);
  v->DoDraw();
  
 
  //auto top = gGeoManager->GetTopVolume()->FindNode("top")->GetVolume();
  //TEveGeoTopNode* trk = new TEveGeoTopNode(gGeoManager, top->FindNode("det_argoneut"));
 
  // name, a, z, rho
   /*
  TGeoMaterial *vacuum=new TGeoMaterial("vacuum",0,0,0);
  TGeoMaterial *Fe=new TGeoMaterial("Fe",55.845,26,7.87);

  // name, numed?, material
  TGeoMedium *Air=new TGeoMedium("Vacuum",0,vacuum);
  TGeoMedium *Iron=new TGeoMedium("Iron",1,Fe);
   */

  // transformation vectors/matrices
  const TVector3 beampos(-31387.58422,
			 -3316.402543,
			 -60100.2414
			 ); // proton on target in uboone coords * cm
  TEveVector beamposeve(beampos.X(), beampos.Y(), beampos.Z());  
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

  /// event script (separate these two)

 

  auto tracksXYZ = new TEveStraightLineSet("StraightLinesXYZ");
  tracksXYZ->SetLineColor(kRed);
  tracksXYZ->SetLineWidth(2);

  Double_t xdev = atan2(0.235,975.); // x deviation
  Double_t ydev = atan2(0.2  ,975.); // y deviation

  
  TFile *f  = new TFile(fstr);
  TTreeReader reader("mCP",f);
  TTreeReaderValue<TLorentzVector> Mom(reader,"Mom");
  TTreeReaderValue<TLorentzVector> Pos(reader,"Pos");
  TTreeReaderValue<Float_t> weight_decay(reader,"WeightDecay");
  TTreeReaderValue<Float_t> weight_meson(reader,"WeightMeson");
  int events = 0;
  while ( reader.Next() ) {
    events++;
    //if ( events < 100000 ) {
      //cout  << Pt << endl;
    //} else break;
    if ( UbooneAcceptanceChecker(Pos->Vect(),Mom->Vect()) > 0) {
      TVector3 Mom3 = Mom->Vect();
      Mom3 = rot*Mom3.Unit();
      TEveVector trackDir(Mom3.X(), Mom3.Y() ,Mom3.Z());
      Float_t length = 70000;
      
      TEveVector trackEnd = beamposeve + trackDir * length;
      tracksXYZ->AddLine(beampos.X(), beampos.Y(), beampos.Z(), trackEnd.fX, trackEnd.fY, trackEnd.fZ );
    }

    //argoneut
    if ( abs(atan2(Mom->Px(),Mom->Pz())) < xdev &&
	 abs(atan2(Mom->Py(),Mom->Pz())) < ydev ) {
      TVector3 Mom3 = Mom->Vect();
      Mom3 = rot*Mom3.Unit();
      TEveVector trackDir(Mom3.X(), Mom3.Y() ,Mom3.Z());
      Float_t length = 70000;
      
      TEveVector trackEnd = beamposeve + trackDir * length;
      //tracksXYZ->AddLine(beampos.X(), beampos.Y(), beampos.Z(), trackEnd.fX, trackEnd.fY, trackEnd.fZ );
      
    }

    
  }

    gEve->AddElement(tracksXYZ);
  
  // -- Load TPC geometry
  //geomGentleTPC();

  //gGeoManager->GetVolume("det_argoneut")->Draw("gl");   
  gEve->Redraw3D(kTRUE);
}
     
