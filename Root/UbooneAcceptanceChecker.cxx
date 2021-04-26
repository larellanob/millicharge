Double_t intersects(const TVector3 &,
		    const TVector3 &,
		    const TVector3 &,
		    const TVector3 &);

Double_t UbooneAcceptanceChecker(TVector3 pos, TVector3 mom )
{
  // not accepted until proven otherwise
  bool accepted = false;
  Double_t distance_traveled = 0;
  //tpcactive dimensions in det coords:
  //X 10, 246.35      --> 
  //Y -105.53, 107.47 --> 
  //Z [10.1, 1026.8]  --> 
  
  // detector half-dimensions in cm
  const TVector3 det_half_dims(0.5*(246.35-10.),0.5*(107.47+105.53),.5*(1026.8-10.1));
  
  // ~= {0.,0.,46536.3525}; // detector centre in det coords
  const TVector3 det_centre = TVector3(10.,-105.53,10.1) + det_half_dims;
  
  /* beam rotation, from /uboone/app/users/zarko/windowRotation/correctWindow.c */
  const TVector3 beampos(-31387.58422,
			 -3316.402543,
			 -60100.2414); // proton on target in uboone coords * cm
  
  
  // defines TRotation c++ 'lambda' named rot
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
  
  pos = (rot*pos)+beampos;
  mom = rot*mom;
  // these should be now in uboone frame of reference
  
  return intersects(pos,mom,det_centre,det_half_dims);
  
  
  //return false;
}



Double_t intersects(const TVector3 &orig, // origin
		const TVector3 &dir, // direction of the mcp?
		const TVector3 &det_centre, // center of the detector
		const TVector3 &det_half_dims//half the size of the detector?
		)
//		  double* lambdas)
{
  const TVector3& unit_dir = dir.Unit();// unitary vector for direction
  int n_intersects = 0; // number of interections?
  unsigned int ilam = 0; // index of the lambda

  TVector3 intersection1(0,0,0);
  TVector3 intersection2(0,0,0);
  
  for(int coord = 0; coord < 3; ++coord) {
    for(int side = -1; side < 2; side += 2) {
      // two sides, -1 and 1
      //(left and right/top and bottom/front and back) respect to centre
      const double plane = det_centre[coord] + side * det_half_dims[coord]; // centre +- half size == wall
      const double lambda = (plane - orig[coord])/unit_dir[coord];
      if(lambda < 0) continue; // no backwards-going scalars
      bool intersects_planes[2] = {false, false};
      unsigned int iplane = 0;
      for(int other_coord = 0; other_coord < 3; ++other_coord) {
	if(other_coord == coord) continue;
	const double oth_plane = lambda * unit_dir[other_coord] + orig[other_coord];
	if(std::abs(oth_plane - det_centre[other_coord]) < det_half_dims[other_coord]) {
	  intersects_planes[iplane]=true;
	}
	iplane++;
      }
      if(intersects_planes[0] && intersects_planes[1]) {
	n_intersects++;
	if ( intersection1 != TVector3(0,0,0) ) {
	  intersection2 = TVector3(
				   orig[0]+lambda*unit_dir[0],
				   orig[1]+lambda*unit_dir[1],
				   orig[2]+lambda*unit_dir[2]
				   );
	} else if ( intersection1 == TVector3(0,0,0) ) {
	  intersection1 = TVector3(
				   orig[0]+lambda*unit_dir[0],
				   orig[1]+lambda*unit_dir[1],
				   orig[2]+lambda*unit_dir[2]
				   );
	}
	
	/*
	  if(ilam < 2) {
	  lambdas[ilam++] = lambda;
	  }
	*/
      }
    }
  }

  // verbose debug
  /*
  if ( n_intersects == 1 ) {
    cout << "WHAT? intersecting only one plane??" << endl;
  } else if ( n_intersects > 2 ) {
    cout << "WHAT? intersecting more than 2 planes??" << endl;
  } else if ( n_intersects == 2 ) {
    cout << "INTERSECTING 2, nothing to see here" << endl;
  }
  */
  /*
  if ( n_intersects == 2 ) {
    
    if ( (intersection1-intersection2).Mag()  > 500. ) {
      intersection1.Print();
      intersection2.Print();
      
      cout << "Travelled distance " << (intersection1-intersection2).Mag() << endl;
    }
  }
  */
  /*
  */

  // returns distance travelled inside the detector
  // if it doesn't reach the detector returns 0
  return (intersection1-intersection2).Mag();
  //return n_intersects >= 2; // returns a bool
  // true if n_int >=2 , false if <
}

