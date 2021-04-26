# Millicharged particles search code

Millicharge (mCP) particle analysis for microboone.

This analysis aims to reproduce the methods and follow the search
strategy outlined in [[1]](https://arxiv.org/abs/1902.03246).

# Quick start guide

0. Requirements:
   * `ROOT 6.xx`.
   
   * Access to the simulation files, at least the text files. They can
     be found in the Manchester cluster or uboonegpvm* machines at
     Fermilab.

1. Use [`txt24v.cxx`](txt24v.cxx) to generate ROOT files from the
   initial text files.

2. Use [`GenerateMcps.cxx`](GenerateMcps.cxx) to generate mCP particles of
   specified mass and charge. A `bash` script is provided in
   [`run_GenerateMcps.sh`](run_GenerateMcps.sh) to generate the values to
   be plotted, so you only need to run
   >bash ./run_GenerateMcps.sh

3. Run [`FilterAccepted.cxx`](FilterAccepted.cxx) to filter the mCPs
   which pass through the uboone detector, or better yet run
   [`run_FilterAccepted.sh`](run_FilterAccepted.sh) to run all mass
   values. Other detector geometries are also possible to validate the
   flux results with different publications.

   ~~ArgoNeuT plot (not matching published results) is obtained running~~
   **Broken for the moment!**

   > root -l -b -q FilterAccepted.cxx(3)
   
   DUNE plot (matching published results) is obtained with 

   > root -l -b -q FilterAccepted.cxx(3,true)

# Detailed info

  **Not quite up to date!**

## Simulation files for pi0 and eta produced in the NuMI target

   The starting point of the analysis are pi0 and eta particles
   produced at the NuMI target. The particles are stored as text files
   containing momentum (px,py,pz,E) and position (x,y,z,t)
   information, importance weights, and an event id. They are stored
   in the Manchester cluster in
   >/afs/hep.man.ac.uk/d/lartpc-RanD/millicharge/meson_flux/

   For convenience these text files were converted into ROOT trees, a
   copy of which can be found in the Fermilab machines in
   >/uboone/data/users/arellano/millicharge/sim/*tree.root

   The ROOT and text files contain the same information. The text to
   ROOT files conversion step is done using [`txt24v.cxx`](txt24v.cxx)

## Meson -> mCP decay

   The mesons from NuMI are then made to decay using the ROOT
   `TGenPhaseSpace` class. The generated decay is Dalitz decay, and
   the mCP pair is stored as TLorentzVector objects in the `sim`
   directory.

   The decay is performed using [`GenerateMcps.cxx`](GenerateMcps.cxx),
   which takes as argument the meson type (pi0 or eta) and horn mode
   (fhc or rhc), and the mass and charge of the mCPs to be
   generated. It calls [`CrossSection.cxx`](CrossSection.cxx) and
   [`DiffCrossSection.cxx`](DiffCrossSection.cxx) to compute the
   branching fraction of the decay based on different published
   estimates; the first one generates a constant value for all events
   while the second one depends on the event and is stored as an event
   variable.

## Detector geometrical acceptance

   The number of mCPs entering the ArgoNeuT detector is calculated
   (and returned) using
   [`AcceptedArgoneut.cxx`](AcceptedArgoneut.cxx), which takes as
   argument a file (as a string to its location, and which itself
   contains the mCP mass, charge and cross sections), a method to be
   used in the calculation, and whether to use the ArgoNeuT or DUNE
   detectors.
   
   The most correct calculations are obtained using `WEIGHT` == 3
   (total cross section) or 4 (differential cross section) in the
   second argument. The factors that enter the calculation are

   > N_{accepted} = N_{entering detector geometry} * POT normalization
     * Simulation weights * Decay fraction

   * N_{entering detector geometry}: considers "box" detector, and
     particles entering through the front face of the box. ArgoNeuT
     and DUNE geometries available.

   * POT normalization: The simulated files contain 500,000 POT, so
     the normalization is either (10^20/500,000) for ArgoNeuT or
     (10^21/500,000) for DUNE.

   * Simulation weights: The NuMI mesons have an "importance weight"
     and `TGenPhaseSpace` also returns a weight, "decay weight". The
     decay weight is normalized using the number of generated events,
     `(events/sum of decay weights)`.

   * Decay fraction: This is the contentious part, different
     publications use different values. The idea is that the phase
     space is supressed at high meson masses. Currently I'm using the
     branching shown in eqs. (1) and (2) in
     [[2]](https://journals.aps.org/prd/pdf/10.1103/PhysRevD.100.015043)
     and this correctly reproduces the DUNE estimate values in that
     same publication. A differential decay fraction is found in
     eq. (2) in [[3]](https://arxiv.org/pdf/2010.07941.pdf)

## Plotting

   Finally, plotting is done using
   [`PlotAcceptedVariables.cxx`](PlotAcceptedVariables.cxx). This can
   either recalculate the decay fraction and save the data points in
   `hist/` or it can just read from these saved data points and just
   plot by entering a negative value in the `WEIGHT` argument. The
   plots are saved in `img/`


# References

1. *Millicharged Particles in Liquid Argon Neutrino Experiments* Roni
Harnik, Zhen Liu, Ornella Palamara, https://arxiv.org/abs/1902.03246

2. *Proton Fixed-Target Scintillation Experiment to Search for
Minicharged Particles* Kevin J. Kelly, Yu-Dai Tsai, Phys. Rev. D 100,
015043 (2019) (NOTE: equation (2) is different in Phys Rev D and arXiv
version)

3. *FORMOSA: Looking Forward to Millicharged Dark Sectors*, Saeid
Foroughi-Abari, Felix Kling, Yu-Dai Tsai,
https://arxiv.org/pdf/2010.07941.pdf