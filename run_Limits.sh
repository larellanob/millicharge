#!/usr/bin/bash

#detector=argoneut_published
#detector=argoneut
detector=uboone

root -l -b -q Root/ClearParameterSpaceTH2F.cxx\(\"$detector\",true\)

for th in 0.1 0.6 0.8 1.0
do
    for meson in pi0 eta
    do
	for mcpmass in 0.010 0.020 0.030 0.050 0.060
	do
	    time root -l -b -q GenerateSensitivityHistograms.cxx\(\"'sim/mCP_uboone_q_0.010_m_'$mcpmass'_fhc_'$meson's.root'\",$th\)
	done
    done
    
    for mcpmass in 0.100 0.200 0.250
    do
	time root -l -b -q GenerateSensitivityHistograms.cxx\(\"'sim/mCP_uboone_q_0.010_m_'$mcpmass'_fhc_etas.root'\",$th\)
    done
done
## at this point call PlotPublishedLimits
