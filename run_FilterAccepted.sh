#!/usr/bin/bash

# pi0 accessible
for meson in pi0 eta
do
    for horn in fhc #rhc
    do
	for mcpmass in 0.010 0.020 0.030 0.050 0.060
	do
	    time root -l -b -q FilterAccepted.cxx\(\"'sim/mCP_q_0.010_m_'$mcpmass'_fhc_'$meson's.root'\"\)
	done
    done
done

# eta accessible
for horn in fhc # rhc
do
    for mcpmass in 0.100 0.200 0.250
    do
	time root -l -b -q FilterAccepted.cxx\(\"'sim/mCP_q_0.010_m_'$mcpmass'_fhc_'$meson's.root'\"\)
    done
done
