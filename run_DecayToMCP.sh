#!/usr/bin/bash

# pi0 accessible
for meson in pi0 eta
do
    for horn in fhc #rhc
    do
	for mcpmass in 0.01 0.02 0.03 0.05 0.06
	do
	    time root -l -b -q DecayToMCP.cxx\(\"$meson\",\"$horn\",$mcpmass\)
	done
    done
done

# eta accessible
for horn in fhc # rhc
do
    for mcpmass in 0.1 0.2 0.25
    do
	time root -l -b -q DecayToMCP.cxx\(\"eta\",\"$horn\",$mcpmass\)
    done
done
