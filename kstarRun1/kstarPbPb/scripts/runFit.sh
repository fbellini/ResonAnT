#!/bin/bash

echo "We are now in" `pwd`
mkdir proj
mkdir roofit
mkdir roofit/img
mkdir rawYields
mkdir rawYields/img

AODN=$1
export KSTARFITMACROPATH=/Users/bellini/alice/macro/kstar/roofit
export KSTARMACROPATH=/Users/bellini/alice/macro/kstar
#root -q -l -b "${KSTARMACROPATH}/projectKStar.C(\"aod0${AODN}_kstar.root\",\"RsnOut_Tof20sigma\", \"proj/kstar\", \"77\", kFALSE, kFALSE, kFALSE)"

for inf in 1.30 #1.35 1.40
do
    #root -q -l -b "${KSTARMACROPATH}/kstarAnalysis.C(\"proj\", \"proj_aod0${AODN}_kstar.root\", ${inf}, 1.50, -1, -1, kFALSE, kFALSE)"

    for func in 'POLY2' 'POLY3' 'LAND'
    do  
	for min in 0.74 # 0.76 0.78 0.8 # 0.76
	do
	    for max in 1.10 # 1.15 1.20  #1.04
	    do 
		# LIKE SIGN
	    	#root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C( \"proj/kstar_centBin0${cbin}.root\", 0, \"${func}\", 2, 14, 4, ${min}, ${max})"
		# MIXING
		root -q -l -b "${KSTARFITMACROPATH}/fitKStar.C(\"sub_EMnorm${inf}-1.50_aod0${AODN}_kstar.root\", 1, \"${func}\", -1, 2, 14, 4, ${min}, ${max},kTRUE)"
		for chi2 in 1.2 1.5 
		do
      		    root -q -l -b "${KSTARMACROPATH}/MakeSpectra.C(\"roofit\", \"fitEM_${func}_${min}-${max}_EMnorm${inf}-1.50_aod0${AODN}_kstar.root\", \"sub_EMnorm${inf}-1.50_aod0${AODN}_kstar.root\", \"${func}_range${min}-${max}\", ${chi2})"
		done
	    done
	done
    done
done
