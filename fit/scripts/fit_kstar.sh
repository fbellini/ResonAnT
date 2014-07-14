#!/bin/bash
startptbin=0

export KSTARFITMACROPATH=/Users/bellini/alice/macro/kstar/fit
for cbin in 0 1 2 3 4 
do
    for min in 0.72 #0.74 0.8 #0.72 0.76 #0.74 0.76
    do
	for max in 1.10 1.14 #1.02 1.04 1.06 1.10 
	do
            for func in  'BW+POLY2' 'BW+POLY3' #'BOLZ+POLY2' 'BOLZ+POLY3'  'BW+POLY1'
	    do
	        # LIKE SIGN
	    	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C( \"proj/tof2s_centBin0${cbin}.root\", kLike, \"${func}\", ${min}, ${max}, kTRUE, kFALSE, kTRUE, 313, 7.0, 1, ${startptbin}, 14, ${cbin}, kFALSE, kTRUE)"
		# MIXING
		#root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin0${cbin}.root\", kMixing, \"${func}\", ${min}, ${max}, kTRUE, kFALSE, kTRUE, 313, 7.0, 1, ${startptbin}, -1, kFALSE, kTRUE)"
	    done
	done
    done
done



