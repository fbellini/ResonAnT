#!/bin/bash
export KSTARFITMACROPATH=/Users/bellini/alice/macro/kstar/fit
for cbin in 0 #1 2 3 #4 
do
    for min in 0.70 0.72 0.75 # 0.76
    do
	for max in 1.1 1.3 1.5 #1.04
	do
            for func in 'BW+POLY1' 'BW+POLY2' 'BW+POLY3'
	    do
	        # LIKE SIGN
	    	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C( \"proj/kstar_centBin0${cbin}.root\", kLike, \"${func}\", ${min}, ${max}, kTRUE, kFALSE, kTRUE, 313, 7.0, 1, 0, -1, kTRUE, kTRUE)"
		# MIXING
		root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/kstar_centBin0${cbin}.root\", kMixing, \"${func}\", ${min}, ${max}, kTRUE, kFALSE, kTRUE, 313, 7.0, 1, 0, -1, kTRUE, kTRUE)"
	    done
	done
    done
done



