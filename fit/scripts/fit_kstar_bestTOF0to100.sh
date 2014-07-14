#!/bin/bash

export KSTARFITMACROPATH=/Users/bellini/alice/macro/kstar/fit
export bg=$1          #1 = kLike, 0=kMixing
export fixgamma=$2    #1 = fixed, 0=free width
export normLike=1     

root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  2, 3, 0, ${fixgamma}, kTRUE)"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  3, 4, 0, ${fixgamma}, kTRUE)"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  4, 5, 0, ${fixgamma}, kTRUE)"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  5, 6, 0, ${fixgamma}, kTRUE)"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  6, 7, 0, ${fixgamma}, kTRUE)"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  7, 8, 0, ${fixgamma}, kTRUE)"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  8, 9, 0, ${fixgamma}, kTRUE)"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  9, 10, 0, ${fixgamma}, kTRUE)"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 10, 11, 0, ${fixgamma}, kTRUE)"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 11, 12, 0, ${fixgamma}, kTRUE)"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 12, 13, 0, ${fixgamma}, kTRUE)"

#different fit range for LS and EM
if [ $bg -eq 1 ]
then   #LS
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 13, 14, 0, ${fixgamma}, kTRUE)"
else   #mixing
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof2s_centBin00.root\", ${bg}, \"BW+POLY2\", 0.76, 1.04, ${normLike}, 0, 1, 313, 7.0, 1, 13, 14, 0, ${fixgamma}, kTRUE)"
fi


ls fit/tof2s_*/*.root > fit/fitResults.txt
root -l -b -q "/Users/bellini/alice/macro/localMergeFiles.C(\"fit/best_fit_poly2.root\",\"fit/fitResults.txt\")"

if [ $bg -eq 1 ] 
then
    if [ $normLike -eq 1 ]
    then
	if [ $fixgamma -eq 1 ] 
	then
#	    mv fit fit_LSnorm13_bestRange_fixedW
	    mv fit fit_bestRange_fixedW
	else
#	    mv fit fit_LSnorm13_bestRange_freeW
	    mv fit fit_bestRange_freeW
	fi
    else
	mv fit fit_LSnotNorm_bestRange_fixedW
    fi
else
    if [ $fixgamma -eq 1 ] 
    then
	mv fit fit_EM_bestRange_fixedW
    else
	mv fit fit_EM_bestRange_freeW
    fi
fi
# if [ ! -d "fit/rawYields" ]; then
#  mkdir fit/rawYields
#  cd fit
#  root -l -b -q "/Users/bellini/alice/macro/kstar/MakeRawSpectra.C(1,\"best_fit_poly2.root\",\"Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1818_tof2s_train215-216.root\",\"_BWPOLY2\",\"fit: BW+poly2\", 10.0))"
# fi
