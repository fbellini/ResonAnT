#!/bin/bash

export bg=$1 #1 = kLike, 0=kMixing 3=bin counting 4=bin counting 2gamma
export fixgamma=$2  #1=fixed 0=free gamma
export resBg=$3
export pid="tpc3s_tof3sveto"    #"comb100"  #"combined" #25s
export norm=1 #$4
sigFit="BW" #"REL" for rel-BW, "BOLZ" for rel-BWxPS, "VOIGT" for voigtian
wmin=0.5
wmax=1.5
#Nmax=0
export KSTARFITMACROPATH="/Users/bellini/alice/macro/kstar/fit/"
export KSTARFITMACRO="FitInvMass.C"

[ "$resBg" == "" ] && { resBg="poly2"; }
[ $bg -eq 3 ] && { fixgamma=1; }

if [ "$pid" == "tpc3s_tof3sveto" ]
then 
    if [ $bg -eq 0 ]  #event mixing
    then
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, 7.0, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.68, 1.06, ${norm}, 0, 1, 313, 7.0, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
#tpc3s_tof3sveto
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, 7.0, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.78, 1.04, ${norm}, 0, 1, 313, 7.0, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
#combinedbest
#root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.73, 1.05, ${norm}, 0, 1, 313, 7.0, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax}"
#root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.15, ${norm}, 0, 1, 313, 7.0, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
#
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.10, ${norm}, 0, 1, 313, 7.0, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.06, ${norm}, 0, 1, 313, 7.0, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, 7.0, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, 7.0, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.10, ${norm}, 0, 1, 313, 7.0, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"

    else #like-sign bg and all the others

	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
#tpc3s_tof3sveto
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, 7.0, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.78, 1.04, ${norm}, 0, 1, 313, 7.0, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
#combinedbest
#root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.73, 1.05, ${norm}, 0, 1, 313, 7.0, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
#root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.15, ${norm}, 0, 1, 313, 7.0, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
#
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, 7.0, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, 7.0, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, 7.0, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, 7.0, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, 7.0, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
    fi
fi
ls fit/${pid}_*/*.root > fit/fitResults.txt
root -l -b -q "/Users/bellini/alice/macro/localMergeFiles.C(\"fit/best_fit_${resBg}.root\",\"fit/fitResults.txt\")"
# ls fit/${pid}_centBin0${cent}*/*.root > fit/fitResults_cent${cent}.txt
# root -l -b -q "/Users/bellini/alice/macro/localMergeFiles.C(\"fit/best_fit_${resBg}_cent${cent}.root\",\"fit/fitResults_cent${cent}.txt\")"

if [ $bg -eq 4 ] 
then
    mv fit bcEM_norm${norm}_${sigFit}${resBg}
else
    if [ $bg -eq 3 ] 
    then
	mv fit bcLS_norm${norm}_${sigFit}${resBg}
    else  
#	[ "$resBg" == "poly2" ] && { resBg=""; } #reset for names
	if [ $bg -eq 1 ] 
	then
	    if [ $fixgamma -eq 1 ] 
	    then
		mv fit fitLS_norm${norm}_${sigFit}${resBg}_fixedW
	    else
		mv fit fitLS_norm${norm}_${sigFit}${resBg}_freeW
	    fi
	else
	    if [ $fixgamma -eq 1 ] 
	    then
		mv fit fitEM_norm${norm}_${sigFit}${resBg}_fixedW
	    else
		mv fit fitEM_norm${norm}_${sigFit}${resBg}_freeW
	    fi
	fi
    fi
#    resBg="poly2"
fi