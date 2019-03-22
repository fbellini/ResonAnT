#!/bin/bash

export bg=$1 #1 = kLike, 0=kMixing 3=bin counting 4=bin counting 2gamma
export fixgamma=$2  #1=fixed 0=free gamma
export resBg=$3
export tpcpid="tpc3s_tof3sveto"   #"comb_binA"  #"combined" #25s
export normLike=1
wmin=0.5
wmax=1.5
#Nmax=0
export KSTARFITMACROPATH="/Users/bellini/alice/macro/kstar/fit/"
export KSTARFITMACRO="FitInvMass.C"

[ "$resBg" == "" ] && { resBg="poly2"; }
[ $bg -eq 3 ] && { fixgamma=1; }

###########################
# Centrality 0 - LIKE SIGN
###########################
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  0,  3, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.05, ${normLike}, 0, 1, 313, 7.0, 1,  3,  4, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.06, ${normLike}, 0, 1, 313, 7.0, 1,  4,  5, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  5, 10, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.06, ${normLike}, 0, 1, 313, 7.0, 1, 10, 12, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 12, 13, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.06, ${normLike}, 0, 1, 313, 7.0, 1, 13, 14, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 14, 19, 0, ${fixgamma},  ${wmin}, ${wmax})"

###########################
# Centrality 1 - LIKE SIGN
###########################
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.06, ${normLike}, 0, 1, 313, 7.0, 1,  0,  1, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  1,  3, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.66, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  3,  4, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.06, ${normLike}, 0, 1, 313, 7.0, 1,  4,  5, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.68, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  5,  6, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  6, 11, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 11, 12, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.06, ${normLike}, 0, 1, 313, 7.0, 1, 12, 13, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 13, 19, 1, ${fixgamma},  ${wmin}, ${wmax})"

###########################
# Centrality 2 - LIKE SIGN
###########################
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  0,  1, 2, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.12, ${normLike}, 0, 1, 313, 7.0, 1,  1,  2, 2, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.68, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  2,  4, 2, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  4, 19, 2, ${fixgamma},  ${wmin}, ${wmax})"

###########################
# Centrality 3 - LIKE SIGN
###########################
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  0,  2, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.68, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  2,  4, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  4,  7, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.12, ${normLike}, 0, 1, 313, 7.0, 1,  7,  8, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  8, 13, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.76, 1.06, ${normLike}, 0, 1, 313, 7.0, 1, 13, 19, 3, ${fixgamma},  ${wmin}, ${wmax})"

###########################
# Centrality 4 - LIKE SIGN
###########################
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  0,  8, 4, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  8,  9, 4, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.75, 1.05, ${normLike}, 0, 1, 313, 7.0, 1,  9, 10, 4, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 10, 11, 4, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/${KSTARFITMACRO}(\"proj/${tpcpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 11, 19, 4, ${fixgamma},  ${wmin}, ${wmax})"

ls fit/${tpcpid}_*/*.root > fit/fitResults.txt
root -l -b -q "/Users/bellini/alice/macro/localMergeFiles.C(\"fit/best_fit_${resBg}.root\",\"fit/fitResults.txt\")"
for cent in 0 1 2 3 4 
do
    ls fit/${tpcpid}_centBin0${cent}*/*.root > fit/fitResults_cent${cent}.txt
    root -l -b -q "/Users/bellini/alice/macro/localMergeFiles.C(\"fit/best_fit_${resBg}_cent${cent}.root\",\"fit/fitResults_cent${cent}.txt\")"
done


if [ $bg -eq 4 ] 
then
    mv fit binCounting_2gamma_${resBg}
else
    if [ $bg -eq 3 ] 
    then
	mv fit binCounting_${resBg}
    else  
	[ "$resBg" == "poly2" ] && { resBg=""; } #reset for names
	if [ $bg -eq 1 ] 
	then
	    if [ $normLike -eq 1 ]
	    then
		if [ $fixgamma -eq 1 ] 
		then
#	    mv fit fit_LSnorm13_bestRange_fixedW
		    mv fit fit_bestRange${resBg}_fixedW
		else
#	    mv fit fit_LSnorm13_bestRange_freeW
		    mv fit fit_bestRange${resBg}_freeW
		fi
	    else
		mv fit fit_LSnotNorm_bestRange${resBg}_fixedW
	    fi
	else
	    if [ $fixgamma -eq 1 ] 
	    then
		mv fit fit_EM_bestRange${resBg}_fixedW
	    else
		mv fit fit_EM_bestRange${resBg}_freeW
	    fi
	fi
    fi
    resBg="poly2"
fi

# if [ ! -d "fit/rawYields" ]; then
#  mkdir fit/rawYields
#  cd fit
#  root -l -b -q "/Users/bellini/alice/macro/kstar/MakeRawSpectra.C(1,\"best_fit_poly2.root\",\"Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TPC_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1717_tpc2s_train215-216.root\",\"_BW${resBg}\",\"fit: BW+poly2\", 10.0))"
# fi
