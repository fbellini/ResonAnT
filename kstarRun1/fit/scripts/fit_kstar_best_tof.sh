#!/bin/bash

export bg=$1 #1 = kLike, 0=kMixing 3=bin count 4=bin count 2 gamma
export fixgamma=$2
export resBg=$3

export tofpid="2s"  #25s
export normLike=1      #0 for syst
wmin=0.5
wmax=1.5

export KSTARFITMACROPATH=/Users/bellini/alice/macro/kstar/fit

[ "$resBg" == "" ] && { resBg="poly2"; }
[ $bg -eq 3 ] && { fixgamma=1; }
#################################
# Centrality 0 
#################################
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  2, 3, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  3, 4, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  4, 5, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  5, 6, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  6, 7, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  7, 8, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  8, 9, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1,  9, 10, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 10, 11, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 11, 12, 0, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 12, 13, 0, ${fixgamma},  ${wmin}, ${wmax})"

#different fit range for LS and EM
if [ $bg -eq 1 ]
then   #LS
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 13, 14, 0, ${fixgamma},  ${wmin}, ${wmax})"
else   #mixing
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin00.root\", ${bg}, \"BW+${resBg}\", 0.76, 1.04, ${normLike}, 0, 1, 313, 7.0, 1, 13, 14, 0, ${fixgamma},  ${wmin}, ${wmax})"
fi


#################################
# Centrality 1 - LIKE SIGN
#################################
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 2, 3, 1, ${fixgamma},  ${wmin}, ${wmax})"

if [ $bg -eq 1 ]
then   #LS
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.20, ${normLike}, 0, 1, 313, 7.0, 1, 3, 4, 1, ${fixgamma},  ${wmin}, ${wmax})"
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.20, ${normLike}, 0, 1, 313, 7.0, 1, 4, 5, 1, ${fixgamma},  ${wmin}, ${wmax})"
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 5, 6, 1, ${fixgamma},  ${wmin}, ${wmax})"
else #Mixing
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.04, ${normLike}, 0, 1, 313, 7.0, 1, 3, 4, 1, ${fixgamma},  ${wmin}, ${wmax})"
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.04, ${normLike}, 0, 1, 313, 7.0, 1, 4, 5, 1, ${fixgamma},  ${wmin}, ${wmax})"
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.06, ${normLike}, 0, 1, 313, 7.0, 1, 5, 6, 1, ${fixgamma},  ${wmin}, ${wmax})"
fi
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 6, 7, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 7, 8, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.06, ${normLike}, 0, 1, 313, 7.0, 1, 8, 9, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1, 9, 10, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1,10, 11, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1,11, 12, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1,12, 13, 1, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin01.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1,13, 14, 1, ${fixgamma},  ${wmin}, ${wmax})"

#################################
# Centrality 2 - LIKE SIGN
#################################
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 2, 3, 2, ${fixgamma},  ${wmin}, ${wmax})"

if [ $bg -eq 1 ]
then   #LS
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 3, 4, 2, ${fixgamma},  ${wmin}, ${wmax})"
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 4, 5, 2, ${fixgamma},  ${wmin}, ${wmax})"
else #mixing
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.06, ${normLike}, 0, 1, 313, 7.0, 1, 3, 4, 2, ${fixgamma},  ${wmin}, ${wmax})"
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.06, ${normLike}, 0, 1, 313, 7.0, 1, 4, 5, 2, ${fixgamma},  ${wmin}, ${wmax})"
fi

root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 5, 6, 2, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 6, 7, 2, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 7, 8, 2, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 8, 9, 2, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 9, 10, 2, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 10,11, 2, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 11,12, 2, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 12,13, 2, ${fixgamma},  ${wmin}, ${wmax})"

if [ $bg -eq 1 ]
then   #LS
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 13,14, 2, ${fixgamma},  ${wmin}, ${wmax})"
else
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin02.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 13,14, 2, ${fixgamma},  ${wmin}, ${wmax})"
fi

#################################
# Centrality 3 - LIKE SIGN
#################################
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.70, 1.06, ${normLike}, 0, 1, 313, 7.0, 1, 2, 3, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 3, 4, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 4, 5, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 5, 6, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1, 6, 7, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 7, 8, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 8, 9, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1, 9, 10, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1, 10,11, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1,11,12, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1,12,13, 3, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1,13,14, 3, ${fixgamma},  ${wmin}, ${wmax})"
#root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin03.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1,13,14, 3, ${fixgamma},  ${wmin}, ${wmax})"

#################################
# Centrality 4 - LIKE SIGN
#################################
if [ $bg -eq 1 ] 
then #ls
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 2, 3, 4, ${fixgamma},  ${wmin}, ${wmax})"
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 3, 4, 4, ${fixgamma},  ${wmin}, ${wmax})"
else #mixing
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.68, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 2, 3, 4, ${fixgamma},  ${wmin}, ${wmax})"
    root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.74, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 3, 4, 4, ${fixgamma},  ${wmin}, ${wmax})"
fi
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 4, 5, 4, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.10, ${normLike}, 0, 1, 313, 7.0, 1, 5, 6, 4, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1, 6, 7, 4, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1, 7, 8, 4, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1, 8, 9, 4, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1, 9, 10, 4, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1, 10,11, 4, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1,11,12, 4, ${fixgamma},  ${wmin}, ${wmax})"
root -q -l -b "${KSTARFITMACROPATH}/FitInvMass.C(\"proj/tof${tofpid}_centBin04.root\", ${bg}, \"BW+${resBg}\", 0.72, 1.14, ${normLike}, 0, 1, 313, 7.0, 1,12,13, 4, ${fixgamma},  ${wmin}, ${wmax})"

ls fit/tof${tofpid}*/*.root > fit/fitResults.txt
root -l -b -q "/Users/bellini/alice/macro/localMergeFiles.C(\"fit/best_fit_${resBg}.root\",\"fit/fitResults.txt\")"

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
# if [ ! -d "fit/rawYields" ]; then
#  mkdir fit/rawYields
#  cd fit
#  root -l -b -q "/Users/bellini/alice/macro/kstar/MakeRawSpectra.C(1,\"best_fit_poly2.root\",\"Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s/_sum_EMnorm1.30-1.50_cut1818_tof2s_train215-216.root\",\"_BW${resBg}\",\"fit: BW+${resBg}\", 10.0))"
# fi
