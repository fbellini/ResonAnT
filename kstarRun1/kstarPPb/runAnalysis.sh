#!/bin/bash

pAworkdir="/Users/fbellini/alice/resonances/kstar_pA5.02TeV/LF_pPb17-21/"
cutCode="2424"
PIDdir="tpc2s_tof3sveto"
fitmacrodir="${ASD}/ResonAnT"
export TRAINOUT=$1
export CUTSET=$2
export newworkdir="cut${CUTSET}_${PIDdir}"

cd ${pAworkdir}
#make projections and subtractions for roofit macros
root -l -b -q ${fitmacrodir}/projectKStar.C"(\"${TRAINOUT}\",\"RsnOut_cut${CUTSET}\",\"cut${CUTSET}_${PIDdir}\",\"${cutCode}\", kTRUE)"
mkdir ${newworkdir}
mv *cut${CUTSET}*.root ${newworkdir}
cd ${newworkdir}

#${pAworkdir}/runAna.sh ${detector}

#run root fit macro - roofit to be added
#run poly2 fits for different bg in best ranges for each 
#and for each variation of the width
#cd ${pAworkdir}/${detector}_ANA/${PIDdir}
#for bgtype in 1 #0 3
#do
#    for fixedwidth in 1 #0
#    do
#	${fitmacrodir}/fit_kstar_best_${detector}.sh ${bgtype} ${fixedwidth} 'poly2'
#    done
#done
