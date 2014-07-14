#!/bin/bash

detector="tpc"
PIDdir="2s"
pAworkdir="/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/"
fitmacrodir="${ASD}/kstar/fit"

cd ${pAworkdir}
ls -lath

#make projections and subtractions for roofit macros
${pAworkdir}/runAna.sh ${detector}

#run root fit macro - roofit to be added
#run poly2 fits for different bg in best ranges for each 
#and for each variation of the width
cd ${pAworkdir}/${detector}_ANA/${PIDdir}
for bgtype in 1 #0 3
do
    for fixedwidth in 1 #0
    do	
	${fitmacrodir}/fit_kstar_best_${detector}.sh ${bgtype} ${fixedwidth} 'poly2'
    done
done
