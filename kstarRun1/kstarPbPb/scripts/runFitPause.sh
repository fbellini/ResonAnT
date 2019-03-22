#!/bin/bash

echo "We are now in" `pwd`

function pause(){
   read -p "$*"
}

AODN=$1
export KSTARFITMACROPATH=/Users/bellini/alice/macro/kstar/roofit
export KSTARMACROPATH=/Users/bellini/alice/macro/kstar
for inf in 1.30 #1.35 1.40
do
    for cent in 3 #0 1 2
    do 
	for func in 'POLY2'
	do  
	    for pt in  2 3 
	    do
		for min in 0.76 0.78 #0.8 # 0.76
		do
		    for max in 1.10 1.30  #1.04
		    do
			root -l "${KSTARFITMACROPATH}/fitKStar.C(\"../sub_EMnorm${inf}-1.50_aod0${AODN}_kstar.root\", 1, 1, \"${func}\", ${cent}, ${pt}, ${pt}+1, 1, ${min}, ${max}, kTRUE,kTRUE)"
			echo "Func=${func} - Min=${min}, Max=${max} - Cent=${cent} - Pt=${pt}"
			pause 'Press [Enter] key to continue...'
			mv roofit/*Pt${pt}*.* roofit/cent${cent}_Pt${pt}/.
		    done
		done
	    done
	done
    done
done
