#!/bin/bash

main()
    {
	export bg=$1 #1 = kLike, 0=kMixing 3=bin counting 4=bin counting 2gamma
	export fixgamma=$2  #1=fixed 0=free gamma
	export resBg=$3
	export pid="tpc2s_tof3sveto"    #"comb100"  #"combined" #25s 
	export norm=$4
	sigFit="BW" #"REL" for rel-BW, "BOLZ" for rel-BWxPS, "VOIGT" for voigtian
	wmin=0.50
	wmax=1.50
	mrange=7.0  #default = 7.0
#Nmax=0
	export FITPATH="/Users/bellini/alice/macro/kstar/fit/FitInvMass.C"
#	export KSTARFITMACRO="FitInvMass.C"
	
	#set default background function to poly2
	[ "$resBg" == "" ] && { resBg="poly2"; }

	#set fixed gamma for bin counting 
	[ $bg -eq 3 ] && { fixgamma=1; }
		
	if [[ "$pid" == "tpc3s_tof3sveto" ]]; then
	  runTpc3sTof3sVeto
	fi
	
	if [[ "$pid" == "tpc2s_tof3sveto" ]]; then
	    if [[ "$resBg" == "sys" ]]; then
    		runTpc2sTof3sVetoFcnSys		
	    else
		runTpc2sTof3sVeto
	    fi
       	fi
	
	if [[ "$pid" == "tpc2s_tof4sveto" ]]; then
	    runTpc2sTof4sVeto
	fi
	
	if [[ "$pid" == "tpc3s_tof4sveto" ]];	then 
	    runTpc3sTof4sVeto
	fi
	if [[ "$pid" == "combinedBest" ]];	then 
	    runCombBest
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
	fi
	
       	# if [[ "$pid" == "tpc2s_tof3sveto" ]]; then
	#     runTpc2sTof3sVetoFitSys
	# fi
    } #end main
    
runTpc3sTof3sVeto(){
    
    if [ $bg -eq 0 ]  #event mixing
    then
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.68, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.02, ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 17, 18, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 18, 19, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 19, 20, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 20, 21, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 21, 22, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 22, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 23, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 24, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"

    else #like-sign bg and all the others

	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.78, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 17, 18, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 18, 19, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 19, 20, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 20, 21, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 21, 22, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 22, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 23, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 24, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"
    fi

}

runTpc2sTof3sVetoFcnSys()
{
    if [ $bg -eq 0 ]  #event mixing
    then
	echo "Event mixing"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly3\", 0.70, 1.20, ${norm}, 0, 1, 313, ${mrange}, 1,  0,  1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly3\", 0.74, 1.20, ${norm}, 0, 1, 313, ${mrange}, 1,  1,  2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly1\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  2,  8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly3\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  8,  9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly3\", 0.74, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  9, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"  
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly3\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 11, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly3\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 24, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"
    else
	echo "Like-sign"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly1\", 0.75, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  0, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly3\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly3\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  6, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly3\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly3\", 0.74, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  9, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly1\", 0.75, 1.05, ${norm}, 0, 1, 313, ${mrange}, 1, 14, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+poly1\", 0.66, 1.20, ${norm}, 0, 1, 313, ${mrange}, 1, 23, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"
    fi
}

runTpc2sTof3sVeto()
{
    if [ $bg -eq 0 ]  #event mixing
    then
       	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	if [ $fixgamma -eq 1 ]
	then
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.68, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	else
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	fi
	
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
	if [ $fixgamma -eq 1 ]
	then
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.02, ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	else
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	fi
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
	if [ $fixgamma -eq 1 ]
	then
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	else
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	fi
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
	if [ $fixgamma -eq 1 ]
	then
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
	else
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.12, ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
	fi
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 17, 18, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 18, 19, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 19, 20, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 20, 21, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 21, 22, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 22, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 23, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 24, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"
    else #like-sign bg and all the others
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.78, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 17, 18, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 18, 19, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 19, 20, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 20, 21, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 21, 22, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 22, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 23, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 24, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"
    fi
}

runTpc2sTof3sVetoFitSys()
{
    deltax=0.02
    for fixgamma in 0 
    do
	for left in -1 0 1
	do
	    for right in -1 0 1
	    do
		if [ $bg -eq 0 ]  #event mixing
		then
       		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.68+(${deltax}*${left}), 1.06+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76+(${deltax}*${left}), 1.02+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74+(${deltax}*${left}), 1.06+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74+(${deltax}*${left}), 1.04+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76+(${deltax}*${left}), 1.04+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74+(${deltax}*${left}), 1.06+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74+(${deltax}*${left}), 1.06+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 17, 18, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 18, 19, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 19, 20, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 20, 21, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 21, 22, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 22, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76+(${deltax}*${left}), 1.04+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 23, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76+(${deltax}*${left}), 1.04+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 24, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    
		    ls fit/${pid}_*/*.root > fit/fitResults.txt
		    root -l -b -q "/Users/bellini/alice/macro/localMergeFiles.C(\"fit/best_fit_${resBg}.root\",\"fit/fitResults.txt\")"
# ls fit/${pid}_centBin0${cent}*/*.root > fit/fitResults_cent${cent}.txt
# root -l -b -q "/Users/bellini/alice/macro/localMergeFiles.C(\"fit/best_fit_${resBg}_cent${cent}.root\",\"fit/fitResults_cent${cent}.txt\")"

		    if [ $fixgamma -eq 1 ] 
		    then
			mv fit fitEM_norm${norm}_${sigFit}${resBg}_fixedW_l${left}_r${right}
		    else
			mv fit fitEM_norm${norm}_${sigFit}${resBg}_freeW_l${left}_r${right}
		    fi
		else #like-sign bg and all the others
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72+(${deltax}*${left}), 1.08+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72+(${deltax}*${left}), 1.08+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72+(${deltax}*${left}), 1.08+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.78+(${deltax}*${left}), 1.04+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76+(${deltax}*${left}), 1.04+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 17, 18, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 18, 19, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 19, 20, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 20, 21, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 21, 22, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70+(${deltax}*${left}), 1.10+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 22, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76+(${deltax}*${left}), 1.04+(${deltax}*${right}), ${norm}, 0, 1, 313, ${mrange}, 1, 23, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
		    ls fit/${pid}_*/*.root > fit/fitResults.txt
		    root -l -b -q "/Users/bellini/alice/macro/localMergeFiles.C(\"fit/best_fit_${resBg}.root\",\"fit/fitResults.txt\")"
# ls fit/${pid}_centBin0${cent}*/*.root > fit/fitResults_cent${cent}.txt
# root -l -b -q "/Users/bellini/alice/macro/localMergeFiles.C(\"fit/best_fit_${resBg}_cent${cent}.root\",\"fit/fitResults_cent${cent}.txt\")"		
		    if [ $fixgamma -eq 1 ] 
		    then
			mv fit fitLS_norm${norm}_${sigFit}${resBg}_fixedW_l${left}_r${right}
		    else
			mv fit fitLS_norm${norm}_${sigFit}${resBg}_freeW_l${left}_r${right}
		    fi		
		fi
	    done 
	done
    done
}


runTpc2sTof4sVeto(){
    
    if [ $bg -eq 0 ]  #event mixing
    then
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.68, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.02, ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 17, 18, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 18, 19, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 19, 20, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 20, 21, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 21, 22, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 22, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 23, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 24, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"

    else #like-sign bg and all the others

	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.78, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 17, 18, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 18, 19, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 19, 20, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 20, 21, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 21, 22, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 22, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 23, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 24, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"
    fi

}



runTpc3sTof4sVeto()
    {
	if [ $bg -eq 0 ]  #event mixing
	then
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.68, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.02, ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 17, 18, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 18, 19, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 19, 20, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 20, 21, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 21, 22, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 22, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 23, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 24, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"

	else #like-sign bg and all the others

	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.78, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 17, 18, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 18, 19, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 19, 20, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 20, 21, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 21, 22, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 22, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 23, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
	    root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 24, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"
	fi
    }

runCombBest(){
    
    if [ $bg -eq 0 ]  #event mixing
    then
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.68, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.73, 1.05, ${norm}, 0, 1, 313, ${mrange}, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax}"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.15, ${norm}, 0, 1, 313, ${mrange}, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.02, ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.06, ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.74, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 17, 18, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 18, 19, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 19, 20, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 20, 21, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 21, 22, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 22, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 23, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 24, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"

    else #like-sign bg and all the others

	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  0, 1, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  1, 2, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1,  2, 3, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.73, 1.05, ${norm}, 0, 1, 313, ${mrange}, 1,  3, 4, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.15, ${norm}, 0, 1, 313, ${mrange}, 1,  4, 5, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  5, 6, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  6, 7, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.72, 1.08, ${norm}, 0, 1, 313, ${mrange}, 1,  7, 8, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.78, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  8, 9, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1,  9, 10, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 10, 11, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 11, 12, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 12, 13, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 13, 14, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 14, 15, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 15, 16, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 16, 17, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 17, 18, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 18, 19, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 19, 20, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 20, 21, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 21, 22, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.70, 1.10, ${norm}, 0, 1, 313, ${mrange}, 1, 22, 23, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 23, 24, 100, ${fixgamma}, ${wmin}, ${wmax})"
	root -q -l -b "${FITPATH}(\"proj/${pid}_centBin00.root\", ${bg}, \"${sigFit}+${resBg}\", 0.76, 1.04, ${norm}, 0, 1, 313, ${mrange}, 1, 24, 25, 100, ${fixgamma}, ${wmin}, ${wmax})"
    fi

}




main $@
