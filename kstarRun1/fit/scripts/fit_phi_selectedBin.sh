#!/bin/bash

#
# FUNCTION
#
#void FitInvMass(
#  EFit        mode = kLike,
#  const char *func = "BW+POLY3",
#   Double_t    viewMin = 0.74,
#   Double_t    viewMax = 1.08,
#   
#   Bool_t      likeNorm = kFALSE,
#   Bool_t      pause = kTRUE,
#   Bool_t      data = kTRUE,
#   Int_t       PDG = 313,
#   const char *filein = "proj/kstar.root",
#   
#   Double_t    nsigmaPeak = 7.0,
#   Int_t       nrebin = 1,
#   Int_t       startBin = 0,
#   Int_t       stopBin = -1,#
#
#   Bool_t      fixGamma  = kFALSE,
#   Bool_t      fixSigma  = kFALSE,
#   Double_t    multSigma = 1.0)
for minSigma in 0.10 0.20
do  
    for min in 0.996 0.998 0.999
    do
	for max in 1.060 1.070 1.080
	do
	    for multSigma in 1.0 #0.9 1.1
	    do
		
         #FUNCTION
		for func in 'VOIGT+POLY3' 'VOIGT+POLY2'
		do
		    root -q -l -b "FitInvMass.C(kFunction, \"${func}\", ${min}, ${max}, kFALSE, kFALSE, kTRUE, 333, \"proj/phi.root\", 7.0, 1, 1, 2, kTRUE, kFALSE, ${multSigma},1-${minSigma},2-${minSigma})"
		done
		
         # LIKE SIGN
		for func in 'VOIGT+POLY1' 'VOIGT+POLY2'
		do
		    root -q -l -b "FitInvMass.C(kLike, \"${func}\", ${min}, ${max}, kFALSE, kFALSE, kTRUE, 333, \"proj/phi.root\", 7.0, 1, 1, 2, kTRUE, kFALSE, ${multSigma},1-${minSigma},2-${minSigma})"
		done
		
         # MIXING
		for func in 'VOIGT+POLY1' 'VOIGT+POLY2'
		do
		    root -q -l -b "FitInvMass.C(kMixing, \"${func}\", ${min}, ${max}, kFALSE, kFALSE, kTRUE, 333, \"proj/phi.root\", 7.0, 1, 1, 2, kTRUE, kFALSE, ${multSigma},1-${minSigma},2-${minSigma})"
		done
	    done
	done
    done
done
