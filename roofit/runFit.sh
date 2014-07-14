#!/bin/bash

#help function
function show_help()
{
    echo "runFit [abchlpPrw] -f <function> -x <x_mix> -X <x_max> -i <file>"   
    echo " Function, fit range and input file are MANDATORY.
           Options: 
                   -a             all centrality bins and pt bins  
                   -b <bg>        set specified background <bg> between EM=1, LS=2 (default)
                   -c <cent>      select centrality <cent>
                   -f <function>  sets fit function to <function>
                   -h             shows this help function
                   -i <file>      sets input file to <file>
                   -x <x_min>     sets minimum of fit range to <x_min>
                   -X <x_max>     sets minimum of fit range to <x_max>
                   -l             enables max likelihood fit (default is chi2)
                   -m             run multiple ranges fits
                   -p <pt_min>    select min pt bin <pt_min> to fit
                   -P <pt_max>    select max pt bin <pt_max> to fit 
                   -r             enable rebin by factor 2
                   -w             enable variable BW width and promts for limiting values"    
}

[ $# -eq 0 ] && { show_help; exit 1; }  		

#default variables
BG=2
ISREBIN=0
ISWFIXED=1
ISCHI2=1
DIROUT="./roofit"
RUNMULTI=0
FUNCTION=
XMIN=
XMAX=
FILEIN=   #"_sum_EMnorm1.30-1.50_cut1818_tof3s_train215-216.root"
selectedcent=-1
startpt=0
stoppt=13
splitc=1
minW=0.0
maxW=1.0
SAVER=
#function to pause
function pause(){
    read "$*"
}


#function to define default settings
function default_fit()
{
    echo "Default settings applied (all pt and centrality, LS bg, chi2 fit, fixed W, no rebin)"
    BG=2
    ISREBIN=0
    ISCHI2=1
    DIROUT="./fit"
    ISWFIXED=1
    fit_all
}

#function to fit all bins - change extremes according to binning
function fit_all()
{
#   selectedcent=-1
#   startpt=1
#   stoppt=13
   splitc=4
   echo "Fit all bins (pt ${startpt}-${stoppt})"
}

#function to read from input limits for W if not fixed
function read_Wlimits()
{ 
    echo "Please insert min value for BW width"
    read minW
    echo "Please insert max value for BW width"
    read maxW
}

#reads options
while getopts ":ab:c:f:hi:x:X:lmp:P:rw" option; do

    case $option in
	a)  fit_all
	    ;;
	b) selectedbg=$OPTARG
	    ;;
	c) selectedcent=$OPTARG
	    ;;
	f) FUNCTION=$OPTARG
	    ;;
	h) show_help
	    ;;
	i) FILEIN=$OPTARG
	    ;;
	x) XMIN=$OPTARG
	    ;;
	X) XMAX=$OPTARG
	    ;;
	l) ISCHI2=0
	    ;;
	m) RUNMULTI=1
	    ;;
	p) startpt=$OPTARG
	    ;;
	P) stoppt=$OPTARG
	    ;;
	r) ISREBIN=1
	    ;;
	w) read_Wlimits
	    ;;
	*) default_fit
	    ;;
    esac
done

#define function to launch fit with different ranges
function run_fit()
{
    root -l -b -q $ASD/kstar/roofit/fitKStarPro.C"(\"${FILEIN}\",\"${DIROUT}\",${ISREBIN}, ${ISCHI2}, ${BG},\"${FUNCTION}\",${selectedcent},${startpt},${stoppt},${splitc},${XMIN},${XMAX},kTRUE,kTRUE,kTRUE,0.0,2.0, ${ISWFIXED}, ${minW}, ${maxW})"
}

#define function to launch fit with different ranges - one bin per output file
function run_multiple_ranges()
{
    for selectedcent in 0 1 2 3 4
    do
	for startpt in 3 4 5 6 7 8 9 10 11 12 13 14
	do
	    for rmin in 0.68 0.70 0.72 0.74
	    do
		for rmax in 1.08 1.10 1.12 1.15
		do
		    root -l $ASD/kstar/roofit/fitKStarPro.C"(\"${FILEIN}\",\"${DIROUT}\",${ISREBIN}, ${ISCHI2}, ${BG},\"${FUNCTION}\",${selectedcent},${startpt},${startpt}+1,1,${rmin},${rmax},kTRUE,kTRUE,kTRUE,0.0,2.0, ${ISWFIXED}, ${minW}, ${maxW})"
		    
		    echo "Save result? Press 'd' to discard, 'q' to quit, any other key to save..."
		    read -n 1 SAVER
		    if [ "$SAVER" = "d" ] 
		    then
			if [ "$BG" -eq "1" ]
			then
			    echo "removing file: ${DIROUT}/*fitEM_${FUNCTION}_${rmin}-${rmax}_*cent${selectedcent}_Pt${startpt}.root"
			    rm ${DIROUT}/*fitEM_${FUNCTION}_${rmin}-${rmax}_*cent${selectedcent}_Pt${startpt}.root 
			else
			    echo "removing file: ${DIROUT}/*fitLS_${FUNCTION}_${rmin}-${rmax}_*cent${selectedcent}_Pt${startpt}.root"
			    rm ${DIROUT}/*fitLS_${FUNCTION}_${rmin}-${rmax}_*cent${selectedcent}_Pt${startpt}.root  
			fi
		    else
			if [ "$SAVER" = "q" ]
			then 
			    exit;
			else
			    echo "YOUR RESULT IS SAVED!"			
			fi
		    fi
		done
	    done
	done
    done
}

#define function to launch fit with different ranges - all bins in one output file
function run_multiple_ranges2()
{
    for rmin in 0.68 0.7 0.72
    do
	for rmax in 1.08 1.1 1.12
	do
	    root -l -b -q $ASD/kstar/roofit/fitKStarPro.C"(\"${FILEIN}\",\"${DIROUT}\",${ISREBIN}, ${ISCHI2}, ${BG},\"${FUNCTION}\",-1 ,1, 13, 4,${rmin},${rmax},kTRUE,kTRUE,kTRUE,0.0,2.0, ${ISWFIXED}, ${minW}, ${maxW})"
	done
    done
}
###############
# MAIN EXEC
###############

#checks parameters
[ ! -e ${DIROUT} ] && { echo "Outpur dir does not exist. Creating..."; mkdir ${DIROUT}; }

[ -z $FUNCTION ] && { echo "Please specify fit function..."; exit 1; }
[ -z $FILEIN ] && { echo "Please specify input file..."; exit 1; }

[ ! -e $FILEIN ] && { echo "File does not exist."; exit 2; }

if [ $RUNMULTI -eq 1 ] 
then 
    run_multiple_ranges
else
    [ -z $XMIN ]  && { echo "Please specify fit range min..."; exit 1; }
    [ -z $XMAX ]  && { echo "Please specify fit range max..."; exit 1; }
    run_fit
fi
