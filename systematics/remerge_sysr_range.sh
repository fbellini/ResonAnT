#workdir="/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TPC_ANA/ana2s"
workdir="/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s"
bg=1 #1 = kLike, 0=kMixing 3=bin counting
normLike=1
resBg="poly2"
pid="2s"

wmin=0.5
wmax=1.5
deltax=0.02
fitoutputfolder=""

for fixgamma in 1 0 
do 
    for extrange in -1 0 1
    do
	for rigthrange in -1 0 1
	do
	    if [ $fixgamma -eq 1 ] 
	    then
		fitoutputfolder="fit_syst_${resBg}_fixedW_l${extrange}_r${rigthrange}"
	    else
		fitoutputfolder="fit_syst_${resBg}_width${wmin}-${wmax}_l${extrange}_r${rigthrange}"
	    fi
	    cd ${workdir}/${fitoutputfolder}
	    echo "*************** SONO IN "`pwd`
	    root -l -b -q "/Users/bellini/alice/macro/localMergeFiles.C(\"best_fit_${resBg}.root\",\"goodfitResults.txt\")"  
	done 
    done 
done
 ls ${workdir}/fit_syst_*/best_fit_*.root > ${workdir}/range_width_systematics_files.lst