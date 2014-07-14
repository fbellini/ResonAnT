#!/bin/bash

systworkdir="/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s/systematics"
rangelist="${systworkdir}/range_width_systematics_files.lst"
filebins="${systworkdir}/../_sum_EMnorm1.30-1.50_cut1818_tof2s_train215-216.root"
centralspectra="/Users/bellini/alice/resonances/kstar_pA5.02TeV/pAexpress215-216/TOF_ANA/ana2s/fit_syst_poly2_fixedW_l1_r0/best_fit_poly2.root"
barlow=true

cd ${systworkdir}
root -l<<EOF
.L $ASD/kstar/systematics/systBinary.C
.L $ASD/kstar/systematics/systematics_pA.C
systBinary_pA_ResBg("TOF", ${barlow})
systBinary_pA_LSnorm("TOF", ${barlow})
systBinary_pA_PID("TOF", ${barlow})
systematics_pA("${rangelist}", 5.0, "${filebins}","range_width", 0,"${centralspectra}")
EOF

# systBinary_pA_EM("TOF", ${barlow})
# systBinary_pA_freeW("TOF", ${barlow})
# systBinary_pA_BinCount("TOF", ${barlow})

ls systUncert/syst_*.root |tee output_syst_files.lst
for isystfile in `cat output_syst_files.lst`
do
    root -l -b -q "$ASD/kstar/systematics/compareSyst2Stat.C(\"${isystfile}\", 1, 1)"
done