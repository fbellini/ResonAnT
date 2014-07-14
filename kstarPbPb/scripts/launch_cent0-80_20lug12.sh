#!bin/bash
cd /Users/bellini/alice/resonances/myKstar/AOD_PbPb2010/train_output/analysis_cent0-80
aliroot -b<<EOF
.L kstar/projectKStar.C
projectKStar("analysisAOD_0-80.root");
.q
EOF
cd /Users/bellini/alice/resonances/myKstar/AOD_PbPb2010/train_output/analysis_cent0-80/roofit
aliroot -b<<EOF
.L kstar/kstarAnalysis.C
kstarAnalysis("/Users/bellini/alice/resonances/myKstar/AOD_PbPb2010/train_output/analysis_cent0-80/proj/");
.q
EOF
for ibg in 1 2 
do
    aliroot -b<<EOF
.L kstar/roofit/fitKStar.C
fitKStar_ByPtBins("sub_proj_analysisAOD_0-80.root",ibg);
.q
EOF
done
ls
pwd 