#!/bin/bash
set -x 


#for idx in 0 1 2 3 4 5 6 7
#do
idx=$1

    comb_card=cards/combined_fit_mbin${idx}.txt
    workspace=workspaces/combined_fit_toy_impacts_${idx}.root

    text2workspace.py $comb_card --keyword-value M_BIN=${idx} -P Analysis.DYAna.my_model:dy_AFB -o $workspace --channel-masks
    combine -M FitDiagnostics $workspace --saveWorkspace --skipBOnlyFit --robustFit 1 -t 1 --toysNoSystematics --setParameters Afb=0.6,A0=0.05 #-v 10
    python scripts/my_diffNuisances.py fitDiagnostics.root -p Afb --skipFitB -g impacts/toy_pulls_mbin${idx}.root
#done



