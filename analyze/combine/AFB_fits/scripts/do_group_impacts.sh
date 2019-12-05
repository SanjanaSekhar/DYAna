#!/bin/bash

#stolen from https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/tree/81x-root606/data/tutorials/groups

#idx=$1
vars=(ee16_fake_shape mumu16_fake_shape  ee17_fake_shape mumu17_fake_shape ee18_fake_shape mumu18_fake_shape  pdfs)
#vars=(pdfs)

set -x 
for idx in 0 1 2 3 4 5
do
    comb_card=cards/combined_fit_mbin${idx}.txt
    workspace=workspaces/combined_fit_impacts_${idx}.root
    text2workspace.py $comb_card --keyword-value M_BIN=${idx} -P Analysis.DYAna.my_model:dy_AFB -o $workspace
    combine -M MultiDimFit $workspace --saveWorkspace  -n _nom  #>& logs/nomfitinit_${idx}
    combine -M FitDiagnostics -d higgsCombine_nom.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _nom1 #&> logs/nomfit_${idx}
    cp impacts/mbin${idx}.json impacts/group_mbin${idx}.json
    for par in ${vars[*]};  do
        combine -M FitDiagnostics --freezeNuisanceGroups $par -d higgsCombine_nom.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _${par} #&> logs/${par}fit_${idx}
        ./add_group_impact.py higgsCombine_${par}.FitDiagnostics.mH120.root higgsCombine_nom1.FitDiagnostics.mH120.root ${par} impacts/group_mbin${idx}.json
    done

    #plotImpacts.py -i impacts/group_mbin${idx}.json -o impacts/impact_plot_group_mbin${idx} --POI Afb #-t impacts/rename.json
    #rm higgsCombine*


done

