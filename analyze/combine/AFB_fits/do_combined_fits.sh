#!/bin/bash
set -x 
year=17

for idx in 0 1 2 3 4 5
do
    workspace=workspaces/combined_fit_y_${idx}.root
    plotdir=postfit_plots/combined_mbin${idx}
    for year in 16 17 18
    do
        card=cards/combined_fit_y${year}_mbin${idx}.txt
        cp combined_fit_template.txt $card
        sed -i 's/YR/${year}' $card
        cat cards/y20${year}_mbin${idx}_bins.txt >> $card
    done
    comb_card=cards/combined_fit_mbin${idx}.txt
    combineCards.py Y16=cards/combined_fit_y16_mbins${idx}.txt Y17=cards/combined_fit_y17_mbins${idx}.txt Y18=cards/combined_fit_y18_mbins${idx}.txt > ${comb_card}
    text2workspace.py $comb_card --keyword-value M_BIN=${idx} -P Analysis.B2GTTrees.my_model:dy_AFB -o $workspace --channel-masks
    combine -M FitDiagnostics $workspace --saveShapes --saveWorkspace --skipBOnlyFit --plot --robustFit 1
    #PostFitShapesFromWorkspace -w higgsCombineTest.FitDiagnostics.mH120.root -f fitDiagnostics.root:fit_s --postfit -o shapes/combined_fit_shapes_mbin${idx}.root
    rm -r $plotdir
    mkdir $plotdir
    mv *.png $plotdir
    echo "fit_s->Print();" > cmd.txt
    echo ".q" >> cmd.txt
    root -l -b fitDiagnostics.root < cmd.txt > fit_results/combined_fit_results_mbin${idx}.txt
    rm combine_logger.out
    rm higgsCombineTest.FitDiagnostics.mH120.root 
    rm fitDiagnostics.root
done


