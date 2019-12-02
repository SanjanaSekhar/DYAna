#!/bin/bash
set -x 
do_plot=${1:-false}
no_sys=${2:-false}

temp_card=combined_fit_template.txt
if $no_sys; then
    temp_card=combined_fit_template_nosys.txt
fi

#mask_params="mask_Y16_emu16=1,mask_Y17_emu17=1,mask_Y18_emu18=1" 

#for idx in 5 6 7
for idx in 0 1 2 3 4 5 6 7
do
#idx=1
    workspace=workspaces/combined_fit_yall_${idx}.root
    plotdir=postfit_plots/combined_yall_mbin${idx}
    for year in 16 17 18
    do
        card=cards/combined_fit_y${year}_mbin${idx}.txt
        cp $temp_card $card
        sed -i "s/YR/${year}/g" $card
        #if [ ! $no_sys ]; then 
            cat cards/y20${year}_mbin${idx}_bins.txt >> $card
        #fi
    done
    comb_card=cards/combined_fit_mbin${idx}.txt
    combineCards.py Y16=cards/combined_fit_y16_mbin${idx}.txt Y17=cards/combined_fit_y17_mbin${idx}.txt Y18=cards/combined_fit_y18_mbin${idx}.txt > ${comb_card}
    text2workspace.py $comb_card --keyword-value M_BIN=${idx} -P Analysis.DYAna.my_model:dy_AFB -o $workspace --channel-masks
    if $do_plot; then
        combine -M FitDiagnostics $workspace --plot --saveWorkspace --skipBOnlyFit  #--setParameters $mask_params #-v 10
        python scripts/my_diffNuisances.py fitDiagnostics.root -p Afb --skipFitB -g impacts/comb_pulls_mbin${idx}.root
        rm -r $plotdir
        mkdir $plotdir
        mv *.png $plotdir
    else
        combine -M FitDiagnostics $workspace --saveWorkspace --skipBOnlyFit -v 1
    fi
    #combine -M MultiDimFit $workspace --saveWorkspace --saveFitResult --robustFit 1 #-v 10
    #PostFitShapesFromWorkspace -w higgsCombineTest.FitDiagnostics.mH120.root -f fitDiagnostics.root:fit_s --postfit -o shapes/combined_fit_shapes_mbin${idx}.root
    echo "fit_s->Print();" > cmd.txt
    echo ".q" >> cmd.txt
    root -l -b fitDiagnostics.root < cmd.txt > fit_results/combined_fit_results_yall_mbin${idx}.txt
    rm combine_logger.out
    rm higgsCombineTest.FitDiagnostics.mH120.root 
    rm fitDiagnostics.root
done


