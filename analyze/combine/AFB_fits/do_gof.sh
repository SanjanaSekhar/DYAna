#!/bin/bash

set -x 
label=combined
label_num=0
gof_file=GoodnessOfFit/${label}_fit.txt
teststat=saturated
mask_params="mask_Y16_ee16_ss=1,mask_Y17_ee17_ss=1,mask_Y18_ee18_ss=1" 
rm $gof_file
idx=2
#for idx in 0 1 2 3 4 5
#do

    comb_card=cards/combined_fit_mbin${idx}.txt
    workspace=workspaces/${label}_fit_gof_${idx}.root
    text2workspace.py $comb_card --keyword-value M_BIN=${idx} -P Analysis.DYAna.my_model:dy_AFB -o $workspace --channel-masks
    combine -M GoodnessOfFit -d $workspace  --algo=${teststat} --setParametersForEval $mask_params
    combine -M GoodnessOfFit -d $workspace  --algo=${teststat} -t 10 -s 123 --setParameters Afb=0.6,A0=0.05 --setParametersForEval ${mask_params}

    echo "mbin ${idx}" >> $gof_file
    echo ".x GoodnessOfFit/gof_helper.C($label_num, $idx)"  > cmds.txt
    echo ".q"  >> cmds.txt
    root -l -b < cmds.txt  >> $gof_file
    mv higgsCombineTest.GoodnessOfFit.mH120.123.root GoodnessOfFit/${label}_bin${idx}_toys.root
    rm higgsCombineTest.GoodnessOfFit.mH120.root

#done
cat $gof_file



