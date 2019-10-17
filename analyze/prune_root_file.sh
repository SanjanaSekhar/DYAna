#!/bin/bash

f_in=$1

f_temp=${f_in}temp.root

save_trees=(T_WJets T_QCD)


iter=0
#copy over trees that I want
for tree in ${save_trees[*]}
do
    if ((${iter}==0));
    then
        rootcp --recreate ${f_in}:${tree} ${f_temp}
    else
        rootcp ${f_in}:${tree} ${f_temp}
    fi
    iter=$iter+1
done

mv ${f_temp} ${f_in}

