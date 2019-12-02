#!/bin/bash
set -x 

cd impacts/
if [ -f file_0.root ]; then

    rm mbin?.json
    rename file mbin *.root
    rename .root .json *.root 
fi

for idx in 0 1 2 3 4 5 6 7
do
    plotImpacts.py -i mbin_${idx}.json -o impact_plot_mbin${idx} -t rename.json --POI Afb

done

cd ..
