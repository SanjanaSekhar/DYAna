#!/bin/bash
set -x 

#for idx in 0 1 2 3 4 5
for idx in 0 1 2 3 4 5 6 7
do
    plotImpacts.py -i impacts/mbin${idx}.json -o impacts/impact_plot_mbin${idx} -t impacts/rename.json --POI Afb

done
