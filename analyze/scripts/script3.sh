#!/bin/bash
set -x

cd LQ_Analysis/DYAna/analyze/combine/AFB_fits
#mkdir limits
#mkdir temp

mkdir sys
python scripts/LQ_check_sys_uncs.py  --chan ee --q d -o sys --ending 052524   --expected true --mLQ 3000 --nToys 1 -s 1619
xrdcp -f sys/* $1 
