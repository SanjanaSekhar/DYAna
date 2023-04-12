#!/bin/bash
set -x

cd LQ_Analysis/DYAna/analyze/combine/AFB_fits
#mkdir limits
#mkdir temp

mkdir imps
python scripts/LQ_do_impacts.py --mLQ 2000 --chan ee --q u -o imps --ending 041123 
xrdcp -f imps/* $1 
