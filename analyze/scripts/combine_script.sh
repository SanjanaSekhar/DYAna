
#!/bin/bash


set -x

cd Analysis/DYAna/analyze/combine/AFB_fits

mkdir temp
python scripts/do_bias_test.py --nToys 50 -o temp/  --mbin $3 --Afb 0.0 --A0 0.0
#python scripts/do_gof.py --nToys 100 -o temp/  --mbin $3 --mask_ee_ss
xrdcp -f temp/* $1
