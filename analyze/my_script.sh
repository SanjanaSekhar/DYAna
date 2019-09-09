#!/bin/bash


set -x

ls -la
pwd
cd Analysis/DYAna/analyze
mkdir output_files
echo ".x FakeRate/SingleElectron_data_measure_fake_rate.C($2,$3);" > cmd.txt
#echo '.x ElEl/ElEl_reco_background_batch.C('$2','$3', "EOS_files/2017/WT_files.txt");' > cmd.txt
#0 is pdfs, 1 is other sys
#echo ".x combine/make_sys_templates.C($2,$3, 1);" > cmd.txt
echo ".q" >> cmd.txt
root -l -b < cmd.txt
xrdcp -f output_files/*.root $1
