#!/bin/bash


set -x

ls -la
pwd
cd Analysis/DYAna/analyze
mkdir output_files
#echo ".x MuMu/MuMu_reco_data_batch.C($2,$3);" > cmd.txt
echo '.x ElEl/ElEl_reco_data_batch.C('$2','$3', "EOS_files/2016/SingleMuon_files.txt", 2016);' > cmd.txt
#0 is pdfs, 1 is other sys
#echo ".x combine/make_sys_templates.C($2,$3, 1);" > cmd.txt
echo ".q" >> cmd.txt
root -l -b < cmd.txt
xrdcp -f output_files/*.root $1
