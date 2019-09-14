#!/bin/bash


set -x

ls -la
pwd
cd Analysis/DYAna/analyze
mkdir output_files
#echo '.x FakeRate/SingleElectron_mc_contam_fake_rate.C('$2','$3', "EOS_files/2018/non_QCD_files.txt", 2018);' > cmd.txt
echo '.x EMu/EMu_reco_background.C('$2','$3', "EOS_files/2018/DY_files.txt", 2016);' > cmd.txt
#0 is pdfs, 1 is other sys
#echo ".x combine/make_sys_templates.C($2,$3, 1);" > cmd.txt
echo ".q" >> cmd.txt
root -l -b < cmd.txt
xrdcp -f output_files/*.root $1
