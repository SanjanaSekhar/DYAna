#!/bin/bash


set -x

ls -la
pwd
cd Analysis/DYAna/analyze
mkdir output_files
echo '.x ElEl/ElEl_reco_mc.C('$2','$3', "EOS_files/2017/DY_files_mlow.txt", 2017);' > cmd.txt
#echo '.x FakeRate/SingleMuon_data_measure_fake_rate.C('$2','$3', "EOS_files/2018/SingleMuon_files.txt", 2018);' > cmd.txt
#echo '.x ElEl/ElEl_reco_data.C('$2','$3', "EOS_files/2018/SingleElectron_files.txt", 2018);' > cmd.txt
#0 is pdfs, 1 is other sys
#echo ".x combine/make_sys_templates.C($2,$3, 2016, 1);" > cmd.txt
echo ".q" >> cmd.txt
root -l -b < cmd.txt
xrdcp -f output_files/*.root $1
