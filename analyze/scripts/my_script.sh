#!/bin/bash


set -x

ls -la
pwd
cd Analysis/DYAna/analyze
mkdir output_files
#echo '.x MuMu/MuMu_reco_mc.C('$2','$3', "EOS_files/2016/DY_files.txt", 2016);' > cmd.txt
#echo '.x ElEl/ElEl_reco_mc.C('$2','$3', "EOS_files/2016/DY_files.txt", 2016);' > cmd.txt
#echo '.x FakeRate/SingleMuon_data_measure_fake_rate.C('$2','$3', "EOS_files/2016/SingleMuon_files.txt", 2016);' > cmd.txt
#echo '.x EMu/EMu_reco_data.C('$2','$3', "EOS_files/2018/SingleMuon_files.txt", 2018);' > cmd.txt

#0 is pdfs, 1 is other sys
xrdcp root://cmseos.fnal.gov//store/user/oamram/Condor_inputs/output_files.tgz .
tar -xvf output_files.tgz
echo ".x combine/make_sys_templates.C($2,$3, 2018, 0);" > cmd.txt


echo ".q" >> cmd.txt
root -l -b < cmd.txt
xrdcp -f output_files/*.root $1
