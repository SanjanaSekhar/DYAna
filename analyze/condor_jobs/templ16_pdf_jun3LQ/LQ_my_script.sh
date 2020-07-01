#!/bin/bash


set -x

ls -la
pwd
cd LQ_Analysis/DYAna/analyze
mkdir output_files

#0 is pdfs, 1 is other sys
xrdcp root://cmseos.fnal.gov//store/user/sasekhar/Condor_inputs/output_files.tgz .
tar -xvf output_files.tgz
echo ".x combine/LQ_make_sys_templates.C($2,$3, 2016, 0);" > cmd.txt


echo ".q" >> cmd.txt
root -l -b < cmd.txt
xrdcp -f output_files/*.root $1