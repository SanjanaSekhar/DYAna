#!/bin/bash


set -x

ls -la
pwd
cd Analysis/DYAna/analyze
mkdir output_files

#0 is pdfs, 1 is other sys
xrdcp root://cmseos.fnal.gov//store/user/oamram/Condor_inputs/output_files.tgz .
tar -xvf output_files.tgz
echo ".x combine/make_sys_templates.C($2,$3, YEAR, TYPE);" > cmd.txt


echo ".q" >> cmd.txt
root -l -b < cmd.txt
xrdcp -f output_files/*.root $1
