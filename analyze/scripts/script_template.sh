#!/bin/bash


set -x

ls -la
pwd
cd Analysis/DYAna/analyze
mkdir output_files
echo '.x PREFIX/EXEC('$2','$3', "EOS_files/YEAR/EOSINPUT", YEAR);' > cmd.txt

echo ".q" >> cmd.txt
root -l -b < cmd.txt
xrdcp -f output_files/*.root $1
