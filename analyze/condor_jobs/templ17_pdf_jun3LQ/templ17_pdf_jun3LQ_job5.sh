#!/bin/bash
# Job Number 5, of 15 
set -x 
source /cvmfs/cms.cern.ch/cmsset_default.sh
pwd
export SCRAM_ARCH=slc6_amd64_gcc530
xrdcp root://cmseos.fnal.gov//store/user/sasekhar/Condor_inputs/DY_CMSSW.tgz CMSSW.tgz 
cat LQ_my_script.sh 
tar -xzf CMSSW.tgz 
mv LQ_my_script.sh DY_analysis/src/ 
cd DY_analysis/src 
eval `scramv1 runtime -sh`
scram b ProjectRename 
xrdcp root://cmseos.fnal.gov//store/user/sasekhar/Condor_inputs/Analysis.tgz tarDir.tgz
tar -xvzf tarDir.tgz 
scram b -j 
eval `scramv1 runtime -sh`
./LQ_my_script.sh root://cmseos.fnal.gov//store/user/sasekhar/Condor_outputs/templ17_pdf_jun3LQ/file_5.root 15 5 
cd ${_CONDOR_SCRATCH_DIR} 
rm -rf CMSSW_10_5_0
