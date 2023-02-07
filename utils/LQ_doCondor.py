#!/usr/bin/env python

# doCondor.py #############################################################################
# Python driver for submitting condor jobs 
# Oz Amram


# ------------------------------------------------------------------------------------

import ROOT as r

import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from numpy import arange
from itertools import product
#from BaconAna.Utils.makeFilelist import *

default_args = []

# Options
parser = OptionParser()
parser = OptionParser(usage="usage: %prog analyzer outputfile [options] \nrun with --help to get list of options")
parser.add_option("-o", "--outdir", default='condor_jobs/',
        help="output for analyzer. This will always be the output for job scripts.")
parser.add_option("-n", "--name", default='', 
        help="Name of job. Will be used for eos output and local directory")
parser.add_option("-v", "--verbose", dest="verbose", default=False, action="store_true", 
        help="Spit out more info")

# Make condor submission scripts options
parser.add_option("--njobs", dest="nJobs", type='int', default=-1,
        help="Split into n jobs, will automatically produce submission scripts")
parser.add_option("-s", "--script", dest="script", default="scripts/LQ_my_script.sh",
        help="sh script to be run by jobs (if splitting, should take eosoutput, nJobs and iJob as args)")
parser.add_option("--dry-run", dest="dryRun", default=False, action="store_true", 
        help="Do nothing, just create jobs if requested")
parser.add_option("--combine", dest='with_combine', default=False, action="store_true", help="include combine in tarball")

# Monitor options (submit,check,resubmit failed)  -- just pass outodir as usual but this time pass --monitor sub --monitor check or --monitor resub
parser.add_option("--sub", default=False, action="store_true", help="Submit jobs")
parser.add_option("--status", default=False, action="store_true", help="Check on submitted jobs")

parser.add_option("-e", "--haddEOS", dest='haddEOS', default = False, action='store_true',  help="Hadd EOS files together and save in output_files/YEAR directory")
parser.add_option("-g", "--getEOS", default = False, action='store_true',  help="Get EOS files and save to out directory")
parser.add_option("-y", "--year", dest='year', type='string', default = 2016,  help="Year for output file location")

parser.add_option("--tar", dest='tar', default = False, action='store_true',  help="Create tarball of current directory")
parser.add_option("--tarname", dest='tarname', default = "LQ_Analysis", help="Name of directory to tar (relative to cmssw_base)")
parser.add_option("--tarexclude", dest='tarexclude', default = '', 
        help="Name of directories to exclude from the tar (relative to cmssw_base), format as comma separated string (eg 'dir1, dir2') ")
parser.add_option("--dy", dest='DY', default = False, action="store_true",  help="Shortcut to create tarball for DY AFB analysis (DYAna directory)")
parser.add_option("--cmssw", default = False, action="store_true",  help="Shortcut to create tarball for CMSSW environment for DY AFB analysis (with combine)")
parser.add_option("--root_files", dest='root_files', default = False, action="store_true",  help="Shortcut to create tarball for root files of AFB analysis")
parser.add_option("--no_rename", default = False, action="store_true",  help="Don't rename files for storing on EOS")



cwd = os.getcwd()
(options, args) = parser.parse_args()
#if len(args) < 1 and (not options.monitor or not options.tar) : sys.exit('Error -- must specify ANALYZER')
cmssw_ver = os.getenv('CMSSW_VERSION', 'CMSSW_10_5_0')
cmssw_base = os.getenv('CMSSW_BASE')
xrd_base = 'root://cmseos.fnal.gov/'
EOS_home = '/store/user/sasekhar/'
EOS_base = xrd_base + EOS_home

# write job
def write_job(out, name, nJobs, iJob, eosout=''):
    #print 'job_i %i nfiles %i subjobi %i'%(i,n,j)
    cwd = os.getcwd()
    eos_an_file = EOS_base + 'Condor_inputs/' + options.tarname + '.tgz'
    eos_cmssw_file = EOS_base + 'Condor_inputs/' + 'LQ_CMSSW' + '.tgz'

    sub_file = open('%s/%s_job%d.sh' % (out, name, iJob), 'w')
    sub_file.write('#!/bin/bash\n')
    sub_file.write('# Job Number %d, of %d \n' % (iJob, nJobs))
    sub_file.write('set -x \n')
    sub_file.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
    sub_file.write('pwd\n')
    sub_file.write('export SCRAM_ARCH=slc6_amd64_gcc530\n')
    if(not options.with_combine):
        sub_file.write('eval `scramv1 project CMSSW %s`\n'% (cmssw_ver))
        sub_file.write('cat LQ_my_script.sh \n')
        sub_file.write('mv LQ_my_script.sh %s/src/ \n'% (cmssw_ver))
        sub_file.write('cd %s/src\n'%(cmssw_ver))
        sub_file.write('eval `scramv1 runtime -sh`\n')
    else:
        sub_file.write('xrdcp %s LQ_CMSSW.tgz \n' % eos_cmssw_file) 
        sub_file.write('cat LQ_my_script.sh \n')
        sub_file.write('tar -xzf LQ_CMSSW.tgz \n')
        sub_file.write('mv LQ_my_script.sh CMSSW_10_5_0/src/ \n')
        sub_file.write('cd CMSSW_10_5_0/src \n')
        sub_file.write('eval `scramv1 runtime -sh`\n')
        sub_file.write('scram b ProjectRename \n')

    sub_file.write('xrdcp %s tarDir.tgz\n' %eos_an_file)
    sub_file.write('tar -xvzf tarDir.tgz \n')
    #sub_file.write('rm -r tarDir.tgz \n')
    sub_file.write('scram b -j \n')
    sub_file.write('eval `scramv1 runtime -sh`\n')
    sub_file.write('./LQ_my_script.sh %s %i %i \n' % (eosout, nJobs,iJob))
    sub_file.write('cd ${_CONDOR_SCRATCH_DIR} \n')
    sub_file.write('rm -rf %s\n' % cmssw_ver)
    sub_file.close()
    os.system('chmod +x %s' % os.path.abspath(sub_file.name))

# write condor submission script
def submit_jobs(lofjobs):
    script_location = os.path.abspath(options.outdir + options.name + "/LQ_my_script.sh")
    for sub_file in lofjobs:
        #os.system('rm -f %s.stdout' % sub_file)
        #os.system('rm -f %s.stderr' % sub_file)
        #os.system('rm -f %s.log' % sub_file)
        #os.system('rm -f %s.jdl'% sub_file)
        condor_file = open('%s.jdl' % sub_file, 'w')
        condor_file.write('universe = vanilla\n')
        condor_file.write('Executable = %s\n'% sub_file)
        condor_file.write('Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )\n')
        #condor_file.write('request_disk = 500000\n') # modify these requirements depending on job
        #if(options.with_combine): 
        condor_file.write('request_memory = 24028\n')
        condor_file.write('request_cpus = 8\n')
        condor_file.write('Should_Transfer_Files = YES\n')
        condor_file.write("Transfer_Input_Files = %s, %s \n" %(script_location, sub_file))
        condor_file.write('WhenToTransferOutput = ON_EXIT \n')
        #condor_file.write('use_x509userproxy = true\n')
        #condor_file.write('x509userproxy = $ENV(X509_USER_PROXY)\n')
        condor_file.write('Output = %s.stdout\n' % os.path.abspath(condor_file.name))
        condor_file.write('Error = %s.stdout\n' % os.path.abspath(condor_file.name))
        condor_file.write('Log = %s.log\n' % os.path.abspath(condor_file.name))
        condor_file.write('Queue 1\n')
        condor_file.close()
        os.system('chmod +x %s'% os.path.abspath(condor_file.name))
        os.system('condor_submit %s'%(os.path.abspath(condor_file.name)))



#When you use -h with -c, it causes tar to archive the files symbolic links point to, instead of the linking themselves.
if options.tar:
    tar_cmd = "tar" 
    excludeList = options.tarexclude.split(',')
    if options.DY:
        print "Using DY tarball options"
        excludeList = ['LQ_Analysis/DYAna/analyze/input_files', 'LQ_Analysis/DYAna/analyze/condor_jobs', 'LQ_Analysis/DYAna/analyze/output_files', 'LQ_Analysis/DYAna/generator_stuff', 
                       'LQ_Analysis/DYAna/test','LQ_Analysis/DYAna/plots/','LQ_Analysis/DYAna/analyze/combine/AFB_fits/fit_results', #'LQ_Analysis/DYAna/analyze/combine/templates',
                       'LQ_Analysis/DYAna/analyze/combine/AFB_fits/postfit_plots' ]
        options.tarname = "LQ_Analysis"
        for item in excludeList:
            #tar_cmd += " --exclude='`%s`' " % ("echo $CMSSW_BASE/src/" + item)
            tar_cmd += " --exclude='%s' " % (item)
        tar_cmd += " --exclude='%s' " %'.git' 
        tar_cmd += " --exclude='%s' " %'*.tgz' 
        tar_cmd += " --exclude='%s' " %'.nfs*' 
        tar_cmd += " --exclude='%s' " %'LQ_Analysis/DYAna/analyze/LQ_my_script.sh'
        tar_cmd += " -zchf %s -C %s %s" % (options.tarname + ".tgz", "$CMSSW_BASE/src/", options.tarname)
        print(tar_cmd)
    if options.root_files:
        print "Tarring root files"
        options.tarname = "output_files"
        include_files = ['2016/MuMu16*dy*','2016/MuMu16*ttbar*', '2016/MuMu16*wt*','2016/MuMu16*diboson*', '2016/MuMu16*phot*', '2016/MuMu16*data*', '2016/MuMu16*fakes*',
                         '2016/ElEl16*dy*','2016/ElEl16*ttbar*', '2016/ElEl16*wt*','2016/ElEl16*diboson*', '2016/ElEl16*phot*', '2016/ElEl16*data*', '2016/ElEl16*fakes*',
                         '2017/MuMu17*dy*','2017/MuMu17*ttbar*', '2017/MuMu17*wt*','2017/MuMu17*diboson*', '2017/MuMu17*phot*', '2017/MuMu17*data*', '2017/MuMu17*fakes*',
                         '2017/ElEl17*dy*','2017/ElEl17*ttbar*', '2017/ElEl17*wt*','2017/ElEl17*diboson*', '2017/ElEl17*phot*', '2017/ElEl17*data*', '2017/ElEl17*fakes*',
                         '2018/MuMu18*dy*','2018/MuMu18*ttbar*', '2018/MuMu18*wt*','2018/MuMu18*diboson*', '2018/MuMu18*phot*', '2018/MuMu18*data*', '2018/MuMu18*fakes*',
                         '2018/ElEl18*dy*','2018/ElEl18*ttbar*', '2018/ElEl18*wt*','2018/ElEl18*diboson*', '2018/ElEl18*phot*', '2018/ElEl18*data*', '2018/ElEl18*fakes*'
                         ]
        tar_cmd = "tar -chf %s" % (options.tarname + ".tgz")
        for item in include_files:
            tar_cmd += " " + "output_files/"+item
    if options.cmssw:
        print("tarring CMSSW")
        tar_cmd += " --exclude='CMSSW_10_5_0/src/LQ_Analysis/*' " 
        tar_cmd += " --exclude='CMSSW_10_5_0/src/Analysis/*' " 
        tar_cmd += " --exclude='%s' " %'*.tgz' 
        tar_cmd += " --exclude='%s' " %'*.git*' 
        tar_cmd += " --exclude='%s' " %'.nfs*' 
        tar_cmd += " -zchf %s -C %s %s" % ("LQ_CMSSW" + ".tgz", "$CMSSW_BASE/../", 'CMSSW_10_5_0')#where to run this from 
        options.tarname = "LQ_CMSSW"


    print "Executing tar command %s \n" % tar_cmd
    os.system(tar_cmd)
    cp_cmd = "xrdcp -f %s %s" %(options.tarname + ".tgz", EOS_base + "Condor_inputs/")
    print cp_cmd
    os.system(cp_cmd)
    rm_cmd = "rm %s" %(options.tarname + ".tgz")
   # os.system(rm_cmd)
    sys.exit("Finished tarring")

elif (options.haddEOS):
    if(options.outdir != "condor_jobs/"): o_dir = options.outdir
    else: o_dir = "output_files/" + str(options.year) + "/" 
    hadd_cmd = "hadd -f " + o_dir + options.name + ".root"
    xrdfsls = "xrdfs root://cmseos.fnal.gov ls"
    hadd_cmd += " `%s -u %s | grep '.root' `" %(xrdfsls, EOS_home + 'Condor_outputs/' + options.name)
    print "Going to execute cmd %s: " % hadd_cmd
    os.system(hadd_cmd)

elif (options.getEOS):
    print("Getting files and outputting to %s" % options.outdir)
    result = subprocess.check_output(["./scripts/get_crab_file_list.sh", EOS_home + 'Condor_outputs/' + options.name])
    print(result)
    for f in result.splitlines():
        cmd = "xrdcp  -f %s %s" % (f, options.outdir)
        #print "Going to execute cmd %s: " % cmd
        os.system(cmd)

    




elif options.nJobs > 0:
# -- MAIN
    if(options.name == ""): sys.exit("ERROR: MUST PROVIDE JOB NAME \n")
    os.system('rm -r %s' % (options.outdir + options.name))
    eos_dir_name = EOS_base + 'Condor_outputs/' + options.name
    #os.system("eosrm -r %s" % eos_dir_name)
    os.system('mkdir -p %s' % (options.outdir + options.name))
    os.system('cp %s %s/LQ_my_script.sh' %(options.script, options.outdir + options.name))

    for iJob in xrange(options.nJobs):
        eos_file_name = EOS_base + 'Condor_outputs/' + options.name + ('/file_%i.root' % iJob)
        if(options.no_rename): eos_file_name = EOS_base + 'Condor_outputs/' + options.name + '/'
        write_job(options.outdir + options.name, options.name, options.nJobs, iJob, eos_file_name)

# submit jobs by looping over job scripts in output dir
if options.sub:
    odir = options.outdir + options.name

    # pick up job scripts in output directory (ends in .sh)
    os.system('xrdfs %s mkdir %s' % (xrd_base, EOS_home + 'Condor_outputs/' + options.name))
    lofjobs = []
    for root, dirs, files in os.walk(odir):
        for f in fnmatch.filter(files, '%s_*.sh' %options.name):
            lofjobs.append('%s/%s' % (os.path.abspath(root), f))
    print 'Submitting %d jobs from directory %s' % (len(lofjobs), odir)
    submit_jobs(lofjobs)
    sys.exit("Finished submitting")
