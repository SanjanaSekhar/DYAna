
import ROOT
from ROOT import *
from utils import *

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--nToys",  default=100, type='int', help="How many toys to run")
parser.add_option("--mbin",  default=0, type='int', help="Which mass bin to run on ")
parser.add_option("--teststat", default="saturated", help="Which gof test to use (saturated,KS,AD)")
parser.add_option("--mask_ee_ss", default=False, action="store_true", help="Mask ee_ss channels from gof calculation (not from fit)")
parser.add_option("--mask_ee", default=False, action="store_true", help="Mask ee_ss channels from gof calculation (not from fit)")
parser.add_option("--mask_mumu", default=False, action="store_true", help="Mask ee_ss channels from gof calculation (not from fit)")
parser.add_option("--prefit", default=False, action="store_true", help="Sample toys from prefit uncs")
parser.add_option("--notfreq", default=False, action="store_true", help="Don't use toysFrequentist option")
parser.add_option("--reuse_fit", default=False, action="store_true", help="Reuse initial fit from previous run to save time")
parser.add_option("-y", "--year", default = -1, type='int', help="Only do fits for this single year (2016,2017, or 2018), default is all years")
parser.add_option("-o", "--odir", default="", help = "output directory")
parser.add_option("--noSymMCStats", default = False, action="store_true",  help="Don't add constraints to mcStat nuisances")

(options, args) = parser.parse_args()

chan = "combined"

toys_freq = "--toysFrequentist" 
if(options.notfreq):
    toys_freq = ""

extra_params = ""
if(options.mask_ee or options.mask_mumu):
    extra_params = "--setParameters "
if(options.mask_ee):
    extra_params+="mask_Y16_ee16=1,mask_Y17_ee17=1,mask_Y18_ee18=1" 
if(options.mask_mumu):
    extra_params+="mask_Y16_mumu16=1,mask_Y17_mumu17=1,mask_Y18_mumu18=1" 

extra_params += "--X-rtd MINIMIZER_no_analytic"
afb = 0.0
a0 = 0.05


workspace = "workspaces/%s_gof_tests_%i.root" % (chan, options.mbin)
if(not options.prefit):
    if(not options.reuse_fit):
        make_workspace(workspace, options.mbin, year = options.year, symMCStats = not (options.noSymMCStats) )
        print_and_do("combine -M MultiDimFit -d %s --saveFit --saveWorkspace --robustFit 1 %s" % (workspace, extra_params))

    fitted_afb, fitted_a0 = setSnapshot(Afb_val = -1., mdf = True)
    print_and_do("combine -M GoodnessOfFit -d %s  --algo=%s %s" % (workspace,options.teststat, extra_params))
    print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root --snapshotName initialFit %s --bypassFrequentistFit --saveToys -t %i  --setParameters Afb=%.2f,A0=%.2f" 
            % (toys_freq, options.nToys, afb,a0))
    print_and_do("combine -M GoodnessOfFit -d %s --algo=%s --toysFile higgsCombineTest.GenerateOnly.mH120.123456.root -t %i %s %s" %(workspace, options.teststat, options.nToys, toys_freq, extra_params))
else:
    make_workspace(workspace, options.mbin, year = options.year, symMCStats = not (options.noSymMCStats) )
    print_and_do("combine -M GoodnessOfFit -d %s  --algo=%s %s" % (workspace,options.teststat, extra_params))
    s = 123456
    #Do in 1 step (equivalent)
    #print_and_do("combine -M GoodnessOfFit -m 121 -d %s --algo=%s -s %i %s --saveToys -t %i --setParameters Afb=%.2f,A0=%.2f %s "
    #        % (workspace, options.teststat, s, toys_freq, options.nToys,  afb, a0, extra_params))
    #Do in 2 steps
    print_and_do("combine -M GenerateOnly -d %s -s %i %s --saveToys -t %i --setParameters Afb=%.2f,A0=%.2f" 
            % (workspace, s, toys_freq, options.nToys, afb, a0))
    print_and_do("combine -M GoodnessOfFit -d %s --algo=%s --toysFile higgsCombineTest.GenerateOnly.mH120.%i.root -s %i %s -t %i %s" 
            %(workspace, options.teststat, s, s, toys_freq, options.nToys, extra_params))

gof_helper(chan, mbin = options.mbin, odir = options.odir, teststat = options.teststat)


#print_and_do("mv higgsCombineTest.GoodnessOfFit.mH120.123456.root %s%s_bin%i_toys.root" % (options.odir, chan, options.mbin))
#print_and_do("rm higgsCombineTest.GoodnessOfFit.mH120.root")



