
import ROOT
from ROOT import *
from LQ_utils import *

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--nToys",  default=1, type='int', help="How many toys to run")
parser.add_option("--mbin",  default=0, type='int', help="Which mass bin to run on ")
parser.add_option("--teststat", default="saturated", help="Which gof test to use (saturated,KS,AD)")
parser.add_option("--mask_ee_ss", default=False, action="store_true", help="Mask ee_ss channels from gof calculation (not from fit)")
parser.add_option("--mask_ee", default=False, action="store_true", help="Mask ee_ss channels from gof calculation (not from fit)")
parser.add_option("--mask_mumu", default=False, action="store_true", help="Mask ee_ss channels from gof calculation (not from fit)")
parser.add_option("--prefit", default=False, action="store_true", help="Sample toys from prefit uncs")
parser.add_option("--notfreq", default=False, action="store_true", help="Don't use toysFrequentist option")
parser.add_option("--reuse_fit", default=False, action="store_true", help="Reuse initial fit from previous run to save time")
parser.add_option("-y", "--year", default = -1, type='int', help="Only do fits for this single year (2016,2017, or 2018), default is all years")
parser.add_option("-o", "--odir", default="gof", help = "output directory")
parser.add_option("--noSymMCStats", default = True, action="store_true",  help="Don't add constraints to mcStat nuisances")
parser.add_option("--chan",  default="ee", type="string", help="What channels to run the fit over (combined, ee, or mumu)")
parser.add_option("--q",  default="u", type="string", help="What channels to run the fit over (combined, u, or d)")
parser.add_option("--mLQ",  default=1000, type='int', help="mLQ")
parser.add_option("--vec",  default=False, help="is vec?")
parser.add_option("--plot",  default=False, help="copy plots from eos to local")
#parser.add_option("-o", "--odir", default="LQ_cards/condor/", help = "output directory")

(options, args) = parser.parse_args()

chan = options.chan
#options.year = 2016
options.gen_level = False
is_vec = options.vec 
mLQ = options.mLQ 
no_sys = False
fake_data = True 
no_LQ = False

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

#extra_params += " --X-rtd MINIMIZER_no_analytic"
A4 = 1.61
yLQ2 = 0.0
seed = 780865

if options.plot:

    for chan in ['mumu']:
        for q in ['u']:
            for mLQ in [2500]:
		#Condor_outputs/gof_ee_u_m4000_y-2001//gof_eu_mLQ4000_-1.png
		# wrong -> Condor_outputs/gof_mumu_d_m5000_y-2001/gof_md_mLQ5000_-2001.png gof/
                print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/sasekhar/Condor_outputs/gof_%s_%s_m%i_y%i/gof_%s_mLQ%i_%i.png gof/gof_%s_mLQ%i_%i_splitrap.png"%(chan, q+('_vec' if is_vec else ''), mLQ, options.year-2000,chan[0]+q+('_vec' if is_vec else ''), mLQ, options.year,chan[0]+q+('_vec' if is_vec else ''), mLQ, options.year))

else:

    workspace = "workspaces/LQ_%s_gof_tests_%i.root" % (chan,options.year)
    if(not options.prefit):
        if(not options.reuse_fit):
            make_workspace(workspace, options.gen_level, options.chan, options.q, is_vec, no_LQ, no_sys, fake_data, mLQ, year = options.year,noSymMCStats = options.noSymMCStats)
            print_and_do("combine %s -M MultiDimFit   --saveWorkspace --saveFitResult --robustFit 1 --trackErrors yLQ2 %s   -s %i   -n _%i " %(workspace, extra_params, seed, seed))

        fitted_yLQ2 = setSnapshot(yLQ2_val = -1., mdf = True, s = seed, freeze = False)
        print_and_do("combine -M GoodnessOfFit -d %s  --algo=%s %s -n _%s" % (workspace,options.teststat, extra_params,options.chan[0]+options.q+('_vec' if is_vec else '')))
        print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root --snapshotName initialFit %s --bypassFrequentistFit --saveToys -t %i  --setParameters yLQ2=%.2f,A4=%.2f" 
            % (toys_freq, options.nToys, yLQ2, A4))
        print_and_do("combine -M GoodnessOfFit -d %s --algo=%s --toysFile higgsCombineTest.GenerateOnly.mH120.123456.root -t %i %s %s -n _%s" %(workspace, options.teststat, options.nToys, toys_freq, extra_params, options.chan[0]+options.q+('_vec' if is_vec else '')))
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

    gof_helper(chan, options.q, mLQ, is_vec, options.year, odir = options.odir, teststat = options.teststat)


#print_and_do("mv higgsCombineTest.GoodnessOfFit.mH120.123456.root %s%s_bin%i_toys.root" % (options.odir, chan, options.mbin))
#print_and_do("rm higgsCombineTest.GoodnessOfFit.mH120.root")



