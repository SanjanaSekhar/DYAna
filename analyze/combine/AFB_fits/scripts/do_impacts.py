
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
parser.add_option("-o", "--odir", default="", help = "output directory")

(options, args) = parser.parse_args()

chan = "combined"

extra_params = ""


pars_corr="alphaS,alphaDen,RENORM,FAC,dy_xsec,bk_xsec,gam_xsec,"
pars16="Pu16,BTAG16,elScaleStat16,elScaleSyst16,elScaleGain16,elSmear16,elHLT16,elID16,elRECO16,muRC16,muID16,muHLT16,lumi16,ee16_qcd,mu16_qcd,R_ee16_os_fakes,"
pars17="Pu17,BTAG17,elScaleStat17,elScaleSyst17,elScaleGain17,elSmear17,elHLT17,elID17,elRECO17,muRC17,muID17,muHLT17,lumi17,ee17_qcd,mu17_qcd,R_ee17_os_fakes,"
pars18="Pu18,BTAG18,elScaleStat18,elScaleSyst18,elScaleGain18,elSmear18,elHLT18,elID18,elRECO18,muRC18,muID18,muHLT18,lumi18,ee18_qcd,mu18_qcd,R_ee18_os_fakes"

pars= pars_corr + pars16+pars17+pars18

workspace = "workspaces/%s_impacts_%i.root" % (chan, options.mbin)
make_workspace(workspace, options.mbin)

print_and_do("combineTool.py -M Impacts -d %s --doInitialFit " % workspace)
print_and_do("combineTool.py -M Impacts -d %s --doFits --named %s --parallel %i" % (workspace, pars, options.nThreads))
print_and_do("combineTool.py -M Impacts -d %s -o %s/impacts_mbin%i.json --named %s" % (workspace, options.odir, options.mbin, pars))
print_and_do("plotImpacts.py -i %i/impacts_mbin%i.json -o %s/impact_plot_mbin%s -t impacts/rename.json --POI Afb" % (options.odir, options.mbin, options.odir, options.mbin))


print_and_do("combine -M MultiDimFit -d %s --saveFit --saveWorkspace --robustFit 1" % (workspace))
fitted_afb, fitted_a0 = setSnapshot(Afb_val = -1., mdf = True)
print_and_do("combine -M GoodnessOfFit -d %s  --algo=%s %s" % (workspace,options.teststat, extra_params))
print("Based on initial fit, injecting Afb = %.3f A0 = %.3f"  %(fitted_afb, fitted_a0))
print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root --snapshotName initialFit --toysFrequentist --bypassFrequentistFit --saveToys -t %i  --setParameters Afb=%.2f,A0=%.2f" 
        % (options.nToys, fitted_afb, fitted_a0))
print_and_do("combine -M GoodnessOfFit -d %s --algo=%s --toysFile higgsCombineTest.GenerateOnly.mH120.123456.root -t %i %s" %(workspace, options.teststat, options.nToys, extra_params))

gof_helper(chan, mbin = options.mbin, odir = options.odir)


print_and_do("mv higgsCombineTest.GoodnessOfFit.mH120.123456.root %s%s_bin%i_toys.root" % (options.odir, chan, options.mbin))
print_and_do("rm higgsCombineTest.GoodnessOfFit.mH120.root")



