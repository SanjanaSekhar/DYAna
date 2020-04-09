
import ROOT
from ROOT import *
from utils import *

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--mbin",  default=0, type='int', help="Which mass bin to run on ")
parser.add_option("--nThreads",  default=1, type='int', help="Which mass bin to run on ")
parser.add_option("-o", "--odir", default="impacts/", help = "output directory")

(options, args) = parser.parse_args()

chan = "combined"

extra_params = ""


pars_corr="alphaS,alphaDen,RENORM,FAC,REFAC,dy_xsec,bk_xsec,gam_xsec,"


pars_year = ["lumi", "Pu", "BTAG", "METJER", "METJEC", "elScaleStat", "elScaleSyst", "elScaleGain", "elSmear", "elHLTBAR", 
"elIDBAR", "elRECOBAR", "elHLTEND", "elIDEND", "elRECOEND", "muRC", "muIDBAR", "muISOBAR", "muHLTBAR", "muIDEND", "muISOEND", "muHLTEND"]

pars16 = [par+"16" for par in pars_year]
pars16.append("prefire16")
pars16.append("mu16_fakes_norm")
pars16.append("el16_fakes_norm")
pars17 = [par+"17" for par in pars_year]
pars17.append("prefire17")
pars17.append("mu17_fakes_norm")
pars17.append("el17_fakes_norm")
pars18 = [par+"18" for par in pars_year]
pars18.append("METHEM18")
pars18.append("mu18_fakes_norm")
pars18.append("el18_fakes_norm")

pars= pars_corr
"""
for par in pars16:
    pars += par +","
for par in pars17:
    pars += par +","
for idx,par in enumerate(pars18):
        pars += par +","
        """

#remove last comma
pars=pars[:-1]


workspace = "workspaces/%s_impacts_%i.root" % (chan, options.mbin)
make_workspace(workspace, options.mbin)

print_and_do("combineTool.py -M Impacts -m 125 -d %s --doInitialFit " % workspace)
print_and_do("combineTool.py -M Impacts -m 125 -d %s --doFits --named %s --parallel %i" % (workspace, pars, options.nThreads))
print_and_do("combineTool.py -M Impacts -m 125 -d %s -o %s/impacts_mbin%i.json --named %s" % (workspace, options.odir, options.mbin, pars))
print_and_do("plotImpacts.py -i %s/impacts_mbin%i.json -o %s/impact_plot_mbin%s --POI Afb --blind --height 800" % (options.odir, options.mbin, options.odir, options.mbin))
print_and_do("rm higgsCombine_*")


