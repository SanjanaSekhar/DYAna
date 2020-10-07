
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


pars_corr="dy_xsec,top_xsec,db_xsec,gam_xsec,"


pars_year = ["alphaS", "A0Den", "lumi", "elFakes", "muFakes", "Pu", "BTAG", "METJER", "METJEC", "elScaleStat", "elScaleSyst", "elScaleGain", "elSmear", 
"elHLTBARPTHIGH", "elIDBARPTHIGH" , "elRECOBARPTHIGH", "elHLTENDPTHIGH", "elIDENDPTHIGH", "elRECOENDPTHIGH", "elHLTBARPTLOW" , "elIDBARPTLOW" , 
"elRECOBARPTLOW" , "elHLTENDPTLOW"  , "elIDENDPTLOW"   , "elRECOENDPTLOW",
"muRC", "muIDBAR", "muISOBAR", "muHLTBAR", "muIDEND", "muISOEND", "muHLTEND",
"emucostrw1b", "emucostrw2b", "emucostrw3b", "emucostrw4b",
"elfakesrw1b", "elfakesrw2b", "elfakesrw3b", "elfakesrw4b",
]

pars_yearc = ["RENORM", "REFAC", "FAC"]
#pars_yearc = []









pars16 = [par+"16" for par in pars_year]
pars16.append("prefire16")
pars17 = [par+"17" for par in pars_year]
pars17.append("prefire17")
pars18 = [par+"18" for par in pars_year]
pars18.append("METHEM18")

for par in pars_yearc:
    pars16.append(par + "16")
    pars17.append(par + "1718")

pars= pars_corr
for par in pars16:
    pars += par +","
for par in pars17:
    pars += par +","
for par in pars18:
    pars += par +","

#remove last comma
pars=pars[:-1]



workspace = "workspaces/%s_impacts_%i.root" % (chan, options.mbin)
make_workspace(workspace, options.mbin)

print("Num pars = %i " % (len(pars16) + len(pars17) + len(pars18)))
print(pars)

print_and_do("combineTool.py -M Impacts -m 125 -d %s --doInitialFit " % workspace)
print_and_do("combineTool.py -M Impacts -m 125 -d %s --doFits --named %s --parallel %i" % (workspace, pars, options.nThreads))
print_and_do("combineTool.py -M Impacts -m 125 -d %s -o %s/impacts_mbin%i.json --named %s" % (workspace, options.odir, options.mbin, pars))
print_and_do("pyton scripts/my_plotImpacts.py -i %s/impacts_mbin%i.json -o %s/impact_plot_mbin%s --POI Afb" % (options.odir, options.mbin, options.odir, options.mbin))
print_and_do("rm higgsCombine_*")


