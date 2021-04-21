
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


        
all_sys = ["METJEC", "elScaleSyst", "elScaleStat","elScaleGain", "elSmear", "muRC", "Pu", "BTAGCOR","BTAGUNCOR", "BTAGLIGHT" ,
            "muHLTBAR", "muIDBAR", "muISOBAR",  "muHLTEND", "muIDEND", "muISOEND",  "muIDSYS", "muISOSYS",  
            "elHLTBARPTHIGH", "elIDBARPTHIGH", "elRECOBARPTHIGH", "elHLTENDPTHIGH", "elIDENDPTHIGH", "elRECOENDPTHIGH",
            "elHLTBARPTLOW", "elIDBARPTLOW", "elRECOBARPTLOW", "elHLTENDPTLOW", "elIDENDPTLOW", "elRECOENDPTLOW",
            "ptrw1b", "ptrw2b", "ptrw3b", "ptrw4b", "ptrw5b", "ptrw6b", "ptrw7b",
            "emucostrw1b", "emucostrw2b", "emucostrw3b", "emucostrw4b",
            "elfakesrw1b", "elfakesrw2b", "elfakesrw3b", "elfakesrw4b",
            "mufakesrw1b", "mufakesrw2b", "mufakesrw3b", "mufakesrw4b",
            "RENORM", "FAC", "REFAC", "alphaS", 
            "dy_xsec","db_xsec","top_xsec","gam_xsec" ,"elFakes", "muFakes", 
            "lumiXY" ,"lumiLS" ,"lumiDB" ,"lumiBC", "lumiGS" ,"lumi", 
            ]

correlate_all = ["elScaleSyst", "elSmear", "Pu", "muIDSYS", "muISOSYS", "elRECOBARPTHIGH", "elRECOENDPTHIGH", "elRECOBARPTLOW", "elRECOENDPTLOW",
                 "elIDBARPTHIGH", "elIDENDPTHIGH", "elIDBARPTLOW", "elIDENDPTLOW", 
                 "dy_xsec","db_xsec"  ,"top_xsec","gam_xsec",
                 "lumiXY" ,"lumiLS" ,"lumiDB" ,"lumiBC" , "lumiGS",
                 "BTAGCOR"
                 ] 

correlate_1718 = ["ptrw1b", "ptrw2b", "ptrw3b", "ptrw4b", "ptrw5b", "ptrw6b", "ptrw7b", 
                    "emucostrw1b", "emucostrw2b", "emucostrw3b", "emucostrw4b",
                    "RENORM", "FAC", "REFAC", "alphaS" ]


pars16 = []
pars17 = []
pars18 = []
pars_comb = []

for i in range(1,61):
    pars_comb.append("pdf" + str(i))

for par in all_sys:
    if(par not in correlate_all): 
        pars16.append(par + "16")
        if(par not in correlate_1718):
            pars17.append(par + "17")
            pars18.append(par + "18")
        else:
            pars_comb.append(par + "1718")

    else:
        pars_comb.append(par)



pars16.append("prefire16")
pars17.append("prefire17")
pars18.append("METHEM18")


par_str = ""

for par in (pars16 + pars17 + pars18 + pars_comb):
    par_str += par +","

#remove last comma
par_str=par_str[:-1]



workspace = "workspaces/%s_impacts_%i.root" % (chan, options.mbin)
make_workspace(workspace, options.mbin)

print("Num pars = %i " % (len(pars16) + len(pars17) + len(pars18) + len(pars_comb)))
print(par_str)

print_and_do("combineTool.py -M Impacts -m 125 -d %s --doInitialFit " % workspace)
print_and_do("combineTool.py -M Impacts -m 125 -d %s --doFits --named %s --parallel %i" % (workspace, par_str, options.nThreads))
print_and_do("combineTool.py -M Impacts -m 125 -d %s -o %s/impacts_mbin%i.json --named %s" % (workspace, options.odir, options.mbin, par_str))
print_and_do("python scripts/my_plotImpacts.py -i %s/impacts_mbin%i.json -o %s/impact_plot_afb_mbin%s --POI Afb --blind" % (options.odir, options.mbin, options.odir, options.mbin))
print_and_do("python scripts/my_plotImpacts.py -i %s/impacts_mbin%i.json -o %s/impact_plot_a0_mbin%s --POI A0 --blind" % (options.odir, options.mbin, options.odir, options.mbin))
print_and_do("rm higgsCombine_*")


