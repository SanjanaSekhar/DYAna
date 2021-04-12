import operator
import ROOT
from ROOT import *
from utils import *
from add_group_impact import *

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--mbin",  default=0, type='int', help="Which mass bin to run on ")
parser.add_option("-o", "--odir", default="sys_uncs/", help = "output directory")

(options, args) = parser.parse_args()

chan = "combined"


individual_pars = [ "dy_xsec", "db_xsec",  "top_xsec", "gam_xsec",  "elFakesYR",  "muFakesYR", "Pu", "prefireYR", "METJECYR", "muRCYR"]
group_pars =["pdfs", "lumisYR","elScalesYR", "elHLTsYR", "elIDs", "elRECOs",   "muIDsYR", "muHLTsYR", "RFscalesYRC","emucostrwsYRC","elfakesrwsYR","mufakesrwsYR", "ptrwsYR",
        "BTAGSYR"] 

sys_name_conv = dict()

sys_name_conv['dy_xsec'] = "DY cross section"
sys_name_conv['db_xsec'] = "Diboson cross section"
sys_name_conv['top_xsec'] = "$\\ttbar$ cross section"
sys_name_conv['gam_xsec'] = "$\\gamma\\gamma$ cross section"
sys_name_conv['elFakesYR'] = "Electron Fakes Normalization"
sys_name_conv['muFakesYR'] = "Muon Fakes Normalization"
sys_name_conv['Pu'] = "Pilup"
sys_name_conv['prefireYR'] = "Prefire"
sys_name_conv['METJECYR'] = "MET Uncertainties"
sys_name_conv['muRCYR'] = "Muon Scale"
sys_name_conv['pdfs'] = "PDFs"
sys_name_conv['lumisYR'] = "Luminosity"
sys_name_conv['elHLTsYR'] = "Electron Trigger"
sys_name_conv['elIDsYR'] = "Electron Identification/Isolation"
sys_name_conv['elRECOs'] = "Electron Reconstruction"
sys_name_conv['muIDsYR'] = "Muon Identification/Isolation"
sys_name_conv['muHLTsYR'] = "Muon Trigger"
sys_name_conv['RFScalesYRC'] = "Renormalization/Factorization Scales"
sys_name_conv['emucostrwYRC'] = "$e\\mu$ Shape Corrections"
sys_name_conv['elfakesrwsYR'] = "Electron Fakes Shape"
sys_name_conv['mufakesrwsYR'] = "Muon Fakes Shape"
sys_name_conv['ptrwsYR'] = "DY $p_{T}$ Correction"




def par_to_freezestr(par):
    if("YRC" in par):
        par16 = par.replace("YR", "16")
        par1718 = par.replace("YRC", "1718")
        return par16 + "," + par1718
    elif("YR" in par):
        if("prefire" not in par):
            par16 = par.replace("YR", "16")
            par17 = par.replace("YR", "17")
            par18 = par.replace("YR", "18")
            return par16 + "," + par17 + "," + par18
        else: #prefire only sys with no 18
            return "prefire16,prefire17"
    else:
        return par

        





workspace = "workspaces/%s_impacts_%i.root" % (chan, options.mbin)
make_workspace(workspace, options.mbin)

print_and_do("combine -M MultiDimFit %s --saveWorkspace  -n _nom" % workspace) 
print_and_do("""combine -M FitDiagnostics -d higgsCombine_nom.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _nom1""") #&> logs/nomfit_${idx}
d = dict()
n = 0
for indi_par in individual_pars:
    n+=1
    freeze_str = par_to_freezestr(indi_par)
    print_and_do("""combine -M FitDiagnostics --freezeParameters %s -d higgsCombine_nom.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _%s """ % (freeze_str, indi_par))
    sys_unc = compute_sys("nom1", indi_par)
    d[indi_par] = sys_unc

for group_par in group_pars:
    freeze_str = par_to_freezestr(group_par)
    print_and_do("""combine -M FitDiagnostics --freezeNuisanceGroups %s -d higgsCombine_nom.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _%s """ % (freeze_str, group_par))
    sys_unc = compute_sys("nom1", group_par)
    d[group_par] = sys_unc

print(d)
os.system("mkdir %s \n" % options.odir)
with open("%s/mbin%i_sys_uncs.txt" % (options.odir, options.mbin), 'w') as f_out:
    sorted_d = sorted(d.items(), key=operator.itemgetter(1))
    f_out.write("Systematic uncertainties (values x1000) for bin %i \n")
    for sys_name, val in sorted_d[::-1]:
        if (sys_name in sys_name_conv.keys()):
                out_name = sys_name_conv[sys_name]
        else:
            out_name = sys_name
        f_out.write("%s & %.3f \n" % (out_name, val*1000.))


print_and_do("rm higgsCombine* fitDiagnostics*")

