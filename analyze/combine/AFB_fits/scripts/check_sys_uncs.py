import operator
import ROOT
from ROOT import *
from utils import *
from add_group_impact import *

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--mbin",  default=0, type='int', help="Which mass bin to run on ")
parser.add_option("--nThreads",  default=1, type='int', help="Which mass bin to run on ")
parser.add_option("-o", "--odir", default="impacts/", help = "output directory")

(options, args) = parser.parse_args()

chan = "combined"


individual_pars = [ "dy_xsec", "db_xsec",  "top_xsec", "gam_xsec",  "elFakesYR",  "muFakesYR", "Pu",  "BTAGYR", "prefireYR", "METJECYR", "muRCYR"]
group_pars =["pdfs", "lumisYR","elScalesYR", "elHLTsYR", "elIDs", "elRECOs",   "muIDsYR", "muHLTsYR", "RFscalesYRC","emucostrwsYRC","elfakesrwsYR","mufakesrwsYR", "ptrwsYR"] 



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
    if(n>3): break

for group_par in group_pars:
    freeze_str = par_to_freezestr(group_par)
    print_and_do("""combine -M FitDiagnostics --freezeNuisanceGroups %s -d higgsCombine_nom.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _%s """ % (freeze_str, group_par))
    sys_unc = compute_sys("nom1", group_par)
    d[group_par] = sys_unc
    break

print(d)
with open("mbin%i_sys_uncs.txt" % options.mbin, 'w') as f_out:
    sorted_d = sorted(d.items(), key=operator.itemgetter(1))
    for sys_name, val in sorted_d[::-1]:
        f_out.write("%s %.3e \n" % (sys_name, val))


print_and_do("rm higgsCombine* fitDiagnostics*")

