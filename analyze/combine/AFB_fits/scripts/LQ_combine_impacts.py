#!/usr/bin/env python
import operator
import ROOT
import string
from ROOT import *
from utils import *
from add_group_impact import *

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--chan",  default="ee", help="Which chan to run ")
parser.add_option("--q",  default="u", help="Which q to run ")
parser.add_option("-o", "--odir", default="sys_uncs/", help = "output directory")
parser.add_option("-i", "--fin", default="test.json", help = "Input json")

(options, args) = parser.parse_args()


if options.chan == "mumu" and options.q = "u": options.fin = "impacts/mumu_u_impacts_mLQ1000.json"
elif options.chan == "mumu" and options.q = "d": options.fin = "impacts/mumu_d_impacts_mLQ1000.json"
elif options.chan == "ee" and options.q = "u": options.fin = "impacts/ee_u_impacts_mLQ1000.json"
if options.chan == "ee" and options.q = "d": options.fin = "impacts/ee_d_impacts_mLQ1000.json"

individual_pars = [ "dy_xsec", "db_xsec",  "top_xsec", "gam_xsec",  "elFakesYR",  "muFakesYR", "Pu", "prefire", "METJEC", "muRC"]
group_pars =[  "RFscales", "emucostrws", "ptrws", "pdfs", "elScales", "elHLTs", "elIDs", "elRECOs",   "muIDs", "muHLTs", 
                "elfakesrws","mufakesrws", "BTAGS"] 
extras = ["MET", "lumi"]

sys_name_conv = dict()

sys_name_conv['dy_xsec'] = "DY cross section"
sys_name_conv['db_xsec'] = "Diboson cross section"
sys_name_conv['top_xsec'] = "$\\ttbar$ cross section"
sys_name_conv['gam_xsec'] = "$\\gamma\\gamma$ cross section"
sys_name_conv['elFakes'] = "Electron Fakes Normalization"
sys_name_conv['muFakes'] = "Muon Fakes Normalization"
sys_name_conv['Pu'] = "Pileup"
sys_name_conv['prefire'] = "Prefire"
sys_name_conv['MET'] = "MET Uncertainties"
sys_name_conv['muRC'] = "Muon Scale"
sys_name_conv['pdfs'] = "PDFs"
sys_name_conv['lumi'] = "Luminosity"
sys_name_conv['elHLTs'] = "Electron Trigger"
sys_name_conv['elIDs'] = "Electron Identification/Isolation"
sys_name_conv['elRECOs'] = "Electron Reconstruction"
sys_name_conv['muIDs'] = "Muon Identification/Isolation"
sys_name_conv['muHLTs'] = "Muon Trigger"
sys_name_conv['RFscales'] = "Renormalization/Factorization Scales"
sys_name_conv['emucostrws'] = "$e\\mu$ Shape Corrections"
sys_name_conv['elfakesrws'] = "Electron Fakes Shape"
sys_name_conv['mufakesrws'] = "Muon Fakes Shape"
sys_name_conv['ptrws'] = "DY $p_{T}$ Correction"
sys_name_conv['BTAGS'] = "b-tagging Uncertainty"
sys_name_conv['elScales'] = "Electron Scale"


def add_entry(d, entry):

    sys_name = str(entry['name'])
    key_name = sys_name.translate(None, string.digits)
    groups = entry['groups']
    if(len(groups) > 0):
        group = groups[0] #extra groups are redundant 

        for g_name in group_pars:
            if(g_name in group):
                key_name = g_name
                break

    for key in extras:
        if(key in sys_name):
            key_name = key
            break


    val = entry['impact_yLQ2']**2

    if(key_name in d.keys()): d[key_name] += val
    else: d[key_name] = val

def sqrt(d):
    d_new = dict()
    for key in d.keys():
        cont = d[key]
        new_cont = cont**(0.5)
        d_new[key] = new_cont
    return d_new


json_file = open(options.fin, 'r')
data = json.load(json_file)
d = dict()

for param in data['params']:
    add_entry(d, param)

d = sqrt(d)


#print(d)
os.system("mkdir %s \n" % options.odir)
with open("%s/%s_%s_sys_uncs.txt" % (options.odir, options.chan, options.q), 'w') as f_out:
    sorted_d = sorted(d.items(), key=operator.itemgetter(1))
    f_out.write("Systematic uncertainties (values x1000) for S_%s%s \n" % (options.chan,options.q))
    for sys_name, val in sorted_d[::-1]:
        if (sys_name in sys_name_conv.keys()):
                out_name = sys_name_conv[sys_name]
        else:
            out_name = sys_name
        #print(sys_name, out_name)
        f_out.write("%s & %.3f \\\\ \n" % (out_name, val*1000.))

