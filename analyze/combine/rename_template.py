import ROOT
from optparse import OptionParser
import sys

parser = OptionParser()
parser.add_option("-i", "--in_name",  default = '', help="Input systematic name")
parser.add_option("-o", "--out_name", default = '', help="Output systematic name")
options, args = parser.parse_args()

fin = sys.argv[1]

f = ROOT.TFile.Open(fin, "UPDATE")
n_bins = 8


if(len(options.in_name) == 0):
    print("missing input systematic name")
    sys.exit(1)
if(len(options.out_name) == 0):
    print("missing output systematic name")
    sys.exit(1)


for i in range(n_bins):
    f.cd("w%i" % i)
    keys = ROOT.gDirectory.GetListOfKeys().Clone()
    for k in keys:
        #print(h.GetName())
        name = k.GetName()
        if(options.in_name in name):
            h = ROOT.gDirectory.Get(name)
            new_name = name.replace(options.in_name,  options.out_name)
            print("copying %s to %s" % (name, new_name))
            h_clone = h.Clone(new_name)
            h_clone.Write()
    #ROOT.gDirectory.ls()
            

