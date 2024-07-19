import ROOT
from optparse import OptionParser
import sys

mLQ = 2500
year = 2016
ext = "071924_splitrap"
fin = "combine/templates/LQm%i_merge_templates%i_%s.root"%(mLQ,year-2000,ext)
f_in = ROOT.TFile.Open(fin, "READ")
fout = "combine/templates/LQm%i_48bins_templates%i_%s.root"%(mLQ,year-2000,ext)
f_out = ROOT.TFile.Open(fout, "RECREATE")
f_in.cd()
f_in.cd("LQ")
keys = ROOT.gDirectory.GetListOfKeys().Clone()

for k in keys:

    f_in.cd()
    f_in.cd("LQ")
    name = k.GetName()
    h = ROOT.gDirectory.Get(name)
    
    f_out.cd()
    f_out.mkdir("LQ")
    f_out.cd("LQ")

    h_48bins = ROOT.TH1F(name,name,48,0,48)
    h_48bins.SetDirectory(0)
    for i in range(1,49):
        h_48bins.SetBinContent(i, h.GetBinContent(i))

    h_48bins.Write()

f_out.Close()
f_in.Close()
