import ROOT
from optparse import OptionParser
import sys

mLQ = 2500
year = 2018
ext = "020923"
fin = "../templates/LQm%i_merge_templates%i_%s.root"%(mLQ,year-2000,ext)
f_in = ROOT.TFile.Open(fin, "READ")
fout = "../templates/LQm%i_48bins_templates%i_%s.root"%(mLQ,year-2000,ext)
f_out = ROOT.TFile.Open(fout, "RECREATE")
f_in.cd()
f_in.cd("LQ")
keys = ROOT.gDirectory.GetListOfKeys().Clone()

for i,k in enumerate(keys):

    f_in.cd()
    f_in.cd("LQ")
    name = k.GetName()
    h = ROOT.gDirectory.Get(name)
    
    f_out.cd()
    if i==0: f_out.mkdir("LQ")
    f_out.cd("LQ")

    h_48bins = ROOT.TH1F(name,name,48,0,48)
    h_48bins.SetDirectory(0)
    for i in range(1,49):
        h_48bins.SetBinContent(i, h.GetBinContent(i))
	h_48bins.SetBinError(i, h.GetBinError(i))
    h_48bins.Write()

f_out.Close()
f_in.Close()
