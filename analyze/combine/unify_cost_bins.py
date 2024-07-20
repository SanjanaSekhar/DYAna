import ROOT
from optparse import OptionParser
import sys

mLQ = 2500
year = 2018
ext = "071924_splitrap"
fin = "../templates/LQm%i_merge_templates%i_%s.root"%(mLQ,year-2000,ext)
f_in = ROOT.TFile.Open(fin, "UPDATE")
f_in.cd()
f_in.cd("LQ")
keys = ROOT.gDirectory.GetListOfKeys().Clone()

for idx,k in enumerate(keys[:3]):

    f_in.cd()
    f_in.cd("LQ")
    name = k.GetName()
    print("Fixing cost bins in ", name)
    h = ROOT.gDirectory.Get(name)

    h_splitrap = ROOT.TH1F(name+"splitrap",name+"splitrap",108,0,108)
    h_splitrap.SetDirectory(0)

    j=1

    for i in range(1,h.GetNbinsX()+1):
        if i%42 <= 24:
            if i%8==1 or i%8==7: # 1st+2nd OR 7th+8th
                binc = h.GetBinContent(i) + h.GetBinContent(i+1)
                print("unifying bins %i and %i"%(i,i+1))
                i+=1
            else: 
                binc = h.GetBinContent(i)
                print("bin ", i)
        else: 
            binc = h.GetBinContent(i)
            print("bin ",i)         
        h_splitrap.SetBinContent(j,binc)
        j+=1  

    h_splitrap.Write()


f_in.Close()
