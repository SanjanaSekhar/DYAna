import ROOT
from optparse import OptionParser
import sys
import numpy as np

ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 6001;")

mLQ = 2500
year = 2018
ext = "071924_splitrap"
fin = "combine/templates/LQm%i_merge_templates%i_%s.root"%(mLQ,year-2000,ext)
f_in = ROOT.TFile.Open(fin, "UPDATE")
f_in.cd()
f_in.cd("LQ")
keys = ROOT.gDirectory.GetListOfKeys().Clone()

for idx,k in enumerate(keys):

    f_in.cd()
    f_in.cd("LQ")
    name = k.GetName()
    #print("Fixing cost bins in ", name)
    h = ROOT.gDirectory.Get(name)

    h_splitrap = ROOT.TH1F(name,name,108,0,108)
    h_splitrap.SetDirectory(0)

    j=1
    flag = 0
    for i in range(1,h.GetNbinsX()+1):
        if i%42 <= 24:
            if i%8==1 or i%8==7: # 1st+2nd OR 7th+8th
                binc = h.GetBinContent(i) + h.GetBinContent(i+1)
		bine = np.sqrt(h.GetBinError(i)**2 + h.GetBinError(i+1)**2)
                #print("unifying bins %i and %i"%(i,i+1))
		flag = 1
            else: 
                if flag==0: 
		    binc = h.GetBinContent(i)
		    bine = h.GetBinError(i)
                    #print("bin ", i)
		else: flag = 0
        else: 
            binc = h.GetBinContent(i)
	    bine = h.GetBinError(i)
            #print("bin ",i)         
        h_splitrap.SetBinContent(j,binc)
	h_splitrap.SetBinError(j,bine)
        j+=1  

    h_splitrap.Write("",ROOT.TObject.kOverwrite)


f_in.Close()
