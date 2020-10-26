from ROOT import *
import os
import sys

def hist_equal(h1,h2):
    for i in range(1, h1.GetNbinsX()+1):
        for j in range(1, h1.GetNbinsY()+1):
            if(abs(h1.GetBinContent(i,j) - h2.GetBinContent(i,j)) > 0.5):
                return False
    return abs(h1.Integral() - h2.Integral()) < 0.1

f1_name = sys.argv[1]
f2_name = sys.argv[2]

f1 = TFile.Open(f1_name)
f2 = TFile.Open(f2_name)

gDirectory.ls()

n_m_bins = 8

my_excludes = ['pdf', 'top', 'db']



for mbin in range(n_m_bins):
    if(mbin < 7): continue
    f1.cd()
    gDirectory.cd("w%i" % mbin)
    keys = gDirectory.GetListOfKeys()
    for key in keys:
        key_name = key.GetName()
        skip = False
        for exc in my_excludes:
            if (exc in key_name):
                skip = True
        if(skip): continue
        #print("Adding key %s" % key_name)
        f1.cd()
        gDirectory.cd("w%i" % mbin)
        h1 = gDirectory.Get(key_name)
        f2.cd()
        gDirectory.cd("w%i" % mbin)
        h2 = gDirectory.Get(key_name)

        if(not hist_equal(h1,h2)):
            print("Mbin %i: Key %s not equal between two templates" % (mbin, key.GetName()))
            h1.Print("range")
            h2.Print("range")

