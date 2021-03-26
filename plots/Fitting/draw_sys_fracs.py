from ROOT import *
import numpy as np
import os
from optparse import OptionParser
from optparse import OptionGroup
from string import digits
import copy

sys_keys = ["pdf", "refac", "lep_eff", "mc_xsec", "fakes", "ptrw", "pileup", "emucostrw", "other"]

color_dict = dict()
color_dict["pdf"] = kGreen +3
color_dict["recac"] = kOrange + 7
color_dict["lep_eff"] = kRed
color_dict["mc_xsec"] = kBlue
color_dict["fakes"] = kRed -7 
color_dict["other"] = kGray
color_dict["pileup"] =  kCyan
color_dict["emucostrw"] = kMagenta + 4



def get_frac_diffs(h_sys, h_nom, h_tot):
    nbins = h_sys.GetNbinsX()
    avg_diff = 0.
    fracs = np.array([0]*nbins)
    for ibin in range(1, nbins + 1):
        c_sys = h_sys.GetBinContent(ibin)
        c_nom = h_nom.GetBinContent(ibin)
        c_tot = h_tot.GetBinContent(ibin)

        diff = abs(c_sys - c_nom)/c_tot
        fracs[ibin -1] = diff 
    return fracs


def setup_dict(nBins):
    d_sys = dict()

    for key in sys_keys:
        d_sys[key] = np.array([0]*nBins)
    
    return d_sys

def dict_to_hists(d, nBins):
    hists = []
    for key in sys_keys:
        h = TH1F("h_"+key,key, nBins, 0.5, 0.5 + nBins)
        h.SetLineColor(color_dict[key])
        h.SetLineWidth(1)
        for i in range(i,nBins+1):
            h.SetBinContent(i, d[key][i])

        hists.append[h]
    return hists





def add_sys(d, fracs, sys_name):
    if("pdf" in sys_name): key_name = 'pdf'
    elif('RENORM' in sys_name or 'FAC' in sys_name or 'REFAC' in sys_name): key_name = 'refac'
    elif('ID' in sys_name or 'RECO' in sys_name or 'HLT' in sys_name or 'ISO' in sys_name ): key_name = "lep_eff"
    elif("xsec" in sys_name and "fakes" not in sys_name): key_name = "mc_xsec"
    elif("fakes" in sys_name): key_name = "fakes"
    elif("ptrw" in sys_name): key_name = "ptrw"
    elif("pu" in sys_name): key_name = "pileup"
    elif("emu" in sys_name): key_name = "emucostrw"
    else: key_name = "other"

    d[key_name] += fracs



def get_sys_dict(mbin, year, chan):

    f_tot_name = "../analyze/combine/AFB_fits/postfit_plots/combined_mbin%i/combined_fit_shapes_mbin%i.root" % (mbin, mbin)
    f_tot = TFile.Open(f_tot_name)
    prefit_dir = "Y%i_%s%i_prefit" % (year, chan,  year)

    f_tot.cd(prefit_dir)
    h_tot= gDirectory.Get('TotalProcs')


    f_in_name = "../analyze/combine/AFB_fits/templates%i.root" % year

    f = TFile.Open(f_in_name)
    gDirectory.cd("w%i" % mbin)
    h = gDirectory.Get("ee%i_fpl" % year)
    nBins = h.GetNBinsX()
    keys = gDirectory.GetListOfKeys()


    sys_dict = setup_sys_dict(nBins)



    ee_base_strs = [ 'ee%i_fpl', 'ee%i_top', 'ee%i_db', 'ee%i_qcd', 'ee%i_gam'  ]
    mumu_base_strs = [ 'mumu%i_fpl', 'mumu%i_top', 'mumu%i_db', 'mumu%i_qcd', 'mumu%i_gam'  ]
    #combine these names into a single systematic
    removes = ['Up', 'Down', 'PTHIGH', 'PTLOW', 'BAR', 'END']
    #plus template gets factor of 3/4s in norm
    xsec_uncs = [0.03 * 0.75, 0.05, 0.04, 0.5, 0.4]
    mumu_xsec_names = ['DY_xsec', 'top_xsec', 'diboson_xsec', 'mumu_fakes_xsec', 'gamgam_xsec']
    ee_xsec_names = ['DY_xsec', 'top_xsec', 'diboson_xsec', 'ee_fakes_xsec', 'gamgam_xsec']
    my_excludes = []

    if(chan == 'ee'):
        base_strs = ee_base_strs
        xsec_names = ee_xsec_names
    else:
        base_strs = mumu_base_strs
        xsec_names = mumu_xsec_names


    for idx,base in enumerate(base_strs):
        if('%i' in base):
            base = base % year
        #print("Doing base %s" % base)
        h_base = gDirectory.Get(base)
        for key in keys:
            key_name = key.GetName()
            if(key_name == base): continue
            skip = False
            for exc in my_excludes:
                if (exc in key_name):
                    skip = True
            if('pdf' in key_name and ('fpl' not in base and 'fmn' not in base)): skip = True
            if(('RENORM' in key_name or 'FAC' in key_name) and ('fpl' not in base and 'fmn' not in base)): skip = True
            if(skip): continue
            if (base in key_name):
                #print("Adding key %s" % key_name)
                h = gDirectory.Get(key_name)
                fracs = get_frac_diffs(h, h_base, h_tot)
                #plus templates get normalization factor of 0.75 in fit
                if('fpl' in base):
                    fracs *= 0.75
                
                #average up and down templates
                fracs *= 0.5

                parse = (key_name.split("_"))
                sys_name = parse[2]
                add_sys(d_fracs, fracs, sys_name)


        #do xsec uncertainties
        #xsec_unc = xsec_uncs[idx]
        #change = xsec_unc * h_base.Integral()/h_tot.Integral()
        ##print(xsec_names[idx], xsec_unc, h_base.Integral(), h_tot.Integral())
        #raw_sys_dict[xsec_names[idx]] = change


    return d_fracs


def add_sys(d, fracs, sys_name):
    if("pdf" in sys_name): key_name = 'pdf'
    elif('RENORM' in sys_name or 'FAC' in sys_name or 'REFAC' in sys_name): key_name = 'refac'
    elif('ID' in sys_name or 'RECO' in sys_name or 'HLT' in sys_name or 'ISO' in sys_name ): key_name = "lep_eff"
    elif("xsec" in sys_name and "fakes" not in sys_name): key_name = "mc_xsec"
    elif("fakes" in sys_name): key_name = "fakes"
    elif("ptrw" in sys_name): key_name = "ptrw"
    elif("pu" in sys_name): key_name = "pileup"
    elif("emu" in sys_name): key_name = "emucostrw"
    else: key_name = "other"

    d[key_name] += fracs


parser = OptionParser()
parser.add_option("-m", "--mbin", type = 'int', default=1, help="mass bin")

(options, args) = parser.parse_args()

chans= ["ee16",  "ee17", "ee18", "mumu16", "mumu17", "mumu18"]
years= [16,17,18, 16, 17,18]

gROOT.SetBatch(1)





for idx,chan in enumerate(chans):
    year = years[idx]

    if('ee' in  chan): sys_dict = get_sys_dict(options.mbin, year, 'ee')
    if('mumu' in  chan): sys_dict = get_sys_dict(options.mbin, year, 'mumu')

    hists = dict_to_hists(sys_dict)

    c = TCanvas(chan, "", 1600, 1000)

    hists[0].GetXaxis.SetTitle("Template Bin")
    hists[0].GetYaxis.SetTitle("Fractional Uncertainty")

    for h in hists:
        h.Draw("hist same")

    c.Print("test.png")
    break





