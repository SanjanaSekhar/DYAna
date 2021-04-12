from ROOT import *
import numpy as np
import os
from optparse import OptionParser
from optparse import OptionGroup
from string import digits
import copy

sys_keys = ["pdf", "refac", "lep_eff", "mc_xsec", "fakes", "ptrw", "pileup", "emucostrw", "other", "btag"]

color_dict = dict()
color_dict["pdf"] = kGreen +3
color_dict["refac"] = kOrange + 7
color_dict["lep_eff"] = kRed + 2
color_dict["mc_xsec"] = kBlue
color_dict["fakes"] = kRed -7 
color_dict["other"] = kGray
color_dict["btag"] = kGreen
color_dict["ptrw"] = kOrange
color_dict["pileup"] =  kCyan
color_dict["emucostrw"] = kMagenta + 4

name_dict = dict()
name_dict["pdf"] = "PDFs"
name_dict["refac"] = "Renorm. & Fac. Scales"
name_dict["lep_eff"] = "Lepton Efficiencies"
name_dict["mc_xsec"] = "MC Cross Sections + Lumi"
name_dict["fakes"] = "Fakes Estimate"
name_dict["other"] = "Other"
name_dict["btag"] = "b-tagging"
name_dict["ptrw"] = "DY p_{T} Reweighting"
name_dict["pileup"] =  "Pileup"
name_dict["emucostrw"] = "e#mu Shape Correction"




def get_frac_diffs(h_sys, h_nom, h_tot):
    nbins = h_sys.GetNbinsX()
    avg_diff = 0.
    fracs = np.array([0.]*nbins, dtype=np.float64)
    for ibin in range(1, nbins + 1):
        c_sys = h_sys.GetBinContent(ibin)
        c_nom = h_nom.GetBinContent(ibin)
        c_tot = h_tot.GetBinContent(ibin)

        diff2 = ((c_sys - c_nom)/c_tot)**2
        fracs[ibin -1] = diff2
    return fracs


def setup_sys_dict(nBins):
    d_sys = dict()

    for key in sys_keys:
        d_sys[key] = np.array([0.]*nBins, dtype=np.float64)
    
    return d_sys

def dict_to_hists(d):
    hists = []
    for key in sys_keys:
        nBins = d[key].shape[0]
        h = TH1F("h_"+key ,key, nBins, 0.5, 0.5 + nBins)
        h.SetLineColor(color_dict[key])
        h.SetLineWidth(4)
        for i in range(1,nBins+1):
            h.SetBinContent(i, d[key][i-1])

        hists.append(h)
    return hists

def avg_sqrt(d):
    d_new = dict()
    for key in sys_keys:
        cont = d[key]
        new_cont = np.sqrt(cont)/2
        d_new[key] = new_cont
    return d_new







def add_sys(d, fracs, sys_name):
    if("pdf" in sys_name): key_name = 'pdf'
    elif('RENORM' in sys_name or 'FAC' in sys_name or 'REFAC' in sys_name): key_name = 'refac'
    elif('ID' in sys_name or 'RECO' in sys_name or 'HLT' in sys_name or 'ISO' in sys_name ): key_name = "lep_eff"
    elif("xsec" in sys_name and "fakes" not in sys_name): key_name = "mc_xsec"
    elif("fakes" in sys_name): key_name = "fakes"
    elif("ptrw" in sys_name): key_name = "ptrw"
    elif("Pu" in sys_name): key_name = "pileup"
    elif("emu" in sys_name): key_name = "emucostrw"
    elif("BTAG" in sys_name): key_name = "btag"
    else: 
        key_name = "other"
        if(np.mean(fracs) > 1e-3):
            print(sys_name, np.mean(fracs))

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
    nBins = h.GetNbinsX()
    keys = gDirectory.GetListOfKeys()


    sys_dict = setup_sys_dict(nBins)



    ee_base_strs = [ 'ee%i_fpl', 'ee%i_top', 'ee%i_db', 'ee%i_qcd', 'ee%i_gam'  ]
    mumu_base_strs = [ 'mumu%i_fpl', 'mumu%i_top', 'mumu%i_db', 'mumu%i_qcd', 'mumu%i_gam'  ]
    #combine these names into a single systematic
    removes = ['Up', 'Down', 'PTHIGH', 'PTLOW', 'BAR', 'END']
    #plus template gets factor of 3/4s in norm
    xsec_uncs = [0.03 * 0.75, 0.05, 0.04, 0.5, 0.4]
    lumi_unc = 0.025
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
                    fracs = fracs * (0.75**2)
                
                parse = (key_name.split("_"))
                sys_name = parse[2]
                add_sys(sys_dict, fracs, sys_name)


        #do xsec uncertainties
        xsec_unc = (xsec_uncs[idx]**2 + lumi_unc**2)**(0.5)
        xsec_frac = np.array([xsec_unc * h_base.Integral()/h_tot.Integral()] *nBins, dtype=np.float64)
        xsec_frac = (2*xsec_frac)**2#factor of 2 to be consistent with up/down averaging
        #print(xsec_names[idx], xsec_unc, h_base.Integral(), h_tot.Integral())
        if("qcd" in base): 
            s_key = "fakes"
        else: 
            s_key = "mc_xsec"
        add_sys(sys_dict, xsec_frac, s_key)  


    d_final = avg_sqrt(sys_dict)
    return d_final



parser = OptionParser()
parser.add_option("-m", "--mbin", type = 'int', default=1, help="mass bin")
parser.add_option("-o", "--plot_dir", default="test/", help="Plotting directory")

(options, args) = parser.parse_args()

chans= ["ee16",  "ee17", "ee18", "mumu16", "mumu17", "mumu18"]
years= [16,17,18, 16, 17,18]

gROOT.SetBatch(1)
gStyle.SetOptStat(0)





for idx,chan in enumerate(chans):
    year = years[idx]

    if('ee' in  chan): sys_dict = get_sys_dict(options.mbin, year, 'ee')
    if('mumu' in  chan): sys_dict = get_sys_dict(options.mbin, year, 'mumu')

    hists = dict_to_hists(sys_dict)

    hmaxs = [0.1, 0.1, 0.1, 0.15, 0.15, 0.2, 0.25, 0.25]

    m_bins = [150, 170, 200,  250, 320, 510, 700, 1000, 14000]
    m_low = m_bins[options.mbin]
    m_high = m_bins[options.mbin+1]

    c = TCanvas(chan, "", 1600, 1000)
    hmax = 0.

    for h in hists:
        hmax = hmaxs[options.mbin]

    h_dummy = hists[0].Clone("h_dummy")
    h_dummy.Reset()
    h_dummy.SetTitle("Channel: %s, Mass %i-%i GeV" % (chan, m_low, m_high)) 
    h_dummy.GetXaxis().SetTitle("Template Bin")
    h_dummy.GetYaxis().SetTitle("Fractional Uncertainty")
    h_dummy.SetMaximum(hmax)
    h_dummy.SetLineColorAlpha(kBlack, 0.)
    h_dummy.Draw("hist")

    for h in hists:
        h.Draw("hist same")

    leg = TLegend(0.25, 0.2)
    for h in hists:
        leg.AddEntry(h, name_dict[h.GetTitle()], "l")

    leg.Draw()
    c.Print("%s/m%i_%s_sysuncs.png" % (options.plot_dir, options.mbin, chan))






