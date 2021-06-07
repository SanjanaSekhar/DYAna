from ROOT import *
import os
from optparse import OptionParser
from optparse import OptionGroup
from string import digits
import copy


def frac_change(h_sys, h_nom, h_tot):
    nbins = h_sys.GetNbinsX()
    avg_diff = 0.
    for ibin in range(1, nbins + 1):
        c_sys = h_sys.GetBinContent(ibin)
        c_nom = h_nom.GetBinContent(ibin)
        c_tot = h_tot.GetBinContent(ibin)

        diff = abs(c_sys - c_nom)/c_tot
        avg_diff += diff
    avg_diff /= nbins
    return avg_diff



def get_sys_dict(mbin, year, chan):

    f_tot_name = "../analyze/combine/AFB_fits/postfit_plots/combined_mbin%i/combined_fit_shapes_mbin%i.root" % (mbin, mbin)
    f_tot = TFile.Open(f_tot_name)
    prefit_dir = "Y%i_%s%i_prefit" % (year, chan,  year)

    f_tot.cd(prefit_dir)
    h_tot= gDirectory.Get('TotalProcs')


    f_in_name = "../analyze/combine/AFB_fits/templates%i.root" % year

    f = TFile.Open(f_in_name)
    gDirectory.cd("w%i" % mbin)
    keys = gDirectory.GetListOfKeys()


    raw_sys_dict = dict()
    sys_dict = dict()

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
                change = frac_change(h, h_base, h_tot)
                #plus templates get normalization factor of 0.75 in fit
                if('fpl' in base):
                    change *= 0.75 
                parse = (key_name.split("_"))
                sys_name = parse[2]

                #if('ptrw' in key_name):
                #    print(key_name, change)



                if(sys_name in raw_sys_dict.keys()):
                    raw_sys_dict[sys_name] += change
                else:
                    raw_sys_dict[sys_name] = change

        #do xsec uncertainties
        xsec_unc = xsec_uncs[idx]
        change = xsec_unc * h_base.Integral()/h_tot.Integral()
        #print(xsec_names[idx], xsec_unc, h_base.Integral(), h_tot.Integral())
        raw_sys_dict[xsec_names[idx]] = change


    for key in raw_sys_dict.keys():
        sys_ = key
        for r in removes:
            sys_ = sys_.replace(r, '')
        sys = sys_.translate(None, digits) # remove numbers
        #halve effect because each sys has up and down versions
        if('xsec' not in key): fac = 0.5
        else: fac = 1.

        if(sys in sys_dict.keys()):
            sys_dict[sys] += fac * raw_sys_dict[key]
        else:
            sys_dict[sys] = fac * raw_sys_dict[key]

    return sys_dict


parser = OptionParser()
parser.add_option("-m", "--mbin", type = 'int', default=1, help="mass bin")

(options, args) = parser.parse_args()

years= [16,17,18]

gROOT.SetBatch(1)




year_sys_dict = [dict() for y in years]


for idx, year in enumerate(years):

    ee_sys_dict = get_sys_dict(options.mbin, year, 'ee')
    mumu_sys_dict = get_sys_dict(options.mbin, year, 'mumu')

    year_sys_dict[idx] = copy.copy(mumu_sys_dict)
    for key in ee_sys_dict.keys():
        if key not in mumu_sys_dict.keys() or mumu_sys_dict[key] < 1e-6:
            year_sys_dict[idx][key] = ee_sys_dict[key]
        else:
            year_sys_dict[idx][key] = (ee_sys_dict[key] + mumu_sys_dict[key])/2.

    print("Year %i:" % year)
    keys = year_sys_dict[idx].keys()
    keys.sort()
    for ent in keys:
        print ent, year_sys_dict[idx][ent]

final_sys_dict = dict()

for key in year_sys_dict[0].keys():
    final_sys_dict[key] = (year_sys_dict[0][key] + year_sys_dict[1][key] + year_sys_dict[2][key]) / 3.

print("\n\nAveraged:\n")
keys = final_sys_dict.keys()
keys.sort()
for ent in keys:
    print ent, final_sys_dict[ent]

