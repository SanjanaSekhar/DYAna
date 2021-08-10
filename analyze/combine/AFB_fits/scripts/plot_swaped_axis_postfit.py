from ROOT import *
from plot_postfit import *

def convert22bin(idx):
#8+8+6 bins -> 6x3 bins
    n_eta_bins = 3
    cos_maxbin=5
    if(idx <= 16): #merge edge cos(theta) bins
        cos_bin = max(0, min(((idx-1) % 8) -1, cos_maxbin))
        eta_bin = (idx-1) / 8
    else:
        cos_bin = (idx - 1 - 16) % 6
        eta_bin = 2 + (idx - 1 -16) /6

    gbin = (n_eta_bins * cos_bin) + eta_bin + 1
    return gbin


def convert28bin(idx):
#8+8+6+6 bins -> 6x4 bins
    n_eta_bins = 4
    cos_maxbin=5
    if(idx <= 16): #merge edge cos(theta) bins
        cos_bin = max(0, min(((idx-1) % 8) -1, cos_maxbin))
        eta_bin = (idx-1) / 8
    else:
        cos_bin = (idx - 1 - 16) % 6
        eta_bin = 2 + (idx - 1 -16) /6


    gbin = (n_eta_bins * cos_bin) + eta_bin + 1
    return gbin



def swap_axes(h):

    nBins = h.GetNbinsX()
    if(nBins == 22):
        coversion_fn = convert22bin
        nBins_new = 18
    elif(nBins == 28):
        coversion_fn = convert28bin
        nBins_new = 24
    else:
        print("Hist has %i bins? Can only handle 22 or 28. Exiting" % nBins)
        sys.exit(1)

    h_new = TH1F(h.GetName() + "swap", "", nBins_new, 0, nBins_new)

    for idx in range(1,nBins+1):
        new_bin = coversion_fn(idx)

        cont = h.GetBinContent(idx)
        err = h.GetBinError(idx)

        old_cont = h_new.GetBinContent(new_bin)
        old_err = h_new.GetBinError(new_bin)

        new_cont = cont + old_cont
        new_err = (err**2 + old_err**2)**(0.5)
        h_new.SetBinContent(new_bin, new_cont)
        h_new.SetBinError(new_bin, new_err)

    return h_new
        





parser = OptionParser()
parser.add_option("--input", "-i", default = "", help="Input file")
parser.add_option("--output", "-o", default = "", help="output directory")
parser.add_option("--mbin", "-m", type = 'int', default = 0, help="Mass bin (for plot label)")
parser.add_option("--year", "-y", type = 'int', default = -1, help="Year (-1 for all) ")
parser.add_option("--ss",   default = False, action='store_true',  help="Fit was done with ee_ss region too")
(options, args) = parser.parse_args()


#fin_ = "combined_fit_shapes_mbin1.root"
#odir = "postfit_plots/combined_fit_mbin1"
#mbin = 1
if(options.year < 0):
    years = [2016, 2017, 2018]
else:
    years = [options.year]

h_names = ["gam", "db", "qcd", "top",  "tautau", "dy"]
h_ss_names = ["bk", "dy", "qcd"]


m_bins = [150, 170, 200,  250, 320, 510, 700, 1000, 14000]


label_color_map = dict()
label_color_map['dy'] = ("DY Signal", kRed + 1)
label_color_map['top'] = ("t#bar{t} + tW ", kBlue)
label_color_map['db'] = ("WW + WZ + ZZ",  kGreen +3)
label_color_map['tautau'] = ("DY #tau#tau Bkg.", kMagenta + 4)
label_color_map['gam'] = ("\\gamma\\gamma \\rightarrow \\ell\\ell ", kOrange)
label_color_map['qcd'] = ("WJets + QCD", kRed - 7)

datastyle = "pe0x0"

dirs = ["Y%i_mumu%i_postfit/", "Y%i_ee%i_postfit/"]
if(options.ss): dirs = dirs = ["Y%i_mumu%i_postfit/", "Y%i_ee%i_postfit/", "Y%i_ee%i_ss_postfit/"]
f_in = TFile.Open(options.input)
gam_frac_avg = 0.
for year in years:
    for idx, dir_name in enumerate(dirs):
        dir_ = dir_name % (year % 2000, year % 2000)
        h_tot = f_in.Get(dir_ + "TotalProcs")
        h_tot = h_tot.Clone("h_tot_c%i_y%i" %(idx, year))
        h_data = f_in.Get(dir_ + "data_obs")
        h_data = h_data.Clone("h_data_c%i_y%i" %(idx, year))

        h_tot_swap = swap_axes(h_tot)
        h_data_swap = swap_axes(h_data)

        h_tot_sig = f_in.Get(dir_ + "TotalSig")
        h_tot_sig = h_tot_sig.Clone("h_tot_sig_c%i_y%i" %(idx, year))

        mbin_low = m_bins[options.mbin]
        mbin_high = m_bins[options.mbin+1]

        title = ""
        if(idx == 0): title = "Muons %i %i-%i GeV" % (year, mbin_low, mbin_high)
        elif(idx == 1): title = "Electrons %i %i-%i GeV" % (year, mbin_low, mbin_high)
        elif(idx == 2): title = "Electrons Samesign %i %i-%i GeV" % (year, mbin_low, mbin_high)
        else:
            print("Idx is %i ?" % idx)
        
        if(idx == 2): name_list = h_ss_names
        else: name_list = h_names
        hist_list = []
        color_list = []
        label_list = []

        for name in name_list:
            if(name == "dy"):
                h = h_tot_sig.Clone("h_%s_c%i_y%i" %(name, idx, year))

            else:
                h = f_in.Get(dir_ + name)
                if(h != None):
                    h = h.Clone("h_%s_c%i_y%i" %(name, idx, year))
            h_swap = swap_axes(h)
            if(h != None):
                #h.Print()
                hist_list.append(h_swap)
                label_list.append(label_color_map[name][0])
                color_list.append(label_color_map[name][1])

        makeCan(dir_[:-1] + "swapped_axis", options.output, [h_data_swap], bkglist=[hist_list], totlist=[h_tot_swap], 
                colors = color_list, bkgNames = label_list, titles = [title], xtitle = "Swapped Template Bin", year = year, mbin = options.mbin ) 



