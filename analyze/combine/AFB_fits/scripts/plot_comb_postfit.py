from ROOT import *
from plot_postfit import *

if (__name__ == "__main__"):
    parser = OptionParser()
    parser.add_option("--input", "-i", default = "", help="Input file")
    parser.add_option("--output", "-o", default = "", help="Input directory")
    parser.add_option("--mbin", "-m", type = 'int', default = 0, help="Mass bin (for plot label)")
    (options, args) = parser.parse_args()


    years = [2016, 2017, 2018]

    h_names = ["gam", "db", "qcd", "top",  "dy"]
    h_ss_names = ["bk", "dy", "qcd"]


    m_bins = [150, 170, 200,  250, 320, 510, 700, 1000, 14000]


    label_color_map = dict()
    label_color_map['dy'] = ("DY Signal", DY_c)
    label_color_map['top'] = ("t#bar{t} + Single Top", ttbar_c)
    label_color_map['db'] = ("WW + WZ + ZZ",  diboson_c)
    #label_color_map['tautau'] = ("DY #tau#tau Bkg.", tautau_c)
    label_color_map['gam'] = ["\\gamma\\gamma \\rightarrow {ell}{ell} ", gamgam_c]
    label_color_map['qcd'] = ("WJets + QCD", qcd_c)

    datastyle = "pe0x0"

    fracs = dict()
    for name in h_names:
        fracs[name] = 0.

    dirs = ["Y%i_mumu%i_postfit/", "Y%i_ee%i_postfit/"]
    tot_names = ["Run2_mumu_postfit", "Run2_ee_postfit"]


    f_in = TFile.Open(options.input)



    for idx, dir_name in enumerate(dirs):
        

        mbin_low = m_bins[options.mbin]
        mbin_high = m_bins[options.mbin+1]

        if(idx == 0): 
            if(mbin_low < 1000): title = "Muons %i-%i GeV" % (mbin_low, mbin_high)
            else: title = "Muons %i+ GeV" % (mbin_low)
            outname = "Run2_mumu_postfit"
        if(idx == 1): 
            if(mbin_low < 1000): title = "Electrons %i-%i GeV" % (mbin_low, mbin_high)
            else: title = "Electrons %i+ GeV" % (mbin_low)
            outname = "Run2_ee_postfit"

        f_in.cd()
        h_tot = f_in.Get(tot_names[idx])

        hist_list = [None] *len(h_names)
        color_list = []
        label_list = []

        for year in years:
            dir_ = dir_name % (year % 2000, year % 2000)
            h_data_ = f_in.Get(dir_ + "data_obs")
            #copy to new hist so poisson error bars work


            h_tot_ = f_in.Get(dir_ + "TotalProcs")
            h_tot_sig_ = f_in.Get(dir_ + "TotalSig")


            if(year == 2016):
                h_data = h_data_.Clone("h_data_c%i_comb" %(idx))
                h_tot_dir = h_tot_.Clone("h_tot_sig_c%i_comb" %(idx))
            else:
                h_data.Add(h_data_)
                h_tot_dir.Add(h_tot_)


            

            for h_idx, name in enumerate(h_names):
                if(name == "dy"):
                    h = h_tot_sig_.Clone("h_%s_c%i_y%i" %(name, idx, year))
                else:
                    h = f_in.Get(dir_ + name)
                    if(h != None):
                        h = h.Clone("h_%s_c%i_y%i" %(name, idx, year))
                if(h != None):
                    #print(hist_list)
                    if(year == 2016):
                        hist_list[h_idx] = h
                        if(name == 'gam'):
                            if(idx == 0): 
                                lab = label_color_map[name][0].format(ell='\\mu')
                            else:
                                lab = label_color_map[name][0].format(ell='e')
                        else:
                            lab = label_color_map[name][0]

                        label_list.append(lab)
                        color_list.append(label_color_map[name][1])
                    else:
                        hist_list[h_idx].Add(h)


        #copy data hist to get poisson errors to work
        h_data_pois = h_data.Clone("h_data_c%i_y%i" %(idx, year))
        h_data_pois.Reset()
        for b in range(h_data.GetXaxis().GetNbins() + 1):
            h_data_pois.SetBinContent(b, h_data.GetBinContent(b))

        #h_data_pois.Print("range")

        #h_tot.Print("range")
        #h_tot_dir.Print("range")

        for b in range(h_tot.GetXaxis().GetNbins() + 1):
            h_tot_dir.SetBinError(b, h_tot.GetBinError(b))

        if(options.mbin <=2):
            ratio_range = (0.91, 1.09)
            NDiv = 205
        elif(options.mbin <=4):
            ratio_range = (0.851, 1.149)
            NDiv = 205
        elif(options.mbin <=6):
            ratio_range = (0.6, 1.4)
            NDiv = 203
        elif(options.mbin == 7):
            ratio_range = (0.2, 1.8)
            NDiv = 203

        makeCan(outname, options.output, [h_data], bkglist=[hist_list], totlist=[h_tot_dir], colors = color_list, bkgNames = label_list, 
                titles = [title], xtitle = "Template Bin", year = -1, datastyle=datastyle, mbin = options.mbin, ratio_range = ratio_range, NDiv = NDiv) 

