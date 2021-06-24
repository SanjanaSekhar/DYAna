from ROOT import *
import os
from optparse import OptionParser
from optparse import OptionGroup
import numpy as np

parser = OptionParser()
parser.add_option("-i", "--fin", default='', help="Input file with templates")
parser.add_option("-y", "--year", type = 'int', default=16, help="Year")
parser.add_option("-o", "--plot_dir", default='../plots/', help="Directory to output plots")


def get_sym_bin(idx, nBins):
    if(idx < 16): #2x8 cos(theta) bins
        n_cos_bins = 8
        cos_bin = idx % n_cos_bins
        eta_bin = idx / n_cos_bins
        opp_cos_bin = (n_cos_bins - 1 - cos_bin) % n_cos_bins
        sym_bin = eta_bin * n_cos_bins + opp_cos_bin
    else: #eta bins have 6 cos bins (one or two eta bins depending on mass)
        n_cos_bins = 6
        diff_idx = idx - 16
        cos_bin = diff_idx % n_cos_bins
        eta_bin = diff_idx / n_cos_bins
        opp_cos_bin = (n_cos_bins - 1 - cos_bin) % n_cos_bins
        sym_bin = 16 + eta_bin*n_cos_bins + opp_cos_bin

    return sym_bin


(options, args) = parser.parse_args()

n_m_bins = 8


chan_names = ["mumu%i_", "ee%i_"]
sig_names = [ 'fpl', 'fmn', 'alpha'] 
sym_bkg_names = ['top',  'qcd' ]
nonsym_bkg_names = ['db', 'tautau',  'gam' ]

gROOT.SetBatch(1)

f = TFile.Open(options.fin)
gDirectory.ls()

AFB = 0.6
A0 = 0.05

if(not os.path.exists(options.plot_dir)):
    print("Making directory %s" % options.plot_dir)
    os.system("mkdir %s" % options.plot_dir)


def add_corr_hist(h_tot, corr_mats, h):
    n_bins = h_tot.GetNbinsX()
    for i in range(n_bins):
        err_i = h.GetBinError(i+1)
        opp_i = get_sym_bin(i, n_bins)
        err_opp_i = h.GetBinError(opp_i+1)
        corr_mats[i][0][0] += err_i**2
        corr_mats[i][0][1] += err_i*err_opp_i
        corr_mats[i][1][0] += err_i*err_opp_i
        corr_mats[i][1][1] += err_opp_i**2

    h_tot.Add(h)

def add_uncorr_hist(h_tot, corr_mats, h):
    n_bins = h_tot.GetNbinsX()
    for i in range(n_bins):
        err_i = h.GetBinError(i+1)
        opp_i = get_sym_bin(i, n_bins)
        err_opp_i = h.GetBinError(opp_i+1)
        corr_mats[i][0][0] += err_i**2
        corr_mats[i][1][1] += err_opp_i**2

    h_tot.Add(h)

def correlation_coeff(mat):
    sigma1 = np.sqrt(mat[0][0])
    sigma2 = np.sqrt(mat[1][1])
    covariance = mat[0][1]
    return covariance/sigma1/sigma2




for chan_ in chan_names:
    chan = chan_ % options.year

    for mbin in range(1, n_m_bins):
        f.cd("w%i" % mbin)
        keys = gDirectory.GetListOfKeys()
        h_pl = gDirectory.Get(chan + "fpl")
        h_mn = gDirectory.Get(chan + "fmn")
        h_alpha = gDirectory.Get(chan + "alpha")

        h_sig = h_pl.Clone("h_sig")
        h_sig.Reset()
        h_tot = h_sig.Clone("h_tot")
        h_tot.Reset()
        n_bins = h_tot.GetNbinsX()
        corr_mats = [np.array([[0.,0.],[0.,0.,]]) for i in range(n_bins)]

        alpha = 2.*A0/(2-A0)
        norm = 3./4./(2+alpha)
        n_pl = (norm + AFB)
        n_mn = (norm - AFB)
        n_alpha = norm * alpha


        h_sig.Add(h_pl, h_mn, n_pl, n_mn)
        h_alpha.Scale(n_alpha)
        h_sig.Add(h_alpha)

        add_corr_hist(h_tot, corr_mats,h_sig)

        for name in sym_bkg_names:
            h_name = chan + name
            h = gDirectory.Get(h_name)
            add_corr_hist(h_tot, corr_mats, h)

        for name in nonsym_bkg_names:
            h_name = chan + name
            h = gDirectory.Get(h_name)
            add_uncorr_hist(h_tot, corr_mats, h)
            h.Print("range")


        h_tot.Print("range")

        h_corr = h_tot.Clone("h_corr")
        h_corr.Reset()

        for i in range(n_bins):
            corr_coeff = correlation_coeff(corr_mats[i])
            val = 1. - corr_coeff
            tot_err = h_tot.GetBinError(i+1)
            h_corr.SetBinContent(i+1, val)
        h_corr.Print("range")

        
        c1 = TCanvas("c1", "", 1600, 1000) 
        if(mbin >= 5): h_corr.SetMaximum(0.6)
        else: h_corr.SetMaximum(0.1)
        h_corr.SetTitle("Correlation Coeffecients: %s Mass Bin %i" % (chan[:-1], mbin))
        h_corr.GetXaxis().SetTitle("Template Bin")
        h_corr.GetYaxis().SetTitle("1 - Corr. Coeff.")
        h_corr.SetLineColor(kBlue)
        h_corr.Draw("hist")
        c1.Print(options.plot_dir + chan + ("mbin%i_" % mbin) + "corr_coeffs.png")


