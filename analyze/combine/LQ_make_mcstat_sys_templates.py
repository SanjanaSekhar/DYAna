import ROOT
from optparse import OptionParser
import sys

ext = "020923"
mLQ_list = [1500]


def get_sym_bin(idx, nBins):
    # 3 mass bins, 8+6+6 cost bins
    print("nBins = ",nBins)
    if nBins < 20:
	n_cos_bins = 6
        cos_bin = idx % n_cos_bins
        eta_bin = idx / n_cos_bins
        opp_cos_bin = (n_cos_bins - 1 - cos_bin) % n_cos_bins
        sym_bin = eta_bin * n_cos_bins + opp_cos_bin
    else:	
        if(idx % 20 >= 0 and idx % 20 < 8):
            n_cos_bins = 8
            cos_bin = idx % n_cos_bins
            eta_bin = idx / n_cos_bins
            opp_cos_bin = (n_cos_bins - 1 - cos_bin) % n_cos_bins
            sym_bin = eta_bin * n_cos_bins + opp_cos_bin
        else:
            n_cos_bins = 6
            cos_bin = idx % n_cos_bins
            eta_bin = idx / n_cos_bins
            opp_cos_bin = (n_cos_bins - 1 - cos_bin) % n_cos_bins
            sym_bin = eta_bin * n_cos_bins + opp_cos_bin
    
    return sym_bin


for mLQ in mLQ_list:
	for year in [2016,2017,2018]:

		fin = "combine/templates/LQm%i_merge_templates%i_%s.root"%(mLQ,year-2000,ext)
		
		for chan in ["mumu","ee"]:
			print("year is %i, chan is %s" % (year, chan))

			f = ROOT.TFile.Open(fin, "UPDATE")
			f.cd("LQ")
			keys = ROOT.gDirectory.GetListOfKeys().Clone()
			
			mcstats_procs = ["fpl","fmn","alpha","tautau","LQpure_u","LQint_u","LQpure_d","LQint_d","LQpure_u_vec","LQint_u_vec","LQpure_d_vec","LQint_d_vec"]

			
			for proc in mcstats_procs:
				
				proc_name = chan + str(year)[2:4] + '_' + mcstats_procs
				print("Creating templates for the MC stat uncs for process: %s"%proc_name)

				h = ROOT.gDirectory.Get(proc_name)
				for idx in range(1,h.GetNBinsX()+1):
					sym_bin = get_sym_bin(idx, h.GetNBinsX())
					err = h.GetBinError(idx)
					
					new_name = proc_name.replace(proc_name,  proc_name + "_MCStatBin%i"%idx)
					
					if(new_name+"Up" in keys or new_name+"Down" in keys): break
					print("Creating %s"%new_name)

					# UP template
					h_clone = h.Clone(new_name+"Up")
					h_clone.SetBinContent(idx, h.GetBinContent() + err)
					h_clone.SetBinContent(sym_bin, h.GetBinContent() + err)
					h_clone.Write()
					# DOWN template
					h_clone = h.Clone(new_name+"Down")
					h_clone.SetBinContent(idx, h.GetBinContent() - err)
					h_clone.SetBinContent(sym_bin, h.GetBinContent() - err)
					h_clone.Write()
