import ROOT
from optparse import OptionParser
import sys

ext = "020923"
mLQ_list = [1500,2000,2500]


def get_sym_bin(idx, nBins):
    # 3 mass bins, 8+6+6 cost bins
    #print("nBins = ",nBins)
    if nBins < 20:
	n_cos_bins = 6
        cos_bin = idx % n_cos_bins
        eta_bin = idx / n_cos_bins
        opp_cos_bin = (n_cos_bins - 1 - cos_bin) % n_cos_bins
        sym_bin = eta_bin * n_cos_bins + opp_cos_bin
    else:	
        if((idx) % 20 >= 0 and (idx) % 20 < 8):
            n_cos_bins = 8
	    idx_new = idx - 20*(idx/20)
            cos_bin = idx_new % n_cos_bins
            eta_bin = idx_new / n_cos_bins
            opp_cos_bin = (n_cos_bins -1 - cos_bin) % n_cos_bins
            sym_bin = 20*(idx/20) + eta_bin * n_cos_bins + opp_cos_bin
	    #else: sym_bin = eta_bin * n_cos_bins + opp_cos_bin
        else:
            n_cos_bins = 6
	    idx_new = idx - 8 - 20*(idx/20)
            cos_bin = idx_new % n_cos_bins
            eta_bin = idx_new / n_cos_bins
            opp_cos_bin = (n_cos_bins -1 - cos_bin) % n_cos_bins
            sym_bin = 8 + 20*(idx/20) + eta_bin * n_cos_bins + opp_cos_bin
    
    return sym_bin+1


for mLQ in mLQ_list:
	for year in [2016,2017,2018]:

		fin = "combine/templates/LQm%i_merge_templates%i_%s.root"%(mLQ,year-2000,ext)
		f = ROOT.TFile.Open(fin, "UPDATE")
		f.cd("LQ")

		for chan in ["ee","mumu"]:

			print("year is %i, chan is %s" % (year, chan))
			
			mcstats_procs = ["fpl","fmn","alpha","tautau","LQpure_u","LQint_u","LQpure_d","LQint_d","LQpure_u_vec","LQint_u_vec","LQpure_d_vec","LQint_d_vec"]

			
			for proc in mcstats_procs:
				
				proc_name = chan + str(year)[2:4] + '_' + proc
				#print("Creating templates for the MC stat uncs for process: %s"%proc_name)

				h = ROOT.gDirectory.Get(proc_name)
				for idx in range(1,h.GetNbinsX()+1):
					keys = ROOT.gDirectory.GetListOfKeys().Clone()
					sym_bin = get_sym_bin(idx-1, h.GetNbinsX())
					err = h.GetBinError(idx)
					
					new_name = proc_name.replace(proc_name,  proc_name + "_MCStatBin%i"%idx)
					new_name2 = proc_name.replace(proc_name,  proc_name + "_MCStatBin%i"%sym_bin)
					#print("New names: %s and %s"%(new_name,new_name2))	
					if((new_name+"Up" in keys) or (new_name+"Down" in keys) or (new_name2+"Up" in keys) or (new_name2+"Down" in keys)): continue
					print("Creating %s"%new_name)

					# UP template
					h_clone = h.Clone(new_name+"Up")
					h_clone.SetBinContent(idx, h.GetBinContent(idx) + err)
					h_clone.SetBinContent(sym_bin, h.GetBinContent(sym_bin) + err)
					h_clone.Write()
					# DOWN template
					h_clone = h.Clone(new_name+"Down")
					h_clone.SetBinContent(idx, h.GetBinContent(idx) - err)
					h_clone.SetBinContent(sym_bin, h.GetBinContent(sym_bin) - err)
					h_clone.Write()
		f.Close()
