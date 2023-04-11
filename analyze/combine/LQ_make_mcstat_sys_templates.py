import ROOT
import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from combine.AFB_fits.scripts.LQ_utils import *

ext = "020923"

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--mLQ",  default=1000, type='int', help="mLQ")
(options, args) = parser.parse_args()

mLQ = options.mLQ


for mLQ in [5500]:
	for year in [2016,2017,2018]:

		fin = "combine/templates/LQm%i_merge_templates%i_%s.root"%(mLQ,year-2000,ext)
		f = ROOT.TFile.Open(fin, "UPDATE")
		f.cd("LQ")

		for chan in ["ee","mumu"]:

			print("making signal temps: year is %i, chan is %s, mLQ is %i" % (year, chan, mLQ))
			
			mcstats_procs = ["fpl","fmn","alpha","LQpure_u","LQint_u","LQpure_d","LQint_d","LQpure_u_vec","LQint_u_vec","LQpure_d_vec","LQint_d_vec"]

			
			for proc in mcstats_procs:
				
				proc_name = chan + str(year)[2:4] + '_' + proc
				#print("Creating templates for the MC stat uncs for process: %s"%proc_name)

				h = ROOT.gDirectory.Get(proc_name)
				for idx in range(1,h.GetNbinsX()+1):
					#keys = ROOT.gDirectory.GetListOfKeys().Clone()

					sym_bin = get_sym_bin(idx-1, h.GetNbinsX())
					if sym_bin < idx : continue
					err = h.GetBinError(idx)
					h.SetBinError(idx,1e-6) #setting nominal stat unc to zero for autmcstat	
					h.SetBinError(sym_bin,1e-6)
					new_name = proc_name.replace(proc_name,  proc_name + "_MCStatBin%i"%idx)
					new_name2 = proc_name.replace(proc_name,  proc_name + "_MCStatBin%i"%sym_bin)
					#print("New names: %s and %s"%(new_name,new_name2))	
					#if((new_name+"Up" in keys) or (new_name+"Down" in keys) or (new_name2+"Up" in keys) or (new_name2+"Down" in keys)): continue
					#print("Creating %s"%new_name)
					#flag = 0
					#for s in keys:
					#	if s.GetName()+"Up"==new_name or s.GetName()+"Down"==new_name or s.GetName()+"Up"==new_name2 or s.GetName()+"Down"==new_name2: flag = 1
					#if flag == 1: continue
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
		del f
	'''			

	for chan in ["ee","mumu"]:

		print("making bkg templates: year is %i, chan is %s, mLQ is %i"  % (year, chan, mLQ))
		
		mcstats_bkg_procs = ["top","db","gam","qcd"]

		
		for proc in mcstats_bkg_procs:
			
			proc_name = chan + str(year)[2:4] + '_' + proc
			#print("Creating templates for the MC stat uncs for process: %s"%proc_name)

			h = ROOT.gDirectory.Get(proc_name)
			for idx in range(1,h.GetNbinsX()+1):
				keys = ROOT.gDirectory.GetListOfKeys().Clone()
				#sym_bin = get_sym_bin(idx-1, h.GetNbinsX())
				err = h.GetBinError(idx)
				
				new_name = proc_name.replace(proc_name,  proc_name + "_MCStatBin%i"%idx)
				#new_name2 = proc_name.replace(proc_name,  proc_name + "_MCStatBin%i"%sym_bin)
				#print("New names: %s and %s"%(new_name,new_name2))	
				#if((new_name+"Up" in keys) or (new_name+"Down" in keys)): continue
				#print("Creating %s"%new_name)
				flag = 0
				for s in keys:
					if s.GetName()+"Up"==new_name or s.GetName()+"Down"==new_name: flag = 1
				if flag == 1: continue
				# UP template
				h_clone = h.Clone(new_name+"Up")
				h_clone.SetBinContent(idx, h.GetBinContent(idx) + err)
				#h_clone.SetBinContent(sym_bin, h.GetBinContent(sym_bin) + err)
				h_clone.Write()
				
				# DOWN template
				h_clone = h.Clone(new_name+"Down")
				h_clone.SetBinContent(idx, h.GetBinContent(idx) - err)
				#h_clone.SetBinContent(sym_bin, h.GetBinContent(sym_bin) - err)
				h_clone.Write()
	'''			
	
