import ROOT
import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
#from AFB_fits.scripts.LQ_utils import *
from math import sqrt

ext = "020923"

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--mLQ",  default=1000, type='int', help="mLQ")
(options, args) = parser.parse_args()

mLQ = options.mLQ
chan = "ee"

for mLQ in [2500]:
	for year in [2016,2017,2018]:

		fin = "combine/templates/LQm%i_merge_templates%i_%s.root"%(mLQ,year-2000,ext)
		f = ROOT.TFile.Open(fin, "UPDATE")
		f.cd("LQ")


		print("making signal temps: year is %i, chan is %s, mLQ is %i" % (year, chan, mLQ))
		
		hilow_procs = ["elHLTBARPT", "elIDBARPT", "elRECOBARPT", "elHLTENDPT", "elIDENDPT", "elRECOENDPT"]

		procs = ['top','db','fpl','fmn','alpha','tautau','gam','LQpure_u','LQint_u','LQpure_d','LQint_d','LQpure_u_vec','LQint_u_vec','LQpure_d_vec','LQint_d_vec']
		
		for proc in hilow_procs:
			for p in procs:
				hi_up = ROOT.gDirectory.Get("%s%s_%s_%sHIGH%sUp"%(chan,year-2000,p,proc,year-2000))
				hi_down = ROOT.gDirectory.Get("%s%s_%s_%sHIGH%sDown"%(chan,year-2000,p,proc,year-2000))
				lo_up = ROOT.gDirectory.Get("%s%s_%s_%sLOW%sUp"%(chan,year-2000,p,proc,year-2000))
				lo_down = ROOT.gDirectory.Get("%s%s_%s_%sLOW%sDown"%(chan,year-2000,p,proc,year-2000))
				h = ROOT.gDirectory.Get("%s%s_%s"%(chan,year-2000,p))
				new_name_up = "%s%s_%s_%s%sUp"%(chan,year-2000,p,proc,year-2000)
				new_name_down = "%s%s_%s_%s%sDown"%(chan,year-2000,p,proc,year-2000)				
				h_clone_up = h.Clone(new_name_up)
				h_clone_down = h.Clone(new_name_down)
				for idx in range(1,h.GetNbinsX()+1):
					#keys = ROOT.gDirectory.GetListOfKeys().Clone()
					nom = h.GetBinContent(idx)
					nom_e = h.GetBinError(idx)
					new_up = hi_up.GetBinContent(idx) + lo_up.GetBinContent(idx) - nom
					new_down =  hi_down.GetBinContent(idx) + lo_down.GetBinContent(idx) - nom
					new_up_e = sqrt(hi_up.GetBinError(idx)**2 + lo_up.GetBinError(idx)**2)
					new_down_e = sqrt(hi_down.GetBinError(idx)**2 + lo_down.GetBinError(idx)**2)
					h_clone_up.SetBinContent(idx, new_up)
					h_clone_up.SetBinError(idx, new_up_e)
					h_clone_down.SetBinContent(idx, new_down)
					h_clone_down.SetBinError(idx, new_down_e)
				h_clone_up.Write()
				h_clone_down.Write()
		del f
