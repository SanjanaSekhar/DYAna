import ROOT
import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
#from AFB_fits.scripts.LQ_utils import *

ext = "020923"

for year in [2016]:
	fin = "combine/templates/LQ_data_templates%i.root"%(year-2000)
	f = ROOT.TFile.Open(fin, "UPDATE")
	f.cd("LQ")
	h_ee_data = ROOT.gDirectory.Get("ee%i_data_obs"%(year-2000))
	h_mumu_data = ROOT.gDirectory.Get("mumu%i_data_obs"%(year-2000))
	
	for mLQ in range(1000,5500,500):

		fin = "combine/templates/LQ_data_templates%i.root"%(year-2000)
		f = ROOT.TFile.Open(fin, "UPDATE")
		f.cd("LQ")
		h_ee_data = ROOT.gDirectory.Get("ee%i_data_obs"%(year-2000))
		h_mumu_data = ROOT.gDirectory.Get("mumu%i_data_obs"%(year-2000))

		#f.Close()
		#del f

		fin2 = "combine/templates/LQm%i_merge_templates%i_%s.root"%(mLQ,year-2000,ext)
		f2 = ROOT.TFile.Open(fin2, "UPDATE")
		f2.cd("LQ")
		h_ee_fakedata = ROOT.gDirectory.Get("ee%i_data_obs"%(year-2000))
		h_mumu_fakedata = ROOT.gDirectory.Get("mumu%i_data_obs"%(year-2000))

		h_ee_fakedata2 = h_ee_fakedata.Clone("ee%i_fakedata_obs"%(year-2000))
		h_mumu_fakedata2 = h_mumu_fakedata.Clone("mumu%i_fakedata_obs"%(year-2000))

		for idx in range(1,h_ee_data.GetNbinsX()+1):

			h_ee_fakedata.SetBinContent(idx, h_ee_data.GetBinContent(idx))
			h_ee_fakedata.SetBinError(idx, h_ee_data.GetBinError(idx))

			h_mumu_fakedata.SetBinContent(idx, h_mumu_data.GetBinContent(idx))
			h_mumu_fakedata.SetBinError(idx, h_mumu_data.GetBinError(idx))

		h_ee_fakedata2.Write()
		h_mumu_fakedata2.Write()
		h_ee_fakedata.Write("",ROOT.TObject.kOverwrite)
		h_mumu_fakedata.Write("",ROOT.TObject.kOverwrite)
		
