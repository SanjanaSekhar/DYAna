#!/usr/bin/env python
import sys
sys.argv.append( '-b-' )
import ROOT as r
import json
r.gROOT.SetBatch(True)
sys.argv.remove( '-b-' )

from ROOT import TFile
from pprint import pprint

def getTripletFromFile(fName):
	f = TFile.Open(fName)
	# f.Print()
	t = f.Get('limit')
	# t.Scan('*')
	vals = []
	for e in t:
		if e.quantileExpected == -1:
			continue
		vals.append( (e.quantileExpected, e.limit) )

	# print vals
	return UncertaintiesFromTriplet(vals)

def UncertaintiesFromTriplet(triplet):
	#This assumes that the values are ordered as central(0),lower(1),upper(2)
	#This could be found from the quantileExpected value, but I was lazy
	#up = triplet[2][1] - triplet[0][1]
	#dn = triplet[1][1] - triplet[0][1]
	# print dn, up
	print(triplet[0],triplet[1],triplet[2])
	up = triplet[2] - triplet[0]
	dn = triplet[1] - triplet[0]
	print("(dn,up) = ",(dn,up))
	return (dn, up)

from math import sqrt
def quadDiffPair(unc1, unc2):
	dn = sqrt( abs(unc1[0]**2 - unc2[0]**2 ))
	up = sqrt( abs(unc1[1]**2 - unc2[1]**2) )
	return (-dn, up)

def quadDiffAvg(unc1, unc2):
	dn = sqrt( abs(unc1[0]**2 - unc2[0]**2 ))
	up = sqrt( abs(unc1[1]**2 - unc2[1]**2) )
	return (dn + up)/2.0

def compute_sys(nom_l, sys_l, seed = -1):
	if(seed <0):
		sys_f = "higgsCombine_%s.FitDiagnostics.mH120.root" % sys_l
		nom_f = "higgsCombine_%s.FitDiagnostics.mH120.root" % nom_l
	else:
		#sys_f = "higgsCombine_%s.FitDiagnostics.mH120.%i.root" % (sys_l,seed)
		#nom_f = "higgsCombine_%s.FitDiagnostics.mH120.%i.root" % (nom_l,seed)
	#sys_f = "higgsCombine_%s.MultiDimFit.mH120.%i.root" % (sys_l,seed)
		#nom_f = "higgsCombine_%s.MultiDimFit.mH120.%i.root" % (nom_l,seed)
		sys_f = "multidimfit_%s.root" % (sys_l)
		nom_f = "multidimfit_%s.root" % (nom_l)
	
	_file0 = TFile.Open(nom_f)
	fit_mdf = _file0.Get("fit_mdf")
	nom = fit_mdf.floatParsFinal()
	nom_index_yLQ2 = nom.index("yLQ2")
	nom_vals = nom.at(nom_index_yLQ2)
	nom_unc = []
	nom_unc.append(nom_vals.getValV())
	nom_unc.append(nom_vals.getErrorLo())
	nom_unc.append(nom_vals.getErrorHi())
	nom_unc = UncertaintiesFromTriplet(nom_unc)
	

	_file1 = TFile.Open(sys_f)
	fit_mdf = _file1.Get("fit_mdf")
	sys = fit_mdf.floatParsFinal()
	sys_index_yLQ2 = sys.index("yLQ2")
	sys_vals = sys.at(sys_index_yLQ2)
	sys_unc = []
	sys_unc.append(sys_vals.getValV())
	sys_unc.append(sys_vals.getErrorLo())
	sys_unc.append(sys_vals.getErrorHi())
	new_unc = UncertaintiesFromTriplet(sys_unc)
	#nom_unc = getTripletFromFile(nom_f)
	#new_unc = getTripletFromFile(sys_f)
	sys_unc = quadDiffAvg(nom_unc, new_unc)
	return sys_unc

##print diff
#fileIn = sys.argv[4]
#
#json_file = open(fileIn, 'r')
#data = json.load(json_file)
#nom_afb = data['POIs'][0]['fit'][1]
#my_param = {}
#my_param['Afb'] = [nom_afb + diff[0], nom_afb, nom_afb + diff[1]]
#my_param['fit'] = [1.0,1.0,1.0]
#my_param['groups'] = []
#my_param['impact_Afb'] = max(abs(diff[0]), abs(diff[1]))
#my_param['name'] = sys.argv[3]
#my_param['prefit'] = [1.0,1.0,1.0]
#my_param["type"]= "Unconstrained"
#data['params'].append(my_param)
#json_file.close()
#out_file = open(fileIn, 'w')
#json.dump(data, out_file, separators=(', \n', ': \n'))


