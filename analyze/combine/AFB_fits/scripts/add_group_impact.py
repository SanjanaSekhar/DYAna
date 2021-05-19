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
    up = triplet[2][1] - triplet[0][1]
    dn = triplet[1][1] - triplet[0][1]
    # print dn, up
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

def compute_sys(sys_l, nom_l, seed = -1):
    if(seed <0):
        sys_f = "higgsCombine_%s.FitDiagnostics.mH120.root" % sys_l
        nom_f = "higgsCombine_%s.FitDiagnostics.mH120.root" % nom_l
    else:
        sys_f = "higgsCombine_%s.FitDiagnostics.mH120.%i.root" % (sys_l,seed)
        nom_f = "higgsCombine_%s.FitDiagnostics.mH120.%i.root" % (nom_l,seed)


    nom_unc = getTripletFromFile(nom_f)
    new_unc = getTripletFromFile(sys_f)
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


