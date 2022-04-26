import ROOT
from optparse import OptionParser
import sys

ext = "042422_ud_nometcut"
mLQ_list = [1000,1500,2000,2500,3000,3500,4000]
#if(len(sys.argv) < 3):
#    print("Need filename and year")
#    sys.exit(1)
for mLQ in mLQ_list:
    for year in [2016,2017,2018]:

        fin = "combine/templates/LQm%i_merge_templates%i_%s.root"%(mLQ,year-2000,ext)
    #    fin = sys.argv[1]
    #    year = int(sys.argv[2])
        print("year is %i" % year)

        f = ROOT.TFile.Open(fin, "UPDATE")
        #n_bins = 8

        """
        {"_METJEC", "_METHEM", "_prefire", "_elScaleSyst", "_elScaleStat","_elScaleGain", "_elSmear", "_muRC", "_Pu", "_BTAGCOR", "_BTAGUNCOR", "_BTAGLIGHT",
        "_muHLTBAR", "_muIDBAR", "_muISOBAR",  "_muHLTEND", "_muIDEND", "_muISOEND",  "_muIDSYS", "_muISOSYS",  
        "_elHLTBARPTHIGH", "_elIDBARPTHIGH", "_elRECOBARPTHIGH", "_elHLTENDPTHIGH", "_elIDENDPTHIGH", "_elRECOENDPTHIGH",
        "_elHLTBARPTLOW", "_elIDBARPTLOW", "_elRECOBARPTLOW", "_elHLTENDPTLOW", "_elIDENDPTLOW", "_elRECOENDPTLOW",
        "_ptrw1b", "_ptrw2b", "_ptrw3b", "_ptrw4b", "_ptrw5b", "_ptrw6b", "_ptrw7b",
        "_emucostrw1b", "_emucostrw2b", "_emucostrw3b", "_emucostrw4b",
        "_elfakesrw1b", "_elfakesrw2b", "_elfakesrw3b", "_elfakesrw4b",
        "_mufakesrw1b", "_mufakesrw2b", "_mufakesrw3b", "_mufakesrw4b",
        "_RENORM", "_FAC", "_REFAC", "_A0Den", "_alphaS",
        """
        correlate_all = ["elScaleSyst", "elSmear", "Pu", "muIDSYS", "muISOSYS", "elRECOBARPTHIGH", "elRECOENDPTHIGH", "elRECOBARPTLOW", "elRECOENDPTLOW",
                         "elIDBARPTHIGH", "elIDENDPTHIGH", "elIDBARPTLOW", "elIDENDPTLOW", "BTAGCOR"] 

        correlate_1718 = ["ptrw1b", "ptrw2b", "ptrw3b", "ptrw4b", "ptrw5b", "ptrw6b", "ptrw7b", 
                            "emucostrw1b", "emucostrw2b", "emucostrw3b", "emucostrw4b",
                            "RENORM", "FAC", "REFAC", "alphaS" ]

        #for i in range(n_bins):
        f.cd("LQ")
        keys = ROOT.gDirectory.GetListOfKeys().Clone()


        
        for k in keys:
            #print(h.GetName())
            name = k.GetName()
            for sys_a in correlate_all:
                if(sys_a in name):
                    h = ROOT.gDirectory.Get(name)
                    new_name = name.replace(sys_a + str(year - 2000),  sys_a)
                    if(new_name in keys): break
                    print("copying %s to %s" % (name, new_name))
                    h_clone = h.Clone(new_name)
                    h_clone.Write()
                    break

            if(year > 2016):
                for sys_1718 in correlate_1718:
                    if(sys_1718 in name):
                        h = ROOT.gDirectory.Get(name)
                        new_name = name.replace(sys_1718 + str(year - 2000),  sys_1718+"1718")
                        if(new_name in keys): break
                        print("copying %s to %s" % (name, new_name))
                        h_clone = h.Clone(new_name)
                        h_clone.Write()
                        break
        keys = ROOT.gDirectory.GetListOfKeys().Clone()
        for k in keys:

            name = k.GetName()
            if ("LQint_d" in name):
                print("Flipping ",name)
                h = ROOT.gDirectory.Get(name)
                #h_clone = h.Clone(name)
                h.Scale(-1)
                h.Write("",ROOT.TObject.kOverwrite)
                
                    

