import os,sys

def print_and_do(s):
    print("%s" %s)
    os.system(s)


files = ["GamGamToEE_ElEl_M-120to200_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",                    "GamGamToMuMu_ElEl_M-120to200_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_ElEl_M-120to200_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",                   "GamGamToMuMu_ElEl_M-1400to2300_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_ElEl_M-1400to2300_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",                  "GamGamToMuMu_ElEl_M-1400to2300_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_ElEl_M-1400to2300_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",                 "GamGamToMuMu_ElEl_M-200to400_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_ElEl_M-200to400_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",                    "GamGamToMuMu_ElEl_M-200to400_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_ElEl_M-200to400_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",                   "GamGamToMuMu_ElEl_M-2300toInf_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_ElEl_M-2300toInf_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",                   "GamGamToMuMu_ElEl_M-2300toInf_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_ElEl_M-2300toInf_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",                  "GamGamToMuMu_ElEl_M-400to800_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_ElEl_M-400to800_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",                    "GamGamToMuMu_ElEl_M-400to800_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_ElEl_M-400to800_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",                   "GamGamToMuMu_ElEl_M-800to1400_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_ElEl_M-800to1400_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",                   "GamGamToMuMu_ElEl_M-800to1400_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_ElEl_M-800to1400_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",                  "GamGamToMuMu_InelEl_SY_M-120to200_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_InelEl_SY_M-120to200_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",       "GamGamToMuMu_InelEl_SY_M-120to200_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_InelEl_SY_M-120to200_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",      "GamGamToMuMu_InelEl_SY_M-1400to2300_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_InelEl_SY_M-1400to2300_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",     "GamGamToMuMu_InelEl_SY_M-1400to2300_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_InelEl_SY_M-1400to2300_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",    "GamGamToMuMu_InelEl_SY_M-200to400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_InelEl_SY_M-200to400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",       "GamGamToMuMu_InelEl_SY_M-200to400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_InelEl_SY_M-200to400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",      "GamGamToMuMu_InelEl_SY_M-2300toInf_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_InelEl_SY_M-2300toInf_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",      "GamGamToMuMu_InelEl_SY_M-2300toInf_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_InelEl_SY_M-2300toInf_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",     "GamGamToMuMu_InelEl_SY_M-400to800_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_InelEl_SY_M-400to800_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",       "GamGamToMuMu_InelEl_SY_M-400to800_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_InelEl_SY_M-400to800_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",      "GamGamToMuMu_InelEl_SY_M-800to1400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_InelEl_SY_M-800to1400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",      "GamGamToMuMu_InelEl_SY_M-800to1400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_InelEl_SY_M-800to1400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",     "GamGamToMuMu_InelInel_SY_M-120to200_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_InelInel_SY_M-120to200_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",     "GamGamToMuMu_InelInel_SY_M-120to200_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_InelInel_SY_M-120to200_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",    "GamGamToMuMu_InelInel_SY_M-1400to2300_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_InelInel_SY_M-1400to2300_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",   "GamGamToMuMu_InelInel_SY_M-1400to2300_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_InelInel_SY_M-1400to2300_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",  "GamGamToMuMu_InelInel_SY_M-200to400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_InelInel_SY_M-200to400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",     "GamGamToMuMu_InelInel_SY_M-200to400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_InelInel_SY_M-200to400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",    "GamGamToMuMu_InelInel_SY_M-2300toInf_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_InelInel_SY_M-2300toInf_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",    "GamGamToMuMu_InelInel_SY_M-2300toInf_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_InelInel_SY_M-2300toInf_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",   "GamGamToMuMu_InelInel_SY_M-400to800_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_InelInel_SY_M-400to800_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",     "GamGamToMuMu_InelInel_SY_M-400to800_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_InelInel_SY_M-400to800_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",    "GamGamToMuMu_InelInel_SY_M-800to1400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",
"GamGamToEE_InelInel_SY_M-800to1400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_first.lhe",    "GamGamToMuMu_InelInel_SY_M-800to1400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",
"GamGamToEE_InelInel_SY_M-800to1400_Pt-15toInf_CepGen-lpair-pythia6_13TeV_v3_tosplit_second.lhe",  
]

n_events = 25000
evts_per_file = 500
nFiles = n_events / evts_per_file

for fname in files:
    if('2300' not in fname): continue
    print_and_do("xrdcp root://cmseos.fnal.gov//store/user/oamram/LPair/lhe/raw/%s ." % (fname))
    f  = open(fname, "r")
    header = []
    line = f.readline()
    while("<event>" not in line):
        header.append(line)
        line = f.readline()
    print("got buffer")

    for iFile in range(nFiles):
        fnew_name = fname[:-4] + ("_%i" % iFile) + ".lhe"
        print("making %s" % fnew_name)
        fout = open(fnew_name, "w")
        evt_counter = 0
        for l in header:
            fout.write(l)
        if(iFile ==0): fout.write(line)
        while(evt_counter < evts_per_file):
            line = f.readline()
            fout.write(line)
            if("</event>" in line): 
                evt_counter +=1
        fout.write("</LesHouchesEvents>")
        fout.close()

        print_and_do("xrdcp -f %s root://cmseos.fnal.gov//store/user/oamram/LPair/lhe/split/%s" % (fnew_name, fnew_name))
        print_and_do("rm %s" % fnew_name)

    f.close()

    print_and_do("rm %s" %fname)
