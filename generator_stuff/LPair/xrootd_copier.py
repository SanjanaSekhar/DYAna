import os

def print_and_do(s):
    print("%s" %s)
    os.system(s)

files = ["GamGamToEE_ElEl_M-120to200_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",           "GamGamToMuMu_ElEl_M-120to200_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe",
                                                                                                    "GamGamToMuMu_ElEl_M-120to200_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_second.lhe",
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
files = ['GamGamToMuMu_ElEl_M-120to200_Pt-15toInf_CepGen-lpair_13TeV_v3_tosplit_first.lhe']


for f in files:
    print_and_do("xrdcp root://cmsxrootd.fnal.gov//store/user/bbilin/lhe/%s ." % f)
    print_and_do("xrdcp %s root://cmseos.fnal.gov//store/user/oamram/LPair/lhe/raw/%s" % (f, f))
    print_and_do("rm %s" %f)

    
