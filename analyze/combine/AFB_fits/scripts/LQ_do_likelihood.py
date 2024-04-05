
from LQ_utils import *
import ROOT
from ROOT import *
import matplotlib.pyplot as plt
from array import array
plt.ioff()


parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--chan",  default="combined", type="string", help="What channels to run the fit over (combined, ee, or mumu)")
parser.add_option("--q",  default="combined", type="string", help="What channels to run the fit over (combined, u, or d)")
parser.add_option("--plot",  default=False, action="store_true", help="make plots")
parser.add_option("--no_sys",  default=False, action="store_true", help="Use fit template without any shape systematics")
parser.add_option("--fake_data",  default=True, action="store_true", help="Use fit template without any shape systematics and no fakes")
parser.add_option("--vec",  default=False, action="store_true", help="vec")
parser.add_option("-y", "--year", default = -1, type='int', help="Only do fits for this single year (2016,2017, or 2018), default is all years")
parser.add_option("--mLQ", default = 2500, type='int', help="LQ mass")
parser.add_option("--poi", default = 'yLQ2', type='string', help="poi to run the likelihood_scan")
parser.add_option("--ending", default = 'yLQ2', type='string', help="ending string")
parser.add_option("--statuncs", default = False,  help="freeze allConstrainedNuisances")
parser.add_option("--noSymMCStats", default = True, action="store_true",  help="Don't add constraints to mcStat nuisances")
parser.add_option("--gen_level",  default=False, action="store_true", help="gen level fits")
#parser.add_option("--fake_data",  default=True, action="store_true", help="Use fit template without any shape systematics and no fakes")
parser.add_option("--no_LQ",  default=False, action="store_true", help="For sanity check purposes remove LQ temps")
parser.add_option("-o", "--odir", default="likelihood_scans/", help = "output directory")
(options, args) = parser.parse_args()

def save_likelihoods(f,label):
    deltaNLL, poi_list = [],[]

    limit_tree = f.Get("limit")
    poi_value = array('f',[0])
    limit_tree.SetBranchAddress("%s"%poi, poi_value)

    for i in range(limit_tree.GetEntries()):

        limit_tree.GetEntry(i)
        deltaNLL.append(limit_tree.deltaNLL)
        poi_list.append(poi_value[0])
        print(poi_value)
    f.Close()
    #print(poi_list)
    #print(deltaNLL)
    idx = np.argsort(np.array(poi_list))
    poi_list = np.array(poi_list)[idx]
    deltaNLL = np.array(deltaNLL)[idx]
    with open('%s/like_scan_%s%s_%s%s_m%i_%s.txt'%(options.odir, label, options.chan,options.q,("_vec" if is_vec else ""), mLQ, poi), 'w') as f:
        for ylq,dnll in zip(poi_list, deltaNLL):
            f.write("%f %f\n" %(ylq,2*dnll))
        #print(np.amax(poi_list),np.amin(poi_list))

is_vec = options.vec
statuncs = options.statuncs
#options.gen_level = False
extra_params=""
mLQ = options.mLQ
poi = options.poi
ending = options.ending
   

#extra_params += " --cminApproxPreFitTolerance 1.0 --cminDefaultMinimizerTolerance 0.5 --cminDefaultMinimizerStrategy 0 "
if statuncs: extra_params += " --freezeParameters allConstrainedNuisances"
    
fit_name = options.chan
if(options.no_sys): 
    fit_name +="_nosys"


if(options.year > 0): fit_name +="_y%i" % (options.year % 2000)
fit_name+="_"+options.q

if(options.no_LQ): fit_name+="_noLQ"

if options.chan=="ee" and options.gen_level : fit_name+="_gen_level_SMdata_nlosys"

if is_vec: fit_name+="_vec"
if statuncs: fit_name += "_statuncs"

print("\n fit_name = ", fit_name)

if options.plot:
    poi_list=["yLQ2"]
    '''
    poi_list = ["MCStatBin1", "MCStatBin2", "MCStatBin3", "MCStatBin4", "MCStatBin9", "MCStatBin10",
     "MCStatBin11", "MCStatBin15", "MCStatBin16", "MCStatBin17", "MCStatBin21", "MCStatBin22", "MCStatBin23",
      "MCStatBin24", "MCStatBin29", "MCStatBin30", "MCStatBin31", "MCStatBin35", "MCStatBin36", "MCStatBin37", 
      "MCStatBin41", "MCStatBin42", "MCStatBin43", "MCStatBin44", "MCStatBin49", "MCStatBin50", "MCStatBin51", 
      "MCStatBin55", "MCStatBin56", "MCStatBin57",
      ]

    for i in range(1,61):
        poi_list.append("pdf" + str(i))


    for i in range(1,61):
        poi_list.append("prop_binY18_bin" + str(i))
    '''
    for poi in poi_list:

        #print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/ssekhar/Condor_outputs/likelihood_%s_%s%s_%s/like_scan_%s_%s%s_m%s_%s.txt %s"
        #            %(options.chan, options.q, ("_vec" if is_vec else ""),  poi, options.chan, options.q, ("_vec" if is_vec else ""), mLQ, poi, options.odir))
        label = ""
        respull = []
        with open('%s/like_scan_%s%s_%s%s_m%i_%s.txt'%(options.odir, label, options.chan, options.q, ("_vec" if is_vec else ""), mLQ, poi), 'r') as f:
            for line in f.readlines():
                respull.append(line.split(' '))

        respull = np.asarray(respull, dtype=float)
        poi_list = respull[:,0].tolist()
        deltaNLL = respull[:,1].tolist()
        label = "expected_"
        respull = []
        with open('%s/like_scan_%s%s_%s%s_m%i_%s.txt'%(options.odir, label, options.chan, options.q, ("_vec" if is_vec else ""), mLQ, poi), 'r') as f:
            for line in f.readlines():
                respull.append(line.split(' '))

        respull = np.asarray(respull, dtype=float)
        poi_list_exp = respull[:,0].tolist()
        deltaNLL_exp = respull[:,1].tolist()

        #plt.xlim(-1,1)
        plt.ylim(0,5)          
        plt.plot(poi_list,deltaNLL,label='data')
        plt.plot(poi_list_exp,deltaNLL_exp,label='asimov dataset')
        plt.plot(poi_list,len(poi_list)*[1],linestyle='dashed',c='g')
        plt.plot(poi_list,len(poi_list)*[2],linestyle='dashed',c='g')
        plt.xlabel("%s"%poi)
        plt.ylabel("-2deltaLL")
        plt.legend()
        plt.title("Likelihood Scan: channel %s %s, mLQ = %.2f TeV"%(options.chan,options.q,mLQ))
        plt.savefig("%s/like_scan_%s_%s_%s.jpg"%(options.odir,options.chan,options.q,poi))
        plt.close()

        
else:

    print(" \n \n Starting fit for LQ m = %i\n\n",mLQ)

    workspace="workspaces/%s_LQ.root" % (options.chan)
    make_workspace(workspace, options.gen_level, options.chan, options.q, is_vec, options.no_LQ, options.no_sys, options.fake_data, mLQ, year = options.year,noSymMCStats = True)
    
    label = "expected_"
    print_and_do("combine %s -M MultiDimFit  --algo grid --points 300 --squareDistPoiStep  --autoRange 2 -P %s --floatOtherPOIs 1   --saveWorkspace --saveFitResult --robustFit 1  %s -t -1" %(workspace, poi,  extra_params))
    print_and_do("cp higgsCombineTest.MultiDimFit.mH120.root higgsCombineTest.MultiDimFit.%s%s_%s_%s.root"%(label,poi,options.chan,options.q))
    f = ROOT.TFile.Open("higgsCombineTest.MultiDimFit.%s%s_%s_%s.root"%(label,poi,options.chan,options.q),"READ")
    save_likelihoods(f,label)
    
    label = ""
    print_and_do("combine %s -M MultiDimFit  --algo grid --points 300 --squareDistPoiStep  --autoRange 2 -P %s --floatOtherPOIs 1   --saveWorkspace --saveFitResult --robustFit 1  %s " %(workspace, poi,  extra_params))

    print_and_do("cp higgsCombineTest.MultiDimFit.mH120.root higgsCombineTest.MultiDimFit.%s_%s_%s.root"%(poi,options.chan,options.q))
    f = ROOT.TFile.Open("higgsCombineTest.MultiDimFit.%s%s_%s_%s.root"%(label,poi,options.chan,options.q),"READ")
    save_likelihoods(f,label)
     
	    

    
