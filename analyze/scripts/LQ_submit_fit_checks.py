import numpy as np
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

def print_and_do(s):
    print(s)
    return os.system(s)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--ending",  default="021223", help="date")
parser.add_option("--bias_tests",  default=False, help="do bias tests")
parser.add_option("--limits",  default=False, help="do limits")
parser.add_option("--combine_review",  default=False, help="do review")
parser.add_option("--impacts",  default=False, help="do impacts")
parser.add_option("--sys_uncs",  default=False, help="do sys uncs")
parser.add_option("--likelihood",  default=False, help="do like scans")
parser.add_option("--gof",  default=False, help="do gof")
(options, args) = parser.parse_args()

n_m_bins = 1
date = options.ending


if options.bias_tests:
    #mLQ_list = [5000,9000]
    mLQ_list = [5500,6500,7000]
    for mLQ in mLQ_list:
        
	cmds = [
        #"python scripts/LQ_do_bias_test.py  --yLQ 2.5 --chan ee --q d --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        #"python scripts/LQ_do_bias_test.py  --yLQ 3.0 --chan ee --q d --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        #"python scripts/LQ_do_bias_test.py  --yLQ 0.35 --chan ee --q d --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
	#"python scripts/LQ_do_bias_test.py  --yLQ 3.5 --chan ee --q d --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        #"python scripts/LQ_do_bias_test.py  --yLQ 2.5 --chan ee --q u --is_vec True --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        #"python scripts/LQ_do_bias_test.py  --yLQ 3.0 --chan ee --q u --is_vec True --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        #"python scripts/LQ_do_bias_test.py  --yLQ 0.35 --chan ee --q u --is_vec True --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
	#"python scripts/LQ_do_bias_test.py  --yLQ 3.5 --chan ee --q u --is_vec True --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        #"python scripts/LQ_do_bias_test.py  --yLQ 2.5 --chan mumu --q u --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        #"python scripts/LQ_do_bias_test.py  --yLQ 3.0 --chan mumu --q u --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        #"python scripts/LQ_do_bias_test.py  --yLQ 0.35 --chan mumu --q u --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
	#"python scripts/LQ_do_bias_test.py  --yLQ 3.5 --chan mumu --q u --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        #"python scripts/LQ_do_bias_test.py  --yLQ 2.5 --chan mumu --q d --is_vec True --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 3.0 --chan mumu --q d --is_vec True --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        #"python scripts/LQ_do_bias_test.py  --yLQ 0.35 --chan mumu --q d --is_vec True --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
	"python scripts/LQ_do_bias_test.py  --yLQ 3.5 --chan mumu --q d --is_vec True --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        
        ]
	
        labels = [
                #"bias_test_yLQ2.5_ee_d_m",  "bias_test_yLQ3.0_ee_d_m",  "bias_test_yLQ3.5_ee_d_m",  #"bias_test_yLQ0.6_ee_d_m", 
		
                #"bias_test_yLQ2.5_ee_u_vec_m", "bias_test_yLQ3.0_ee_u_vec_m","bias_test_yLQ3.5_ee_u_vec_m", #"bias_test_yLQ0.6_ee_u_vec_m",
        	#"bias_test_yLQ2.5_mumu_u_m", "bias_test_yLQ3.0_mumu_u_m", "bias_test_yLQ3.5_mumu_u_m",#"bias_test_yLQ0.6_mumu_u_m",  
		#"bias_test_yLQ2.5_mumu_d_vec_m", 
		"bias_test_yLQ3.0_mumu_d_vec_m", "bias_test_yLQ3.5_mumu_d_vec_m",# "bias_test_yLQ0.6_mumu_d_vec_m"
                
        ]
	'''
	cmds = [
        "python scripts/LQ_do_bias_test.py  --yLQ 0.25 --chan mumu --q u  --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
	"python scripts/LQ_do_bias_test.py  --yLQ 0.5 --chan ee --q d  --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 0.25 --chan ee --q d  --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
	"python scripts/LQ_do_bias_test.py  --yLQ 0.5 --chan mumu --q d --is_vec True --nToys 10 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        ]
	labels = [ "bias_test_yLQ0.25_mumu_u_m",
	"bias_test_yLQ0.5_ee_d_m","bias_test_yLQ0.25_ee_d_m",
	"bias_test_yLQ0.5_mumu_d_vec_m"
	]
	'''
        #sys = ["MCStatBin16,MCStatBin17,MCStatBin18,autoMCStats","elScales16,elScales17,elScales18","RFscales16,RFscales1718","elHLTs16,elHLTs17,elHLTs18"]
        #sys = ["xsecs","MCStatBin16,MCStatBin17,MCStatBin18,autoMCStats","RFscales16,RFscales1718","muIDs16,muIDs17,muIDs18","muHLTs16,muHLTs17,muHLTs18","pdfs"]
        sys = [""]
        cpy_cmd = "xrdcp -f temp/* $1 \n"

        total_jobs = 500

        for i,cmd in enumerate(cmds):

            for job_idx in range(total_jobs/10):
                for s in sys:
                    #regular templates
                    script_name = "scripts/script3.sh"
                    print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
                    script_file = open(script_name, 'a+')
                    script_file.write("mkdir temp\n")
                    script_file.write(cmd + " --job %i \n"%(job_idx))
                    script_file.write(cpy_cmd)
                    script_file.close()
                    #print_and_do("cat %s" % script_name)
                    print_and_do("chmod +x %s" % script_name)
                    print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s%i_no%s_%i_%s"  % (n_m_bins, script_name, labels[i], mLQ, s.replace(",",""), job_idx, date))
                    print_and_do("rm scripts/script3.sh")

if options.limits:
    cmds = [
  
    #"python scripts/LQ_get_limits.py --chan ee --q u  -o limits/ --ending %s  "%date,
    #"python scripts/LQ_get_limits.py --chan ee --q d  -o limits/ --ending %s  "%date,
    #"python scripts/LQ_get_limits.py --chan mumu --q u -o limits/ --ending %s "%date,
    #"python scripts/LQ_get_limits.py --chan mumu --q d -o limits/ --ending %s "%date,
    #"python scripts/LQ_get_limits.py --chan ee --q u --vec True -o limits/ --ending %s "%date,
    #"python scripts/LQ_get_limits.py --chan ee --q d --vec True -o limits/ --ending %s  "%date,
    "python scripts/LQ_get_limits.py --chan mumu --q u --vec True -o limits/ --ending %s "%date,
    #"python scripts/LQ_get_limits.py --chan mumu --q d --vec True -o limits/ --ending %s  "%date,
 
    ]

    labels = [
        #"limits_ee_u","limits_ee_d","limits_mumu_u","limits_mumu_d",
        #"limits_ee_u_vec","limits_ee_d_vec","limits_mumu_u_vec",
        #"limits_mumu_d_vec"
        #"limits_ee_s","limits_mumu_s",
        #"limits_ee_d_vec","limits_mumu_d_vec"
	"limits_mumu_u_vec_HybridNew"
    ]
    
    year = -1
    cpy_cmd = "xrdcp -f limits_HybridNew/* $1 \n"

    for i,cmd in enumerate(cmds):
	   for m in [2500]:
	   	for point in np.arange(0.01,0.05,0.005):
	    		for ite in range(0,10):
            			#regular templates
            			script_name = "scripts/script3.sh"
            			print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
            			script_file = open(script_name, 'a+')
            			script_file.write("mkdir limits_HybridNew\n")
            			script_file.write(cmd+" --HybridNew true --mLQ %i --year %i --inject_yLQ2 %.5f --ntoys 20 --iterations 1 \n"%(m,year,point))
            			script_file.write(cpy_cmd)
            			script_file.close()
            			#print_and_do("cat %s" % script_name)
            			print_and_do("chmod +x %s" % script_name)
            			#print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s_m%i_y%i_%s"  % (n_m_bins, script_name, labels[i], m, year-2000, date))
            			
				print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename -s %s -n %s_m%i_p%.5f_%i"% (n_m_bins, script_name, labels[i], m, point, ite))
				print_and_do("rm scripts/script3.sh")

if options.gof:

    nToys = 300
    year = -1    
    for mLQ in [2500]:


        cmds = [
      
        #"python scripts/LQ_do_gof.py --chan ee --q u  ",
        #"python scripts/LQ_do_gof.py --chan ee --q d  ",
        #"python scripts/LQ_do_gof.py --chan mumu --q u ",
        #"python scripts/LQ_do_gof.py --chan mumu --q d ",
        #"python scripts/LQ_do_gof.py --chan ee --q u --vec True ",
        #"python scripts/LQ_do_gof.py --chan ee --q d --vec True ",
        "python scripts/LQ_do_gof.py --chan mumu --q u --vec True ",
        #"python scripts/LQ_do_gof.py --chan mumu --q d --vec True ",
     
        ]

        labels = [
            #"gof_ee_u","gof_ee_d",#"gof_mumu_u","gof_mumu_d",
            #"gof_ee_u_vec","gof_ee_d_vec",
	    "gof_mumu_u_vec"#,"gof_mumu_d_vec"
            #"gof_ee_s","gof_mumu_s",
            #"gof_ee_d_vec","gof_mumu_d_vec"
        ]


        cpy_cmd = "xrdcp -f gof/* $1 \n"

        for i,cmd in enumerate(cmds):
       
        #for point in np.arange(0.28,1.5,0.005):
        #for q in [0.025,0.16,0.5,0.84,0.975]:
            #regular templates
            script_name = "scripts/script3.sh"
            print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
            script_file = open(script_name, 'a+')
            script_file.write("mkdir gof_b_only\n")
            script_file.write(cmd+" -o gof/ --mLQ %i --nToys %i --year %i\n"%(mLQ, nToys, year))
            script_file.write(cpy_cmd)
            script_file.close()
            #print_and_do("cat %s" % script_name)
            print_and_do("chmod +x %s" % script_name)
            print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s_m%i_y%i"  % (n_m_bins, script_name, labels[i], mLQ, year-2000))
            print_and_do("rm scripts/script3.sh")

if options.combine_review:

    cmds = [
  
    "python scripts/for_combine_review.py --chan ee --q u -o CR\n",
    "python scripts/for_combine_review.py --chan ee --q d -o CR\n",
    "python scripts/for_combine_review.py --chan mumu --q u -o CR\n",
    "python scripts/for_combine_review.py --chan mumu --q d -o CR\n",
    "python scripts/for_combine_review.py --chan ee --q u --vec True -o CR\n",
    "python scripts/for_combine_review.py --chan ee --q d --vec True -o CR\n",
    "python scripts/for_combine_review.py --chan mumu --q u --vec True -o CR\n",
    "python scripts/for_combine_review.py --chan mumu --q d --vec True -o CR\n",
 
    ]

    labels = [
        "CR_ee_u","CR_ee_d","CR_mumu_u","CR_mumu_d",
        "CR_ee_u_vec","CR_ee_d_vec","CR_mumu_u_vec","CR_mumu_d_vec"
        #"limits_ee_s","limits_mumu_s",
        #"limits_ee_s_vec","limits_mumu_s_vec"
    ]


    cpy_cmd = "xrdcp -f CR/* $1 \n"

    for i,cmd in enumerate(cmds):
    #for m in range(1000,9500,500):
    #for point in np.arange(0.28,1.5,0.005):
    #for q in [0.025,0.16,0.5,0.84,0.975]:
        #regular templates
        script_name = "scripts/script3.sh"
        print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
        script_file = open(script_name, 'a+')
        script_file.write("mkdir CR\n")
        script_file.write(cmd)
        script_file.write(cpy_cmd)
        script_file.close()
        #print_and_do("cat %s" % script_name)
        print_and_do("chmod +x %s" % script_name)
        print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s"  % (n_m_bins, script_name, labels[i]))
        print_and_do("rm scripts/script3.sh")


if options.impacts:

    cmds = [
  
    "python scripts/LQ_do_impacts.py --mLQ 2500 --chan ee --q u -o imps --ending %s \n"%date,
    "python scripts/LQ_do_impacts.py --mLQ 2500 --chan ee --q d -o imps --ending %s  \n"%date,
    "python scripts/LQ_do_impacts.py --mLQ 2500 --chan mumu --q u -o imps --ending %s  \n"%date,
    "python scripts/LQ_do_impacts.py --mLQ 2500 --chan mumu --q d -o imps --ending %s \n "%date,
    "python scripts/LQ_do_impacts.py --mLQ 2500 --chan ee --q u --vec True -o imps --ending %s \n "%date,
    "python scripts/LQ_do_impacts.py --mLQ 2500 --chan ee --q d --vec True -o imps --ending %s  \n"%date,
    "python scripts/LQ_do_impacts.py --mLQ 2500 --chan mumu --q u --vec True -o imps --ending %s  \n"%date,
    "python scripts/LQ_do_impacts.py --mLQ 2500 --chan mumu --q d --vec True -o imps --ending %s  \n"%date,
 
    ]

    labels = [
        "imps_ee_u","imps_ee_d","imps_mumu_u","imps_mumu_d",
        "imps_ee_u_vec","imps_ee_d_vec","imps_mumu_u_vec", "imps_mumu_d_vec"
        #"limits_ee_s","limits_mumu_s",
        #"limits_ee_s_vec","limits_mumu_s_vec"
    ]


    cpy_cmd = "xrdcp -f imps/* $1 \n"

    for i,cmd in enumerate(cmds):
    #for m in range(1000,9500,500):
    #for point in np.arange(0.28,1.5,0.005):
    #for q in [0.025,0.16,0.5,0.84,0.975]:
        #regular templates
        script_name = "scripts/script3.sh"
        print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
        script_file = open(script_name, 'a+')
        script_file.write("mkdir imps\n")
        script_file.write(cmd)
        script_file.write(cpy_cmd)
        script_file.close()
        #print_and_do("cat %s" % script_name)
        print_and_do("chmod +x %s" % script_name)
        print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s"  % (n_m_bins, script_name, labels[i]))
        print_and_do("rm scripts/script3.sh")

if options.sys_uncs:

    cmds = [
  
    "python scripts/LQ_check_sys_uncs.py  --chan ee --q u -o sys --ending %s "%date,
    #"python scripts/LQ_check_sys_uncs.py  --chan ee --q d -o sys --ending %s  "%date,
    "python scripts/LQ_check_sys_uncs.py  --chan mumu --q u -o sys --ending %s  "%date,
    #"python scripts/LQ_check_sys_uncs.py  --chan mumu --q d -o sys --ending %s  "%date,
    "python scripts/LQ_check_sys_uncs.py  --chan ee --q u --vec True -o sys --ending %s  "%date,
    #"python scripts/LQ_check_sys_uncs.py  --chan ee --q d --vec True -o sys --ending %s  "%date,
    "python scripts/LQ_check_sys_uncs.py  --chan mumu --q u --vec True -o sys --ending %s  "%date,
    #"python scripts/LQ_check_sys_uncs.py  --chan mumu --q d --vec True -o sys --ending %s  "%date,
 
    ]

    labels = [
        "sys_ee_u","sys_mumu_u",
        "sys_ee_u_vec","sys_mumu_u_vec"
        #"limits_ee_s","limits_mumu_s",
        #"limits_ee_s_vec","limits_mumu_s_vec"
    ]

    #ntoys = 10
    seed = 3456
    cpy_cmd = "xrdcp -f sys/* $1 \n"
    for mass in [2500]:
	for toy in range(1):
            for i,cmd in enumerate(cmds):
            #for m in range(1000,9500,500):
            #for point in np.arange(0.28,1.5,0.005):
            #for q in [0.025,0.16,0.5,0.84,0.975]:
            #regular templates
                script_name = "scripts/script3.sh"
                print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
                script_file = open(script_name, 'a+')
                script_file.write("mkdir sys\n")
                script_file.write(cmd+"  --mLQ %i \n"%(mass))
                script_file.write(cpy_cmd)
                script_file.close()
                #print_and_do("cat %s" % script_name)
                print_and_do("chmod +x %s" % script_name)
                print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s_%i"  % (n_m_bins, script_name, labels[i], mass))
                print_and_do("rm scripts/script3.sh")
                

if options.likelihood:

    cmds = [
  
    
    "python scripts/LQ_do_likelihood.py --mLQ 2500 --chan ee --q u  -o likelihood_scans ",
    "python scripts/LQ_do_likelihood.py --mLQ 2500 --chan ee --q d  -o likelihood_scans ",
    "python scripts/LQ_do_likelihood.py --mLQ 2500 --chan ee --q u --vec true -o likelihood_scans ",
    "python scripts/LQ_do_likelihood.py --mLQ 2500 --chan ee --q d --vec true -o likelihood_scans ",
    "python scripts/LQ_do_likelihood.py --mLQ 2500 --chan mumu --q u  -o likelihood_scans ",
    "python scripts/LQ_do_likelihood.py --mLQ 2500 --chan mumu --q d  -o likelihood_scans ",
    "python scripts/LQ_do_likelihood.py --mLQ 2500 --chan mumu --q u --vec true -o likelihood_scans ",
    "python scripts/LQ_do_likelihood.py --mLQ 2500 --chan mumu --q d --vec true -o likelihood_scans ",
 
    ]

    labels = [
     "likelihood_ee_u","likelihood_ee_d",
     "likelihood_ee_u_vec","likelihood_ee_d_vec",
     "likelihood_mumu_u","likelihood_mumu_d","likelihood_mumu_u_vec","likelihood_mumu_d_vec"
    ]


    cpy_cmd = "xrdcp -f likelihood_scans/* $1 \n"
    
    #poi = "RFscales16,RFscales1718"
 
    #["muISOBAR17", "muISOEND17","muIDEND17","muIDBAR17", "mufakesrw1b17", "mufakesrw2b17", "mufakesrw3b17", "mufakesrw4b17", "muISOBAR16", "muISOEND16","muIDEND16","muIDBAR16", "mufakesrw1b18", "mufakesrw2b18", "mufakesrw3b18", "mufakesrw4b18"]
    #poi_list = ["elScaleStat16", "elSmear", "elScaleGain16", "elScaleSyst","elScaleStat17","elScaleGain17","elScaleGain18", "elScaleStat17"]
    #poi_list = ["elIDENDPT", "elIDBARPT"]
    poi_list=["RENORM16", "alphaS16", "REFAC16", "FAC16","REFAC1718", "RENORM1718", "FAC1718", "alphaS1718"]
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
        for year in [-1]:
           for i,cmd in enumerate(cmds):
               #for m in range(1000,9500,500):
               #for point in np.arange(0.28,1.5,0.005):
               #for q in [0.025,0.16,0.5,0.84,0.975]:
               #regular templates
               script_name = "scripts/script3.sh"
               print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
               script_file = open(script_name, 'a+')
               script_file.write("mkdir likelihood_scans\n")
               script_file.write(cmd + " --poi %s --year %i\n "%(poi,year))
               script_file.write(cpy_cmd)
               script_file.close()
               #print_and_do("cat %s" % script_name)
               print_and_do("chmod +x %s" % script_name)
               print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s_%s_%s"  % (n_m_bins, script_name, labels[i], year,poi))
               print_and_do("rm scripts/script3.sh")
