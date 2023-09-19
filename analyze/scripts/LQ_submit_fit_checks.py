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
  
    "python scripts/LQ_get_limits.py --chan ee --q u  -o limits/ --ending %s  "%date,
    "python scripts/LQ_get_limits.py --chan ee --q d  -o limits/ --ending %s  "%date,
    "python scripts/LQ_get_limits.py --chan mumu --q u -o limits/ --ending %s "%date,
    "python scripts/LQ_get_limits.py --chan mumu --q d -o limits/ --ending %s "%date,
    "python scripts/LQ_get_limits.py --chan ee --q u --vec True -o limits/ --ending %s "%date,
    "python scripts/LQ_get_limits.py --chan ee --q d --vec True -o limits/ --ending %s  "%date,
    "python scripts/LQ_get_limits.py --chan mumu --q u --vec True -o limits/ --ending %s "%date,
    "python scripts/LQ_get_limits.py --chan mumu --q d --vec True -o limits/ --ending %s  "%date,
 
    ]

    labels = [
        "limits_ee_u","limits_ee_d","limits_mumu_u","limits_mumu_d",
        "limits_ee_u_vec","limits_ee_d_vec","limits_mumu_u_vec","limits_mumu_d_vec"
        #"limits_ee_s","limits_mumu_s",
        #"limits_ee_d_vec","limits_mumu_d_vec"
    ]


    cpy_cmd = "xrdcp -f limits/* $1 \n"

    for i,cmd in enumerate(cmds):
	   for m in range(5000,5500,500):
	    #for point in np.arange(0.28,1.5,0.005):
	    #for q in [0.025,0.16,0.5,0.84,0.975]:
            #regular templates
            script_name = "scripts/script3.sh"
            print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
            script_file = open(script_name, 'a+')
            script_file.write("mkdir limits\n")
            script_file.write(cmd+" --mLQ %i \n"%(m))
            script_file.write(cpy_cmd)
            script_file.close()
            #print_and_do("cat %s" % script_name)
            print_and_do("chmod +x %s" % script_name)
            print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s_m%i_%s"  % (n_m_bins, script_name, labels[i], m, date))
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
  
    "python scripts/LQ_do_impacts.py --mLQ 2000 --chan ee --q u -o imps --ending %s \n"%date,
    "python scripts/LQ_do_impacts.py --mLQ 2000 --chan ee --q d -o imps --ending %s  \n"%date,
    "python scripts/LQ_do_impacts.py --mLQ 2000 --chan mumu --q u -o imps --ending %s  \n"%date,
    "python scripts/LQ_do_impacts.py --mLQ 2000 --chan mumu --q d -o imps --ending %s \n "%date,
    "python scripts/LQ_do_impacts.py --mLQ 2000 --chan ee --q u --vec True -o imps --ending %s \n "%date,
    "python scripts/LQ_do_impacts.py --mLQ 2000 --chan ee --q d --vec True -o imps --ending %s  \n"%date,
    "python scripts/LQ_do_impacts.py --mLQ 2000 --chan mumu --q u --vec True -o imps --ending %s  \n"%date,
    "python scripts/LQ_do_impacts.py --mLQ 2000 --chan mumu --q d --vec True -o imps --ending %s  \n"%date,
 
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
  
    "python scripts/LQ_check_sys_uncs.py --mLQ 2000 --chan ee --q u -o sys --ending %s \n"%date,
    "python scripts/LQ_check_sys_uncs.py --mLQ 2000 --chan ee --q d -o sys --ending %s  \n"%date,
    "python scripts/LQ_check_sys_uncs.py --mLQ 2000 --chan mumu --q u -o sys --ending %s  \n"%date,
    "python scripts/LQ_check_sys_uncs.py --mLQ 2000 --chan mumu --q d -o sys --ending %s \n "%date,
    "python scripts/LQ_check_sys_uncs.py --mLQ 2000 --chan ee --q u --vec True -o sys --ending %s \n "%date,
    "python scripts/LQ_check_sys_uncs.py --mLQ 2000 --chan ee --q d --vec True -o sys --ending %s  \n"%date,
    "python scripts/LQ_check_sys_uncs.py --mLQ 2000 --chan mumu --q u --vec True -o sys --ending %s  \n"%date,
    "python scripts/LQ_check_sys_uncs.py --mLQ 2000 --chan mumu --q d --vec True -o sys --ending %s  \n"%date,
 
    ]

    labels = [
        "sys_ee_u","sys_ee_d","sys_mumu_u","sys_mumu_d",
        "sys_ee_u_vec","sys_ee_d_vec","sys_mumu_u_vec","sys_mumu_d_vec"
        #"limits_ee_s","limits_mumu_s",
        #"limits_ee_s_vec","limits_mumu_s_vec"
    ]


    cpy_cmd = "xrdcp -f sys/* $1 \n"

    for i,cmd in enumerate(cmds):
    #for m in range(1000,9500,500):
    #for point in np.arange(0.28,1.5,0.005):
    #for q in [0.025,0.16,0.5,0.84,0.975]:
        #regular templates
        script_name = "scripts/script3.sh"
        print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
        script_file = open(script_name, 'a+')
        script_file.write("mkdir sys\n")
        script_file.write(cmd)
        script_file.write(cpy_cmd)
        script_file.close()
        #print_and_do("cat %s" % script_name)
        print_and_do("chmod +x %s" % script_name)
        print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s"  % (n_m_bins, script_name, labels[i]))
        print_and_do("rm scripts/script3.sh")


if options.likelihood:

    cmds = [
  

    "python scripts/LQ_do_likelihood.py --mLQ 2000 --chan mumu --q d -o likelihood_scans ",

 
    ]

    labels = [
        "likelihood_mumu_d"
    ]


    cpy_cmd = "xrdcp -f likelihood_scans/* $1 \n"
    poi_list = ["muISOBAR17", "muISOEND17","muIDEND17","muIDBAR17", "mufakesrw1b17", "mufakesrw2b17", "mufakesrw3b17", "mufakesrw4b17", "muISOBAR16", "muISOEND16","muIDEND16","muIDBAR16", "mufakesrw1b18", "mufakesrw2b18", "mufakesrw3b18", "mufakesrw4b18"]
    '''
    poi_list = ["dy_xsec","nlo_sys", "muISOBAR18", "muISOEND18", "muIDSYS", "muIDEND18", "muISOSYS", "muIDBAR18", "mufakesrw1b16", "mufakesrw2b16", "mufakesrw3b16", "mufakesrw4b16"]
    
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
        for i,cmd in enumerate(cmds):
        #for m in range(1000,9500,500):
        #for point in np.arange(0.28,1.5,0.005):
        #for q in [0.025,0.16,0.5,0.84,0.975]:
            #regular templates
            script_name = "scripts/script3.sh"
            print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
            script_file = open(script_name, 'a+')
            script_file.write("mkdir likelihood_scans\n")
            script_file.write(cmd + " --poi %s \n "%poi)
            script_file.write(cpy_cmd)
            script_file.close()
            #print_and_do("cat %s" % script_name)
            print_and_do("chmod +x %s" % script_name)
            print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s_%s"  % (n_m_bins, script_name, labels[i], poi))
            print_and_do("rm scripts/script3.sh")
