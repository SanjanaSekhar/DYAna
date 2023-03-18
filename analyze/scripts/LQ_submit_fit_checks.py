import numpy as np
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

def print_and_do(s):
    print(s)
    return os.system(s)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--ending",  default="021223", help="date")
parser.add_option("--limits",  default=False, help="do limits")
(options, args) = parser.parse_args()

n_m_bins = 1
date = options.ending


if not options.limits:

    mLQ_list = [5000]
    for mLQ in mLQ_list:
        cmds = [
        "python scripts/LQ_do_bias_test.py  --yLQ 0.0 --chan ee --q d --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 0.2 --chan ee --q d --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 0.4 --chan ee --q d --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
	"python scripts/LQ_do_bias_test.py  --yLQ 0.6 --chan ee --q d --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 0.0 --chan ee --q u --is_vec True --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 0.2 --chan ee --q u --is_vec True --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 0.4 --chan ee --q u --is_vec True --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
	"python scripts/LQ_do_bias_test.py  --yLQ 0.6 --chan ee --q u --is_vec True --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 0.0 --chan mumu --q u --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 0.2 --chan mumu --q u --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 0.4 --chan mumu --q u --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
	"python scripts/LQ_do_bias_test.py  --yLQ 0.6 --chan mumu --q u --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 0.0 --chan mumu --q d --is_vec True --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 0.2 --chan mumu --q d --is_vec True --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        "python scripts/LQ_do_bias_test.py  --yLQ 0.4 --chan mumu --q d --is_vec True --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
	"python scripts/LQ_do_bias_test.py  --yLQ 0.6 --chan mumu --q d --is_vec True --nToys 50 -o temp/ --mLQ %i --ending %s"%(mLQ,date),
        ]

        labels = [
                "bias_test_yLQ0.0_ee_d_m",  "bias_test_yLQ0.2_ee_d_m",  "bias_test_yLQ0.4_ee_d_m",  "bias_test_yLQ0.6_ee_d_m", 
		"bias_test_yLQ0.0_ee_u_vec_m", "bias_test_yLQ0.2_ee_u_vec_m","bias_test_yLQ0.4_ee_u_vec_m", "bias_test_yLQ0.8_ee_u_vec_m",
        	"bias_test_yLQ0.0_mumu_u_m", "bias_test_yLQ0.2_mumu_u_m", "bias_test_yLQ0.4_mumu_u_m","bias_test_yLQ0.6_mumu_u_m",  
		"bias_test_yLQ0.0_mumu_d_vec_m", "bias_test_yLQ0.2_mumu_d_vec_m", "bias_test_yLQ0.4_mumu_d_vec_m", "bias_test_yLQ0.6_mumu_d_vec_m"
                
        ]

        cpy_cmd = "xrdcp -f temp/* $1 \n"

        total_jobs = 300

        for i,cmd in enumerate(cmds):

            for job_idx in range(total_jobs/50):

                #regular templates
                script_name = "scripts/script3.sh"
                print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
                script_file = open(script_name, 'a+')
                script_file.write("mkdir temp\n")
                script_file.write(cmd + " --job %i\n"%job_idx)
                script_file.write(cpy_cmd)
                script_file.close()
                #print_and_do("cat %s" % script_name)
                print_and_do("chmod +x %s" % script_name)
                print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s%i_%i_%s"  % (n_m_bins, script_name, labels[i], mLQ, job_idx, date))
                print_and_do("rm scripts/script3.sh")

else:
    cmds = [
  
    "python scripts/LQ_get_limits.py --chan ee --q u  -o limits/ --ending %s --ntoys 500 --iterations 10 "%date,
    "python scripts/LQ_get_limits.py --chan ee --q d  -o limits/ --ending %s  --ntoys 500 --iterations 10 "%date,
    "python scripts/LQ_get_limits.py --chan mumu --q u -o limits/ --ending %s  --ntoys 500 --iterations 10 "%date,
    "python scripts/LQ_get_limits.py --chan mumu --q d -o limits/ --ending %s  --ntoys 500 --iterations 10 "%date,
    "python scripts/LQ_get_limits.py --chan ee --q u --vec True -o limits/ --ending %s  --ntoys 500 --iterations 10 "%date,
    "python scripts/LQ_get_limits.py --chan ee --q d --vec True -o limits/ --ending %s  --ntoys 500 --iterations 10 "%date,
    "python scripts/LQ_get_limits.py --chan mumu --q u --vec True -o limits/ --ending %s  --ntoys 500 --iterations 10 "%date,
    "python scripts/LQ_get_limits.py --chan mumu --q d --vec True -o limits/ --ending %s  --ntoys 500 --iterations 10 "%date,
 
    ]

    labels = [
        "limits_ee_u","limits_ee_d","limits_mumu_u","limits_mumu_d",
        "limits_ee_u_vec","limits_ee_d_vec","limits_mumu_u_vec","limits_mumu_d_vec"
        #"limits_ee_s","limits_mumu_s",
        #"limits_ee_s_vec","limits_mumu_s_vec"
    ]


    cpy_cmd = "xrdcp -f limits/* $1 \n"

    for i,cmd in enumerate(cmds):
	for point in np.arange(0.2,1.0,0.05):
	    for q in [0.025,0.16,0.5,0.84,0.975]:
                #regular templates
                script_name = "scripts/script3.sh"
                print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
                script_file = open(script_name, 'a+')
                script_file.write("mkdir limits\n")
                script_file.write(cmd+" --inject_yLQ2 %f --quantile %f\n"%(point,q))
                script_file.write(cpy_cmd)
                script_file.close()
                #print_and_do("cat %s" % script_name)
                print_and_do("chmod +x %s" % script_name)
                print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s_yLQ2%.2f_q%.3f_%s"  % (n_m_bins, script_name, labels[i], point, q, date))
                print_and_do("rm scripts/script3.sh")

