
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

def print_and_do(s):
    print(s)
    return os.system(s)

cmds = [
#"python scripts/check_sys_uncs.py -o temp/  --mbin $3 --diff \n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit\n",
#"python scripts/check_sys_uncs.py -o temp/  --mbin $3 \n",
#"python scripts/do_impacts.py -o temp/  --mbin $3  --nThreads 4\n",
#"python scripts/do_impacts.py -o temp/  --mbin $3  --nThreads 4 --diff\n",
#"python scripts/do_gof.py --nToys 1000 -o temp/  --mbin $3 --teststat saturated\n",
#"python scripts/do_gof.py --nToys 500 -o temp/  --mbin $3 --teststat saturated --prefit --mask_ee \n",
#"python scripts/do_gof.py --nToys 500 -o temp/  --mbin $3 --teststat saturated --prefit --mask_mumu \n",
#"python scripts/do_impacts.py -o temp/  --mbin $3  --nThreads 4\n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit --noSymMCStats \n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit --mask_ee -y 16\n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit --mask_ee -y 17\n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit --mask_ee -y 18\n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit --mask_mumu -y 16\n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit --mask_mumu -y 17\n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit --mask_mumu -y 18\n",
#"python scripts/do_impacts.py -o temp/  --mbin $3 --nThreads 4 --expected \n",
"python scripts/LQ_do_bias_test.py  --yLQ 0.0 --chan ee --q u -o temp/\n",
"python scripts/LQ_do_bias_test.py  --yLQ 0.5 --chan ee --q u -o temp/\n",
"python scripts/LQ_do_bias_test.py  --yLQ 1.0 --chan ee --q u -o temp/\n",
"python scripts/LQ_do_bias_test.py  --yLQ 0.0 --chan mumu --q d -o temp/\n",
"python scripts/LQ_do_bias_test.py  --yLQ 0.5 --chan mumu --q d -o temp/\n",
"python scripts/LQ_do_bias_test.py  --yLQ 1.0 --chan mumu --q d -o temp/\n",
#"python scripts/LQ_do_bias_test.py --nToys 100  --yLQ 0.0 --chan ee --q u \n",
#"python scripts/do_bias_test.py --nToys 100 -o temp/  --mbin $3 --Afb 0.6 --A0 0.1 \n",
#"python scripts/do_bias_test.py --nToys 100 -o temp/  --mbin $3 --Afb 0.0 --A0 0.00 \n",
#"python scripts/do_bias_test.py --nToys 100 -o temp/  --mbin $3 --Afb 0.0 --A0 0.1 \n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated -y 2016 \n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated -y 2017 \n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated -y 2018 \n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated -y 2016 --mask_mumu\n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated -y 2016 --mask_ee\n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit -y 2016 \n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit -y 2017 \n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit -y 2018 \n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit -y 2018 --mask_mumu \n",
#"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated --prefit -y 2018 --mask_ee \n"
]

labels = [
        #"sys_uncs_diff",
        #"gof_prefit",
        #"sys_uncs",
        #"impacts",  
        #"impacts_diff",  
        #'gof_prefit_mumu',
        #'gof_prefit_ee',
        #'gof_prefit_1000toys',
        #'gof_postfit_1000toys',
        #'gof_noSym',
        #'gof_mumu16',
        #'gof_mumu17',
        #'gof_mumu18',
        #'gof_ee16',
        #'gof_ee17',
        #'gof_ee18',
        #"gof_prefit",
        #"gof_postfit", 
        #"expected_impacts",  
        "bias_test_yLQ0.0_ee_u",  "bias_test_yLQ0.5_ee_u",  "bias_test_yLQ1.0_ee_u",  "bias_test_yLQ0.0_mumu_d", "bias_test_yLQ0.5_mumu_d", "bias_test_yLQ1.0_mumu_d"
]


cpy_cmd = "xrdcp -f temp/* $1 \n"

date = "062422"
n_m_bins = 1


for i,cmd in enumerate(cmds):

    #regular templates
    script_name = "scripts/script3.sh"
    print_and_do("cp scripts/LQ_combine_template.sh %s" % script_name)
    script_file = open(script_name, 'a+')
    script_file.write(cmd)
    script_file.write(cpy_cmd)
    script_file.close()
    #print_and_do("cat %s" % script_name)
    print_and_do("chmod +x %s" % script_name)
    print_and_do("python LQ_doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s_%s"  % (n_m_bins, script_name, labels[i], date))
    print_and_do("rm scripts/script3.sh")


