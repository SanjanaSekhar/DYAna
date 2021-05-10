
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

def print_and_do(s):
    print(s)
    return os.system(s)

cmds = [
"python scripts/do_gof.py --nToys 200 -o temp/  --mbin $3 --teststat saturated\n",
"python scripts/do_impacts.py -o temp/  --mbin $3 \n",
#"python scripts/do_bias_test.py --nToys 100 -o temp/  --mbin $3 --Afb 0.6 --A0 0.00 \n",
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

labels = ["gof_postfit", 
        "impacts",  
        #"bias_test_afb6_a00", "bias_test_afb6_a01", "bias_test_afb0_a00", "bias_test_afb0_a01",
]


cpy_cmd = "xrdcp -f temp/* $1 \n"

date = "may8"
n_m_bins = 8


for i,cmd in enumerate(cmds):

    #regular templates
    script_name = "scripts/script3.sh"
    print_and_do("cp scripts/combine_template.sh %s" % script_name)
    script_file = open(script_name, 'a+')
    script_file.write(cmd)
    script_file.write(cpy_cmd)
    script_file.close()
    #print_and_do("cat %s" % script_name)
    print_and_do("chmod +x %s" % script_name)
    print_and_do("python doCondor.py --njobs %i --combine --sub --no_rename  -s %s -n %s_%s"  % (n_m_bins, script_name, labels[i], date))
    print_and_do("rm scripts/script3.sh")


