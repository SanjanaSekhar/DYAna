
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

def print_and_do(s):
    print(s)
    return os.system(s)

masses=[1500.]
years = [2016, 2017,2018]
#types = [1]
#labels = ["sys"]
njobs = 1
ending = "aug20"
for mass in masses:
    for year in years:
            script_name = "scripts/script2.sh"
            print_and_do("cp scripts/LQ_nonsys_template.sh %s" % script_name)
            print_and_do("sed -i  s/YEAR/%i/g %s" % (year, script_name))
            print_and_do("sed -i  s/MASS/%i/g %s" % (mass, script_name))
            print_and_do("chmod +x %s" % script_name)
            print_and_do("python LQ_doCondor.py --njobs %i --combine --sub  -s %s -n templ%i_nonsys_%s_m%i"  % (njobs, script_name,  year % 2000,  ending, mass))
            print_and_do("rm scripts/script2.sh")

