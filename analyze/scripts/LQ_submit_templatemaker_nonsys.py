
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

def print_and_do(s):
    print(s)
    return os.system(s)

masses=[2000.,2500.]#,3000.,3500.,4000.,4500.,5000.]#,5500.,6000.,6500.,7000.,7500.,8000.,8500.,9000.]
#masses=[1000.]
years = [2016,2017,2018]
#types = [1]
#labels = ["sys"]
njobs = 1
ending = "122122_cs"
for mass in masses:
    for year in years:
            script_name = "scripts/script2.sh"
            print_and_do("cp scripts/LQ_nonsys_template.sh %s" % script_name)
            print_and_do("sed -i  s/YEAR/%i/g %s" % (year, script_name))
            print_and_do("sed -i  s/MASS/%i/g %s" % (mass, script_name))
            print_and_do("chmod +x %s" % script_name)
            print_and_do("python LQ_doCondor.py --njobs %i --combine --sub  -s %s -n templ%i_nonsys_%s_m%i"  % (njobs, script_name,  year % 2000,  ending, mass))
            print_and_do("rm scripts/script2.sh")

