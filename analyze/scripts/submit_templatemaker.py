
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

def print_and_do(s):
    print(s)
    return os.system(s)

years = [2016, 2017,2018]
#years = [2016]
#types = []
types = [0, 1]
labels = ["pdf", "sys"]
#types = [1]
#labels = ["sys"]
njobs = 25
ending = "sep1_unblinded"

for year in years:

    #regular templates
    script_name = "scripts/script2.sh"
    print_and_do("cp scripts/templ_template.sh %s" % script_name)
    print_and_do("sed -i  s/YEAR/%i/g %s" % (year, script_name))
    print_and_do("chmod +x %s" % script_name)
    print_and_do("python doCondor.py --njobs %i --combine --sub  -s %s -n templ%i_%s"  % (1, script_name,  year % 2000, ending))
    print_and_do("rm scripts/script2.sh")


    #systematics
    for j,ttype in enumerate(types):
        script_name = "scripts/script2.sh"
        print_and_do("cp scripts/templ_sys_template.sh %s" % script_name)
        print_and_do("sed -i  s/YEAR/%i/g %s" % (year, script_name))
        print_and_do("sed -i  s/TYPE/%i/g %s" % (ttype, script_name))
        print_and_do("chmod +x %s" % script_name)
        print_and_do("python doCondor.py --njobs %i --combine --sub  -s %s -n templ%i_%s_%s"  % (njobs, script_name,  year % 2000, labels[j], ending))
        print_and_do("rm scripts/script2.sh")


