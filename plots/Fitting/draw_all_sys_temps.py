from ROOT import *
import os
from optparse import OptionParser
from optparse import OptionGroup

parser = OptionParser()
parser.add_option("-i", "--fin", default='', help="Input file with templates")
parser.add_option("-o", "--plot_dir", default='../plots/', help="Directory to output plots")
parser.add_option("-y", "--year", type = 'int', default=16, help="Year")



(options, args) = parser.parse_args()

if(options.year > 2000):
    options.year = options.year % 2000

n_m_bins = 8


base_strs = [ 'ee%i_fpl', 'ee%i_fmn', 'mumu%i_fmn', 'mumu%i_fpl']
#my_excludes = ['RENORM', 'REFAC', 'FAC']
my_excludes = []

gROOT.SetBatch(1)

f = TFile.Open(options.fin)

if(not os.path.exists(options.plot_dir)):
    print("Making directory %s" % options.plot_dir)
    os.system("mkdir %s" % options.plot_dir)


l_base = TLine(0,0,1,1)
l_base.SetLineColor(kBlack)
l_base.SetLineWidth(7)

l_rf = TLine(0,0,1,1)
l_rf.SetLineColor(kRed)
l_rf.SetLineWidth(7)

l_pdf = TLine(0,0,1,1)
l_pdf.SetLineColor(kGreen+4)
l_pdf.SetLineWidth(7)

l_id = TLine(0,0,1,1)
l_id.SetLineColor(kBlue)
l_id.SetLineWidth(7)

l_other = TLine(0,0,1,1)
l_other.SetLineColor(kMagenta)
l_other.SetLineWidth(7)

for mbin in range(n_m_bins):
    gDirectory.cd("w%i" % mbin)
    keys = gDirectory.GetListOfKeys()
    for base in base_strs:
        if('%i' in base):
            base = base % options.year
        print("Doing plot for sys %s " % base)
        c1 = TCanvas("c1", "", 1600, 1000) 
        for key in keys:
            key_name = key.GetName()
            skip = False
            for exc in my_excludes:
                if (exc in key_name):
                    skip = True
            if(skip): continue
            if (base in key_name):
                #print("Adding key %s" % key_name)
                h = gDirectory.Get(key_name)
                color = kMagenta
                if('pdf' in key_name):
                    color = kGreen+4
                if('RENORM' in key_name or 'FAC' in key_name or 'REFAC' in key_name):
                    color = kRed
                if('ID' in key_name or 'RECO' in key_name or 'HLT' in key_name or 'ISO' in key_name ):
                    color = kBlue
                h.SetLineColor(color)
                h.SetLineWidth(1)
                h.Draw("hist same")
        h_base = gDirectory.Get(base)
        h_base.SetLineColor(kBlack)
        h_base.SetLineWidth(3)
        h_base.SetTitle(base)
        h_base.Draw("hist same")


        leg = TLegend(0.2,0.2)
        leg.SetNColumns(2)
        leg.AddEntry(l_base, "Nominal", "l")
        leg.AddEntry(l_rf, "Renorm. + Fac. ", "l")
        leg.AddEntry(l_id, "Lepton Efficiencies", "l")
        leg.AddEntry(l_pdf, "pdf", "l")
        leg.AddEntry(l_other, "Other", "l")
        leg.Draw()


        c1.Print(options.plot_dir + ("Y%i_" % options.year) + ("mbin%i_" % mbin) + base + "_all_sys.png")

    gDirectory.cd("..")

