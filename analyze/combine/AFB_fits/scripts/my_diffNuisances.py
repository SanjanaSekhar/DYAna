#!/usr/bin/env python
import re
from sys import argv, stdout, stderr, exit
import datetime
from optparse import OptionParser
import HiggsAnalysis.CombinedLimit.calculate_pulls as CP 

# tool to compare fitted nuisance parameters to prefit values.
#
# Also used to check for potential problems in RooFit workspaces to be used with combine
# (see https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsPAGPreapprovalChecks)

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
hasHelp = False
for X in ("-h", "-?", "--help"):
    if X in argv:
        hasHelp = True
        argv.remove(X)
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetLabelSize(.03, "XY")
#ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
argv.remove( '-b-' )
if hasHelp: argv.append("-h")


parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--vtol", "--val-tolerance", dest="vtol", default=0.30, type="float", help="Report nuisances whose value changes by more than this amount of sigmas")
parser.add_option("--stol", "--sig-tolerance", dest="stol", default=0.10, type="float", help="Report nuisances whose sigma changes by more than this amount")
parser.add_option("--vtol2", "--val-tolerance2", dest="vtol2", default=2.0, type="float", help="Report severely nuisances whose value changes by more than this amount of sigmas")
parser.add_option("--stol2", "--sig-tolerance2", dest="stol2", default=0.50, type="float", help="Report severely nuisances whose sigma changes by more than this amount")
parser.add_option("-a", "--all",      dest="show_all_parameters",    default=False,  action="store_true", help="Print all nuisances, even the ones which are unchanged w.r.t. pre-fit values.")
parser.add_option("-A", "--abs",      dest="absolute_values",    default=False,  action="store_true", help="Report also absolute values of nuisance values and errors, not only the ones normalized to the input sigma")
parser.add_option("-p", "--poi",      dest="poi",    default="r",    type="string",  help="Name of signal strength parameter (default is 'r' as per text2workspace.py)")
parser.add_option("-f", "--format",   dest="format", default="text", type="string",  help="Output format ('text', 'latex', 'twiki'")
parser.add_option("-g", "--histogram", dest="plotfile", default=None, type="string", help="If true, plot the pulls of the nuisances to the given location")
parser.add_option("", "--pullDef",  dest="pullDef", default="", type="string", help="Choose the definition of the pull, see python/calculate_pulls.py for options")
parser.add_option("", "--skipFitS", dest="skipFitS", default=False, action='store_true', help="skip the S+B fit, instead the B-only fit will be repeated")
parser.add_option("", "--skipFitB", dest="skipFitB", default=False, action='store_true', help="skip the B-only fit, instead the S+B fit will be repeated")
parser.add_option("", "--multidim", dest="multidim", default=False, action='store_true', help="Do procedure for MultDimFit")
parser.add_option("", "--prefit", dest="prefit", default="", type="string", help="Look up prefit nuisances in this separate file")
parser.add_option("--mLQ", "-m", type = 'int', default = 0, help="Mass bin (for plot label)")
(options, args) = parser.parse_args()

mLQ = options.mLQ

if len(args) == 0:
    parser.print_usage()
    exit(1)

if options.pullDef!="" and options.pullDef not in CP.allowed_methods(): exit("Method %s not allowed, choose one of [%s]"%(options.pullDef,",".join(CP.allowed_methods())))

if options.pullDef and options.absolute_values : 
  print "Pulls are always defined as absolute, will modify --absolute_values to False for you"
  options.absolute_values = False 

if options.pullDef : options.show_all_parameters=True

setUpString = "diffNuisances run on %s, at %s with the following options ... "%(args[0],datetime.datetime.utcnow())+str(options)

file_prefit = None
file = ROOT.TFile(args[0])
if(options.prefit != ""): 
    file_prefit = ROOT.TFile(options.prefit)
    if(file_prefit == None): raise RuntimeError, "Cannot open prefit file %s" % args[0]
if file == None: raise RuntimeError, "Cannot open file %s" % args[0]
if(not options.multidim):
    fit_s  = file.Get("fit_s") if not options.skipFitS  else file.Get("fit_b")
    fit_b  = file.Get("fit_b") if not options.skipFitB  else file.Get("fit_s")
else:
    fit_s = fit_b = file.Get("fit_mdf")
    options.skipFitB = True
if(file_prefit == None):
    prefit = file.Get("nuisances_prefit")
else:
    prefit = file_prefit.Get("nuisances_prefit")

if fit_s == None or fit_s.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the signal fit 'fit_s'"     % args[0]
if fit_b == None or fit_b.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the background fit 'fit_b'" % args[0]
if prefit == None or prefit.ClassName() != "RooArgSet":    raise RuntimeError, "File %s does not contain the prefit nuisances 'nuisances_prefit'"  % args[0]


#hardcoded list of nuissances we want to plot
#pars = ["alphaS","alphaDen","RENORM","FAC","REFAC", "dy_xsec","bk_xsec","gam_xsec",
#"Pu16","BTAG16","elScaleStat16","elScaleSyst16","elScaleGain16","elSmear16","elHLT16","elID16","elRECO16","muRC16","muID16","muHLT16","lumi16","ee16_fakes_norm","mu16_fakes_norm","R_ee16_os_fakes",
#"Pu17","BTAG17","elScaleStat17","elScaleSyst17","elScaleGain17","elSmear17","elHLT17","elID17","elRECO17","muRC17","muID17","muHLT17","lumi17","ee17_fakes_norm","mu17_fakes_norm","R_ee17_os_fakes",
#"Pu18","BTAG18","elScaleStat18","elScaleSyst18","elScaleGain18","elSmear18","elHLT18","elID18","elRECO18","muRC18","muID18","muHLT18","lumi18","ee18_fakes_norm","mu18_fakes_norm","R_ee18_os_fakes"]

exclude_pars = ["prop", "pdf", "A0", "Afb"]


isFlagged = {}

# maps from nuisance parameter name to the row to be printed in the table
table = {}

# get the fitted parameters
fpf_b = fit_b.floatParsFinal()
fpf_s = fit_s.floatParsFinal()

pulls = []

nuis_p_i=0
title = "pull" if options.pullDef else "#theta"

"""
def getGraph(hist,shift): 
   gr = ROOT.TGraphAsymErrors()
   gr.SetName(hist.GetName())
   for i in range(hist.GetNbinsX()):
     x = hist.GetBinCenter(i+1)+shift
     y = hist.GetBinContent(i+1)
     e = hist.GetBinError(i+1)
     gr.SetPoint(i,x,y)
     gr.SetPointError(i,float(abs(shift))*0.8,e)
   return gr
"""

a_nuis_s = []
a_nuis_b = []
a_nuis_p = []
a_name = []

for i in range(fpf_s.getSize()):
    exclude = False
    nuis_s = fpf_s.at(i)
    name   = nuis_s.GetName();
    for label in exclude_pars:
        if(label in name): exclude = True
    if(not exclude):
        nuis_b = fpf_b.find(name)
        nuis_p = prefit.find(name)
        a_name.append(name)
        a_nuis_s.append(nuis_s)
        a_nuis_b.append(nuis_b)
        a_nuis_p.append(nuis_p)




# Also make histograms for pull distributions:
n_bins = len(a_name)
hist_fit_b  = ROOT.TH1F("fit_b"   ,"B-only fit Nuisances;;%s "%title,n_bins,0, n_bins)
hist_fit_s  = ROOT.TH1F("fit_s"   ,"S+B fit Nuisances   ;;%s "%title,n_bins,0,n_bins)
hist_prefit = ROOT.TH1F("prefit_nuisancs","Prefit Nuisances    ;;%s "%title,n_bins,0,n_bins)
# Store also the *asymmetric* uncertainties
gr_fit_b    = ROOT.TGraphAsymmErrors(); gr_fit_b.SetTitle("fit_b_g")
gr_fit_s    = ROOT.TGraphAsymmErrors(); gr_fit_s.SetTitle("fit_b_s")


# loop over all fitted parameters
for i in range(len(a_name)):

    name = a_name[i]
    nuis_s = a_nuis_s[i]
    nuis_b = a_nuis_b[i]
    nuis_p = a_nuis_p[i]
    #print(name)

    # keeps information to be printed about the nuisance parameter
    row = []

    flag = False;
    mean_p, sigma_p, sigma_pu, sigma_pd = 0,0,0,0

    if nuis_p == None:
        # nuisance parameter NOT present in the prefit result
        if not options.absolute_values and not (options.pullDef=="unconstPullAsym"): continue
        row += [ "[%.2f, %.2f]" % (nuis_s.getMin(), nuis_s.getMax()) ]

    else:
        # get best-fit value and uncertainty at prefit for this 
        # nuisance parameter
 	if nuis_p.getErrorLo()==0 : nuis_p.setError(nuis_p.getErrorHi())
        mean_p, sigma_p, sigma_pu,sigma_pd = (nuis_p.getVal(), nuis_p.getError(),nuis_p.getErrorHi(),nuis_p.getErrorLo())

	if not sigma_p > 0: sigma_p = (nuis_p.getMax()-nuis_p.getMin())/2
	nuisIsSymm = abs(abs(nuis_p.getErrorLo())-abs(nuis_p.getErrorHi()))<0.01 or nuis_p.getErrorLo() == 0
        if options.absolute_values: 
		if nuisIsSymm : row += [ "%.6f +/- %.6f" % (nuis_p.getVal(), nuis_p.getError()) ]
		else: row += [ "%.6f +%.6f %.6f" % (nuis_p.getVal(), nuis_p.getErrorHi(), nuis_p.getErrorLo()) ]

    for fit_name, nuis_x in [('b', nuis_b), ('s',nuis_s)]:
        if nuis_x == None:
            row += [ " n/a " ]
        else:
	    nuisIsSymm = abs(abs(nuis_x.getErrorLo())-abs(nuis_x.getErrorHi()))<0.01 or nuis_x.getErrorLo() == 0
 	    if nuisIsSymm : nuis_x.setError(nuis_x.getErrorHi())
	    nuiselo = abs(nuis_x.getErrorLo()) if nuis_x.getErrorLo()>0 else nuis_x.getError()
	    nuisehi = nuis_x.getErrorHi()
	    if options.pullDef and nuis_p!=None: 
	    	nx,ned,neu = CP.returnPullAsym(options.pullDef,nuis_x.getVal(),mean_p,nuisehi,sigma_pu,abs(nuiselo),abs(sigma_pd))
	    else: 
	        nx,ned,neu = nuis_x.getVal(), nuiselo, nuisehi

            if nuisIsSymm : row += [ "%+.2f +/- %.2f" % (nx, (abs(ned)+abs(neu))/2) ]
	    else: row += [ "%+.2f +%.2f %.2f" % (nx, neu, ned) ]

            if nuis_p != None:
	        if options.plotfile: 
	          if fit_name=='b':
	    	    nuis_p_i+=1
		    if options.pullDef and nuis_p!=None:
		      #nx,ned,neu = CP.returnPullAsym(options.pullDef,nuis_x.getVal(),mean_p,nuis_x.getErrorHi(),sigma_pu,abs(nuis_x.getErrorLo()),abs(sigma_pd))
		      gr_fit_b.SetPoint(nuis_p_i-1,nuis_p_i-0.5+0.1,nx)
		      gr_fit_b.SetPointError(nuis_p_i-1,0,0,ned,neu)
		    else:
		      gr_fit_b.SetPoint(nuis_p_i-1,nuis_p_i-0.5+0.1,nuis_x.getVal())
		      gr_fit_b.SetPointError(nuis_p_i-1,0,0,abs(nuis_x.getErrorLo()),nuis_x.getErrorHi())
	      	    hist_fit_b.SetBinContent(nuis_p_i,nuis_x.getVal())
	      	    hist_fit_b.SetBinError(nuis_p_i,nuis_x.getError())
	      	    hist_fit_b.GetXaxis().SetBinLabel(nuis_p_i,name)
	      	    gr_fit_b.GetXaxis().SetBinLabel(nuis_p_i,name)
	          if fit_name=='s':
		    if options.pullDef and nuis_p!=None:
		      #nx,ned,neu = CP.returnPullAsym(options.pullDef,nuis_x.getVal(),mean_p,nuis_x.getErrorHi(),sigma_pu,abs(nuis_x.getErrorLo()),abs(sigma_pd))
		      gr_fit_s.SetPoint(nuis_p_i-1,nuis_p_i-0.5-0.1,nx)
		      gr_fit_s.SetPointError(nuis_p_i-1,0,0,ned,neu)
		    else:
		      gr_fit_s.SetPoint(nuis_p_i-1,nuis_p_i-0.5-0.1,nuis_x.getVal())
		      gr_fit_s.SetPointError(nuis_p_i-1,0,0,abs(nuis_x.getErrorLo()),nuis_x.getErrorHi())
	      	    hist_fit_s.SetBinContent(nuis_p_i,nuis_x.getVal())
	      	    hist_fit_s.SetBinError(nuis_p_i,nuis_x.getError())
	      	    hist_fit_s.GetXaxis().SetBinLabel(nuis_p_i,name)
	      	    gr_fit_s.GetXaxis().SetBinLabel(nuis_p_i,name)
		  hist_prefit.SetBinContent(nuis_p_i,mean_p)
		  hist_prefit.SetBinError(nuis_p_i,sigma_p)
	      	  hist_prefit.GetXaxis().SetBinLabel(nuis_p_i,name)

                if sigma_p>0: 
                    if options.pullDef:
			valShift = nx 
			sigShift = 1
                    else: 
		        # calculate the difference of the nuisance parameter
                        # w.r.t to the prefit value in terms of the uncertainty
                        # on the prefit value
			valShift = (nuis_x.getVal() - mean_p)/sigma_p

                        # ratio of the nuisance parameter's uncertainty
                        # w.r.t the prefit uncertainty
                	sigShift = nuis_x.getError()/sigma_p

		else :
			#print "No definition for prefit uncertainty %s. Printing absolute shifts"%(nuis_p.GetName())
			valShift = (nuis_x.getVal() - mean_p)
                	sigShift = nuis_x.getError()

                if options.pullDef:
                    row[-1] += ""
		elif options.absolute_values:
                    row[-1] += " (%+4.2fsig, %4.2f)" % (valShift, sigShift)
                else:
                    row[-1] = " %+4.2f, %4.2f" % (valShift, sigShift)
                
		if fit_name == 'b':
                      pulls.append(valShift)

                if (abs(valShift) > options.vtol2 or abs(sigShift-1) > options.stol2):

                    # severely report this nuisance:
                    # 
                    # the best fit moved by more than 2.0 sigma or the uncertainty (sigma)
                    # changed by more than 50% (default thresholds) w.r.t the prefit values

                    isFlagged[(name,fit_name)] = 2

                    flag = True

                elif (abs(valShift) > options.vtol  or abs(sigShift-1) > options.stol):

                    # report this nuisance:
                    # 
                    # the best fit moved by more than 0.3 sigma or the uncertainty (sigma)
                    # changed by more than 10% (default thresholds) w.r.t the prefit values

                    if options.show_all_parameters: isFlagged[(name,fit_name)] = 1

                    flag = True

                elif options.show_all_parameters:
                    flag = True
	
    # end of loop over s and b

    row += [ "%+4.2f"  % fit_s.correlation(name, options.poi) ]
    if flag or options.show_all_parameters: table[name] = row

#end of loop over all fitted parameters

#----------
# print the results
#----------

#print details
print setUpString
print 

fmtstring = "%-40s     %15s    %15s  %10s"
highlight = "*%s*"
morelight = "!%s!"
pmsub, sigsub = None, None
if options.format == 'text':
    if options.skipFitS: print " option '--skipFitS' set true. s+b Fit is just a copy of the b-only fit"
    if options.skipFitB: print " option '--skipFitB' set true. b-only Fit is just a copy of the s+b fit"
    if options.pullDef:
        fmtstring = "%-40s       %30s    %30s  %10s"
        print fmtstring % ('name',  'b-only fit pull', 's+b fit pull', 'rho')
    elif options.absolute_values:
        fmtstring = "%-40s     %15s    %30s    %30s  %10s"
        print fmtstring % ('name', 'pre fit', 'b-only fit', 's+b fit', 'rho')
    else:
        print fmtstring % ('name', 'b-only fit', 's+b fit', 'rho')
elif options.format == 'latex':
    pmsub  = (r"(\S+) \+/- (\S+)", r"$\1 \\pm \2$")
    sigsub = ("sig", r"$\\sigma$")
    highlight = "\\textbf{%s}"
    morelight = "{{\\color{red}\\textbf{%s}}}"
    if options.skipFitS: print " option '--skipFitS' set true. $s+b$ Fit is just a copy of the $b$-only fit"
    if options.skipFitB: print " option '--skipFitB' set true. $b$-only Fit is just a copy of the $s+b$ fit"
    if options.pullDef:
        fmtstring = "%-40s & %30s & %30s & %6s \\\\"
        print "\\begin{tabular}{|l|r|r|r|} \\hline ";
        print (fmtstring % ('name', '$b$-only fit pull', '$s+b$ fit pull', r'$\rho(\theta, \mu)$')), " \\hline"
    elif options.absolute_values:
        fmtstring = "%-40s &  %15s & %30s & %30s & %6s \\\\"
        print "\\begin{tabular}{|l|r|r|r|r|} \\hline ";
        print (fmtstring % ('name', 'pre fit', '$b$-only fit', '$s+b$ fit', r'$\rho(\theta, \mu)$')), " \\hline"
    else:
        fmtstring = "%-40s &  %15s & %15s & %6s \\\\"
        print "\\begin{tabular}{|l|r|r|r|} \\hline ";
        #what = r"$(x_\text{out} - x_\text{in})/\sigma_{\text{in}}$, $\sigma_{\text{out}}/\sigma_{\text{in}}$"
        what = r"\Delta x/\sigma_{\text{in}}$, $\sigma_{\text{out}}/\sigma_{\text{in}}$"
        print  fmtstring % ('',     '$b$-only fit', '$s+b$ fit', '')
        print (fmtstring % ('name', what, what, r'$\rho(\theta, \mu)$')), " \\hline"
elif options.format == 'twiki':
    pmsub  = (r"(\S+) \+/- (\S+)", r"\1 &plusmn; \2")
    sigsub = ("sig", r"&sigma;")
    highlight = "<b>%s</b>"
    morelight = "<b style='color:red;'>%s</b>"
    if options.skipFitS: print " option '--skipFitS' set true. $s+b$ Fit is just a copy of the $b$-only fit"
    if options.skipFitB: print " option '--skipFitB' set true. $b$-only Fit is just a copy of the $s+b$ fit"
    if options.pullDef:
        fmtstring = "| <verbatim>%-40s</verbatim>  | %-30s  | %-30s   | %-15s  |"
        print "| *name* | *b-only fit pull* | *s+b fit pull* | "
    elif options.absolute_values:
        fmtstring = "| <verbatim>%-40s</verbatim>  | %-15s  | %-30s  | %-30s   | %-15s  |"
        print "| *name* | *pre fit* | *b-only fit* | *s+b fit* | "
    else:
        fmtstring = "| <verbatim>%-40s</verbatim>  | %-15s  | %-15s | %-15s  |"
        print "| *name* | *b-only fit* | *s+b fit* | *corr.* |"
elif options.format == 'html':
    pmsub  = (r"(\S+) \+/- (\S+)", r"\1 &plusmn; \2")
    sigsub = ("sig", r"&sigma;")
    highlight = "<b>%s</b>"
    morelight = "<strong>%s</strong>"
    print """
<html><head><title>Comparison of nuisances</title>
<style type="text/css">
    td, th { border-bottom: 1px solid black; padding: 1px 1em; }
    td { font-family: 'Consolas', 'Courier New', courier, monospace; }
    strong { color: red; font-weight: bolder; }
</style>
</head><body style="font-family: 'Verdana', sans-serif; font-size: 10pt;"><h1>Comparison of nuisances</h1>
<table>
"""
    if options.pullDef:
        print "<tr><th>nuisance</th><th>background fit pull </th><th>signal fit pull</th><th>correlation</th></tr>"
        fmtstring = "<tr><td><tt>%-40s</tt> </td><td> %-30s </td><td> %-30s </td><td> %-15s </td></tr>"
    elif options.absolute_values:
        print "<tr><th>nuisance</th><th>pre fit</th><th>background fit </th><th>signal fit</th><th>correlation</th></tr>"
        fmtstring = "<tr><td><tt>%-40s</tt> </td><td> %-15s </td><td> %-30s </td><td> %-30s </td><td> %-15s </td></tr>"
    else:
        what = "&Delta;x/&sigma;<sub>in</sub>, &sigma;<sub>out</sub>/&sigma;<sub>in</sub>";
        print "<tr><th>nuisance</th><th>background fit<br/>%s </th><th>signal fit<br/>%s</th><th>&rho;(&mu;, &theta;)</tr>" % (what,what)
        fmtstring = "<tr><td><tt>%-40s</tt> </td><td> %-15s </td><td> %-15s </td><td> %-15s </td></tr>"

names = table.keys()
names.sort()
highlighters = { 1:highlight, 2:morelight };
for n in names:
    v = table[n]
    if pmsub  != None: v = [ re.sub(pmsub[0],  pmsub[1],  i) for i in v ]
    if sigsub != None: v = [ re.sub(sigsub[0], sigsub[1], i) for i in v ]
    if (n,'b') in isFlagged: v[-3] = highlighters[isFlagged[(n,'b')]] % v[-3]
    if (n,'s') in isFlagged: v[-2] = highlighters[isFlagged[(n,'s')]] % v[-2]
    if options.format == "latex": n = n.replace(r"_", r"\_")
    if options.absolute_values:
       print fmtstring % (n, v[0], v[1], v[2], v[3])
    else:
       print fmtstring % (n, v[0], v[1], v[2])

if options.format == "latex":
    print " \\hline\n\end{tabular}"
elif options.format == "html":
    print "</table></body></html>"


if options.plotfile:
    import ROOT
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptFit(1)

    canvas_nuis = ROOT.TCanvas("nuisances", "nuisances", 1600, 1000)
    canvas_nuis.SetBottomMargin(0.18)
    hist_fit_e_s = hist_fit_s.Clone("errors_s")
    hist_fit_e_b = hist_fit_b.Clone("errors_b")
    #gr_fit_s = getGraph(hist_fit_s,-0.1)
    #gr_fit_b = getGraph(hist_fit_b, 0.1)
    gr_fit_s.SetLineColor(ROOT.kBlue)
    gr_fit_s.SetMarkerColor(ROOT.kBlue)
    gr_fit_b.SetLineColor(ROOT.kRed)
    gr_fit_b.SetMarkerColor(ROOT.kRed)
    gr_fit_b.SetMarkerStyle(20)
    gr_fit_s.SetMarkerStyle(20)
    gr_fit_b.SetMarkerSize(1.0)
    gr_fit_s.SetMarkerSize(1.0)
    gr_fit_b.SetLineWidth(2)
    gr_fit_s.SetLineWidth(2)
    hist_prefit.SetLineWidth(2)
    hist_prefit.SetTitle("Nuisance Paramaeters")
    hist_prefit.SetLineColor(ROOT.kBlack)
    hist_prefit.SetFillColor(ROOT.kGray)
    hist_prefit.SetMaximum(3)
    hist_prefit.SetMinimum(-3)
    hist_prefit.Draw("E2")
    hist_prefit.Draw("histsame")
    if not options.skipFitB: gr_fit_b.Draw("EPsame")
    if not options.skipFitS: gr_fit_s.Draw("EPsame")
    canvas_nuis.SetGridx()
    canvas_nuis.RedrawAxis()
    canvas_nuis.RedrawAxis('g')
    leg=ROOT.TLegend(0.7,0.8,0.89,0.89)
    leg.SetFillColor(0)
    leg.SetTextFont(42)
    leg.AddEntry(hist_prefit,"Prefit","FL")
    if not options.skipFitB:leg.AddEntry(gr_fit_b,"B-only fit","EPL")
    if not options.skipFitS:leg.AddEntry(gr_fit_s,"S+B fit"   ,"EPL")
    leg.Draw()
    canvas_nuis.Print(options.plotfile + "/" +  "postfit_pulls_m"+str(mLQ)+".png")



    canvas_pferrs = ROOT.TCanvas("post_fit_errs", "post_fit_errs", 1600, 1000)
    canvas_pferrs.SetBottomMargin(0.18)
    for b in range(1,hist_fit_e_s.GetNbinsX()+1): 
      e_prefit = max(1e-6, hist_prefit.GetBinError(b))
      hist_fit_e_s.SetBinContent(b,hist_fit_s.GetBinError(b)/e_prefit)
      hist_fit_e_b.SetBinContent(b,hist_fit_b.GetBinError(b)/e_prefit)
      hist_fit_e_s.SetBinError(b,0)
      hist_fit_e_b.SetBinError(b,0)
    hist_fit_e_s.SetFillColor(ROOT.kBlue)
    hist_fit_e_b.SetFillColor(ROOT.kRed)
    hist_fit_e_s.SetBarWidth(0.4)
    hist_fit_e_b.SetBarWidth(0.4)
    hist_fit_e_b.SetBarOffset(0.45)
    hist_fit_e_b.GetYaxis().SetTitle("#sigma_{#theta}/(#sigma_{#theta} prefit)")
    hist_fit_e_b.SetTitle("Nuisance Parameter Uncertainty Reduction")
    hist_fit_e_b.SetMaximum(1.5)
    hist_fit_e_b.SetMinimum(0)
    if(not options.skipFitB): hist_fit_e_b.Draw("bar")
    hist_fit_e_s.Draw("barsame")
    leg_rat=ROOT.TLegend(0.7,0.8,0.89,0.89)
    leg_rat.SetFillColor(0)
    leg_rat.SetTextFont(42)
    if(not options.skipFitB): leg_rat.AddEntry(hist_fit_e_b,"B-only fit","F")
    leg_rat.AddEntry(hist_fit_e_s,"S+B fit"   ,"F")
    leg_rat.Draw()
    line_one = ROOT.TLine(0,1,hist_fit_e_s.GetXaxis().GetXmax(),1)
    line_one.SetLineColor(1); line_one.SetLineStyle(2); line_one.SetLineWidth(2)
    line_one.Draw()
    canvas_pferrs.RedrawAxis()

    canvas_pferrs.Print(options.plotfile + "/" + "postfit_errs_m"+str(mLQ)+".png")

   

