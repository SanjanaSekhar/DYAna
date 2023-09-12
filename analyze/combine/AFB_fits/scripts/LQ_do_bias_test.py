
import ROOT
from ROOT import *
from LQ_utils import *
import numpy as np
from itertools import *
from array import array
gROOT.SetBatch(True)

gStyle.SetOptFit(1) 




parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--yLQ",  default=0.0, type='float', help="yLQ value to inject")
parser.add_option("--chan",  default="ee", help="channel ee or mumu ")
parser.add_option("--nToys",  default=100, type='int', help="How many toys to run")
parser.add_option("--q",  default="u", help=" channel u,d,c,s ")
parser.add_option("-o", "--odir", default="signal_injection/condor/", help = "output directory")
parser.add_option("--reuse_fit", default=False, action="store_true", help="Reuse initial fit from previous run to save time")
parser.add_option("--prefit", default=False, action="store_true", help="Sample toys from prefit uncs")
parser.add_option("--no_sys",  default=False, action="store_true", help="Use fit template without any shape systematics")
parser.add_option("--freeze",  default="", help="freezeNuisanceGroups")
parser.add_option("--freezeGroups",  default="", help="freezeNuisanceGroups")
parser.add_option("--mLQ",  default=2000, type='int', help="mLQ")
parser.add_option("--is_vec", default=False, action="store_true", help="is vec?")
parser.add_option("--ending",  default="102022", help="png ext")
parser.add_option("--job",  default=1, type='int',help="job index")
parser.add_option("--plot",  default=False, help="make residuals and pulls")
(options, args) = parser.parse_args()



#options.nToys = 200
is_vec = options.is_vec
mLQ = options.mLQ
gen_level = False
no_LQ = False
fake_data = True
no_sys = options.no_sys
year = -1
ending = options.ending
if is_vec: ending = "vec_"+ending
if options.no_sys: ending  = "nosys_"+ending
if options.freeze: ending = "freeze_%s_"%((options.freeze).replace(",","")) + ending
if options.freezeGroups: ending = "freeze_%s_"%((options.freezeGroups).replace(",","")) + ending
#if options.chan=="mumu" and options.q=="d": is_vec = True
yLQ2 = options.yLQ**2
yLQ = options.yLQ
print(options.chan,options.q)
sys = ["xsecs", "pdfs","elRECOs","elIDs"]
if not options.plot:
    workspace = "workspaces/%s_%s_%i_fit_bias_tests.root" % (options.chan, options.q, options.job)
    make_workspace(workspace, gen_level, options.chan, options.q, is_vec, no_LQ , no_sys, fake_data, mLQ, year, True, False)

    #extra_params = " --freezeParameters A0,A4 "
    extra_params = ""
    if options.freeze!="": extra_params = " --freezeParameters %s " %options.freeze
    if options.freezeGroups!="": extra_params = " --freezeNuisanceGroups %s " %options.freezeGroups
    if options.no_sys: extra_params = " --freezeParameters allConstrainedNuisances"
    
    print("Will inject options.yLQ %.2f for options.chan %s%s for all toys " %(options.yLQ,options.chan,options.q))

    res_yLQ2 = []
    pull_yLQ2 = []

    seed = i = 3457 * (options.job+1)
    if(not options.prefit):
        print("Sampling toys based on postfit")
        if(not options.reuse_fit):
            print_and_do("combine -M MultiDimFit -d %s --saveFit --saveWorkspace --robustFit 1 %s -s %i   -n _%i" % (workspace, extra_params,seed,seed))
    nToys_generated = 0 
    #for i in range(options.nToys):
    while nToys_generated < options.nToys:
        i += 1
	if(not options.prefit):
            fitted_yLQ2 = setSnapshot(yLQ2_val = -1., mdf = True, s = seed)
            if(options.no_sys):
                print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root -s %i  --snapshotName initialFit --toysNoSystematics --bypassFrequentistFit --saveToys -t 1  --setParameters yLQ2=%f,A4=1.6,A0=0.05 " % (i,yLQ2))

            else:
                print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root -s %i --snapshotName initialFit --toysFrequentist --bypassFrequentistFit --saveToys -t 1  --setParameters yLQ2=%f,A4=1.6,A0=0.05 " % (i,yLQ2))
        else:

            print_and_do("combine -M GenerateOnly -d %s -s %i  --saveToys -t 1 --toysFrequentist --setParameters yLQ2=%f,A4=1.6,A0=0.05 "% (workspace, i, yLQ2))

        #print_and_do("combine -M MultiDimFit -d %s --saveWorkspace --saveFitResult --toysFile higgsCombineTest.GenerateOnly.mH120.%i.root --toysFrequentist  -t 1 --robustFit 1 --forceRecreateNLL %s -n _%i --freezeParameters A4,A0 --robustHesse=1" %(workspace,  i, extra_params, i))
        if no_sys: print_and_do("combine -M MultiDimFit -d %s --saveWorkspace --saveFitResult -t 1 --toysNoSystematics --toysFile higgsCombineTest.GenerateOnly.mH120.%i.root   --robustFit 1  --forceRecreateNLL %s -n _%i " %(workspace,  i, extra_params, i))
        else: print_and_do("combine -M MultiDimFit -d %s --saveWorkspace --saveFitResult -t 1 --toysFrequentist --toysFile higgsCombineTest.GenerateOnly.mH120.%i.root   --robustFit 1  --forceRecreateNLL %s -n _%i " %(workspace,  i, extra_params, i))
        f_fit = TFile.Open("multidimfit_%i.root"%i)
        if f_fit:
            print_and_do("""echo "fit_mdf->Print();" > cmd.txt""")
            print_and_do(""" cat cmd.txt """)
            print_and_do("""echo ".q" >> cmd.txt """)
            #print_and_do("root -l -b multidimfit.root < cmd.txt > fit_results/%s_m%i.txt" % (fit_name,mLQ))
            print_and_do("root -l -b multidimfit_%i.root < cmd.txt"%i)
            fr = f_fit.Get('fit_mdf')
            myargs = RooArgSet(fr.floatParsFinal())

            yLQ2_fit = myargs.find("yLQ2").getVal()
            yLQ2_err_hi = myargs.find("yLQ2").getErrorHi()
            yLQ2_err_lo = myargs.find("yLQ2").getErrorLo()
            #A0 = myargs.find("A0").getVal()
            #A0_err = myargs.find("A0").getError()
            if yLQ2_err_hi == 0. or yLQ2_err_lo == 0. :
                print("FIT FAILED, SKIPPING TOY ")
                continue
            print("yLQ2 %.3f err %.3f %.3f" % (yLQ2_fit, yLQ2_err_hi, yLQ2_err_lo))
            res_yLQ2.append(yLQ2_fit - yLQ2) 
            nToys_generated += 1
            if (yLQ2_fit - yLQ2) < 0 and yLQ2_err_hi > 0.: pull_yLQ2.append((yLQ2_fit-yLQ2)/ yLQ2_err_hi)
            elif (yLQ2_fit - yLQ2) > 0 and yLQ2_err_lo < 0.: pull_yLQ2.append((yLQ2_fit-yLQ2)/abs(yLQ2_err_lo))
            

            #h_res_afb.Fill(Afb - options.Afb)
            #h_res_a0.Fill(A0 - options.A0)
            #if(Afb_err > 0.): h_pull_afb.Fill((Afb- options.Afb)/ Afb_err)
            #if(A0_err > 0.):  h_pull_a0.Fill((A0- options.A0)/ A0_err)
            f_fit.Close()
            #print_and_do("rm higgsCombineTest.GenerateOnly.mH120.%i.root" % i)
        else:
            print("Can't open multidimfit. Looks like the fit failed.")



    #f_ws = TFile.Open("higgsCombineTest.FitDiagnostics.mH120.123456.root")
    #f_toys = TFile.Open("higgsCombineTest.GenerateOnly.mH120.123456.root")
    #ws = f_ws.Get('w')
    #toy = f_toys.Get('toys/toy_1')
    #getattr(ws, 'import')(toy)
    #ws.writeToFile("toy_ws.root")
    #print_and_do("PostFitShapesFromWorkspace -w toy_ws.root --dataset model_sData  -f fitDiagnostics.root:fit_s -o toy_shapes.root --sampling --samples 100")
    ##print_and_do("python scripts/plot_postfit.py -i toy_ws.root -o test/ -m %i" % (mbin))
    with open('%srespull_%s_%s_%i_yLQ%.1f_%s.txt'%(options.odir,options.chan,options.q,options.job,options.yLQ,ending), 'w') as f:
        for res,pull in zip(res_yLQ2, pull_yLQ2):
            f.write("%f %f\n" %(res,pull))

else:
    
    c4 = TCanvas("c4", "", 900, 900)
    leg = TLegend()
    leg.SetHeader("Channel","C")
    l1 = TLine(0,-1,9500,-1)
    l2 = TLine(0,-0.5,9500,-0.5)
    l3 = TLine(0,0,9500,0)
    l4 = TLine(0,0.5,9500,0.5)
    l5 = TLine(0,1,9500,1)
    l1.SetLineStyle(9)
    l2.SetLineStyle(9)
    l3.SetLineStyle(9)
    l4.SetLineStyle(9)
    l5.SetLineStyle(9)

    for options.chan in ["ee","mumu"]:
        for options.q in ["u","d"]:

	    pull_mean, pull_sigma, m_list = array("d"),array("d"),array("d")
            
	    for options.mLQ in [1000,2500,3500,5000]:
                for options.yLQ in [0.0]:

                    if (options.chan == "ee" and options.q == "u") or (options.chan == "mumu" and options.q == "d"): is_vec = True
                    else: is_vec = False


                    respull = []
		    #pull_mean, pull_sigma, m_list, q_list, chan_list, yLQ_list = [],[],[],[],[],[]
                    for job_idx in range(0,50):
                        if options.yLQ == 0.25: print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/sasekhar/Condor_outputs/bias_test_yLQ%.2f_%s_%s%s_m%s_no%s_%i_%s/respull_%s_%s_%i_yLQ%.1f%s_%s.txt %s%s" % (options.yLQ, options.chan, options.q, ("_vec" if is_vec else ""), options.mLQ, (options.freezeGroups).replace(",",""),job_idx, ending[-6:],options.chan,options.q,job_idx,options.yLQ,("_vec" if is_vec else ""),ending, options.odir,options.mLQ))
                        else: print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/sasekhar/Condor_outputs/bias_test_yLQ%.1f_%s_%s%s_m%s_no%s_%i_%s/respull_%s_%s_%i_yLQ%.1f%s_%s.txt %s%s" % (options.yLQ, options.chan, options.q,("_vec" if is_vec else ""), options.mLQ, (options.freezeGroups).replace(",",""),job_idx, ending[-6:],options.chan,options.q,job_idx,options.yLQ,("_vec" if is_vec else ""), ending, options.odir,options.mLQ))
                        filename = '%s%s/respull_%s_%s_%i_yLQ%.1f%s_%s.txt'%(options.odir,options.mLQ,options.chan,options.q,job_idx,options.yLQ,("_vec" if is_vec else ""), ending)
                        
                        if os.path.isfile(filename):
                            with open(filename, 'r') as f:
                                for line in f.readlines():
                        	    respull.append(line.split(' '))

                    respull = np.asarray(respull, dtype=float)
                    res_yLQ2 = respull[:,0]
                    pull_yLQ2 = respull[:,1]
                    #pull_yLQ2[res_yLQ2<0] = -1.*pull_yLQ2[res_yLQ2<0]
                    res_yLQ2 = res_yLQ2.tolist()
                    pull_yLQ2 = pull_yLQ2.tolist()
                    #print("No. of toys in residuals: ", len(res_yLQ2))
		    #print(res_yLQ2)
                    #print("No. of toys in pulls: ", len(pull_yLQ2))
		    #print(pull_yLQ2)
                    if len(res_yLQ2) != len(pull_yLQ2): print("REDO TESTS for channel %s %s - yLQ = %.2f - mLQ = %s ",options.chan, options.q, options.yLQ, options.mLQ)

                    n_bins = 25
                    h_pull_yLQ2 = TH1F("h_pull_yLQ2", "", n_bins, -4, 4)

                    res_yLQ2_range = max(3*np.std(res_yLQ2), 0.15)

                    h_res_yLQ2 = TH1F("h_res_yLQ2", "", n_bins, -res_yLQ2_range, res_yLQ2_range)



                    def fill_h(arr, h):
                        #print(arr)
                        for arg in arr:
                            h.Fill(arg)

                    fill_h(pull_yLQ2, h_pull_yLQ2)
                    #fill_h(pull_a0, h_pull_a0)
                    fill_h(res_yLQ2, h_res_yLQ2)
                    #fill_h(res_a0, h_res_a0)

                    if options.chan=="mumu": chan_label = "#mu"
                    else : chan_label = "e"

                    c1 = TCanvas("c1", "", 900, 900)
                    c1.cd()
		    h_pull_yLQ2.Fit("gaus")
                    fit_yLQ2= h_pull_yLQ2.GetFunction("gaus")
                    if(fit_yLQ2): fit_yLQ2.SetLineColor(kBlue)
                    h_pull_yLQ2.Draw()
		    pull_mean.append(fit_yLQ2.GetParameter(1))
		    pull_sigma.append(fit_yLQ2.GetParameter(2))
                    if is_vec: 
                        h_pull_yLQ2.SetTitle("Signal Inject Test : Inject g_{%s %s} = %.1f (M_{LQ} = %.1f TeV); freeze: %s" % (chan_label,options.q,options.yLQ,options.mLQ/1000.,options.freezeGroups))
                        h_pull_yLQ2.GetXaxis().SetTitle("Pull g_{%s %s}^2"%(chan_label,options.q))
                        c1.Print("%s%s/bias_test_pull_yLQ%.1f_%s_%s_vec_%s.png" %(options.odir, options.mLQ,options.yLQ, options.chan, options.q,ending))
                    else: 
                        h_pull_yLQ2.SetTitle("Signal Inject Test : Inject y_{%s %s} = %.1f (M_{LQ} = %.1f TeV); freeze %s" % (chan_label,options.q,options.yLQ,options.mLQ/1000.,options.freezeGroups))
                        h_pull_yLQ2.GetXaxis().SetTitle("Pull y_{%s %s}^2"%(chan_label,options.q))
                        c1.Print("%s%s/bias_test_pull_yLQ%.1f_%s_%s_%s.png" %(options.odir, options.mLQ,options.yLQ, options.chan, options.q,ending))


                    # c2 = TCanvas("c1", "", 900, 900)
                    # h_pull_a0.Fit("gaus")
                    # fit_a0= h_pull_a0.GetFunction("gaus")
                    # if(fit_a0): fit_a0.SetLineColor(kBlue)
                    # h_pull_a0.Draw()
                    # h_pull_a0.SetTitle("Signal Inject Test Mass bin %i, Inject A0 %.2f" % (options.mbin, options.A0))
                    # h_pull_a0.GetXaxis().SetTitle("Pull A0")
                    # c2.Print("%sbias_test_mbin%i_Az%.0f.png" %(options.odir, options.mbin, 100.* options.A0))


                    c3 = TCanvas("c3", "", 900, 900)
                    c3.cd()
		    h_res_yLQ2.Fit("gaus")
                    fit_yLQ2= h_res_yLQ2.GetFunction("gaus")
                    if(fit_yLQ2): fit_yLQ2.SetLineColor(kBlue)
                    h_res_yLQ2.Draw()
                    if is_vec: 
                        h_res_yLQ2.SetTitle("Signal Inject Test : Inject g_{%s %s} = %.1f (M_{LQ} = %.1f TeV); freeze: %s" % (chan_label,options.q,options.yLQ,options.mLQ/1000.,options.freezeGroups))
                        h_res_yLQ2.GetXaxis().SetTitle("#Delta g_{%s %s}^2"%(chan_label,options.q))
                        c3.Print("%s%s/bias_test_res_yLQ%.1f_%s_%s_vec_%s.png" %(options.odir, options.mLQ, options.yLQ, options.chan, options.q,ending))
                    else: 
                        h_res_yLQ2.SetTitle("Signal Inject Test : Inject y_{%s %s} = %.1f (M_{LQ} = %.1f TeV); freeze: %s" % (chan_label,options.q,options.yLQ,options.mLQ/1000.,options.freezeGroups))
                        h_res_yLQ2.GetXaxis().SetTitle("#Delta y_{%s %s}^2"%(chan_label,options.q))
                        c3.Print("%s%s/bias_test_res_yLQ%.1f_%s_%s_%s.png" %(options.odir,options.mLQ, options.yLQ, options.chan, options.q,ending))



                    # c4 = TCanvas("c1", "", 900, 900)
                    # h_res_a0.Fit("gaus")
                    # fit_a0= h_res_a0.GetFunction("gaus")
                    # if(fit_a0): fit_a0.SetLineColor(kBlue)
                    # h_res_a0.Draw()
                    # h_res_a0.SetTitle("Signal Inject Test Mass bin %i, Inject A0 %.2f" % (options.mbin, options.A0))
                    # h_res_a0.GetXaxis().SetTitle("#Delta A0")
                    # c4.Print("%sbias_test_res_mbin%i_Az%.0f.png" %(options.odir, options.mbin, 100.* options.A0))
		    m_list.append(options.mLQ)
		    #chan_list.append(options.chan)
		    #q_list.append(options.q)
		    #yLQ_list.append(options.yLQ)
	
	    c4.cd()
	    #draw the list of pulls for 1 channel
            #print(m_list,pull_mean,[0,0,0,0],pull_sigma)
	    gr = TGraphErrors(4, m_list, pull_mean, array("d",[0,0,0,0]), pull_sigma)
            gr.SetTitle("Pulls: Injection y_{LQ} (g_{LQ}) = %.2f"%options.yLQ)
	    gr.SetName('gr')
	    leg.AddEntry('gr', chan_label+"-"+options.q+"%s"%("-vec" if is_vec else ""), 'lep')
	    gr.Draw("AP same")
	    #del gr

    leg.Draw("same")
    l1.Draw("same")
    l2.Draw("same")
    l3.Draw("same")
    l4.Draw("same")
    l5.Draw("same")
    c4.Print("bias_summary_yLQ%.2f_%s.png"%(options.yLQ,ending))
    c4.Close()
	
