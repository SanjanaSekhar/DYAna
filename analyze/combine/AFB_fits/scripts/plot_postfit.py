import ROOT
from ROOT import *
from contextlib import contextmanager
import os, pickle, subprocess, time,random
import math
from math import sqrt
from array import array
from optparse import OptionParser
import CMS_lumi, tdrstyle


gStyle.SetOptStat(0)
gROOT.SetBatch(1)



# colors

DY_c = 2011;
ttbar_c = 2022;
wt_c = 2033;
diboson_c = 2044;
qcd_c = 2055;
tautau_c = 2066;
gamgam_c = 2077;

DY_co = ROOT.TColor(DY_c, 213./255.,94./255.,0., "DY_co", 1);
ttbar_co = ROOT.TColor(ttbar_c,  123./255., 202./255., 255./255., "ttbar_co", 1);
wt_co = ROOT.TColor(wt_c,  86./255., 180./255., 233./255., "wt_co", 1);
diboson_co = ROOT.TColor(diboson_c,  0., 158./255., 115./255., "diboson_co", 1);
qcd_co = ROOT.TColor(qcd_c,  243./255., 168./255., 87./255., "qcd_co", 1);
tautau_co = ROOT.TColor(tautau_c,  19./255., 58./255., 54./255., "tautau_co", 1);
gamgam_co = ROOT.TColor(gamgam_c,  240./255., 228./255., 66./255., "gamgam_co", 1);


def makeCan(name, tag, histlist, bkglist=[],signals=[],totlist = [], colors=[],titles=[],dataName='Data',bkgNames=[],signalNames=[],
        logy=False,rootfile=False,xtitle='',ytitle='',dataOff=False,datastyle='pe',year=1, mbin = 1, ratio_range = None, NDiv = 205, prelim = False):  

    # histlist is just the generic list but if bkglist is specified (non-empty)
    # then this function will stack the backgrounds and compare against histlist as if 
    # it is data. The imporant bit is that bkglist is a list of lists. The first index
    # of bkglist corresponds to the index in histlist (the corresponding data). 
    # For example you could have:
    #   histlist = [data1, data2]
    #   bkglist = [[bkg1_1,bkg2_1],[bkg1_2,bkg2_2]]


    if len(histlist) == 1:
        width = 800
        height = 700
        padx = 1
        pady = 1
    elif len(histlist) == 2:
        width = 1200
        height = 700
        padx = 2
        pady = 1
    elif len(histlist) == 3:
        width = 1600
        height = 700
        padx = 3
        pady = 1
    elif len(histlist) == 4:
        width = 1200
        height = 1000
        padx = 2
        pady = 2
    elif len(histlist) == 6 or len(histlist) == 5:
        width = 1600
        height = 1000
        padx = 3
        pady = 2
    else:
        print 'histlist of size ' + str(len(histlist)) + ' not currently supported'
        print histlist
        return 0

    tdrstyle.setTDRStyle()

    myCan = TCanvas(name,name,width,height)
    myCan.Divide(padx,pady)

    # Just some colors that I think work well together and a bunch of empty lists for storage if needed
    default_colors = [kRed,kMagenta,kGreen,kCyan,kBlue]
    if len(colors) == 0:   
        colors = default_colors
    stacks = []
    legends = []
    legends_list = []
    mains = []
    subs = []
    pulls = []
    logString = ''
    leg_align_right = True
    CMS_align_right = False

    # For each hist/data distribution
    for hist_index, hist in enumerate(histlist):
        # Grab the pad we want to draw in
        myCan.cd(hist_index+1)
        # if len(histlist) > 1:
        thisPad = myCan.GetPrimitive(name+'_'+str(hist_index+1))
        thisPad.cd()        

        # If this is a TH2, just draw the lego
        if hist.ClassName().find('TH2') != -1:
            print 'ERROR: It seems you are trying to plot backgrounds with data on a 2D plot. This is not supported since there is no good way to view this type of distribution.'
        
        # Otherwise it's a TH1 hopefully
        else:
            titleSize = 0.09
            alpha = 1
            if dataOff:
                alpha = 0
            hist.SetLineColorAlpha(kBlack,alpha)
            if 'pe' in datastyle.lower():
                hist.SetMarkerColorAlpha(kBlack,alpha)
                hist.SetMarkerStyle(8)
            if 'hist' in datastyle.lower():
                hist.SetFillColorAlpha(0,0)

            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle(ytitle)
            
            # If there are no backgrounds, only plot the data (semilog if desired)
            if len(bkglist) == 0:
                hist.SetMaximum(1.13*hist.GetMaximum())
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                    hist.SetTitleOffset(1.1)
                hist.Draw(datastyle)
            
            # Otherwise...
            else:
                # Create some subpads, a legend, a stack, and a total bkg hist that we'll use for the error bars
                if not dataOff:
                    mains.append(TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.3, 1, 1))
                    subs.append(TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 1, 0.3))

                else:
                    mains.append(TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.1, 1, 1))
                    subs.append(TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 0, 0))

                leg_align_right = True
                CMS_align_right = False
                x_max = totlist[hist_index].GetMaximumBin()
                nbins = totlist[hist_index].GetXaxis().GetNbins()
                #if(2 *x_max > nbins):
                #    print("Found max val in bin %i, aligning legend on the left" % x_max)
                #    leg_align_right = False
                #    CMS_align_right = True
                if not logy: 
                    y_end = 0.88
                    y_size = 0.2 + 0.03*(len(bkglist[0])+len(signals))
                    x_size = 0.5
                    if(leg_align_right):
                        x_start = 0.42
                    else:
                        x_start = 0.2

                    legends.append(TLegend(x_start,y_end - y_size,x_start + x_size,y_end))
                else: 
                    legends.append(TLegend(0.2,0.11,0.45,0.2+0.02*(len(bkglist[0])+len(signals))))

                stacks.append(THStack(hist.GetName()+'_stack',hist.GetName()+'_stack'))
                legends_list.append([])


                # Set margins and make these two pads primitives of the division, thisPad
                mains[hist_index].SetBottomMargin(0.04)
                mains[hist_index].SetLeftMargin(0.2)
                mains[hist_index].SetRightMargin(0.05)
                mains[hist_index].SetTopMargin(0.08)

                subs[hist_index].SetLeftMargin(0.2)
                subs[hist_index].SetRightMargin(0.05)
                subs[hist_index].SetTopMargin(0.01)
                subs[hist_index].SetBottomMargin(0.5)
                mains[hist_index].SetFillColorAlpha(0,1)
                subs[hist_index].SetFillColorAlpha(0,1)
                mains[hist_index].SetLineColorAlpha(0,1)
                subs[hist_index].SetLineColorAlpha(0,1)
                mains[hist_index].Draw()
                subs[hist_index].Draw()

                # Build the stack
                for bkg_index,bkg in enumerate(bkglist[hist_index]):     # Won't loop if bkglist is empty
                    # bkg.Sumw2()
                    bkg.SetLineColor(kBlack)
                    if logy:
                        bkg.SetMinimum(1e-3)

                    if colors[bkg_index] != None:
                        bkg.SetFillColor(colors[bkg_index])
                        bkg.SetLineColor(colors[bkg_index])
                    else:
                        bkg.SetFillColor(default_colors[bkg_index])
                        bkg.Print()

                    stacks[hist_index].Add(bkg)
                    if bkgNames == []: this_bkg_name = bkg.GetName().split('_')[0]
                    elif type(bkgNames[0]) != list: this_bkg_name = bkgNames[bkg_index]
                    else: this_bkg_name = bkgNames[hist_index][bkg_index]
                    legends_list[hist_index].append((bkg,this_bkg_name,'f'))
                    
                # Go to main pad, set logy if needed
                mains[hist_index].cd()


                # Set y max of all hists to be the same to accomodate the tallest
                max_scaling = 3.0
                histList = [stacks[hist_index],totlist[hist_index],hist]

                yMax = histList[0].GetMaximum()
                maxHist = histList[0]
                for h in range(1,len(histList)):
                    if histList[h].GetMaximum() > yMax:
                        yMax = histList[h].GetMaximum()
                        maxHist = histList[h]
                for h in histList:
                    h.SetMaximum(yMax*max_scaling)
                    if logy == True:
                        h.SetMaximum(yMax*10)
                    else:
                        h.SetMinimum(0.)

                
                mLS = 0.07
                mTS = 0.1
                TOffset = 1.
                if(mbin > 5):
                    TOffset = 0.9
                # Now draw the main pad
                data_leg_title = hist.GetTitle()
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                hist.GetYaxis().SetTitleOffset(TOffset)
                hist.GetXaxis().SetTitleOffset(1.2)
                hist.GetYaxis().SetTitle('Events / bin')
                hist.GetYaxis().SetLabelSize(mLS)
                hist.GetYaxis().SetTitleSize(mTS)
                hist.GetYaxis().SetNdivisions(505)
                hist.GetXaxis().SetLabelOffset(999)
                hist.SetLineWidth(2)


                if logy == True:
                    hist.SetMinimum(1e-3)

                hist.SetBinErrorOption(ROOT.TH1.kPoisson)
                hist.Draw(datastyle)
                #print("Drawing %s %s \n" hist.GetName(), datastyle)

                stacks[hist_index].Draw('same hist')
                #print("Drawing %s same hist \n" stacks[hist_index].GetName())

                # Do the signals
                if len(signals) > 0: 
                    signals[hist_index].SetLineColor(kBlue)
                    signals[hist_index].SetLineWidth(2)
                    if logy == True:
                        signals[hist_index].SetMinimum(1e-3)
                    if signalNames == []: this_sig_name = signals[hist_index].GetName().split('_')[0]
                    legends_list[hist_index].append((signals[hist_index],this_sig_name,'L'))
                    signals[hist_index].Draw('hist same')

                gStyle.SetHatchesLineWidth(2)
                totlist[hist_index].SetLineColor(kWhite)
                totlist[hist_index].SetFillColor(kBlack)
                totlist[hist_index].SetFillStyle(3354)
                totlist[hist_index].SetMarkerStyle(20)
                totlist[hist_index].SetMarkerSize(0.01)

                totlist[hist_index].Draw('e2 same')

                if not dataOff:
                    legends_list[hist_index].append((hist,dataName,datastyle))
                    hist.Draw(datastyle+' same')


                #Draw helpful lines
                if(mbin <=5):
                    line_vals = [8, 16, 22]
                    text_center_bins = [4, 12, 19, 25]
                    text_strs = ["#bf{|y| #epsilon [0, 0.6]}", "#bf{|y| #epsilon [0.6, 1.0]}", "   #splitline{   #bf{|y| #epsilon}}{#bf{[1.0, 1.5]}}", 
                            "#splitline{   #bf{|y| #epsilon}}{#bf{[1.5, 2.4]}}"]
                else:
                    line_vals = [8, 16]
                    text_center_bins = [4, 12, 19]
                    text_strs = ["#bf{|y| #epsilon [0, 0.6]}", "#bf{|y| #epsilon [0.6, 1.0]}", "#splitline{   #bf{|y| #epsilon}}{#bf{[1.0, 2.4]}}"]


                lstyle = 7
                lwidth = 1
                line_eps = 0.05

                #line_max = gPad.GetY2()
                line_max = yMax * 1.5
                lines = []
                texts = []

                for idx in range(len(line_vals)):
                    line_x = line_vals[idx] + line_eps
                    l = TLine(line_x, 0, line_x, line_max)
                    l.SetLineColor(ROOT.kBlack)
                    l.SetLineStyle(lstyle)
                    l.SetLineWidth(lwidth)
                    l.Draw()
                    lines.append(l)



                #text labels
                latext = TLatex()
                latext.SetNDC();
                latext.SetTextColor(kBlack);
                latext.SetTextAlign(22); #center
                latext.SetTextFont(42);
                if(mbin <= 5):
                    latext.SetTextSize(0.050);    
                else:
                    latext.SetTextSize(0.06);    
                text_y = 0.43

                l_margin = gPad.GetLeftMargin();
                r_margin = gPad.GetRightMargin();
                nbins = float(hist.GetNbinsX())

                for idx,text_str in enumerate(text_strs):
                    text_center = l_margin + (text_center_bins[idx] / nbins) * (1.-l_margin - r_margin)
                    latext.DrawLatex(text_center, text_y, text_str)


                legends[hist_index].SetHeader(titles[0], "c")
                legends[hist_index].SetNColumns(2)
                legends[hist_index].SetTextSize(0.05)
                
                for entry in legends_list[hist_index][::-1]:
                    legends[hist_index].AddEntry(entry[0], entry[1], entry[2])

                legends[hist_index].AddEntry(totlist[hist_index], "Sys. unc.", "f")


                ratio, ratio_sys_unc = makeRatio(hist,totlist[hist_index])
                chi2 = 0.
                #for i in range(1, pull.GetNbinsX()+1):
                    #chi2 += pull.GetBinContent(i)**2;
                #print("Chi2/nbin for chan %s is %.1f/%i" % (titles[hist_index], chi2, pull.GetNbinsX()))

                legends[hist_index].AddEntry(ratio_sys_unc, "Total fit. unc.", "f")



                legends[hist_index].SetBorderSize(0)
                legends[hist_index].Draw()
                gPad.RedrawAxis()

                # Draw the pull
                subs[hist_index].cd()
                # Build the pull

                LS = mLS * 0.7/0.3
                #title size given as fraction of pad width, scale up to have same size as main pad
                YTS =  mTS * 0.7/0.3
                XTS =  mTS * 0.7/0.3
                lTOffset = TOffset * 0.27 / 0.7


                if(ratio_range == None):
                    if(mbin <= 4):
                        ratio_range = (0.86, 1.14)
                        NDiv = 205
                    elif(mbin ==5):
                        ratio_range = (0.45, 1.55)
                        NDiv = 203
                    elif(mbin ==6 or mbin == 7):
                        ratio_range = (0.1, 1.9)
                        NDiv = 303




                ratio_sys_unc.GetYaxis().SetRangeUser(ratio_range[0], ratio_range[1])
                ratio_sys_unc.GetYaxis().SetTitleOffset(lTOffset)
                ratio_sys_unc.GetYaxis().SetTickLength(0.04)
                             
                ratio_sys_unc.GetYaxis().SetLabelSize(LS)
                ratio_sys_unc.GetYaxis().SetTitleSize(YTS)
                ratio_sys_unc.GetYaxis().SetNdivisions(NDiv)
                ratio_sys_unc.GetYaxis().SetTitle("Data / fit")

                ratio_sys_unc.GetXaxis().SetRangeUser(0., hist.GetNbinsX()-.08)
                ratio_sys_unc.GetXaxis().SetTitleOffset(1.)
                ratio_sys_unc.GetXaxis().SetLabelOffset(0.05)
                ratio_sys_unc.GetXaxis().SetLabelSize(LS)
                ratio_sys_unc.GetXaxis().SetTitleSize(XTS)
                ratio_sys_unc.GetXaxis().SetTitle(xtitle)
                ratio_sys_unc.GetXaxis().SetTickLength(0.06)


                #ratio_sys_unc.SetFillColor(ROOT.kBlack)
                #ratio_sys_unc.SetFillStyle(3015)
                ratio_sys_unc.SetLineColor(ROOT.kGray)
                ratio_sys_unc.SetFillColor(ROOT.kGray)
                #ratio_sys_unc.SetFillStyle(3015)


                ratio.SetLineStyle(1)
                ratio.SetLineWidth(2)
                ratio_sys_unc.Draw("A3 same")
                ratio.Draw('p0e0Z same')
                
                line = TLine(0, 1.0, hist.GetNbinsX() - 0.08, 1.0)
                line.SetLineStyle(9)
                line.Draw()

                for idx in range(len(line_vals)):
                    line_x = line_vals[idx] + line_eps
                    l = TLine(line_x, ratio_range[0], line_x, ratio_range[1])
                    l.SetLineColor(ROOT.kBlack)
                    l.SetLineStyle(lstyle)
                    l.SetLineWidth(lwidth)
                    l.Draw()
                    lines.append(l)

                if logy == True:
                    mains[hist_index].SetLogy()

                if(prelim): 
                    print("Prelim")
                    CMS_lumi.writeExtraText = True
                else: CMS_lumi.writeExtraText = False
                if(CMS_align_right): CMS_loc = 33
                else: CMS_loc = 11
                CMS_lumi.CMS_lumi(mains[hist_index], year, CMS_loc)

    if rootfile:
        myCan.Print(tag+'/'+name+'.root','root')
    else:
        myCan.Print(tag+'/'+name+'.pdf','pdf')
        myCan.Print(tag+'/'+name+'.png','png')


def reducedCorrMatrixHist(fit_result,varsOfInterest=[]):
    ROOT.gStyle.SetOptStat(0)
    # ROOT.gStyle.SetPaintTextFormat('.3f')
    CM = fit_result.correlationMatrix()
    finalPars = fit_result.floatParsFinal()

    nParams = CM.GetNcols()
    finalParamsDict = {}
    for cm_index in range(nParams):
        if varsOfInterest == []:
            if 'Fail_' not in finalPars.at(cm_index).GetName():
                finalParamsDict[finalPars.at(cm_index).GetName()] = cm_index
        else:
            if finalPars.at(cm_index).GetName() in varsOfInterest:
                finalParamsDict[finalPars.at(cm_index).GetName()] = cm_index

    nFinalParams = len(finalParamsDict.keys())
    out = TH2D('correlation_matrix','correlation_matrix',nFinalParams,0,nFinalParams,nFinalParams,0,nFinalParams)
    out_txt = open('correlation_matrix.txt','w')

    for out_x_index, paramXName in enumerate(sorted(finalParamsDict.keys())):
        cm_index_x = finalParamsDict[paramXName]
        for out_y_index, paramYName in enumerate(sorted(finalParamsDict.keys())):
            cm_index_y = finalParamsDict[paramYName]
            if cm_index_x > cm_index_y:
                out_txt.write('%s:%s = %s\n'%(paramXName,paramYName,CM[cm_index_x][cm_index_y]))
            out.Fill(out_x_index+0.5,out_y_index+0.5,CM[cm_index_x][cm_index_y])

        out.GetXaxis().SetBinLabel(out_x_index+1,finalPars.at(cm_index_x).GetName())
        out.GetYaxis().SetBinLabel(out_x_index+1,finalPars.at(cm_index_x).GetName())
    out.SetMinimum(-1)
    out.SetMaximum(+1)

    return out

def FindCommonString(string_list):
    to_match = ''   # initialize the string we're looking for/building
    for s in string_list[0]:    # for each character in the first string
        passed = True
        for istring in range(1,len(string_list)):   # compare to_match+s against strings in string_list
            string = string_list[istring]
            if to_match not in string:                  # if in the string, add more
                passed = False
            
        if passed == True:
            to_match+=s

    if to_match[-2] == '_':
        return to_match[:-2] 
    else:
        return to_match[:-1]                # if not, return to_match minus final character

    return to_match[:-2]
        
def makeRatio( DATA,BKG):

    nbins = DATA.GetNbinsX()
    x = array('d')
    ratio = array('d')
    y_err_up = array('d')
    y_err_down = array('d')
    sys_err_up = array('d')
    sys_err_down = array('d')
    x_err_low = array('d')
    x_err_high = array('d')
    x_err_zero = array('d', [0.]* nbins)
    y_val1 = array('d', [1.]* nbins)
    
    for ibin in range(1,DATA.GetNbinsX()+1):
        DATAcont = DATA.GetBinContent(ibin)
        BKGcont = BKG.GetBinContent(ibin)
        BKG_errup = BKG.GetBinErrorUp(ibin)
        BKG_errdown = BKG.GetBinErrorLow(ibin)
        DATA_errup = DATA.GetBinErrorUp(ibin)
        DATA_errdown = DATA.GetBinErrorLow(ibin)


        x.append(DATA.GetXaxis().GetBinCenter(ibin))
        x_err_low.append(DATA.GetXaxis().GetBinCenter(ibin) - DATA.GetXaxis().GetBinLowEdge(ibin))
        x_err_high.append(-DATA.GetXaxis().GetBinCenter(ibin) + DATA.GetXaxis().GetBinUpEdge(ibin))

    

        ratio.append(DATAcont/BKGcont)
        y_err_up.append(DATA_errup / BKGcont)
        y_err_down.append(DATA_errdown / BKGcont)
        sys_err_up.append(BKG_errup/BKGcont)
        sys_err_down.append(BKG_errdown/BKGcont)
        
    pull = ROOT.TGraphAsymmErrors(nbins,x,ratio, x_err_zero, x_err_zero, y_err_down, y_err_up)
    #add extra at high-x for ratio plot
    x.append(x[-1] + x_err_high[-1])
    y_val1.append(1)
    x_err_low.append(x_err_low[-1])
    x_err_high.append(0)
    sys_err_up.append(sys_err_up[-1])
    sys_err_down.append(sys_err_down[-1])
    #add extra at low-x for ratio plot
    x.insert(0, x[0] - x_err_low[0])
    y_val1.insert(0,1)
    x_err_low.insert(0,0)
    x_err_high.insert(0, x_err_high[0])
    sys_err_up.insert(0,sys_err_up[0])
    sys_err_down.insert(0,sys_err_down[0])



    sys_unc = ROOT.TGraphAsymmErrors(nbins+2, x, y_val1, x_err_low, x_err_high, sys_err_down,  sys_err_up)
    return pull, sys_unc

def Make_up_down(hist):
    hist_up = hist.Clone(hist.GetName()+'_up')
    hist_down = hist.Clone(hist.GetName()+'_down')

    for xbin in range(1,hist.GetNbinsX()+1):
        errup = hist.GetBinErrorUp(xbin)
        errdown = hist.GetBinErrorLow(xbin)
        nom = hist.GetBinContent(xbin)

        hist_up.SetBinContent(xbin,nom+errup)
        hist_down.SetBinContent(xbin,nom-errdown)

    return hist_up,hist_down

if (__name__ == "__main__"):
    parser = OptionParser()
    parser.add_option("--input", "-i", default = "", help="Input file")
    parser.add_option("--output", "-o", default = "", help="Input directory")
    parser.add_option("--mbin", "-m", type = 'int', default = 0, help="Mass bin (for plot label)")
    parser.add_option("--year", "-y", type = 'int', default = -1, help="Year (-1 for all) ")
    parser.add_option("--ss",   default = False, action='store_true',  help="Fit was done with ee_ss region too")
    (options, args) = parser.parse_args()


#fin_ = "combined_fit_shapes_mbin1.root"
#odir = "postfit_plots/combined_fit_mbin1"
#mbin = 1
    if(options.year < 0):
        years = [2016, 2017, 2018]
    else:
        years = [options.year]
    h_names = ["gam", "db", "qcd", "top",  "dy"]
    h_ss_names = ["bk", "dy", "qcd"]


    m_bins = [150, 170, 200,  250, 320, 510, 700, 1000, 14000]


    label_color_map = dict()
    label_color_map['dy'] = ("DY Signal", DY_c)
    label_color_map['top'] = ("t#bar{t} + Single Top", ttbar_c)
    label_color_map['db'] = ("WW + WZ + ZZ",  diboson_c)
    #label_color_map['tautau'] = ("DY #tau#tau Bkg.", tautau_c)
    label_color_map['gam'] = ("\\gamma\\gamma \\rightarrow {\\ell}{\\ell} ", gamgam_c)
    label_color_map['qcd'] = ("WJets + QCD", qcd_c)

    datastyle = "pe0x0"

    fracs = dict()
    for name in h_names:
        fracs[name] = 0.

    dirs = ["Y%i_mumu%i_postfit/", "Y%i_ee%i_postfit/"]
    if(options.ss): dirs = dirs = ["Y%i_mumu%i_postfit/", "Y%i_ee%i_postfit/", "Y%i_ee%i_ss_postfit/"]
    f_in = TFile.Open(options.input)
    for year in years:
        for idx, dir_name in enumerate(dirs):
            dir_ = dir_name % (year % 2000, year % 2000)
            h_tot = f_in.Get(dir_ + "TotalProcs")
            h_tot = h_tot.Clone("h_tot_c%i_y%i" %(idx, year))
            h_data_orig = f_in.Get(dir_ + "data_obs")
            #copy to new hist so poisson error bars work
            h_data = h_data_orig.Clone("h_data_c%i_y%i" %(idx, year))
            h_data.Reset()
            for b in range(h_data_orig.GetXaxis().GetNbins() + 1):
                h_data.SetBinContent(b, h_data_orig.GetBinContent(b))

            h_tot_sig = f_in.Get(dir_ + "TotalSig")
            h_tot_sig = h_tot_sig.Clone("h_tot_sig_c%i_y%i" %(idx, year))

            mbin_low = m_bins[options.mbin]
            mbin_high = m_bins[options.mbin+1]

            if(idx == 0): title = "Muons %i-%i GeV" % (mbin_low, mbin_high)
            if(idx == 1): title = "Electrons %i-%i GeV" % (mbin_low, mbin_high)
            if(idx == 2): title = "Electrons Samesign %i-%i GeV" % (mbin_low, mbin_high)
            
            if(idx == 2): name_list = h_ss_names
            else: name_list = h_names
            hist_list = []
            color_list = []
            label_list = []

            for name in name_list:
                if(name == "dy"):
                    h = h_tot_sig.Clone("h_%s_c%i_y%i" %(name, idx, year))
                else:
                    h = f_in.Get(dir_ + name)
                    if(h != None):
                        h = h.Clone("h_%s_c%i_y%i" %(name, idx, year))
                if(h != None):
                    #h.Print()
                    hist_list.append(h)
                    label_list.append(label_color_map[name][0])
                    color_list.append(label_color_map[name][1])

                    if("gam" in name):
                        this_frac = h.Integral()/(h_tot_sig.Integral() + h.Integral())
                        print("Chan %i Year %i Name %s frac %.3f \n" % (idx, year, name, this_frac))
                        fracs[name] += this_frac

            makeCan(dir_[:-1], options.output, [h_data], bkglist=[hist_list], totlist=[h_tot], colors = color_list, bkgNames = label_list, 
                    titles = [title], xtitle = "Template bin", year = year, datastyle=datastyle, mbin = options.mbin ) 

    for key in fracs.keys():
        fracs[key] /= (len(years)*len(dirs))

    key_ = "gam"
    print("Average fraction for %s  is  %.3f \n" % (key_, fracs[key_]))


