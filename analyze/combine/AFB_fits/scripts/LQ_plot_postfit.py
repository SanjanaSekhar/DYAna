import ROOT
from ROOT import *
from contextlib import contextmanager
import os, pickle, subprocess, time,random
import math
from math import sqrt
import array
from optparse import OptionParser
import CMS_lumi, tdrstyle

gStyle.SetOptStat(0)
gROOT.SetBatch(1)

def makeCan(name, tag, histlist, bkglist=[],signals=[],totlist = [], colors=[],titles=[],dataName='data',bkgNames=[],signalNames=[],
        logy=False,rootfile=False,xtitle='',ytitle='',dataOff=False,datastyle='pe',year=1):  
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
            if logy == True:
                gPad.SetLogy()
            gPad.SetLeftMargin(0.2)
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle(ytitle)
            hist.GetXaxis().SetTitleOffset(1.5)
            hist.GetYaxis().SetTitleOffset(1.0)
            hist.GetZaxis().SetTitleOffset(1.8)
            if len(titles) > 0:
                hist.SetTitle(titles[hist_index])

            hist.Draw('lego')
            if len(bkglist) > 0:
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
                    y_size = 0.2 + 0.02*(len(bkglist[0])+len(signals))
                    x_size = 0.35
                    if(leg_align_right):
                        x_start = 0.55
                    else:
                        x_start = 0.2

                    legends.append(TLegend(x_start,y_end - y_size,x_start + x_size,y_end))
                else: 
                    legends.append(TLegend(0.2,0.11,0.45,0.2+0.02*(len(bkglist[0])+len(signals))))

                stacks.append(THStack(hist.GetName()+'_stack',hist.GetName()+'_stack'))
                legends_list.append([])


                # Set margins and make these two pads primitives of the division, thisPad
                mains[hist_index].SetBottomMargin(0.04)
                mains[hist_index].SetLeftMargin(0.17)
                mains[hist_index].SetRightMargin(0.05)
                mains[hist_index].SetTopMargin(0.08)

                subs[hist_index].SetLeftMargin(0.17)
                subs[hist_index].SetRightMargin(0.05)
                subs[hist_index].SetTopMargin(0)
                subs[hist_index].SetBottomMargin(0.35)
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
                    else:
                        bkg.SetFillColor(default_colors[bkg_index])

                    stacks[hist_index].Add(bkg)
                    if bkgNames == []: this_bkg_name = bkg.GetName().split('_')[0]
                    elif type(bkgNames[0]) != list: this_bkg_name = bkgNames[bkg_index]
                    else: this_bkg_name = bkgNames[hist_index][bkg_index]
                    legends_list[hist_index].append((bkg,this_bkg_name,'f'))
                    
                # Go to main pad, set logy if needed
                mains[hist_index].cd()


                # Set y max of all hists to be the same to accomodate the tallest
                max_scaling = 2.0
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

                
                mLS = 0.08
                mTS = 0.06
                # Now draw the main pad
                data_leg_title = hist.GetTitle()
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                hist.GetYaxis().SetTitleOffset(1.4)
                hist.GetXaxis().SetTitleOffset(1.2)
                hist.GetYaxis().SetTitle('Events / bin')
                hist.GetYaxis().SetLabelSize(mLS)
                hist.GetYaxis().SetTitleSize(mTS)
                hist.GetXaxis().SetLabelOffset(999)
                if logy == True:
                    hist.SetMinimum(1e-3)
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

                totlist[hist_index].SetFillColor(kBlack)
                totlist[hist_index].SetFillStyle(3354)
                totlist[hist_index].SetMarkerStyle(20)
                totlist[hist_index].SetMarkerSize(0.01)

                totlist[hist_index].Draw('e2 same')

                if not dataOff:
                    legends_list[hist_index].append((hist,dataName,datastyle))
                    hist.Draw(datastyle+' same')


            
                legends[hist_index].SetHeader(title, "c")
                legends[hist_index].SetNColumns(2)
                
                for entry in legends_list[hist_index][::-1]:
                    legends[hist_index].AddEntry(entry[0], entry[1], entry[2])



                legends[hist_index].SetBorderSize(0)
                legends[hist_index].Draw()
                gPad.RedrawAxis()

                # Draw the pull
                subs[hist_index].cd()
                # Build the pull
                pulls.append(Make_Pull_plot(hist,totlist[hist_index]))
                pulls[hist_index].SetFillColor(kGray)
                pulls[hist_index].SetTitle(";"+hist.GetXaxis().GetTitle()+";(Data-Bkg)/#sigma")
                pulls[hist_index].SetStats(0)
                chi2 = 0.
                for i in range(1, pulls[hist_index].GetNbinsX()+1):
                    chi2 += pulls[hist_index].GetBinContent(i)**2;
                print("Chi2/nbin for chan %s is %.1f/%i" % (titles[hist_index], chi2, pulls[hist_index].GetNbinsX()))

                LS = .13
                #title size given as fraction of pad width, scale up to have same size as main pad
                TS =  0.06 * 0.7/0.3

                pulls[hist_index].GetYaxis().SetRangeUser(-2.9,2.9)
                pulls[hist_index].GetYaxis().SetTitleOffset(0.3)
                pulls[hist_index].GetXaxis().SetTitleOffset(1.2)
                             
                pulls[hist_index].GetYaxis().SetLabelSize(LS)
                pulls[hist_index].GetYaxis().SetTitleSize(TS)
                pulls[hist_index].GetYaxis().SetNdivisions(306)
                pulls[hist_index].GetYaxis().SetTitle("(Data-Fit)/#sigma")

                pulls[hist_index].GetXaxis().SetLabelOffset(0.05)
                pulls[hist_index].GetXaxis().SetLabelSize(LS)
                pulls[hist_index].GetXaxis().SetTitleSize(TS)
                pulls[hist_index].GetXaxis().SetTitle(xtitle)


                pulls[hist_index].Draw('hist')

                if logy == True:
                    mains[hist_index].SetLogy()

                if(CMS_align_right): CMS_loc = 33
                else: CMS_loc = 11
                CMS_lumi.CMS_lumi(thisPad, year, CMS_loc)

    if rootfile:
        myCan.Print(tag+'/'+name+'.root','root')
    else:
        myCan.Print(tag+'/'+name+'_m'+str(mLQ)+'.png','png')


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
        
def Make_Pull_plot( DATA,BKG):
    BKGUP, BKGDOWN = Make_up_down(BKG)
    pull = DATA.Clone(DATA.GetName()+"_pull")
    pull.Add(BKG,-1)
    sigma = 0.0
    FScont = 0.0
    BKGcont = 0.0
    for ibin in range(1,pull.GetNbinsX()+1):
        FScont = DATA.GetBinContent(ibin)
        BKGcont = BKG.GetBinContent(ibin)
        
        if FScont>=BKGcont:
            FSerr = DATA.GetBinErrorLow(ibin)
            BKGerr = abs(BKGUP.GetBinContent(ibin)-BKG.GetBinContent(ibin))
        if FScont<BKGcont:
            FSerr = DATA.GetBinErrorUp(ibin)
            BKGerr = abs(BKGDOWN.GetBinContent(ibin)-BKG.GetBinContent(ibin))
        if FSerr != None:
            sigma = sqrt(FSerr*FSerr + BKGerr*BKGerr)
        else:
            sigma = sqrt(BKGerr*BKGerr)
        if FScont == 0.0:
            pull.SetBinContent(ibin, 0.0 )  
        else:
            if sigma != 0 :
                pullcont = (pull.GetBinContent(ibin))/sigma
                pull.SetBinContent(ibin, pullcont)
            else :
                pull.SetBinContent(ibin, 0.0 )
    return pull

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


parser = OptionParser()
parser.add_option("--input", "-i", default = "", help="Input file")
parser.add_option("--output", "-o", default = "", help="Input directory")
parser.add_option("--mLQ", "-m", type = 'int', default = 0, help="mLQ (for plot label)")
parser.add_option("--q", "-q", type = 'string', default = "", help="q = u or d")
parser.add_option("--chan", "-c", type = 'string', default = "", help="ee or mumu")
parser.add_option("--year", "-y", type = 'int', default = -1, help="Year (-1 for all) ")
parser.add_option("--ss",   default = False, action='store_true',  help="Fit was done with ee_ss region too")
parser.add_option("--gen_level", default = False, action='store_true', help="generator level fits")
(options, args) = parser.parse_args()

mLQ = options.mLQ
#fin_ = "combined_fit_shapes_mbin1.root"
#odir = "postfit_plots/combined_fit_mbin1"
#mbin = 1
if(options.year < 0):
    years = [2016, 2017, 2018]
else:
    years = [options.year]
if options.q == "u":
    h_names = ["gam", "LQint_u", "LQpure_u" ,"qcd", "top", "db", "tautau", "dy"]
elif options.q == "d":
    h_names = ["gam", "LQint_d", "LQpure_d" ,"qcd", "top", "db", "tautau", "dy"]
if options.gen_level and options.q == "u":
    h_names = ["LQ", "dy"]
#if options.gen_level and options.q == "d":
#    h_names = ["LQint_d", "LQpure_d", "dy"]
#h_ss_names = ["bk", "dy", "qcd"]


#m_bins = [150, 171, 200,  250, 320, 510, 700, 1000, 14000]


label_color_map = dict()
label_color_map['dy'] = ("DY", kRed + 1)
label_color_map['top'] = ("t#bar{t} + tW ", kBlue)
label_color_map['db'] = ("WW + WZ + ZZ",  kGreen +3)
label_color_map['tautau'] = ("DY #tau#tau Bkg.", kMagenta + 4)
label_color_map['gam'] = ("\\gamma\\gamma \\rightarrow \\ell\\ell ", kOrange)
label_color_map['qcd'] = ("WJets + QCD", kRed - 7)
label_color_map['LQ'] = ("LQpure+int", kBlue)
label_color_map['alpha'] = ("alpha", kGreen)
label_color_map['fmn'] = ("DY minus", kBlue)
label_color_map['fpl'] = ("DY plus", kRed)
if options.q=="u":
    label_color_map['LQint_u'] = ("LQint_u", kGreen)
    label_color_map['LQpure_u'] = ("LQpure_u", kBlue)
if options.q=="d":
    label_color_map['LQint_d'] = ("LQint_d", kGreen)
    label_color_map['LQpure_d'] = ("LQpure_d", kBlue)


datastyle = "pe0x0"

fracs = dict()
for name in h_names:
    fracs[name] = 0.

dirs = ["Y%i_postfit/"]
#if(options.ss): dirs =  ["Y%i_mumu%i_postfit/", "Y%i_ee%i_postfit/", "Y%i_ee%i_ss_postfit/"]
print ("\n file from where the h_tot is being pulled= ",options.input)
f_in = TFile.Open(options.input)
for year in years:
    for idx, dir_name in enumerate(dirs):
        dir_ = dir_name % ( year % 2000)
        print ("\n dir_ = ", dir_)
        h_tot = f_in.Get(dir_ + "TotalProcs")
        h_tot = h_tot.Clone("h_tot_c%i_y%i" %(idx, year))
        h_data = f_in.Get(dir_ + "data_obs")
        h_data = h_data.Clone("h_data_c%i_y%i" %(idx, year))

        h_tot_sig = f_in.Get(dir_ + "TotalSig")
        h_tot_sig = h_tot_sig.Clone("h_tot_sig_c%i_y%i" %(idx, year))

        #mbin_low = m_bins[options.mbin]
        #mbin_high = m_bins[options.mbin+1]

        if(options.chan=="mumu"): title = "Muons %i GeV mLQ = %i" % (year, mLQ)
        if(options.chan=="ee"): title = "Electrons %i  GeV mLQ = %i" % (year,mLQ)
        #if(idx == 2): title = "Electrons Samesign %i  GeV" % (year)
        
        #if(idx == 2): name_list = h_ss_names
        #else: 
        name_list = h_names
        hist_list = []
        color_list = []
        label_list = []

        if not options.gen_level:
            for name in name_list:
                if(name == "dy"):
                    h = h_tot_sig.Clone("h_%s_c%i_y%i" %(name, idx, year))
                else:
                    h = f_in.Get(dir_ + name)
                    if(h != None):
                        h = h.Clone("h_%s_c%i_y%i" %(name, idx, year))
                if(h != None):
                    h.Print()
                    hist_list.append(h)
                    label_list.append(label_color_map[name][0])
                    color_list.append(label_color_map[name][1])

                    if("gam" in name):
                        this_frac = h.Integral()/(h_tot_sig.Integral() + h.Integral())
                        print("Chan %i Year %i Name %s frac %.3f \n" % (idx, year, name, this_frac))
                        fracs[name] += this_frac
        
        else:
            
            h_sig = h_tot_sig.Clone("h_LQ_c%i_y%i" %(idx, year))
            h_dy = f_in.Get(dir_ + "fpl")
            if(h_dy != None):
                h_dy = h_dy.Clone("h_fpl_c%i_y%i" %(idx, year))
            h = f_in.Get(dir_ + "fmn")
            if(h != None):
                h = h.Clone("h_fmn_c%i_y%i" %(idx, year))
                h_dy.Add(h)
            h = f_in.Get(dir_ + "alpha")
            if(h != None):
                h = h.Clone("h_alpha_c%i_y%i" %(idx, year))
                h_dy.Add(h)

            h_sig.Print()
            hist_list.append(h_sig)
            label_list.append(label_color_map["LQ"][0])
            color_list.append(label_color_map["LQ"][1])
            h_dy.Print()
            hist_list.append(h_dy)
            label_list.append(label_color_map["dy"][0])
            color_list.append(label_color_map["dy"][1])

                
        makeCan(dir_[:-1], options.output, [h_data], bkglist=[hist_list], totlist=[h_tot], colors = color_list, bkgNames = label_list, titles = [title], xtitle = "Template Bins" ,year = year, datastyle=datastyle) 

for key in fracs.keys():
    fracs[key] /= (len(years)*len(dirs))

#key_ = "gam"
#print("Average fraction for %s  is  %.3f \n" % (key_, fracs[key_]))
