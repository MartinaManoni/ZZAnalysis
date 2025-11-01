### Draw and decorate plots produced with H4l_fill.py.
# This is just a quick example, for illustration purposes only!!!
# It lacks all the plots and features of the full miniAOD plotter.
#
# run 
# python3 H4l_draw_mZZ_periods2022.py

from __future__ import print_function
import glob
import optparse
import os
import os.path as osp
import sys
from datetime import date
import math
import ctypes
import ROOT
ROOT.gROOT.SetBatch(True)
import CMSGraphics, CMS_lumi
import numpy as np
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True

inFilenameMC2022     = 'H4l_MC_2022.root'
inFilenameMC2022EE   = 'H4l_MC_2022EE.root'
inFilenameData2022   = "H4l_Data_2022.root"
inFilenameData2022EE = "H4l_Data_2022EE.root"

inFilenameMC2023preBPix     = 'H4l_MC_2023preBPix.root'
inFilenameMC2023postBPix   = 'H4l_MC_2023postBPix.root'
inFilenameData2023preBPix   = "H4l_Data_2023preBPix.root"
inFilenameData2023postBPix = "H4l_Data_2023postBPix.root"

outFilename = "Plots_22_23_Eras.root"

## output directory
today = date.today()
print('Creating output dir...')
out_dir = str(today)+'_plots_22_23_Eras'
os.makedirs(out_dir, exist_ok=True) #check if output dir exist

# lumi
lumi_2022 = 34.65 # 1/fb
lumi_CD   = 7.98 # 1/fb
lumi_EFG  = 26.67 # 1/fb

# lumi
lumi_2023 = 27.24 # 1/fb 
lumi_preBPix   = 17.79 # 1/fb
lumi_postBPix  = 9.45 # 1/fb

lumi_Run3 = 62


# plots options
blindPlots = True
blindHLow = 105.
blindHHi  = 140.
blindHM   = 500.
epsilon=0.1
addEmptyBins = True


# Set style matching the one used for HZZ plots
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetErrorX(0)
ROOT.gStyle.SetPadTopMargin(0.05)  
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadLeftMargin(0.16) 
ROOT.gStyle.SetPadRightMargin(0.03)
ROOT.gStyle.SetLabelOffset(0.008, "XYZ")
ROOT.gStyle.SetLabelSize(0.04, "XYZ")
ROOT.gStyle.SetAxisColor(1, "XYZ")
ROOT.gStyle.SetStripDecimals(True)
ROOT.gStyle.SetTickLength(0.03, "XYZ")
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetTitleSize(0.05, "XYZ")
ROOT.gStyle.SetTitleOffset(1.00, "X")
ROOT.gStyle.SetTitleOffset(1.25, "Y")
ROOT.gStyle.SetLabelOffset(0.008, "XYZ")
ROOT.gStyle.SetLabelSize(0.04, "XYZ")

canvasSizeX=910
canvasSizeY=700



#ZX estaimation parameters - taken from 2018 data - approx. normalization, just for visualization purposes
def getZX(h_model, finalState) :
    print("doing XS")
    n_entries = 10000
    bin_down  = 70.
    bin_up    = 3000.
    #lumi2018  = 59.7*1000. # to normalize
    lumiRun3full= 62*1000
   
    f_4e_comb    = ROOT.TF1("f_4e_comb", "TMath::Landau(x, [0], [1])", bin_down, bin_up)
    f_4mu_comb   = ROOT.TF1("f_4mu_comb","TMath::Landau(x, [0], [1])", bin_down, bin_up)
    f_2e2mu_comb = ROOT.TF1("f_2e2mu_comb","[0]*TMath::Landau(x, [1], [2]) + [3]*TMath::Landau(x, [4], [5])", bin_down, bin_up)

    f_4e_comb.SetParameters(141.9, 21.3)
    f_4mu_comb.SetParameters(130.4, 15.6)
    f_2e2mu_comb.SetParameters(0.45,131.1,18.1, 0.55,133.8,18.9)

    yield_Comb_4e_2018    = 24.256/lumiRun3full
    yield_Comb_4mu_2018   = 73.199/lumiRun3full
    yield_Comb_2e2mu_2018 = 86.14/lumiRun3full

    h_4e=h_model.Clone("ZX_4e")
    h_4e.Reset()
#    h_4e.SetFillColor(ROOT.TColor.GetColor("#0331B9"))
    h_4mu=h_4e.Clone("ZX_4mu")
    h_2e2mu=h_4e.Clone("ZX_2e2mu")
    
    h_4e.FillRandom("f_4e_comb"   , n_entries)
    h_4mu.FillRandom("f_4mu_comb"  , n_entries)
    h_2e2mu.FillRandom("f_2e2mu_comb", n_entries)

    h_4e.Scale(yield_Comb_4e_2018/h_4e.Integral())
    h_4mu.Scale(yield_Comb_4mu_2018/h_4mu.Integral())
    h_2e2mu.Scale(yield_Comb_2e2mu_2018/h_2e2mu.Integral())
    

    if(finalState == 'fs_4e'):
        h_total = h_4e
    elif(finalState == 'fs_4mu'):
        h_total = h_4mu
    elif(finalState == 'fs_2e2mu'):
        h_total = h_2e2mu
    elif(finalState == 'fs_4l'):
        h_total=h_4e.Clone("ZX_tot")
        h_total.Add(h_4mu)
        h_total.Add(h_2e2mu)
    else:
        raise ValueError('Error: wrong final state!')

    print('Final State:', finalState)
    print("Z+X integral", h_total.Integral())
    return h_total


#####################
def printCanvases(type="png", path=".") :
    canvases = ROOT.gROOT.GetListOfCanvases()
    for c in canvases :
        c.Print(path+"/"+c.GetTitle()+"."+type)

def printCanvas(c, type="png", name=None, path="." ) :
    if name == None : name = c.GetTitle()
    name=name.replace(">","")
    name=name.replace("<","")
    name=name.replace(" ","_")
    c.Print(path+"/"+name+"."+type)


######################
def Stack(f2022, lumi, version = "_4GeV_", finalState = 'fs_4l'):
    # final state
    if(finalState == 'fs_4e'):
        fs_string = '4e_'
        x_label = "m_{4e} (GeV)"
    elif(finalState == 'fs_4mu'):
        fs_string = '4mu_'
        x_label = "m_{4#mu} (GeV)"
    elif(finalState == 'fs_2e2mu'):
        fs_string = '2e2mu_'
        x_label = "m_{2e2#mu} (GeV)"
    elif(finalState == 'fs_4l'):
        fs_string = ''
        x_label = "m_{4l} (GeV)"
    else:
        raise ValueError('Error: wrong final state!')

    # define histo name
    name = "ZZMass" + version + fs_string
    print('hist name: ', name)

    
    #------------EW------------------#
    WWZ  = f2022.Get(name+"WWZ")
    print("WWZ", WWZ)
    WZZ  = f2022.Get(name+"WZZ")
    print("WZZ", WZZ)
    ZZZ  = f2022.Get(name+"ZZZ")
    print("ZZZ", ZZZ)
    #TTWW = f2022.Get(name+"TTWW")
    #TTZZ = f2022.Get(name+"TTZZ")
    print("qui")
    EWSamples = [WZZ, ZZZ]
    EW = WWZ.Clone("h_EW")
    print("qui 2")
    for i in EWSamples:
        print("DEBUG: Type of i:", type(i), " | Value:", i)
        print("i", i)
        EW.Add(i,1.)
    EW.Scale(lumi*1000.)
    print("qui 3")
    EW.SetLineColor(ROOT.TColor.GetColor("#000099"))
    EW.SetFillColor(ROOT.TColor.GetColor("#0331B9"))

    
    #-----------qqZZ---------------#
    ZZTo4l = f2022.Get(name+"ZZTo4l")
    ZZTo4l.Scale(lumi*1000.) 
       
    ZZTo4l.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ZZTo4l.SetFillColor(ROOT.TColor.GetColor("#99ccff"))
    print("qui 4")
    #-----------signal------------#
    VBF125     = f2022.Get(name+"VBF125")
    ggH125     = f2022.Get(name+"ggH")
    WplusH125  = f2022.Get(name+"WplusH125")
    WminusH125 = f2022.Get(name+"WminusH125")
    ZH125      = f2022.Get(name+"ZH125")
    ttH125     = f2022.Get(name+"ttH125")
    #bbH125     = f2022.Get(name+"bbH125")
    print("qui 5")
    signalSamples = [ggH125, WplusH125, WminusH125, ZH125, ttH125]
    print("signalSamples", signalSamples)
    signal = VBF125.Clone("h_signal")
    for i in signalSamples:
        signal.Add(i, 1.)
    signal.Scale(lumi*1000.) 
    print("qui 7")
    signal.SetLineColor(ROOT.TColor.GetColor("#cc0000"))
    signal.SetFillColor(ROOT.TColor.GetColor("#ff9b9b"))

    
    #------------ggTo-----------------#
    ggTo4mu     = f2022.Get(name+"ggTo4mu") 
    ggTo4e      = f2022.Get(name+"ggTo4e")
    ggTo4tau    = f2022.Get(name+"ggTo4tau")
    ggTo2e2mu   = f2022.Get(name+"ggTo2e2mu")
    ggTo2e2tau  = f2022.Get(name+"ggTo2e2tau")
    ggTo2mu2tau = f2022.Get(name+"ggTo2mu2tau")

    ggZZSamples = [ ggTo4e, ggTo4tau, ggTo2e2mu, ggTo2e2tau, ggTo2mu2tau]
    ggToZZ = ggTo4mu.Clone("h_ggTo")
    for i in ggZZSamples:
        ggToZZ.Add(i,1.)
    ggToZZ.Scale(lumi*1000.)
    ggToZZ.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ggToZZ.SetFillColor(ROOT.TColor.GetColor("#4b78ff"))  
    
    ### ZX
    # from 2018 for now
    hzx=getZX(signal, finalState)
    hzx.Scale(lumi*1000.)
    hzx.SetLineColor(ROOT.TColor.GetColor("#003300"))
    hzx.SetFillColor(ROOT.TColor.GetColor("#669966"))
    
    
    #------------------Stack----------#
    # Set the bin width label depending on version
    if version=="_4GeV_" :
        bin_width = "4 GeV"
    elif version=="_10GeV_" :
        bin_width = "10 GeV"
    else:
        bin_width = "2 GeV"

    # Create the stack with dynamic x-axis label:
    hs = ROOT.THStack("Stack"+version, f"; {x_label} ; Events / {bin_width}")


    print("QUI", hs.GetName(), type(hs), hs.ClassName())
    hs.Add(hzx,"HISTO")
    print(hzx.GetName(), type(hzx), hzx.ClassName())
    hs.Add(EW,"HISTO")
    hs.Add(ggToZZ,"HISTO")
    hs.Add(ZZTo4l,"HISTO")
    hs.Add(signal,"HISTO")
    
    return hs, [hzx, EW, ggToZZ, ZZTo4l, signal]


### Get a TGraph for data, blinded if required
def dataGraph (f, version = "_4GeV_", finalState = 'fs_4l', blind = True):

    # final state
    if(finalState == 'fs_4e'):
        fs_string = '4e_'
    elif(finalState == 'fs_4mu'):
        fs_string = '4mu_'
    elif(finalState == 'fs_2e2mu'):
        fs_string = '2e2mu_'
    elif(finalState == 'fs_4l'):
        fs_string = ''
    else:
        raise ValueError('Error: wrong final state!')

    # define histo name
    name = "ZZMass"+ version + fs_string
    print(name)
    hd = f.Get(name+"Data")
    
    nbinsIn = hd.GetNbinsX()
    nbins = 0
    
    x = np.array([0.]*nbinsIn, dtype='double')
    y = np.array([0.]*nbinsIn, dtype='double')
    errX = np.array([0.]*nbinsIn, dtype='double')  
    UpErr = np.array([0.]*nbinsIn, dtype='double')
    LowErr = np.array([0.]*nbinsIn, dtype='double')

    for i in range (1, nbinsIn):
        if blind and ((hd.GetBinCenter(i)>=blindHLow and hd.GetBinCenter(i)<=blindHHi) or hd.GetBinCenter(i)>=blindHM) : continue
        if not addEmptyBins and hd.GetBinContent(i) == 0 : continue 
        x[nbins]      = hd.GetBinCenter(i)
        y[nbins]      = hd.GetBinContent(i)
        UpErr[nbins]  = hd.GetBinErrorUp(i)
        LowErr[nbins] = hd.GetBinErrorLow(i)
        nbins += 1

    Data = ROOT.TGraphAsymmErrors(nbins,x,y,errX,errX,LowErr,UpErr)
    Data.SetMarkerStyle(20)
    Data.SetLineColor(ROOT.kBlack)
    Data.SetMarkerSize(0.9)
    return Data



## --------------------------------
if __name__ == "__main__" :
    

    fMC2022     = ROOT.TFile.Open(inFilenameMC2022,"READ")
    fMC2022EE   = ROOT.TFile.Open(inFilenameMC2022EE,"READ")
    fData2022   = ROOT.TFile.Open(inFilenameData2022,"READ")
    fData2022EE = ROOT.TFile.Open(inFilenameData2022EE,"READ")

    fMC2023preBPix     = ROOT.TFile.Open(inFilenameMC2023preBPix,"READ")
    fMC2023postBPix   = ROOT.TFile.Open(inFilenameMC2023postBPix,"READ")
    fData2023preBPix   = ROOT.TFile.Open(inFilenameData2023preBPix,"READ")
    fData2023postBPix = ROOT.TFile.Open(inFilenameData2023postBPix,"READ")


    # Labels for log plots
    xlabelsv = [80, 100, 200, 300, 400, 500]
    label_margin = -0.1
    xlabels=[None]*len(xlabelsv)
    for i, label in enumerate(xlabelsv): 
        xlabels[i] = ROOT.TLatex(label, label_margin , str(label))
        xlabels[i].SetTextAlign(23)
        xlabels[i].SetTextFont(42)
        xlabels[i].SetTextSize(0.04)

    ## --- plots     
    periods = ['2022CD', '2022EFG', '2023preBPix', '2023postBPix']
    finalStates = ['fs_4e', 'fs_4mu', 'fs_2e2mu', 'fs_4l']
    for p in periods:
        for fs in finalStates:
            print(p, fs)
            
            if(p == '2022CD'):
                fMC_2 = fMC2022
                fData = fData2022
                lumi = lumi_CD
                lumiText = r'7.98~\mathrm{fb}^{-1}'
            elif(p == '2022EFG'):
                fMC_2 = fMC2022EE
                fData = fData2022EE
                lumi = lumi_EFG
                lumiText = r'26.67~\mathrm{fb}^{-1}'
            elif(p == '2023preBPix'):
                fMC_2 = fMC2023preBPix
                fData = fData2023preBPix
                lumi = lumi_preBPix
                lumiText = r'17.79~\mathrm{fb}^{-1}'
            elif(p == '2023postBPix'):
                fMC_2 = fMC2023postBPix
                fData = fData2023postBPix
                lumi = lumi_postBPix
                lumiText = r'9.45~\mathrm{fb}^{-1}'
            else:
                raise ValueError('Error: wrong data-taking period!')

            ### m4l plot - full range
            HStack, h_list = Stack(fMC_2, lumi, '_4GeV_', fs)
            HData = dataGraph(fData, '_4GeV_', fs, blind=blindPlots)
            HStack_hm = HStack.Clone()
            HData_hm = HData.Clone()
          
            Canvas= ROOT.TCanvas('M4l_'+p+'_'+fs,'M4l_'+p+'_'+fs, canvasSizeX,canvasSizeY)
            Canvas.SetTicks()
            Canvas.SetLogx()
            #ymaxd=HData.GetMaximum()
            xmin=ctypes.c_double(0.)
            ymin=ctypes.c_double(0.)
            xmax=ctypes.c_double(0.)
            ymax=ctypes.c_double(0.)
            HData.ComputeRange(xmin,ymin,xmax,ymax)
            yhmax=math.ceil(max(HStack.GetMaximum(), ymax.value))
            HStack.SetMaximum(yhmax)
            HStack.Draw("histo")
            HStack.GetXaxis().SetRangeUser(70., 300.)
            if blindPlots:
                ROOT.gPad.GetRangeAxis(xmin,ymin,xmax,ymax)
                bblind = ROOT.TBox(blindHLow, 0, blindHHi, ymax.value-epsilon)
                bblind.SetFillColor(ROOT.kGray)
                bblind.SetFillStyle(3002)
                bblind.Draw()
            HData.Draw("samePE1")
            # Hide labels and rewrite them
            HStack.GetXaxis().SetLabelSize(0)
            for label in xlabels :
                label.Draw()
            ROOT.gPad.RedrawAxis()
        
            legend = ROOT.TLegend(0.72,0.70,0.94,0.92)
            legend.AddEntry(HData,"Data", "p")
            legend.AddEntry(h_list[4],"H(125)","f")
            legend.AddEntry(h_list[3],"q#bar{q}#rightarrow ZZ","f")
            legend.AddEntry(h_list[2],"gg#rightarrow ZZ","f")
            legend.AddEntry(h_list[1],"EW","f")
            legend.AddEntry(h_list[0],"Z+X","f")
            legend.SetFillColor(ROOT.kWhite)
            legend.SetLineColor(ROOT.kWhite)
            legend.SetTextFont(43)
            legend.SetTextSize(20)
            legend.Draw()
        
            #draw CMS and lumi text
            CMS_lumi.writeExtraText = True
            CMS_lumi.extraText      = "Preliminary"
            CMS_lumi.lumi_sqrtS     = lumiText + " (13.6 TeV)"
            CMS_lumi.cmsTextSize    = 1 #0.6
            CMS_lumi.lumiTextSize   = 0.7 #0.46
            CMS_lumi.extraOverCmsTextSize = 0.75
            CMS_lumi.relPosX = 0.12
            CMS_lumi.CMS_lumi(Canvas, 0, 0)
            
            Canvas.Update() #very important!!!
            Canvas.Draw()
            Canvas.Write()


            ### Zoomed m4l
            HStack_z, h_list = Stack(fMC_2, lumi, "_2GeV_", fs)
            HData_z = dataGraph(fData, "_2GeV_", fs, blind=blindPlots)
            Canvas_z = ROOT.TCanvas('M4l_z_'+p+'_'+fs,'M4l_z_'+p+'_'+fs,canvasSizeX,canvasSizeY)
            Canvas_z.SetTicks()
            HData_z.ComputeRange(xmin,ymin,xmax,ymax)
            yhmax=math.ceil(max(HStack_z.GetMaximum(), ymax.value))
            HStack_z.SetMaximum(yhmax)
            HStack_z.Draw("histo")
            HStack_z.GetXaxis().SetRangeUser(70., 170.)
            if blindPlots:
                 ROOT.gPad.GetRangeAxis(xmin,ymin,xmax,ymax)
                 bblind_z = ROOT.TBox(blindHLow, 0, blindHHi, ymax.value-epsilon)
                 bblind_z.SetFillColor(ROOT.kGray)
                 bblind_z.SetFillStyle(3002)
                 bblind_z.Draw()
            HData_z.Draw("samePE1")
            ROOT.gPad.RedrawAxis()
            
            legend_z = ROOT.TLegend(0.72,0.70,0.94,0.92)
            legend_z.AddEntry(HData,"Data", "p")
            legend_z.AddEntry(h_list[4],"H(125)","f")
            legend_z.AddEntry(h_list[3],"q#bar{q}#rightarrow ZZ","f")
            legend_z.AddEntry(h_list[2],"gg#rightarrow ZZ","f")
            legend_z.AddEntry(h_list[1],"EW","f")
            legend_z.AddEntry(h_list[0],"Z+X","f")
            legend_z.SetFillColor(ROOT.kWhite)
            legend_z.SetLineColor(ROOT.kWhite)
            legend_z.SetTextFont(43)
            legend_z.SetTextSize(20)
            legend_z.Draw()
            
            #draw CMS and lumi text
            CMS_lumi.writeExtraText = True
            CMS_lumi.extraText      = "Preliminary"
            CMS_lumi.lumi_sqrtS     = lumiText + " (13.6 TeV)"
            CMS_lumi.cmsTextSize    = 1 #0.6
            CMS_lumi.lumiTextSize   = 0.7 #0.46
            CMS_lumi.extraOverCmsTextSize = 0.75
            CMS_lumi.relPosX = 0.12
            CMS_lumi.CMS_lumi(Canvas_z, 0, 0)
            
            Canvas_z.Update()
            Canvas_z.Write()
    
    
            printCanvas(Canvas, "png", path=out_dir)
            printCanvas(Canvas_z, "png", path=out_dir)
    
    
    
    
    
    
    