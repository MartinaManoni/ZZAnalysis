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

inFilenameMC2022     = 'H4l_MC_2022_RDF.root'
inFilenameMC2022EE   = 'H4l_MC_2022EE_RDF.root'
inFilenameData2022   = "H4l_Data_2022_RDF.root"
inFilenameData2022EE = "H4l_Data_2022EE_RDF.root"
inFilenameZX2022   = "H4l_ZX_2022_RDF.root"
inFilenameZX2022EE = "H4l_ZX_2022EE_RDF.root"
outFilename = "Plots_full2022_RDF.root"

## output directory
today = date.today()
print('Creating output dir...')
out_dir = str(today)+'_plots_full2022'
os.makedirs(out_dir, exist_ok=True) #check if output dir exist

# lumi
lumi_2022 = 34.65 # 1/fb 2022 C-G  (= 35.181930231/fb of full 355100_362760 Golden json - 0.097685694 of eraB that we don't use
lumi_CD   = 7.98 # 1/fb
lumi_EFG  = 26.67 # 1/fb


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
def getZX(h_model, finalState):
    print("doing XS")
    n_entries = 10000
    bin_down = 70.
    bin_up = 3000.
    lumiRun3full= 62*1000
   
    f_4e_comb    = ROOT.TF1("f_4e_comb", "TMath::Landau(x, [0], [1])", bin_down, bin_up)
    f_4mu_comb   = ROOT.TF1("f_4mu_comb","TMath::Landau(x, [0], [1])", bin_down, bin_up)
    f_2e2mu_comb = ROOT.TF1("f_2e2mu_comb","[0]*TMath::Landau(x, [1], [2]) + [3]*TMath::Landau(x, [4], [5])", bin_down, bin_up)

    f_4e_comb.SetParameters(141.9, 21.3)
    f_4mu_comb.SetParameters(130.4, 15.6)
    f_2e2mu_comb.SetParameters(0.45,131.1,18.1, 0.55,133.8,18.9)

    yield_Comb_4e_2018    = 24.256/lumiRun3full #7.052
    print("yield_Comb_4e_2018 NOSR", yield_Comb_4e_2018)
    yield_Comb_4mu_2018   = 73.199/lumiRun3full #33.32
    yield_Comb_2e2mu_2018 = 86.14/lumiRun3full #8.905

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


#ZX estaimation parameters - taken from 2018 data - approx. normalization, just for visualization purposes
def getZX_SR(h_model, finalState):
    print("doing ZX for SR")
    n_entries = 10000
    bin_down = 105.
    bin_up = 160.
    lumiRun3full= 62*1000
   
    f_4e_comb    = ROOT.TF1("f_4e_comb", "TMath::Landau(x, [0], [1])", bin_down, bin_up)
    f_4mu_comb   = ROOT.TF1("f_4mu_comb","TMath::Landau(x, [0], [1])", bin_down, bin_up)
    f_2e2mu_comb = ROOT.TF1("f_2e2mu_comb","[0]*TMath::Landau(x, [1], [2]) + [3]*TMath::Landau(x, [4], [5])", bin_down, bin_up)

    f_4e_comb.SetParameters(141.9, 21.3)
    f_4mu_comb.SetParameters(130.4, 15.6)
    f_2e2mu_comb.SetParameters(0.45,131.1,18.1, 0.55,133.8,18.9)

    yield_Comb_4e_2018    = 7.052/lumiRun3full #7.052
    print("yield_Comb_4e_2018 SISR", yield_Comb_4e_2018)
    yield_Comb_4mu_2018   = 33.32/lumiRun3full #33.32
    yield_Comb_2e2mu_2018 = 33.45/lumiRun3full #8.905

    h_4e=h_model.Clone("ZX_4e")
    h_4e.Reset()
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

    print('Final State SR:', finalState)
    print("Z+X integral SR", h_total.Integral())
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
def Stack(f2022, f2022EE, f2022ZX, f2022EEZX, variable="ZZmass_", version = "", finalState = 'fs_4l'):
    print("Stacking")
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
    name = variable + (version if variable == "ZZMass_" else "") + fs_string
    print("variable", variable)
    print("version", version)
    print("fs_string", fs_string)
    print(f"[Stack] Using histogram base name: {name}")

    
    #------------EW------------------#
    #---------------------#
    # 2022
    WWZ  = f2022.Get(name+"WWZ")
    print("WWZ",name)
    WZZ  = f2022.Get(name+"WZZ")
    ZZZ  = f2022.Get(name+"ZZZ")
    EWSamples = [WWZ, WZZ, ZZZ]
    EW = WWZ.Clone("h_EW")
    for i in EWSamples:
        EW.Add(i,1.)
    EW.Scale(lumi_CD*1000.)

    #---------------------#
    # 2022EE
    WWZee  = f2022EE.Get(name+"WWZ")
    WZZee  = f2022EE.Get(name+"WZZ")
    ZZZee  = f2022EE.Get(name+"ZZZ")
    EWSamplesee = [WZZee, ZZZee]
    EWee = WWZee.Clone("h_EWee")
    for i in EWSamplesee:
        EWee.Add(i,1.)
    EWee.Scale(lumi_EFG*1000.)

    # full 2022 histo
    EW.Add(EWee,1.) #full 2022

    EW.SetLineColor(ROOT.TColor.GetColor("#000099"))
    EW.SetFillColor(ROOT.TColor.GetColor("#0331B9"))

    
    #-----------qqZZ---------------#
    # 2022
    ZZTo4l = f2022.Get(name+"ZZTo4l")
    ZZTo4l.Scale(lumi_CD*1000.) 
    # 2022EE
    ZZTo4lee = f2022EE.Get(name+"ZZTo4l")
    ZZTo4lee.Scale(lumi_EFG*1000.) 
    # full 2022 histo
    ZZTo4l.Add(ZZTo4lee,1.) #full 2022

    ZZTo4l.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ZZTo4l.SetFillColor(ROOT.TColor.GetColor("#99ccff"))
    
    #-----------signal------------#
    #2022
    VBF125     = f2022.Get(name+"VBFH125")
    ggH125     = f2022.Get(name+"ggH125") 
    WplusH125  = f2022.Get(name+"WplusH125")
    WminusH125 = f2022.Get(name+"WminusH125")
    ZH125      = f2022.Get(name+"ZH125")
    ttH125     = f2022.Get(name+"ttH125")
    
    signalSamples = [ggH125, WplusH125, WminusH125, ZH125, ttH125]
    print("signalSamples", signalSamples)
    signal = VBF125.Clone("h_signal")
    for i in signalSamples:
        signal.Add(i, 1.)
    signal.Scale(lumi_CD*1000.) 

    # 2022EE (E-G)
    VBF125ee     = f2022EE.Get(name+"VBFH125")
    ggH125ee     = f2022EE.Get(name+"ggH125")
    WplusH125ee  = f2022EE.Get(name+"WplusH125")
    WminusH125ee = f2022EE.Get(name+"WminusH125")
    ZH125ee      = f2022EE.Get(name+"ZH125")
    ttH125ee     = f2022EE.Get(name+"ttH125")
    
    signalSamplesee = [ggH125ee, WplusH125ee, WminusH125ee, ZH125ee, ttH125ee]
    signalee = VBF125ee.Clone("h_signalee")
    for i in signalSamplesee:
        signalee.Add(i, 1.)
    signalee.Scale(lumi_EFG*1000.)

    # full 2022 histo
    signal.Add(signalee,1.) # full 2022
      
    signal.SetLineColor(ROOT.TColor.GetColor("#cc0000"))
    signal.SetFillColor(ROOT.TColor.GetColor("#ff9b9b"))

    
    #------------ggTo-----------------#
    # 2022
    ggTo4mu     = f2022.Get(name+"ggTo4mu_Contin_MCFM701") 
    ggTo4e      = f2022.Get(name+"ggTo4e_Contin_MCFM701")
    ggTo4tau    = f2022.Get(name+"ggTo4tau_Contin_MCFM701")
    ggTo2e2mu   = f2022.Get(name+"ggTo2e2mu_Contin_MCFM701")
    ggTo2e2tau  = f2022.Get(name+"ggTo2e2tau_Contin_MCFM701")
    ggTo2mu2tau = f2022.Get(name+"ggTo2mu2tau_Contin_MCFM701")

    ggZZSamples = [ ggTo4e, ggTo4tau, ggTo2e2mu, ggTo2e2tau, ggTo2mu2tau]
    ggToZZ = ggTo4mu.Clone("h_ggTo")
    for i in ggZZSamples:
        ggToZZ.Add(i,1.)
    ggToZZ.Scale(lumi_CD*1000.)

     # 2022
    ggTo4muee     = f2022EE.Get(name+"ggTo4mu_Contin_MCFM701") 
    ggTo4eee     = f2022EE.Get(name+"ggTo4e_Contin_MCFM701")
    ggTo4tauee    = f2022EE.Get(name+"ggTo4tau_Contin_MCFM701")
    ggTo2e2muee  = f2022EE.Get(name+"ggTo2e2mu_Contin_MCFM701")
    ggTo2e2tauee  = f2022EE.Get(name+"ggTo2e2tau_Contin_MCFM701")
    ggTo2mu2tauee = f2022EE.Get(name+"ggTo2mu2tau_Contin_MCFM701")

    ggZZSamplesee= [ ggTo4eee, ggTo4tauee, ggTo2e2muee, ggTo2e2tauee, ggTo2mu2tauee]
    ggToZZee = ggTo4muee.Clone("h_ggToee")
    for i in ggZZSamplesee:
        ggToZZee.Add(i,1.)
    ggToZZee.Scale(lumi_EFG*1000.)

    ggToZZ.Add(ggToZZee, 1.)

    ggToZZ.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ggToZZ.SetFillColor(ROOT.TColor.GetColor("#4b78ff"))  
    
    ##############################
    #ZX 2022
    hzx= f2022ZX.Get((name+"ZX"))
    #hzx.Scale(lumi_CD*1000.)
    hzxee= f2022EEZX.Get((name+"ZX"))
    #hzxee.Scale(lumi_EFG*1000.)
    hzx.Add(hzxee,1.) #full 2022
    
    hzx.SetLineColor(ROOT.TColor.GetColor("#003300"))
    hzx.SetFillColor(ROOT.TColor.GetColor("#669966"))
    

    ###########################
    ### ZX
    '''if "Z1Mass_SR_" in name or "Z2Mass_SR_" in name:
        print("hzx FOR SR")
        hzx = getZX_SR(signal, finalState)
    else:
        print("hzx NOT FOR SR")
        hzx = getZX(signal, finalState)
    #hzx=getZX(signal, finalState)
    hzx.Scale(lumi_2022*1000.)
    hzx.SetLineColor(ROOT.TColor.GetColor("#003300"))
    hzx.SetFillColor(ROOT.TColor.GetColor("#669966"))'''
    
    
    #------------------Stack----------#
    if variable == "ZZmass_":
        if version=="4GeV_" :
            hs = ROOT.THStack("Stack_4GeV", "; m_{#it{4l}} (GeV) ; Events / 4 GeV" )
        elif version=="10GeV_" :
            hs = ROOT.THStack("Stack_10GeV", "; m_{#it{4l}} (GeV) ; Events / 10 GeV" )
        else:
            hs = ROOT.THStack("Stack_2GeV", "; m_{#it{4l}} (GeV) ; Events / 2 GeV" )
    else:
        hs = ROOT.THStack(f"Stack_{variable}", f"; {variable} (GeV); Events")

    hs.Add(hzx,"HISTO")
    hs.Add(EW,"HISTO")
    hs.Add(ggToZZ,"HISTO")
    hs.Add(ZZTo4l,"HISTO")
    hs.Add(signal,"HISTO")
    
    return hs, [hzx, EW, ggToZZ, ZZTo4l, signal]


### Get a TGraph for data, blinded if required
def dataGraph (f1, f2, variable="ZZMass_", version = "", finalState = 'fs_4l', blind = True):

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
    name = variable + (version if variable == "ZZMass_" else "") + fs_string
    print(name)

    hd1 = f1.Get(name+"Data")
    hd2 = f2.Get(name+"Data")
    hd = hd1.Clone('h_data') # full 2022
    hd.Add(hd2,1.)
    
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
    fZX2022   = ROOT.TFile.Open(inFilenameZX2022,"READ")
    fZX2022EE = ROOT.TFile.Open(inFilenameZX2022EE,"READ")
    of = ROOT.TFile.Open(outFilename,"recreate")


    # Labels for log plots
    xlabelsv = [80, 100, 200, 300, 400, 500]
    label_margin = -0.1
    xlabels=[None]*len(xlabelsv)
    for i, label in enumerate(xlabelsv): 
        xlabels[i] = ROOT.TLatex(label, label_margin , str(label));
        xlabels[i].SetTextAlign(23)
        xlabels[i].SetTextFont(42)
        xlabels[i].SetTextSize(0.04)

    ## --- plots     
    finalStates = ['fs_4e', 'fs_4mu', 'fs_2e2mu', 'fs_4l']
    variables = ['ZZMass_', 'Z1Mass_', 'Z2Mass_', 'Z1Mass_SR_', 'Z2Mass_SR_']

    for fs in finalStates: 
        print(f"\nProcessing final state: {fs}")
    
        for variable in variables:
            # Only loop over versions for ZZMass
            versions = ["2GeV_", "4GeV_"] if variable == "ZZMass_" else [""]

            for version in versions:
                print(f"  Variable: {variable}, version: {version if version else 'no version'}")

                # Load histograms
                HStack_var, h_list_var = Stack(fMC2022, fMC2022EE, fZX2022, fZX2022EE, variable=variable, version=version, finalState=fs)
                HData_var = dataGraph(fData2022, fData2022EE, variable=variable, version=version, finalState=fs, blind=(blindPlots and variable=="ZZMass_"))

                # Create canvas
                canvas_name = f"{variable}{version}_full2022_{fs}".replace("__", "_")
                Canvas = ROOT.TCanvas(canvas_name, canvas_name, canvasSizeX, canvasSizeY)
                Canvas.SetTicks()
                #if variable == "ZZMass":
                    #Canvas.SetLogx()

                xmin = ctypes.c_double(0.)
                ymin = ctypes.c_double(0.)
                xmax = ctypes.c_double(0.)
                ymax = ctypes.c_double(0.)
                HData_var.ComputeRange(xmin, ymin, xmax, ymax)
                HStack_var.SetMaximum(math.ceil(max(HStack_var.GetMaximum(), ymax.value)))
                HStack_var.Draw("histo")

                # Axis range (adjust if needed per variable)
                if variable == "ZZMass":
                    HStack_var.GetXaxis().SetRangeUser(70., 300.)
                elif variable in ["Z1Mass", "Z2Mass", "Z1Mass_SR", "Z2Mass_SR"]:
                    HStack_var.GetXaxis().SetRangeUser(0., 120.)

                # Blinding only for ZZMass
                if blindPlots and variable == "ZZMass":
                    ROOT.gPad.GetRangeAxis(xmin, ymin, xmax, ymax)
                    bblind = ROOT.TBox(blindHLow, 0, blindHHi, ymax.value - epsilon)
                    bblind.SetFillColor(ROOT.kGray)
                    bblind.SetFillStyle(3002)
                    bblind.Draw()

                HData_var.Draw("samePE1")

                # Axis label tweaks (only for ZZMass)
                if variable == "ZZMass":
                    HStack_var.GetXaxis().SetLabelSize(0)
                    for label in xlabels:
                        label.Draw()

                ROOT.gPad.RedrawAxis()

                # Legend
                legend = ROOT.TLegend(0.72, 0.70, 0.94, 0.92)
                legend.AddEntry(HData_var, "Data", "p")
                legend.AddEntry(h_list_var[4], "H(125)", "f")
                legend.AddEntry(h_list_var[3], "q#bar{q}→ZZ", "f")
                legend.AddEntry(h_list_var[2], "gg→ZZ", "f")
                legend.AddEntry(h_list_var[1], "EW", "f")
                legend.AddEntry(h_list_var[0], "Z+X", "f")
                legend.SetFillColor(ROOT.kWhite)
                legend.SetLineColor(ROOT.kWhite)
                legend.SetTextFont(43)
                legend.SetTextSize(20)
                legend.Draw()

                # CMS / Lumi
                CMS_lumi.writeExtraText = True
                CMS_lumi.extraText      = "Preliminary"
                CMS_lumi.lumi_sqrtS     = r"34.65~\mathrm{fb}^{-1} (13.6 TeV)"
                CMS_lumi.cmsTextSize    = 1
                CMS_lumi.lumiTextSize   = 0.7
                CMS_lumi.extraOverCmsTextSize = 0.75
                CMS_lumi.relPosX = 0.12
                CMS_lumi.CMS_lumi(Canvas, 0, 0)

                Canvas.Update()
                Canvas.Draw()
                Canvas.Write()
                printCanvas(Canvas, "png", path=out_dir)


        '''### Zoomed m4l
        HStack_z, h_list = Stack(fMC2022, fMC2022EE, '_2GeV_', fs)
        HData_z = dataGraph(fData2022, fData2022EE,'_2GeV_', fs, blind=blindPlots)
        Canvas_z = ROOT.TCanvas('M4l_z_full2022_'+fs,'M4l_z_full2022_'+fs,canvasSizeX,canvasSizeY)
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
        legend_z.AddEntry(h_list[3],"q#bar{q}#rightarrow ZZ,Z#gamma*","f")
        legend_z.AddEntry(h_list[2],"gg#rightarrow ZZ,Z#gamma*","f")
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
        CMS_lumi.lumi_sqrtS     = r"34.65~\mathrm{fb}^{-1} (13.6 TeV)"
        CMS_lumi.cmsTextSize    = 1 #0.6
        CMS_lumi.lumiTextSize   = 0.7 #0.46
        CMS_lumi.extraOverCmsTextSize = 0.75
        CMS_lumi.relPosX = 0.12
        CMS_lumi.CMS_lumi(Canvas_z, 0, 0)
            
        Canvas_z.Update() #very important!!!
        Canvas_z.Draw()
        Canvas_z.Write()'''