#!/bin/env python3
from __future__ import print_function
import math
import argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw


# python H4l_fill_RDF_ZX.py --mode both --period 2022
# python H4l_fill_Martina.py --mode mc --period 2023preBPix
# python3 H4l_fill_Martina.py --mode both --period 2022
# python3 H4l_fill_Martina.py --mode both --period 2022EE
# python3 H4l_fill_Martina.py --mode both --period 2023preBPix
# python3 H4l_fill_Martina.py --mode both --period 2023postBPix


# Define paths for each period
MC_PATHS = {
    "2022": "/eos/user/m/mmanoni/HZZ_prod_170625/MC/2022/",
    "2022EE": "/eos/user/m/mmanoni/HZZ_prod_170625/MC/2022EE/",
    "2023preBPix": "/eos/user/m/mmanoni/HZZ_prod_170625/MC/2023preBPix/",
    "2023postBPix": "/eos/user/m/mmanoni/HZZ_prod_170625/MC/2023postBPix/",
}

DATA_PATHS = {
    "2022": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2022/Data_eraCD_preEE.root",
    "2022EE": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2022/Data_eraEFG_postEE.root",  # adjust as needed
    "2023preBPix": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2023/Data_eraC_preBPix.root",
    "2023postBPix": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2023/Data_eraD_postBPix.root",
}

ZmassValue = 91.1876
maxEntriesPerSample = 1e12
ROOT.TH1.SetDefaultSumw2()


def fill_histograms(samplename, filename):
    ### ---------------------    
    ## ZZMass

    h_ZZMass2 = ROOT.TH1F("ZZMass_2GeV_" + samplename, "ZZMass_2GeV_" + samplename, 65, 70., 200.)
    h_ZZMass2.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass2.GetYaxis().SetTitle("Events / 2 GeV")

    h_ZZMass4 = ROOT.TH1F("ZZMass_4GeV_" + samplename, "ZZMass_4GeV_" + samplename, 233, 70., 1002.)
    h_ZZMass4.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass4.GetYaxis().SetTitle("Events / 4 GeV")

    h_ZZMass10 = ROOT.TH1F("ZZMass_10GeV_" + samplename, "ZZMass_10GeV_" + samplename, 93, 70., 1000.)
    h_ZZMass10.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass10.GetYaxis().SetTitle("Events / 10 GeV")


    h_ZZMass2_4mu = ROOT.TH1F("ZZMass_2GeV_4mu_" + samplename, "", 65, 70., 200.)
    h_ZZMass2_4e = ROOT.TH1F("ZZMass_2GeV_4e_" + samplename, "", 65, 70., 200.)
    h_ZZMass2_2e2mu = ROOT.TH1F("ZZMass_2GeV_2e2mu_" + samplename, "", 65, 70., 200.)

    h_ZZMass4_4mu = ROOT.TH1F("ZZMass_4GeV_4mu_" + samplename, "", 233, 70., 1002.)
    h_ZZMass4_4e = ROOT.TH1F("ZZMass_4GeV_4e_" + samplename, "", 233, 70., 1002.)
    h_ZZMass4_2e2mu = ROOT.TH1F("ZZMass_4GeV_2e2mu_" + samplename, "", 233, 70., 1002.)

    h_ZZMass10_4mu = ROOT.TH1F("ZZMass_10GeV_4mu_" + samplename, "", 93, 70., 1000.)
    h_ZZMass10_4e = ROOT.TH1F("ZZMass_10GeV_4e_" + samplename, "", 93, 70., 1000.)
    h_ZZMass10_2e2mu = ROOT.TH1F("ZZMass_10GeV_2e2mu_" + samplename, "", 93, 70., 1000.)

    ### ---------------------    
    ## Z1 and Z2 masses
    # Z1
    h_Z1Mass = ROOT.TH1F("Z1Mass_"+samplename,
                         "Z1Mass_"+samplename,40,40.,120.)
    h_Z1Mass.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass.GetYaxis().SetTitle("Events / 2 GeV")
    #4mu
    h_Z1Mass_4mu = ROOT.TH1F("Z1Mass_4mu_"+samplename,
                             "Z1Mass_4mu_"+samplename,40,40.,120.)
    h_Z1Mass_4mu.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass_4mu.GetYaxis().SetTitle("Events / 2 GeV")
    #4e
    h_Z1Mass_4e = ROOT.TH1F("Z1Mass_4e_"+samplename,
                            "Z1Mass_4e_"+samplename,40,40.,120.)
    h_Z1Mass_4e.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass_4e.GetYaxis().SetTitle("Events / 2 GeV")
    #2e2mu
    h_Z1Mass_2e2mu = ROOT.TH1F("Z1Mass_2e2mu_"+samplename,
                               "Z1Mass_2e2mu_"+samplename,40,40.,120.)
    h_Z1Mass_2e2mu.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass_2e2mu.GetYaxis().SetTitle("Events / 2 GeV")

    ### ---------------------    
    # Z2
    h_Z2Mass = ROOT.TH1F("Z2Mass_"+samplename,
                         "Z2Mass_"+samplename,40,0.,80.)
    h_Z2Mass.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass.GetYaxis().SetTitle("Events / 2 GeV")
    #4mu
    h_Z2Mass_4mu = ROOT.TH1F("Z2Mass_4mu_"+samplename,
                             "Z2Mass_4mu_"+samplename,40,0.,80.)
    h_Z2Mass_4mu.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass_4mu.GetYaxis().SetTitle("Events / 2 GeV")
    #4e
    h_Z2Mass_4e = ROOT.TH1F("Z2Mass_4e_"+samplename,
                            "Z2Mass_4e_"+samplename,40,0.,80.)
    h_Z2Mass_4e.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass_4e.GetYaxis().SetTitle("Events / 2 GeV")
    #2e2mu
    h_Z2Mass_2e2mu = ROOT.TH1F("Z2Mass_2e2mu_"+samplename,
                               "Z2Mass_2e2mu_"+samplename,40,0.,80.)
    h_Z2Mass_2e2mu.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass_2e2mu.GetYaxis().SetTitle("Events / 2 GeV")

    #-------------------------------------------------------------------------
    f = ROOT.TFile.Open(filename)
    event = f.Events
    event.SetBranchStatus("*", 0)
    event.SetBranchStatus("run", 1)
    event.SetBranchStatus("luminosityBlock", 1)
    event.SetBranchStatus("*Muon*", 1)
    event.SetBranchStatus("*Electron*", 1)
    event.SetBranchStatus("*ZZCand*", 1)
    event.SetBranchStatus("bestCandIdx", 1)
    event.SetBranchStatus("HLT_passZZ4l", 1)

    nEntries = event.GetEntries()
    isMC = samplename != "Data"
    if isMC:
        event.SetBranchStatus("overallEventWeight", 1)
        genEventSumw = get_genEventSumw(f, maxEntriesPerSample)
    else:
        print(f"Running on data ({nEntries} entries)")

    printEntries = max(5000, nEntries // 10)
    for iEntry in range(nEntries):
        if not event.GetEntry(iEntry):
            continue
        if iEntry % printEntries == 0:
            print("Processing entry", iEntry)

        if event.bestCandIdx == -1 or not event.HLT_passZZ4l:
            continue

        weight = 1.0
        ZZs = Collection(event, 'ZZCand')
        theZZ = ZZs[event.bestCandIdx]
        if isMC:
            weight = (event.overallEventWeight * theZZ.dataMCWeight / genEventSumw)

        m4l = theZZ.mass
        h_ZZMass2.Fill(m4l, weight)
        h_ZZMass4.Fill(m4l, weight)
        h_ZZMass10.Fill(m4l, weight)

        ## Z1Mass
        mZ1=theZZ.Z1mass
        h_Z1Mass.Fill(mZ1,weight)
        ## Z2Mass
        mZ2=theZZ.Z2mass
        h_Z2Mass.Fill(mZ2,weight)

        Z1flav, Z2flav = theZZ.Z1flav, theZZ.Z2flav
        #print("Z1flav", Z1flav)
        if Z1flav == -169 and Z2flav == -169:
            #print("4mu")
            h_ZZMass2_4mu.Fill(m4l, weight)
            h_ZZMass4_4mu.Fill(m4l, weight)
            h_ZZMass10_4mu.Fill(m4l, weight)
            h_Z1Mass_4mu.Fill(mZ1,weight)
            h_Z2Mass_4mu.Fill(mZ2,weight)
        elif Z1flav == -121 and Z2flav == -121:
            #print("4e")
            h_ZZMass2_4e.Fill(m4l, weight)
            h_ZZMass4_4e.Fill(m4l, weight)
            h_ZZMass10_4e.Fill(m4l, weight)
            h_Z1Mass_4e.Fill(mZ1,weight)
            h_Z2Mass_4e.Fill(mZ2,weight)
        elif (Z1flav == -169 and Z2flav == -121) or (Z1flav == -121 and Z2flav == -169):
            #print("2e2mu")
            h_ZZMass2_2e2mu.Fill(m4l, weight)
            h_ZZMass4_2e2mu.Fill(m4l, weight)
            h_ZZMass10_2e2mu.Fill(m4l, weight)
            h_Z1Mass_2e2mu.Fill(mZ1,weight)
            h_Z2Mass_2e2mu.Fill(mZ2,weight)
        else:
            print(f"Warning: unexpected Z flavors {Z1flav}, {Z2flav}")

    f.Close()
    return h_ZZMass2, h_ZZMass4, h_ZZMass10, h_ZZMass2_4mu, h_ZZMass4_4mu, h_ZZMass10_4mu, h_ZZMass2_4e, h_ZZMass4_4e, h_ZZMass10_4e, h_ZZMass2_2e2mu, h_ZZMass4_2e2mu, h_ZZMass10_2e2mu, h_Z1Mass, h_Z1Mass_4mu, h_Z1Mass_4e, h_Z1Mass_2e2mu, h_Z2Mass, h_Z2Mass_4mu, h_Z2Mass_4e, h_Z2Mass_2e2mu


def run_data(period):
    path = DATA_PATHS.get(period)
    if not path:
        raise ValueError(f"Unknown data period: {period}")
    
    data_file = path
    output = ROOT.TFile.Open(f"H4l_Data_{period}.root", "recreate")
    histos = fill_histograms("Data", data_file)
    for h in histos:
        h.SetBinErrorOption(ROOT.TH1.kPoisson)
        output.WriteObject(h, h.GetName())
    output.Close()


def run_mc(period):
    path = MC_PATHS.get(period)
    if not path:
        raise ValueError(f"Unknown MC period: {period}")

    samples = [
        {"name": "WWZ", "filename": path + "WWZ/ZZ4lAnalysis.root"},
        {"name": "WZZ", "filename": path + "WZZ/ZZ4lAnalysis.root"},
        {"name": "ZZZ", "filename": path + "ZZZ/ZZ4lAnalysis.root"},
        {"name": "ggTo4mu", "filename": path + "ggTo4mu_Contin_MCFM701/ZZ4lAnalysis.root"},
        {"name": "ggTo4e", "filename": path + "ggTo4e_Contin_MCFM701/ZZ4lAnalysis.root"},
        {"name": "ggTo4tau", "filename": path + "ggTo4tau_Contin_MCFM701/ZZ4lAnalysis.root"},
        {"name": "ggTo2e2mu", "filename": path + "ggTo2e2mu_Contin_MCFM701/ZZ4lAnalysis.root"},
        {"name": "ggTo2e2tau", "filename": path + "ggTo2e2tau_Contin_MCFM701/ZZ4lAnalysis.root"},
        {"name": "ggTo2mu2tau", "filename": path + "ggTo2mu2tau_Contin_MCFM701/ZZ4lAnalysis.root"},
        {"name": "ZZTo4l", "filename": path + "ZZTo4l/ZZ4lAnalysis.root"},
        {"name": "VBF125", "filename": path + "VBFH125/ZZ4lAnalysis.root"},
        {"name": "ggH", "filename": path + "ggH125/ZZ4lAnalysis.root"},
        {"name": "WplusH125", "filename": path + "WplusH125/ZZ4lAnalysis.root"},
        {"name": "WminusH125", "filename": path + "WminusH125/ZZ4lAnalysis.root"},
        {"name": "ZH125", "filename": path + "ZH125/ZZ4lAnalysis.root"},
        {"name": "ttH125", "filename": path + "ttH125/ZZ4lAnalysis.root"},
    ]

    output = ROOT.TFile.Open(f"H4l_MC_{period}.root", "recreate")
    for sample in samples:
        histos = fill_histograms(sample["name"], sample["filename"])
        for h in histos:
            output.WriteObject(h, h.GetName())
    output.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run H4l histogram filling.")
    parser.add_argument(
        "--mode", choices=["data", "mc", "both"], default="data",
        help="Which type of sample to run on: data, mc, or both"
    )
    parser.add_argument(
        "--period", choices=["2022", "2022EE", "2023preBPix", "2023postBPix"], default="2022",
        help="Data-taking period to process"
    )
    args = parser.parse_args()

    if args.mode == "data":
        run_data(args.period)
    elif args.mode == "mc":
        run_mc(args.period)
    elif args.mode == "both":
        run_mc(args.period)
        run_data(args.period)
