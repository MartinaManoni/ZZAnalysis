#!/bin/env python3
import ROOT
import argparse
import os

ROOT.EnableImplicitMT()

# Import utility functions
from ZZAnalysis.NanoAnalysis.tools import get_genEventSumw
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

# Define paths
MC_PATHS = {
    "2022": "/eos/user/m/mmanoni/HZZ_prod_170625/MC/2022/",
    "2022EE": "/eos/user/m/mmanoni/HZZ_prod_170625/MC/2022EE/",
    "2023preBPix": "/eos/user/m/mmanoni/HZZ_prod_170625/MC/2023preBPix/",
    "2023postBPix": "/eos/user/m/mmanoni/HZZ_prod_170625/MC/2023postBPix/",
}

DATA_PATHS = {
    "2022": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2022/Data_eraCD_preEE.root",
    "2022EE": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2022/Data_eraEFG_postEE.root",
    "2023preBPix": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2023/Data_eraC_preBPix.root",
    "2023postBPix": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2023/Data_eraD_postBPix.root",
}

ZX_PATHS = {
    "2022": "/afs/cern.ch/user/m/mmanoni/HZZ_CMSSSW14_new/CMSSW_14_1_6/src/ZZAnalysis/NanoAnalysis/plotters/ZX_results_2022.root",
    "2022EE": "/afs/cern.ch/user/m/mmanoni/HZZ_CMSSSW14_new/CMSSW_14_1_6/src/ZZAnalysis/NanoAnalysis/plotters/ZX_results_2022EE.root",
    "2023preBPix": "/afs/cern.ch/user/m/mmanoni/HZZ_CMSSSW14_new/CMSSW_14_1_6/src/ZZAnalysis/NanoAnalysis/plotters/ZX_results_2023preBPix.root",
    "2023postBPix": "/afs/cern.ch/user/m/mmanoni/HZZ_CMSSSW14_new/CMSSW_14_1_6/src/ZZAnalysis/NanoAnalysis/plotters/ZX_results_2023postBPix.root",
}

Z_FLAVORS = {
    "4mu": (-169, -169),
    "4e": (-121, -121),
    "2e2mu": [(-169, -121), (-121, -169)],
}



def define_histograms(df, df_SR, samplename, isMC):
    histos = {}
    weight_col = "weight" if isMC else None

    def book(df_input, hname, model, var):
        return df_input.Histo1D(model, var, weight_col) if weight_col else df_input.Histo1D(model, var)

    # General histograms (full range)
    histos[f"ZZMass_2GeV_{samplename}"] = book(df,
        f"ZZMass_2GeV_{samplename}",
        ROOT.RDF.TH1DModel("ZZMass_2GeV_" + samplename, "ZZMass_2GeV_" + samplename, 65, 70., 200.),
        "m4l")

    histos[f"ZZMass_4GeV_{samplename}"] = book(df,
        f"ZZMass_4GeV_{samplename}",
        ROOT.RDF.TH1DModel("ZZMass_4GeV_" + samplename, "ZZMass_4GeV_" + samplename, 233, 70., 1002.),
        "m4l")

    histos[f"ZZMass_10GeV_{samplename}"] = book(df,
        f"ZZMass_10GeV_{samplename}",
        ROOT.RDF.TH1DModel("ZZMass_10GeV_" + samplename, "ZZMass_10GeV_" + samplename, 93, 70., 1000.),
        "m4l")

    histos[f"Z1Mass_{samplename}"] = book(df,
        f"Z1Mass_{samplename}",
        ROOT.RDF.TH1DModel("Z1Mass_" + samplename, "Z1Mass_" + samplename, 40, 40., 120.),
        "Z1mass")

    histos[f"Z2Mass_{samplename}"] = book(df,
        f"Z2Mass_{samplename}",
        ROOT.RDF.TH1DModel("Z2Mass_" + samplename, "Z2Mass_" + samplename, 54, 12., 120.),
        "Z2mass")

    # Signal region only
    histos[f"Z1Mass_SR_{samplename}"] = book(df_SR,
        f"Z1Mass_SR_{samplename}",
        ROOT.RDF.TH1DModel("Z1Mass_SR_" + samplename, "Z1Mass_SR_" + samplename, 40, 40., 120.),
        "Z1mass")

    histos[f"Z2Mass_SR_{samplename}"] = book(df_SR,
        f"Z2Mass_SR_{samplename}",
        ROOT.RDF.TH1DModel("Z2Mass_SR_" + samplename, "Z2Mass_SR_" + samplename, 40, 0., 80.),
        "Z2mass")


    CHANNEL_TO_FINSTATE = {
        "4e": 0,
        "4mu": 1,
        "2e2mu": [2, 3],
    }

    for ch in ["4mu", "4e", "2e2mu"]:
        if samplename == "ZX":
            # Use FinState instead of Z1/Z2 flav
            finstates = CHANNEL_TO_FINSTATE[ch]
            if isinstance(finstates, list):
                selection = " || ".join([f"(finState == {fs})" for fs in finstates])
            else:
                selection = f"(finState == {finstates})"
        else:
            flav = Z_FLAVORS[ch]
            if isinstance(flav, tuple):
                selection = f"(Z1flav == {flav[0]} && Z2flav == {flav[1]})"
            else:
                selection = " || ".join([f"(Z1flav == {f[0]} && Z2flav == {f[1]})" for f in flav])

        
        df_ch = df.Filter(selection)
        df_SR_ch = df_SR.Filter(selection)

        histos[f"ZZMass_2GeV_{ch}_{samplename}"] = book(df_ch,
            f"ZZMass_2GeV_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"ZZMass_2GeV_{ch}_{samplename}", "", 65, 70., 200.),
            "m4l")

        histos[f"ZZMass_4GeV_{ch}_{samplename}"] = book(df_ch,
            f"ZZMass_4GeV_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"ZZMass_4GeV_{ch}_{samplename}", "", 233, 70., 1002.),
            "m4l")
        
        histos[f"ZZMass_10GeV_{ch}_{samplename}"] = book(df_ch,
            f"ZZMass_10GeV_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"ZZMass_10GeV_{ch}_{samplename}", "",93, 70., 1000.),
            "m4l")

        histos[f"Z1Mass_{ch}_{samplename}"] = book(df_ch,
            f"Z1Mass_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"Z1Mass_{ch}_{samplename}", "", 40, 40., 120.),
            "Z1mass")

        histos[f"Z2Mass_{ch}_{samplename}"] = book(df_ch,
            f"Z2Mass_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"Z2Mass_{ch}_{samplename}", "", 54, 12., 120.),
            "Z2mass")

        histos[f"Z1Mass_SR_{ch}_{samplename}"] = book(df_SR_ch,
            f"Z1Mass_SR_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"Z1Mass_SR_{ch}_{samplename}", "", 40, 40., 120.),
            "Z1mass")

        histos[f"Z2Mass_SR_{ch}_{samplename}"] = book(df_SR_ch,
            f"Z2Mass_SR_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"Z2Mass_SR_{ch}_{samplename}", "", 40, 0., 80.),
            "Z2mass")

    return histos

def run_sample(samplename, filepath, output_file, isMC):
    print(f"[INFO] Processing {samplename} from {filepath}")
    if not os.path.exists(filepath):
        print(f"[WARNING] File not found: {filepath}")
        return

    df = ROOT.RDataFrame("Events", filepath)

    # Basic cuts
    df = df.Filter("bestCandIdx != -1").Filter("HLT_passZZ4l")

    if isMC:
        genEventSumw = get_genEventSumw(ROOT.TFile.Open(filepath), 1e12)
        df = df.Define("genEventSumw", str(genEventSumw))
        df = df.Define("weight", "overallEventWeight * ZZCand_dataMCWeight / genEventSumw")

    df = df.Define("m4l", "ZZCand_mass[bestCandIdx]") \
           .Define("Z1mass", "ZZCand_Z1mass[bestCandIdx]") \
           .Define("Z2mass", "ZZCand_Z2mass[bestCandIdx]") \
           .Define("Z1flav", "ZZCand_Z1flav[bestCandIdx]") \
           .Define("Z2flav", "ZZCand_Z2flav[bestCandIdx]")

    df_SR = df.Filter("m4l >= 105 && m4l <= 160")

    histos = define_histograms(df, df_SR, samplename, isMC)

    output_file.cd()
    for hname, h in histos.items():
        h_clone = ROOT.TH1F(hname, hname, h.GetValue().GetNbinsX(), h.GetValue().GetXaxis().GetXmin(), h.GetValue().GetXaxis().GetXmax())
        for i in range(1, h_clone.GetNbinsX() + 1):
            h_clone.SetBinContent(i, h.GetValue().GetBinContent(i))
            h_clone.SetBinError(i, h.GetValue().GetBinError(i))
        h_clone.Write()

def run_zx(period):

    path = ZX_PATHS.get(period)
    if not path:
        raise ValueError(f"Unknown ZX period: {period}")

    if not os.path.exists(path):
        print(f"[WARNING] ZX file not found: {path}")
        return

    df = ROOT.RDataFrame("candTree", path)

    # Define necessary branches
    df = df.Define("m4l", "ZZMass") \
           .Define("Z1mass", "Z1Mass") \
           .Define("Z2mass", "Z2Mass") \
           .Define("Z1flav", "Z1Flav") \
           .Define("Z2flav", "Z2Flav") \
           .Define("finState", "FinState") \
           .Define("weight", "weight1")

    df_SR = df.Filter("m4l >= 105 && m4l <= 160")

    fout = ROOT.TFile.Open(f"H4l_ZX_{period}_RDF.root", "RECREATE")
    histos = define_histograms(df, df_SR, "ZX", isMC=True)

    for hname, h in histos.items():
        h_clone = ROOT.TH1F(hname, hname, h.GetValue().GetNbinsX(), h.GetValue().GetXaxis().GetXmin(), h.GetValue().GetXaxis().GetXmax())
        for i in range(1, h_clone.GetNbinsX() + 1):
            h_clone.SetBinContent(i, h.GetValue().GetBinContent(i))
            h_clone.SetBinError(i, h.GetValue().GetBinError(i))
        h_clone.Write()
    fout.Close()

def run_data(period):
    path = DATA_PATHS.get(period)
    if not path:
        raise ValueError(f"Unknown data period: {period}")
    fout = ROOT.TFile.Open(f"H4l_Data_{period}_RDF.root", "RECREATE")
    run_sample("Data", path, fout, isMC=False)
    fout.Close()

def run_mc(period):
    path = MC_PATHS.get(period)
    if not path:
        raise ValueError(f"Unknown MC period: {period}")

    samples = [
        "WWZ", "WZZ", "ZZZ", "ggTo4mu_Contin_MCFM701", "ggTo4e_Contin_MCFM701", "ggTo4tau_Contin_MCFM701",
        "ggTo2e2mu_Contin_MCFM701", "ggTo2e2tau_Contin_MCFM701", "ggTo2mu2tau_Contin_MCFM701", "ZZTo4l",
        "VBFH125", "ggH125", "WplusH125", "WminusH125", "ZH125", "ttH125"
    ]

    fout = ROOT.TFile.Open(f"H4l_MC_{period}_RDF.root", "RECREATE")
    for sample in samples:
        filename = os.path.join(path, sample, "ZZ4lAnalysis.root")
        run_sample(sample, filename, fout, isMC=True)
    fout.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["data", "mc", "zx", "both"], default="data")
    parser.add_argument("--period", choices=["2022", "2022EE", "2023preBPix", "2023postBPix"], default="2022")
    args = parser.parse_args()

    if args.mode == "data":
        run_data(args.period)
    elif args.mode == "mc":
        run_mc(args.period)
    elif args.mode == "zx":
        run_zx(args.period)
    elif args.mode == "both":
        run_mc(args.period)
        run_data(args.period)