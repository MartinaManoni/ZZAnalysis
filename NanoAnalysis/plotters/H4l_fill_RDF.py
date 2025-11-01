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
    "2024": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/062025/2024_MC/",
}

DATA_PATHS = {
    "2022": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2022/Data_eraCD_preEE.root",
    "2022EE": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2022/Data_eraEFG_postEE.root",
    "2023preBPix": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2023/Data_eraC_preBPix.root",
    "2023postBPix": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2023/Data_eraD_postBPix.root",
    "2024": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/062025/2024_Data/2024_Data.root",
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
        ROOT.RDF.TH1DModel("Z2Mass_SR_" + samplename, "Z2Mass_SR_" + samplename, 54, 12., 120.),
        "Z2mass")

    # Channel-separated histograms
    for ch, flav in Z_FLAVORS.items():
        selection = f"(Z1flav == {flav[0]} && Z2flav == {flav[1]})" if isinstance(flav, tuple) \
                    else " || ".join([f"(Z1flav == {f[0]} && Z2flav == {f[1]})" for f in flav])
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
            ROOT.RDF.TH1DModel(f"Z2Mass_SR_{ch}_{samplename}", "", 54, 12., 120.),
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

    if int(period) < 2024:
        samples = [
            "WWZ", "WZZ", "ZZZ", "ggTo4mu_Contin_MCFM701", "ggTo4e_Contin_MCFM701", "ggTo4tau_Contin_MCFM701",
            "ggTo2e2mu_Contin_MCFM701", "ggTo2e2tau_Contin_MCFM701", "ggTo2mu2tau_Contin_MCFM701", "ZZTo4l",
            "VBFH125", "ggH125", "WplusH125", "WminusH125", "ZH125", "ttH125"
        ]
    else:
        samples = [
            "WWZ", "WZZ", "ZZZ", "ZZTo4l",
            "VBFH125", "ggH125", "WplusH125", "WminusH125", "ttH125"
        ]


    fout = ROOT.TFile.Open(f"H4l_MC_{period}_RDF.root", "RECREATE")
    for sample in samples:
        filename = os.path.join(path, sample, "ZZ4lAnalysis.root")
        run_sample(sample, filename, fout, isMC=True)
    fout.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["data", "mc", "both"], default="data")
    parser.add_argument("--period", choices=["2022", "2022EE", "2023preBPix", "2023postBPix", "2024"], default="2022")
    args = parser.parse_args()

    if args.mode == "data":
        run_data(args.period)
    elif args.mode == "mc":
        run_mc(args.period)
    elif args.mode == "both":
        run_mc(args.period)
        run_data(args.period)
