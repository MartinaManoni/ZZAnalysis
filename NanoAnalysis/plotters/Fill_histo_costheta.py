#!/bin/env python3
from __future__ import print_function
import math
import argparse
from array import array
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

# ---------------------------------------------
# Histogram‑filling script for H→ZZ→4ℓ analysis
# Now fills cosθ₁ (polar angle of leading Z) instead of m₄ℓ.
# Binning: [-1.0,‑0.75,‑0.50,‑0.25,0.0,0.25,0.50,0.75,1.0]
# ---------------------------------------------

# Define paths for each period
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

maxEntriesPerSample = 1e12
ROOT.TH1.SetDefaultSumw2()

# --------------------------------------------------
# Core histogram‑filling routine
# --------------------------------------------------

def fill_histograms(samplename: str, filename: str):
    """Build cosθ₁ histograms for one sample."""

    # costheta1 histograms (inclusive + by final state)
    bin_edges = array('f', [-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0])

    h_costheta1       = ROOT.TH1F(f"costheta1_{samplename}",       f"cos(#theta_{{1}}) {samplename}",       len(bin_edges) - 1, bin_edges)
    h_costheta1_4mu   = ROOT.TH1F(f"costheta1_4mu_{samplename}",   "", len(bin_edges) - 1, bin_edges)
    h_costheta1_4e    = ROOT.TH1F(f"costheta1_4e_{samplename}",    "", len(bin_edges) - 1, bin_edges)
    h_costheta1_2e2mu = ROOT.TH1F(f"costheta1_2e2mu_{samplename}", "", len(bin_edges) - 1, bin_edges)

    for h in [h_costheta1, h_costheta1_4mu, h_costheta1_4e, h_costheta1_2e2mu]:
        h.GetXaxis().SetTitle("cos(#theta_{1})")
        h.GetYaxis().SetTitle("Events")

    # --- open file and configure branches ---
    f = ROOT.TFile.Open(filename)
    event = f.Events
    event.SetBranchStatus("*", 0)
    for br in [
        "run", "luminosityBlock", "bestCandIdx", "HLT_passZZ4l",
        "*Muon*", "*Electron*", "*ZZCand*",
    ]:
        event.SetBranchStatus(br, 1)

    nEntries = event.GetEntries()
    isMC = samplename != "Data"
    if isMC:
        event.SetBranchStatus("overallEventWeight", 1)
        genEventSumw = get_genEventSumw(f, maxEntriesPerSample)
    else:
        print(f"Running on data ({nEntries} entries)")

    printEvery = max(5000, nEntries // 10)

    for iEntry in range(nEntries):
        if not event.GetEntry(iEntry):
            continue
        if iEntry % printEvery == 0:
            print(f"Processing entry {iEntry}/{nEntries}")

        # --- event selection ---
        if event.bestCandIdx == -1 or not event.HLT_passZZ4l:
            continue

        weight = 1.0
        ZZs = Collection(event, "ZZCand")
        theZZ = ZZs[event.bestCandIdx]

        if isMC:
            weight = event.overallEventWeight * theZZ.dataMCWeight / genEventSumw

        # --- fill inclusive histogram ---
        h_costheta1.Fill(theZZ.costheta1, weight)

        # --- fill by final state ---
        Z1flav, Z2flav = theZZ.Z1flav, theZZ.Z2flav
        if Z1flav == -169 and Z2flav == -169:            # 4μ
            h_costheta1_4mu.Fill(theZZ.costheta1, weight)
        elif Z1flav == -121 and Z2flav == -121:          # 4e
            h_costheta1_4e.Fill(theZZ.costheta1, weight)
        elif (Z1flav == -169 and Z2flav == -121) or (Z1flav == -121 and Z2flav == -169):  # 2e2μ
            h_costheta1_2e2mu.Fill(theZZ.costheta1, weight)
        else:
            print(f"Warning: unexpected Z flavors {Z1flav}, {Z2flav}")

    f.Close()

    # Return histogram tuple for writing
    return (
        h_costheta1,
        h_costheta1_4mu,
        h_costheta1_4e,
        h_costheta1_2e2mu,
    )

# --------------------------------------------------
# Data and MC wrappers
# --------------------------------------------------

def run_data(period: str):
    path = DATA_PATHS.get(period)
    if not path:
        raise ValueError(f"Unknown data period: {period}")

    output = ROOT.TFile.Open(f"H4l_Data_{period}_costheta1.root", "recreate")
    for h in fill_histograms("Data", path):
        h.SetBinErrorOption(ROOT.TH1.kPoisson)
        output.WriteObject(h, h.GetName())
    output.Close()


def run_mc(period: str):
    path = MC_PATHS.get(period)
    if not path:
        raise ValueError(f"Unknown MC period: {period}")

    samples = [
        {"name": "WWZ",         "filename": path + "WWZ/ZZ4lAnalysis.root"},
        {"name": "WZZ",         "filename": path + "WZZ/ZZ4lAnalysis.root"},
        {"name": "ZZZ",         "filename": path + "ZZZ/ZZ4lAnalysis.root"},
        {"name": "ggTo4mu",     "filename": path + "ggTo4mu_Contin_MCFM701/ZZ4lAnalysis.root"},
        {"name": "ggTo4e",      "filename": path + "ggTo4e_Contin_MCFM701/ZZ4lAnalysis.root"},
        {"name": "ggTo4tau",    "filename": path + "ggTo4tau_Contin_MCFM701/ZZ4lAnalysis.root"},
        {"name": "ggTo2e2mu",   "filename": path + "ggTo2e2mu_Contin_MCFM701/ZZ4lAnalysis.root"},
        {"name": "ggTo2e2tau",  "filename": path + "ggTo2e2tau_Contin_MCFM701/ZZ4lAnalysis.root"},
        {"name": "ggTo2mu2tau", "filename": path + "ggTo2mu2tau_Contin_MCFM701/ZZ4lAnalysis.root"},
        {"name": "ZZTo4l",      "filename": path + "ZZTo4l/ZZ4lAnalysis.root"},
        {"name": "VBFH125",     "filename": path + "VBFH125/ZZ4lAnalysis.root"},
        {"name": "ggH125",      "filename": path + "ggH125/ZZ4lAnalysis.root"},
        {"name": "WplusH125",   "filename": path + "WplusH125/ZZ4lAnalysis.root"},
        {"name": "WminusH125",  "filename": path + "WminusH125/ZZ4lAnalysis.root"},
        {"name": "ZH125",       "filename": path + "ZH125/ZZ4lAnalysis.root"},
        {"name": "ttH125",      "filename": path + "ttH125/ZZ4lAnalysis.root"},
    ]

    output = ROOT.TFile.Open(f"H4l_MC_{period}_costheta1.root", "recreate")
    for sample in samples:
        for h in fill_histograms(sample["name"], sample["filename"]):
            output.WriteObject(h, h.GetName())
    output.Close()

# --------------------------------------------------
# Main entry
# --------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run H4l cos(theta1) histogram filling.")
    parser.add_argument(
        "--mode", choices=["data", "mc", "both"], default="data",
        help="Which type of sample to run on: data, mc, or both",
    )
    parser.add_argument(
        "--period", choices=["2022", "2022EE", "2023preBPix", "2023postBPix"], default="2022",
        help="Data‑taking period to process",
    )
    args = parser.parse_args()

    if args.mode == "data":
        run_data(args.period)
    elif args.mode == "mc":
        run_mc(args.period)
    elif args.mode == "both":
        run_mc(args.period)
        run_data(args.period)