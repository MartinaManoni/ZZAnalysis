import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
# =========================================================
# YEAR SELECTION
# =========================================================
year = "2023preBPix"

# =========================================================
# FILE PATH BUILDER
# =========================================================
#base_path = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES"
base_path = "root://eosuser.cern.ch//eos/user/m/mmanoni/ZX_studies"
# /eos/user/m/mmanoni/ZX_studies/2022_MC/DYJetsToLL/ZZ4lAnalysis_SKIMMED.root 
file_map = {
    "2022": {
        "DY": f"{base_path}/2022_MC/DYJetsToLL/ZZ4lAnalysis_SKIMMED.root",
        "TT": f"{base_path}/2022_MC/TTto2L2Nu/ZZ4lAnalysis_SKIMMED.root",
    },
    "2022EE": {
        "DY": f"{base_path}/2022EE_MC/DYJetsToLL/ZZ4lAnalysis_SKIMMED.root",
        "TT": f"{base_path}/2022EE_MC/TTto2L2Nu/ZZ4lAnalysis_SKIMMED.root",
    },
    "2023preBPix": {
        "DY": f"{base_path}/2023preBPix_MC/DYJetsToLL/ZZ4lAnalysis_SKIMMED.root",
        "TT": f"{base_path}/2023preBPix_MC/TTto2L2Nu/ZZ4lAnalysis_SKIMMED.root",
    },
    "2023postBPix": {
        "DY": f"{base_path}/2023postBPix_MC/DYJetsToLL/ZZ4lAnalysis_SKIMMED.root",
        "TT": f"{base_path}/2023postBPix_MC/TTto2L2Nu/ZZ4lAnalysis_SKIMMED.root",
    },
}

file_DY = file_map[year]["DY"]
file_TT = file_map[year]["TT"]

# =========================================================
# Helper functions
# =========================================================
def style(h, color, dashed=False):
    h.SetLineColor(color)
    h.SetLineWidth(2)
    if dashed:
        h.SetLineStyle(2)

def normalize(h):
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())

def set_range(hlist):
    ymax = max(h.GetMaximum() for h in hlist)
    for h in hlist:
        h.SetMaximum(ymax * 1.2)
        h.SetMinimum(0)

def final_state_filters():
    return {
        "4e":   "abs(Z1Flav)==121 && abs(Z2Flav)==121",
        "4mu":  "abs(Z1Flav)==169 && abs(Z2Flav)==169",
        "2e2mu": "(abs(Z1Flav)==121 && abs(Z2Flav)==169)",
        "2mu2e": "(abs(Z1Flav)==169 && abs(Z2Flav)==121)",  # optional if you want ordering
    }

# =========================================================
# NEW: flavour histogram helper
# =========================================================
def get_flavour_hist(df, tag, region):
    return df.Histo1D(
        (f"flav_{tag}_{region}", "Jet flavour;Flavour;Jets", 22, 0, 22),
        "jet_partonFlavour",
        "w"
    )

def get_flavour_hist_nj(df, tag, region, nj_cut):
    return df.Filter(nj_cut).Histo1D(
        (f"flav_{tag}_{region}_{nj_cut}", "Jet flavour;Flavour;Jets", 22, 0, 22),
        "jet_partonFlavour",
        "w"
    )

# =========================================================
# Build DataFrames
# =========================================================
def build_dfs(file):
    df_CR = ROOT.RDataFrame("CRZLLTree/candTree", file)
    df_SR = ROOT.RDataFrame("ZZTree/candTree", file)

    df_CR = df_CR.Define("w", "overallEventWeight")
    df_SR = df_SR.Define("w", "overallEventWeight")

    regions_CR = {
        "3P1F": df_CR.Filter("CRflag == 8388608"),
        "2P2F": df_CR.Filter("CRflag == 4194304"),
        "SS":   df_CR.Filter("CRflag == 2097152"),
        "SIP":  df_CR.Filter("CRflag == 21"),
    }

    return df_SR, regions_CR


df_SR_DY, regions_DY = build_dfs(file_DY)
df_SR_TT, regions_TT = build_dfs(file_TT)

# =========================================================
# Build Histograms (UNCHANGED)
# =========================================================
def make_histos(df_SR, regions, tag):
    h_Nj = {}
    h_M  = {}

    h_Nj["SR"] = df_SR.Histo1D((f"h_SR_Nj_{tag}", "Nj;N_{jets};Events", 8, 0, 8), "Nj", "w")
    h_M["SR"]  = df_SR.Histo1D((f"h_SR_M_{tag}",  "m4l;m_{4l} [GeV];Events", 40, 70, 180), "ZZMass", "w")

    for name, df in regions.items():
        h_Nj[name] = df.Histo1D((f"h_{name}_Nj_{tag}", "Nj;N_{jets};Events", 8, 0, 8), "Nj", "w")
        h_M[name]  = df.Histo1D((f"h_{name}_M_{tag}",  "m4l;m_{4l} [GeV];Events", 40, 70, 180), "ZZMass", "w")

    return h_Nj, h_M


h_Nj_DY, h_M_DY = make_histos(df_SR_DY, regions_DY, "DY")
h_Nj_TT, h_M_TT = make_histos(df_SR_TT, regions_TT, "TT")

# =========================================================
# Styling (UNCHANGED)
# =========================================================
colors = {
    "SR": ROOT.kBlack,
    "3P1F": ROOT.kRed,
    "2P2F": ROOT.kBlue,
    "SS": ROOT.kGreen + 2,
    "SIP": ROOT.kMagenta,
}

for key in h_Nj_DY:
    style(h_Nj_DY[key], colors[key])
    style(h_M_DY[key],  colors[key])

    style(h_Nj_TT[key], colors[key], dashed=True)
    style(h_M_TT[key],  colors[key], dashed=True)

# =========================================================
# Plotting functions (UNCHANGED)
# =========================================================
def plot_process(h_Nj, h_M, label):

    for h in h_Nj.values():
        normalize(h)
    set_range(h_Nj.values())

    c1 = ROOT.TCanvas(f"c_Nj_{label}", "", 800, 700)

    first = True
    for key in ["SR","3P1F","2P2F","SS","SIP"]:
        h_Nj[key].Draw("hist" if first else "hist same")
        first = False

    leg = ROOT.TLegend(0.65,0.65,0.88,0.88)
    for key in ["SR","3P1F","2P2F","SS","SIP"]:
        leg.AddEntry(h_Nj[key].GetPtr(), key, "l")
    leg.Draw()

    c1.SaveAs(f"Nj_{label}_{year}_comparison.png")


    set_range(h_M.values())

    c2 = ROOT.TCanvas(f"c_M_{label}", "", 800, 700)

    first = True
    for key in ["SR","3P1F","2P2F","SS","SIP"]:
        h_M[key].Draw("hist" if first else "hist same")
        first = False

    leg = ROOT.TLegend(0.65,0.65,0.88,0.88)
    for key in ["SR","3P1F","2P2F","SS","SIP"]:
        leg.AddEntry(h_M[key].GetPtr(), key, "l")
    leg.Draw()

    c2.SaveAs(f"m4l_{label}_{year}_comparison.png")

# =========================================================
# Plot overlay (UNCHANGED)
# =========================================================
def plot_overlay(h_DY, h_TT, var):

    for region in ["SR","3P1F","2P2F","SS","SIP"]:

        h1 = h_DY[region]
        h2 = h_TT[region]

        normalize(h1)
        normalize(h2)

        set_range([h1, h2])

        c = ROOT.TCanvas(f"c_{var}_{region}", "", 800, 700)

        h1.Draw("hist")
        h2.Draw("hist same")

        leg = ROOT.TLegend(0.65,0.75,0.88,0.88)
        leg.AddEntry(h1.GetPtr(), f"DY {region}", "l")
        leg.AddEntry(h2.GetPtr(), f"TT {region}", "l")
        leg.Draw()

        c.SaveAs(f"{var}_{region}_DY_vs_TT_{year}.png")

# =========================================================
# NEW: FLAVOUR PERCENTAGE COMPUTATION
# =========================================================
def compute_flavour(fr_dict, tag):

    out = {}

    for region, df in fr_dict.items():

        h = get_flavour_hist(df, tag, region).GetValue()

        total = h.Integral()

        if total > 0:
            for i in range(1, h.GetNbinsX()+1):
                h.SetBinContent(i, 100.0 * h.GetBinContent(i) / total)

        out[region] = h

    return out

# =========================================================
# NEW: RUN FLAVOUR STUDY
# =========================================================
flav_DY = compute_flavour(regions_DY, "DY")
flav_TT = compute_flavour(regions_TT, "TT")

# =========================================================
# OPTIONAL DRAWING (NEW)
# =========================================================
def draw_flav(flav, tag):

    for region, h in flav.items():

        c = ROOT.TCanvas(f"c_flav_{tag}_{region}", "", 700, 600)

        h.GetXaxis().SetTitle("Jet flavour (PDG)")
        h.GetYaxis().SetTitle("% jets")

        h.SetLineColor(ROOT.kBlack)
        h.SetFillColor(ROOT.kAzure+1)

        h.Draw("hist")

        c.SaveAs(f"flavour_{tag}_{region}_{year}.png")


draw_flav(flav_DY, "DY")
draw_flav(flav_TT, "TTbar")


def plot_flavour_CR_nj_splitFS(regions_DY, regions_TT, nj_cut, label):

    fs_filters = final_state_filters()

    for region in regions_DY.keys():
        for fs_name, fs_cut in fs_filters.items():

            full_cut = f"({nj_cut}) && ({fs_cut})"

            h_DY = get_flavour_hist_nj(regions_DY[region], "DY", f"{region}_{fs_name}", full_cut).GetValue()
            h_TT = get_flavour_hist_nj(regions_TT[region], "TT", f"{region}_{fs_name}", full_cut).GetValue()

            # normalize to %
            if h_DY.Integral() > 0:
                h_DY.Scale(1.0 / h_DY.Integral())

            if h_TT.Integral() > 0:
                h_TT.Scale(1.0 / h_TT.Integral())

            ymax = max(h_DY.GetMaximum(), h_TT.GetMaximum())

            h_DY.SetMaximum(ymax * 1.25)
            h_TT.SetMaximum(ymax * 1.25)
            h_DY.SetMinimum(0)
            h_TT.SetMinimum(0)

            c = ROOT.TCanvas(f"c_flav_{region}_{fs_name}_{label}", "", 700, 600)

            h_DY.SetLineColor(ROOT.kBlue)
            h_TT.SetLineColor(ROOT.kRed)

            h_DY.SetLineWidth(2)
            h_TT.SetLineWidth(2)

            h_DY.Draw("hist")
            h_TT.Draw("hist same")

            leg = ROOT.TLegend(0.62, 0.62, 0.95, 0.88) 
            leg.SetTextSize(0.035) 
            leg.SetBorderSize(0) 
            leg.SetFillStyle(0) 
            leg.SetHeader(f"{region} | {label}", "C") 
            leg.AddEntry(h_DY, "DY", "l") 
            leg.AddEntry(h_TT, "TT", "l") 
            leg.AddEntry("", "Jet flavour:", "") 
            leg.AddEntry("", "1 d, 2 u, 3 s", "") 
            leg.AddEntry("", "4 c, 5 b, 21 g", "")

            label_tex = label.replace("Nj_leq_1", "N_{j} \\leq 1").replace("Nj_geq_2", "N_{j} \\geq 2")

            leg.SetHeader(f"{region} | {fs_name} | {label_tex}")

            leg.Draw()

            c.SaveAs(f"flavour_CR_{region}_{fs_name}_{label}_{year}.png")

# =========================================================
# Run everything (UNCHANGED)
# =========================================================
plot_process(h_Nj_DY, h_M_DY, "DY")
plot_process(h_Nj_TT, h_M_TT, "TTbar")

plot_overlay(h_Nj_DY, h_Nj_TT, "Nj")
plot_overlay(h_M_DY,  h_M_TT,  "m4l")


# Nj <= 1
plot_flavour_CR_nj_splitFS(regions_DY, regions_TT, "Nj<=1", "Nj_leq_1")

# Nj >= 2
plot_flavour_CR_nj_splitFS(regions_DY, regions_TT, "Nj>=2", "Nj_geq_2")