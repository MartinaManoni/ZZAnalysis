import ROOT

ROOT.gStyle.SetOptStat(0)

# =========================================================
# YEAR SELECTION
# =========================================================
year = "2022"

# =========================================================
# FILE PATH BUILDER
# =========================================================
base_path = "root://eosuser.cern.ch//eos/user/m/mmanoni/ZX_studies"

file_map = {
    "2022": {
        "DY": f"{base_path}/2022_MC/DYJetsToLL/ZZ4lAnalysis_SKIMMED.root",
    },
}

file_DY = file_map[year]["DY"]

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

# =========================================================
# FLAVOUR HISTOGRAM
# =========================================================
def get_flavour_hist(df, tag, region, cut):
    return df.Filter(cut).Histo1D(
        (f"flav_{tag}_{region}_{cut}", "Jet flavour;Flavour;Jets", 10, 0, 10),
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

# =========================================================
# BUILD NJ (unchanged)
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

# =========================================================
# FLAVOUR STUDY (NEW LOGIC)
# =========================================================
def compute_flavour(fr_dict, tag, cut):

    out = {}

    for region, df in fr_dict.items():

        h = df.Filter(cut).Histo1D(
            (f"flav_{tag}_{region}_{cut}", "Jet flavour;Flavour;Jets", 10, 0, 10),
            "jet_partonFlavour",
            "w"
        ).GetValue()

        total = h.Integral()

        if total > 0:
            for i in range(1, h.GetNbinsX()+1):
                h.SetBinContent(i, 100.0 * h.GetBinContent(i) / total)

        out[region] = h

    return out

# =========================================================
# RUN FLAVOUR STUDY (SPLIT IN 2 NJ REGIONS)
# =========================================================
flav_DY_le1 = compute_flavour(regions_DY, "DY", "Nj <= 1")
flav_DY_ge2 = compute_flavour(regions_DY, "DY", "Nj >= 2")

# =========================================================
# DRAW FUNCTION
# =========================================================
def draw_flav(flav, tag):

    for region, h in flav.items():

        c = ROOT.TCanvas(f"c_flav_{tag}_{region}", "", 700, 600)

        h.GetXaxis().SetTitle("Jet flavour (PDG)")
        h.GetYaxis().SetTitle("% jets")

        h.SetLineColor(ROOT.kBlack)
        h.SetFillColor(ROOT.kAzure + 1)

        h.Draw("hist")

        c.SaveAs(f"flavour_{tag}_{region}_{year}.png")

# =========================================================
# DRAW BOTH REGIONS
# =========================================================
draw_flav(flav_DY_le1, "DY_Nj_le1")
draw_flav(flav_DY_ge2, "DY_Nj_ge2")

# =========================================================
# EVERYTHING ELSE UNCHANGED
# =========================================================
plot_process(h_Nj_DY, h_M_DY, "DY")