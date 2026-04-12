import ROOT

ROOT.gStyle.SetOptStat(0)

file = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_MC/DYJetsToLL/ZZ4lAnalysis_SKIMMED.root"

file = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_MC/TTto2L2Nu/ZZ4lAnalysis_SKIMMED.root"
# --- Load dataframes ---
df_CR = ROOT.RDataFrame("CRZLLTree/candTree", file)
df_SR = ROOT.RDataFrame("ZZTree/candTree", file)

# --- Define weights (you can refine later) ---
df_CR = df_CR.Define("w", "overallEventWeight")
df_SR = df_SR.Define("w", "overallEventWeight")

# --- Bit selections ---
df_3P1F = df_CR.Filter("CRflag == 8388608")
df_2P2F = df_CR.Filter("CRflag == 4194304")
df_SS   = df_CR.Filter("CRflag == 2097152")
df_SIP  = df_CR.Filter("CRflag == 21")

# --- Histograms: Nj ---
h_SR_Nj   = df_SR.Histo1D(("h_SR_Nj",   "Nj;N_{jets};Events", 8, 0, 8), "Nj", "w")
h_3P1F_Nj = df_3P1F.Histo1D(("h_3P1F_Nj","Nj;N_{jets};Events", 8, 0, 8), "Nj", "w")
h_2P2F_Nj = df_2P2F.Histo1D(("h_2P2F_Nj","Nj;N_{jets};Events", 8, 0, 8), "Nj", "w")
h_SS_Nj   = df_SS.Histo1D(("h_SS_Nj",   "Nj;N_{jets};Events", 8, 0, 8), "Nj", "w")
h_SIP_Nj  = df_SIP.Histo1D(("h_SIP_Nj", "Nj;N_{jets};Events", 8, 0, 8), "Nj", "w")

# --- Histograms: m4l ---
h_SR_M   = df_SR.Histo1D(("h_SR_M",   "m4l;m_{4l} [GeV];Events", 40, 70, 180), "ZZMass", "w")
h_3P1F_M = df_3P1F.Histo1D(("h_3P1F_M","m4l;m_{4l} [GeV];Events", 40, 70, 180), "ZZMass", "w")
h_2P2F_M = df_2P2F.Histo1D(("h_2P2F_M","m4l;m_{4l} [GeV];Events", 40, 70, 180), "ZZMass", "w")
h_SS_M   = df_SS.Histo1D(("h_SS_M",   "m4l;m_{4l} [GeV];Events", 40, 70, 180), "ZZMass", "w")
h_SIP_M  = df_SIP.Histo1D(("h_SIP_M", "m4l;m_{4l} [GeV];Events", 40, 70, 180), "ZZMass", "w")

# --- Styling function ---
def style(h, color):
    h.SetLineColor(color)
    h.SetLineWidth(2)

style(h_SR_Nj, ROOT.kBlack)
style(h_3P1F_Nj, ROOT.kRed)
style(h_2P2F_Nj, ROOT.kBlue)
style(h_SS_Nj, ROOT.kGreen+2)
style(h_SIP_Nj, ROOT.kMagenta)

style(h_SR_M, ROOT.kBlack)
style(h_3P1F_M, ROOT.kRed)
style(h_2P2F_M, ROOT.kBlue)
style(h_SS_M, ROOT.kGreen+2)
style(h_SIP_M, ROOT.kMagenta)

# --- Canvas Nj ---
c1 = ROOT.TCanvas("c1","Nj comparison",800,700)

for h in [h_SR_Nj, h_3P1F_Nj, h_2P2F_Nj, h_SS_Nj, h_SIP_Nj]:
    integral = h.Integral()
    if integral > 0:
        h.Scale(1.0/integral)

ymax = max(h_SR_Nj.GetMaximum(),
           h_3P1F_Nj.GetMaximum(),
           h_2P2F_Nj.GetMaximum(),
           h_SS_Nj.GetMaximum(),
           h_SIP_Nj.GetMaximum())

for h in [h_SR_Nj, h_3P1F_Nj, h_2P2F_Nj, h_SS_Nj, h_SIP_Nj]:
    h.SetMaximum(ymax*1.2)
    h.SetMinimum(0)


h_SR_Nj.Draw("hist")
h_3P1F_Nj.Draw("hist same")
h_2P2F_Nj.Draw("hist same")
h_SS_Nj.Draw("hist same")
h_SIP_Nj.Draw("hist same")

leg1 = ROOT.TLegend(0.65,0.65,0.88,0.88)
leg1.AddEntry(h_SR_Nj.GetPtr(),   "SR", "l")
leg1.AddEntry(h_3P1F_Nj.GetPtr(), "3P1F", "l")
leg1.AddEntry(h_2P2F_Nj.GetPtr(), "2P2F", "l")
leg1.AddEntry(h_SS_Nj.GetPtr(),   "SS", "l")
leg1.AddEntry(h_SIP_Nj.GetPtr(),  "SIPCR", "l")
leg1.Draw()

c1.SaveAs("Nj_CR_comparison.png")

# --- Canvas m4l ---
c2 = ROOT.TCanvas("c2","m4l comparison",800,700)

# Compute max over all histograms
ymax = max(h_SR_M.GetMaximum(),
           h_3P1F_M.GetMaximum(),
           h_2P2F_M.GetMaximum(),
           h_SS_M.GetMaximum(),
           h_SIP_M.GetMaximum())

# Set same y-range for all histograms
for h in [h_SR_M, h_3P1F_M, h_2P2F_M, h_SS_M, h_SIP_M]:
    h.SetMaximum(ymax*1.2)  # add 20% headroom
    h.SetMinimum(0)

h_SR_M.Draw("hist")
h_3P1F_M.Draw("hist same")
h_2P2F_M.Draw("hist same")
h_SS_M.Draw("hist same")
h_SIP_M.Draw("hist same")

leg2 = ROOT.TLegend(0.65,0.65,0.88,0.88)
leg2.AddEntry(h_SR_M.GetPtr(),   "SR", "l")
leg2.AddEntry(h_3P1F_M.GetPtr(), "3P1F", "l")
leg2.AddEntry(h_2P2F_M.GetPtr(), "2P2F", "l")
leg2.AddEntry(h_SS_M.GetPtr(),   "SS", "l")
leg2.AddEntry(h_SIP_M.GetPtr(),  "SIPCR", "l")
leg2.Draw()

c2.SaveAs("m4l_CR_comparison.png")