#by Matteo Bonanomi
import uproot
import numpy as np
import awkward as ak
from tqdm import tqdm
from collections import defaultdict
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.patches import Patch
import matplotlib.patches as patches
from scipy.stats import chi2
import ROOT

plt.style.use(hep.style.CMS)

def poisson_interval(k, alpha=0.317): 
    a = alpha
    low, high = (chi2.ppf(a/2, 2*k) / 2, chi2.ppf(1-a/2, 2*k + 2) / 2)
    return k-low, high-k

def get_genEventSumw(input_file, maxEntriesPerSample=None):
    f = input_file
    runs = f.Runs
    event = f.Events
    nRuns = runs.GetEntries()
    nEntries = event.GetEntries()

    iRun = 0
    genEventCount = 0
    genEventSumw = 0.

    while iRun < nRuns and runs.GetEntry(iRun):
        genEventCount += runs.genEventCount
        genEventSumw += runs.genEventSumw
        iRun += 1
    
    print("gen=", genEventCount, "sumw=", genEventSumw)

    if maxEntriesPerSample is not None:
        print(f"Scaling to {maxEntriesPerSample} entries")
        if nEntries > maxEntriesPerSample:
            genEventSumw = genEventSumw * maxEntriesPerSample / nEntries
            nEntries = maxEntriesPerSample
        print("    scaled to:", nEntries, "sumw=", genEventSumw)

    return genEventSumw

# File paths
#MATTEOS
#fname_dy = "/eos/user/m/mbonanom/nanoProd_Run3_2022EE_ControlPlots/samplesNanoDY_2022EE_MC_87a15506/DYJetsToLL/ZZ4lAnalysis.root"
#fname_tt = "/eos/user/m/mbonanom/nanoProd_Run3_2022EE_ControlPlots/samplesNanoTT_2022EE_MC_87a15506/TTto2L2Nu/ZZ4lAnalysis.root"
#fname_wz = "/eos/user/m/mbonanom/nanoProd_Run3_2022EE_ControlPlots/samplesNanoWZ_2022EE_MC_87a15506/WZto3LNu/ZZ4lAnalysis.root"
#fname_dt = "/eos/user/m/mbonanom/nanoProd_Run3_Data_ControlPlots/samplesNano_2022_Data_70a8ade8/post_EE/ZZ4lAnalysis.root"

#MARTINAS 2023 preBPix (ERA C) 
fname_dy = "/eos/user/m/mmanoni/HZZ_samples_2023/MC/PROD_samplesNano_2023preBPix_MC_DY_TT_WZ_EleSF/DYJetsToLL/ZZ4lAnalysis.root"
fname_tt = "/eos/user/m/mmanoni/HZZ_samples_2023/MC/PROD_samplesNano_2023preBPix_MC_DY_TT_WZ_EleSF/TTto2L2Nu/ZZ4lAnalysis.root"
fname_wz = "/eos/user/m/mmanoni/HZZ_samples_2023/MC/PROD_samplesNano_2023preBPix_MC_DY_TT_WZ_EleSF/WZto3LNu/ZZ4lAnalysis.root"
fname_dt = "/eos/user/m/mmanoni/HZZ_samples_2023/Data/PROD_samplesNano_2023_Data_eraC_preBPix/Data_eraC.root"

#MARTINAS 2023 postBPix (ERA D) 
'''fname_dy = "/eos/user/m/mmanoni/HZZ_samples_2023/MC/PROD_samplesNano_2023postBPix_MC_DY_TT_WZ_EleSF/DYJetsToLL/ZZ4lAnalysis.root"
fname_tt = "/eos/user/m/mmanoni/HZZ_samples_2023/MC/PROD_samplesNano_2023postBPix_MC_DY_TT_WZ_EleSF/TTto2L2Nu/ZZ4lAnalysis.root"
fname_wz = "/eos/user/m/mmanoni/HZZ_samples_2023/MC/PROD_samplesNano_2023postBPix_MC_DY_TT_WZ_EleSF/WZto3LNu/ZZ4lAnalysis.root"
fname_dt = "/eos/user/m/mmanoni/HZZ_samples_2023/Data/PROD_samplesNano_2023_Data_eraD_postBPix/Data_eraD.root"'''

# Determine luminosity
#17.794	Run3Summer23 era C pre BPix
#9.451	Run3Summer23BPix era D post BPix 
if "preBPix" in fname_dt:
    lumi = 17.794
else:
    lumi = 9.451

# Read files and get data
with uproot.open(fname_dy) as f:
    dy_zc = f['Events/bestZIdx'].array()
    dy_zm = f['Events/ZCand_mass'].array()
    dy_z1f = f['Events/ZCand_flav'].array()
    dy_ow = f['Events/overallEventWeight'].array()
    dy_cnt = get_genEventSumw(ROOT.TFile.Open(fname_dy))
    
with uproot.open(fname_tt) as f:
    tt_zc = f['Events/bestZIdx'].array()
    tt_zm = f['Events/ZCand_mass'].array()
    tt_z1f = f['Events/ZCand_flav'].array()
    tt_ow = f['Events/overallEventWeight'].array()
    tt_cnt = get_genEventSumw(ROOT.TFile.Open(fname_tt))
    
with uproot.open(fname_wz) as f:
    wz_zc = f['Events/bestZIdx'].array()
    wz_zm = f['Events/ZCand_mass'].array()
    wz_z1f = f['Events/ZCand_flav'].array()
    wz_ow = f['Events/overallEventWeight'].array()
    wz_cnt = get_genEventSumw(ROOT.TFile.Open(fname_wz))

with uproot.open(fname_dt) as f:
    dt_zc = f['Events/bestZIdx'].array()
    dt_zm = f['Events/ZCand_mass'].array()
    dt_z1f = f['Events/ZCand_flav'].array()

dt_z1mass = dt_zm[(ak.singletons(dt_zc))]
dt_z1flav = dt_z1f[(ak.singletons(dt_zc))]

m_dy = dy_zm[(ak.singletons(dy_zc))]
zf_dy = dy_z1f[(ak.singletons(dy_zc))]
w_dy = dy_ow

m_tt = tt_zm[(ak.singletons(tt_zc))]
zf_tt = tt_z1f[(ak.singletons(tt_zc))]
w_tt = tt_ow

m_wz = wz_zm[(ak.singletons(wz_zc))]
zf_wz = wz_z1f[(ak.singletons(wz_zc))]
w_wz = wz_ow

# Selections
#169 muons
#121 ele
#(abs(zf_dy) == 169) | (abs(zf_dy) == 121) inclusive
sel_dy = (abs(zf_dy) == 169) | (abs(zf_dy) == 121)
sel_dy = ak.flatten(sel_dy)
m_dy = m_dy[sel_dy]
w_dy = w_dy[sel_dy]

sel_tt = (abs(zf_tt) == 169) | (abs(zf_tt) == 121)
sel_tt = ak.flatten(sel_tt)
m_tt = m_tt[sel_tt]
w_tt = w_tt[sel_tt]

sel_wz = (abs(zf_wz) == 169) | (abs(zf_wz) == 121)
sel_wz = ak.flatten(sel_wz)
m_wz = m_wz[sel_wz]
w_wz = w_wz[sel_wz]

sel_dt = (abs(dt_z1flav) == 169) | (abs(dt_z1flav) == 121)
sel_dt = ak.flatten(sel_dt)
dt_z1mass = dt_z1mass[sel_dt]

m_tot = np.concatenate([ak.flatten(m_dy), ak.flatten(m_tt), ak.flatten(m_wz)])

w_tot_dy = w_dy * 1000 * lumi / dy_cnt
w_tot_tt = w_tt * 1000 * lumi / tt_cnt
w_tot_wz = w_wz * 1000 * lumi / wz_cnt

w_tot = np.concatenate([w_tot_dy, w_tot_tt, w_tot_wz])

CP_BINNING = np.linspace(60, 120, 100)

# The scaling of DY and TT xsecs are there because we spotted some wrong values in the csv used for the prod
# Doesn't really matter if then we normalize MC to Data. OTOH TT is off by one order of magnitude, so better to correct

n_dy, bins = np.histogram(ak.flatten(m_dy), weights=w_dy * 1000 * lumi / dy_cnt * (6225.4 / 5558.0), bins=CP_BINNING)
n_tt, bins = np.histogram(ak.flatten(m_tt), weights=w_tt * 1000 * lumi / tt_cnt * (87.3 / 762.1), bins=CP_BINNING)
n_wz, bins = np.histogram(ak.flatten(m_wz), weights=w_wz * 1000 * lumi / wz_cnt, bins=CP_BINNING)

n_dt, bins_dt = np.histogram(ak.flatten(dt_z1mass), bins=CP_BINNING)
mc_tot = n_dy + n_tt + n_wz

# Plot and normalize to Data
binsc = 0.5 * (bins[1:] + bins[:-1])
fig, ax = plt.subplots(
    nrows=2, 
    ncols=1,
    figsize=(10, 8),
    gridspec_kw={"height_ratios": (3, 1)},
    sharex=True
)

binsc_dt = 0.5 * (bins_dt[1:] + bins_dt[:-1])

ax[0].step(binsc, mc_tot * (sum(n_dt) / sum(mc_tot)), where='mid', label='MC')
ax[0].errorbar(binsc_dt, n_dt, poisson_interval(n_dt), marker='.', color='k', linestyle='None', label='Data')
ax[0].set_ylabel('Events / bin width')

ax[1].errorbar(binsc, mc_tot * (sum(n_dt) / sum(mc_tot)) / n_dt, marker='.', color='k', linestyle='None')
#ax[1].set_xlabel(r'$m_{ee}$')
#ax[1].set_xlabel(r'$m_{\mu\mu}$')
ax[1].set_xlabel(r'$m_{ll}$')
ax[1].set_ylabel('MC/Data')
ax[1].set_ylim(0.5, 2)
# Add legend
ax[0].legend()

hep.cms.label(
    llabel="Preliminary",
    data=True,
    lumi=round(lumi, 2),
    ax=ax[0]
)

plt.savefig("plot_2023_preBPix_SF_INCLUSIVE.png")
