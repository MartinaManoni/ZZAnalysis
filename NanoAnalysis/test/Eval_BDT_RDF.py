import ROOT
import numpy as np
import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT()

ROOT_FILE = "/eos/user/m/mmanoni/test_BDT/PROD_samplesNano_2023preBPix_MC/ggH125/ZZ4lAnalysis.root"
TREE_NAME = "Events"
PT_BINS = np.linspace(0, 200, 21)
PLOT_FILENAME = "bdt_efficiency_comparison_2023preBPix.png"

BDT_VARS_TO_PLOT = {
    "Electron_passBDT": True,
    "Electron_mvaIso_WP80": True,
    "Electron_mvaIso_WP90": True,
    "Electron_mvaNoIso_WP80": True,
    "Electron_mvaNoIso_WP90": True
}

# === INIT DATAFRAME ===
df = ROOT.RDataFrame(TREE_NAME, ROOT_FILE)

# === Redefine nElectron if needed ===
try:
    df = df.Redefine("nElectron", "Electron_pt.size()")
except RuntimeError:
    pass  # If not needed or errors, continue

# === Define Electron_pt_fromZ by manual filtering ===
df = df.Define("Electron_pt_fromZ", """
    std::vector<float> filtered;
    for (unsigned int i = 0; i < Electron_pt.size(); ++i) {
        int genIdx = Electron_genPartIdx[i];
        if (genIdx < 0 || genIdx >= GenPart_pdgId.size()) continue;
        int motherIdx = GenPart_genPartIdxMother[genIdx];
        int grandmotherIdx = (motherIdx >= 0) ? GenPart_genPartIdxMother[motherIdx] : -1;
        int motherPdg = (motherIdx >= 0) ? GenPart_pdgId[motherIdx] : 0;
        int grandmotherPdg = (grandmotherIdx >= 0) ? GenPart_pdgId[grandmotherIdx] : 0;
        if (abs(motherPdg) == 23 || abs(grandmotherPdg) == 23)
            filtered.push_back(Electron_pt[i]);
    }
    return filtered;
""")

# === Total histogram for electrons from Z ===
h_total = df.Histo1D(("total", "Electrons from Z;Reco Electron p_{T} [GeV];Entries", len(PT_BINS)-1, PT_BINS),
                     "Electron_pt_fromZ")

# === Loop over BDT vars, define filtered passing pt vectors ===
efficiencies = {}
for bdt_var, enabled in BDT_VARS_TO_PLOT.items():
    if not enabled:
        continue

    df = df.Define(f"{bdt_var}_fromZ_passed", f"""
        std::vector<float> filtered;
        for (unsigned int i = 0; i < Electron_pt.size(); ++i) {{
            int genIdx = Electron_genPartIdx[i];
            if (genIdx < 0 || genIdx >= GenPart_pdgId.size()) continue;
            int motherIdx = GenPart_genPartIdxMother[genIdx];
            int grandmotherIdx = (motherIdx >= 0) ? GenPart_genPartIdxMother[motherIdx] : -1;
            int motherPdg = (motherIdx >= 0) ? GenPart_pdgId[motherIdx] : 0;
            int grandmotherPdg = (grandmotherIdx >= 0) ? GenPart_pdgId[grandmotherIdx] : 0;
            bool fromZ = (abs(motherPdg) == 23 || abs(grandmotherPdg) == 23);
            if (fromZ && {bdt_var}[i]) {{
                filtered.push_back(Electron_pt[i]);
            }}
        }}
        return filtered;
    """)

    h_pass = df.Histo1D((bdt_var, bdt_var + ";Reco Electron p_{T} [GeV];Entries", len(PT_BINS)-1, PT_BINS),
                       f"{bdt_var}_fromZ_passed")
    efficiencies[bdt_var] = h_pass

# === Extract bin contents properly (skip underflow/overflow) ===
h_total_root = h_total.GetValue()
nbins = h_total_root.GetNbinsX()
total_counts = np.array([h_total_root.GetBinContent(i) for i in range(1, nbins+1)], dtype=float)

plt.figure(figsize=(10, 6))
bin_centers = 0.5 * (PT_BINS[:-1] + PT_BINS[1:])

for bdt_var, hist in efficiencies.items():
    h_root = hist.GetValue()
    passed_counts = np.array([h_root.GetBinContent(i) for i in range(1, nbins+1)], dtype=float)

    with np.errstate(divide='ignore', invalid='ignore'):
        eff = np.divide(passed_counts, total_counts, out=np.zeros_like(passed_counts), where=total_counts != 0)
        err = np.sqrt(eff * (1 - eff) / total_counts)

    total_matched = int(np.sum(total_counts))
    passed_total = int(np.sum(passed_counts))

    plt.errorbar(
        bin_centers, eff, yerr=err, fmt='o', capsize=3,
        label=f"{bdt_var} (eff: {passed_total / total_matched:.3f})"
    )

plt.xlabel("Reco Electron $p_T$ [GeV]")
plt.ylabel("Efficiency")
plt.title("Electron BDT Efficiency vs Reco $p_T$ (Era 2022)")
plt.grid(True)
plt.ylim(0, 1.05)
plt.legend()
plt.tight_layout()
plt.savefig(PLOT_FILENAME)
print(f"\nPlot saved as: {PLOT_FILENAME}")

# === Summary ===
print("\n=== Summary (Reco-Based, From Z) ===")
print(f"Total matched reco electrons from Z: {int(np.sum(total_counts))}")
for bdt_var, hist in efficiencies.items():
    h_root = hist.GetValue()
    passed_counts = np.array([h_root.GetBinContent(i) for i in range(1, nbins+1)], dtype=float)
    passed_total = int(np.sum(passed_counts))
    eff = passed_total / int(np.sum(total_counts)) if total_matched > 0 else 0
    print(f"{bdt_var:25s} => Passed: {passed_total:5d} | Efficiency: {eff:.4f}")
