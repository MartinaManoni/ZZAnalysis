import ROOT
import numpy as np
import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT()

ROOT_FILE = "/eos/user/m/mmanoni/test_BDT/PROD_samplesNano_2022EE_MC_fd6e1fef/ggH125/ZZ4lAnalysis.root"
TREE_NAME = "Events"
PT_BINS = np.linspace(0, 200, 21)
PLOT_FILENAME = "bdt_efficiency_comparison_2022EE_FINAL.png"

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

# === Define Electron_pt_fromZ with baseline filtering ===
df = df.Define("Electron_pt_fromZ", """
    std::vector<float> filtered;
    for (unsigned int i = 0; i < Electron_pt.size(); ++i) {
        // Apply baseline selections
        if (!(Electron_pt[i] > 7 &&
              fabs(Electron_eta[i]) < 2.5 &&
              fabs(Electron_dxy[i]) < 0.5 &&
              fabs(Electron_dz[i]) < 1.0 &&
              fabs(Electron_sip3d[i]) < 4.0))
            continue;

        int genIdx = Electron_genPartIdx[i];
        if (genIdx < 0 || genIdx >= GenPart_pdgId.size()) continue;

        int motherIdx = GenPart_genPartIdxMother[genIdx];
        if (motherIdx < 0 || motherIdx >= GenPart_pdgId.size()) continue;

        int motherPdg = GenPart_pdgId[motherIdx];

        bool fromZ = false;
        if (abs(motherPdg) == 23) {
            fromZ = true;
        } else if (abs(motherPdg) == 11) {
            int grandmotherIdx = GenPart_genPartIdxMother[motherIdx];
            if (grandmotherIdx >= 0 && grandmotherIdx < GenPart_pdgId.size()) {
                int grandmotherPdg = GenPart_pdgId[grandmotherIdx];
                if (abs(grandmotherPdg) == 23)
                    fromZ = true;
            }
        }

        if (fromZ)
            filtered.push_back(Electron_pt[i]);
    }
    return filtered;
""")

# === Total histogram for electrons from Z passing baseline ===
h_total = df.Histo1D(("total", "Electrons from Z (baseline);Reco Electron p_{T} [GeV];Entries", len(PT_BINS)-1, PT_BINS),
                     "Electron_pt_fromZ")

# === Loop over BDT vars, define filtered passing pt vectors with baseline ===
efficiencies = {}
for bdt_var, enabled in BDT_VARS_TO_PLOT.items():
    if not enabled:
        continue

    df = df.Define(f"{bdt_var}_fromZ_passed", f"""
        std::vector<float> filtered;
        for (unsigned int i = 0; i < Electron_pt.size(); ++i) {{
            // Apply baseline selections
            if (!(Electron_pt[i] > 7 &&
                fabs(Electron_eta[i]) < 2.5 &&
                fabs(Electron_dxy[i]) < 0.5 &&
                fabs(Electron_dz[i]) < 1.0 &&
                fabs(Electron_sip3d[i]) < 4.0))
                continue;

            int genIdx = Electron_genPartIdx[i];
            if (genIdx < 0 || genIdx >= GenPart_pdgId.size()) continue;

            int motherIdx = GenPart_genPartIdxMother[genIdx];
            if (motherIdx < 0 || motherIdx >= GenPart_pdgId.size()) continue;

            int motherPdg = GenPart_pdgId[motherIdx];

            bool fromZ = false;
            if (abs(motherPdg) == 23) {{
                fromZ = true;
            }} else if (abs(motherPdg) == 11) {{
                int grandmotherIdx = GenPart_genPartIdxMother[motherIdx];
                if (grandmotherIdx >= 0 && grandmotherIdx < GenPart_pdgId.size()) {{
                    int grandmotherPdg = GenPart_pdgId[grandmotherIdx];
                    if (abs(grandmotherPdg) == 23)
                        fromZ = true;
                }}
            }}

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
plt.title("Electron BDT Efficiency vs Reco $p_T$ (Era 2022EE)")
plt.grid(True)
plt.ylim(0, 1.05)
plt.legend()
plt.tight_layout()
plt.savefig(PLOT_FILENAME)
print(f"\nPlot saved as: {PLOT_FILENAME}")

# === Summary ===
print("\n=== Summary (Reco-Based, From Z, Baseline applied) ===")
print(f"Total matched reco electrons from Z (baseline): {int(np.sum(total_counts))}")
for bdt_var, hist in efficiencies.items():
    h_root = hist.GetValue()
    passed_counts = np.array([h_root.GetBinContent(i) for i in range(1, nbins+1)], dtype=float)
    passed_total = int(np.sum(passed_counts))
    eff = passed_total / int(np.sum(total_counts)) if total_matched > 0 else 0
    print(f"{bdt_var:25s} => Passed: {passed_total:5d} | Efficiency: {eff:.4f}")