import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

# === CONFIG ===
ROOT_FILE = "/eos/user/m/mmanoni/test_BDT/PROD_samplesNano_2022EE_MC/ggH125/ZZ4lAnalysis.root"
TREE_NAME = "Events"
PT_BINS = np.linspace(0, 200, 21)
PLOT_FILENAME = "bdt_efficiency_comparison_2022EE.png"

# Choose which BDT variables to plot (set True/False)
BDT_VARS_TO_PLOT = {
    "Electron_passBDT": True,
    "Electron_mvaIso_WP80": True,
    "Electron_mvaIso_WP90": True,
    "Electron_mvaNoIso_WP80": True,
    "Electron_mvaNoIso_WP90": True
}

# === Open ROOT file and read branches ===
base_branches = [
    "GenPart_eta", "GenPart_phi", "GenPart_pt", "GenPart_pdgId", "GenPart_genPartIdxMother",
    "Electron_eta", "Electron_phi", "Electron_pt", "Electron_genPartIdx"
]
bdt_branches = [key for key, val in BDT_VARS_TO_PLOT.items() if val]
branches = base_branches + bdt_branches

print(f"Loading file: {ROOT_FILE}")
file = uproot.open(ROOT_FILE)
tree = file[TREE_NAME]
events = tree.arrays(branches)
print("File loaded successfully.\n")

# === Initialize counters ===
bdt_efficiencies = {bdt: {"pt_total": [], "pt_pass": [], "count_passed": 0} for bdt in bdt_branches}
total_from_Z = 0

# === Loop over events ===
for event in events:
    gen_eta     = event["GenPart_eta"]
    gen_phi     = event["GenPart_phi"]
    gen_pt      = event["GenPart_pt"]
    gen_pdg     = event["GenPart_pdgId"]
    gen_mother  = event["GenPart_genPartIdxMother"]

    ele_eta     = event["Electron_eta"]
    ele_phi     = event["Electron_phi"]
    ele_pt      = event["Electron_pt"]
    ele_genidx  = event["Electron_genPartIdx"]

    for i in range(len(ele_eta)):
        gen_idx = ele_genidx[i]
        if gen_idx < 0 or gen_idx >= len(gen_pdg):
            continue

        mother_idx = gen_mother[gen_idx] if gen_idx >= 0 else -1
        grandmother_idx = gen_mother[mother_idx] if mother_idx >= 0 and mother_idx < len(gen_mother) else -1

        comes_from_Z = False
        if mother_idx >= 0 and abs(gen_pdg[mother_idx]) == 23:
            comes_from_Z = True
        elif grandmother_idx >= 0 and abs(gen_pdg[grandmother_idx]) == 23:
            comes_from_Z = True

        if not comes_from_Z:
            continue

        total_from_Z += 1

        for bdt_var in bdt_branches:
            passed = event[bdt_var][i]
            bdt_efficiencies[bdt_var]["pt_total"].append(ele_pt[i])
            if passed:
                bdt_efficiencies[bdt_var]["pt_pass"].append(ele_pt[i])
                bdt_efficiencies[bdt_var]["count_passed"] += 1

# === Plotting with Binomial Error Bars ===
plt.figure(figsize=(10, 6))

for bdt_var, data in bdt_efficiencies.items():
    total_hist, _ = np.histogram(data["pt_total"], bins=PT_BINS)
    pass_hist, _ = np.histogram(data["pt_pass"], bins=PT_BINS)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        eff = np.divide(pass_hist, total_hist, out=np.zeros_like(pass_hist, dtype=float), where=total_hist != 0)
        err = np.sqrt(eff * (1 - eff) / total_hist)
    
    bin_centers = 0.5 * (PT_BINS[:-1] + PT_BINS[1:])
    
    plt.errorbar(
        bin_centers, eff, yerr=err,
        fmt='o', capsize=3, label=f"{bdt_var} (eff: {data['count_passed'] / total_from_Z:.3f})"
    )

# === Plot styling ===
plt.xlabel("Reco Electron $p_T$ [GeV]")
plt.ylabel("Efficiency")
plt.title("Electron BDT Efficiency vs Reco $p_T$")
plt.grid(True)
plt.ylim(0, 1.05)
plt.legend()
plt.tight_layout()
plt.savefig(PLOT_FILENAME)
print(f"Plot saved as: {PLOT_FILENAME}\n")

# === Summary Printout ===
print("=== Summary (Reco-Based, From Z) ===")
print(f"Total matched reco electrons from Z: {total_from_Z}")
for bdt_var, data in bdt_efficiencies.items():
    eff = data["count_passed"] / total_from_Z if total_from_Z > 0 else 0
    print(f"{bdt_var:25s} => Passed: {data['count_passed']:5d} | Efficiency: {eff:.4f}")
