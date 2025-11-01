import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

# === CONFIG ===
ROOT_FILE = "BDTWP80.root"
TREE_NAME = "Events"
PT_BINS = np.linspace(0, 200, 21)
PLOT_FILENAME = "bdt_efficiency_reco_WP80_NEW.png"

# === Open ROOT file and read branches ===
branches = [
    "GenPart_eta", "GenPart_phi", "GenPart_pt", "GenPart_pdgId", "GenPart_genPartIdxMother",
    "Electron_eta", "Electron_phi", "Electron_pt", "Electron_passBDT", "Electron_genPartIdx"
]

print(f"Loading file: {ROOT_FILE}")
file = uproot.open(ROOT_FILE)
tree = file[TREE_NAME]
events = tree.arrays(branches)
print("File loaded successfully.\n")

# === Initialize counters ===
pt_total = []
pt_pass = []

total_from_Z = 0
matched_electrons = 0
passed_bdt = 0

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
    ele_pass    = event["Electron_passBDT"]
    ele_genidx  = event["Electron_genPartIdx"]
    #Electron_mvaIso_WP80
    #Electron_mvaIso_WP90
    #Electron_mvaNoIso_WP80
    #Electron_mvaNoIso_WP90

    for i in range(len(ele_eta)):
        gen_idx = ele_genidx[i]
        if gen_idx < 0 or gen_idx >= len(gen_pdg):
            continue

        # Walk up to mother and grandmother
        mother_idx = gen_mother[gen_idx] if gen_idx >= 0 else -1
        grandmother_idx = gen_mother[mother_idx] if mother_idx >= 0 and mother_idx < len(gen_mother) else -1

        # Check if mother or grandmother is a Z boson
        comes_from_Z = False
        if mother_idx >= 0 and abs(gen_pdg[mother_idx]) == 23:
            comes_from_Z = True
        elif grandmother_idx >= 0 and abs(gen_pdg[grandmother_idx]) == 23:
            comes_from_Z = True

        if not comes_from_Z:
            continue

        total_from_Z += 1
        pt_total.append(ele_pt[i])
        if ele_pass[i]:
            pt_pass.append(ele_pt[i])
            passed_bdt += 1

# === Efficiency histogram ===
total_hist, _ = np.histogram(pt_total, bins=PT_BINS)
pass_hist, _ = np.histogram(pt_pass, bins=PT_BINS)

eff = np.divide(pass_hist, total_hist, out=np.zeros_like(pass_hist, dtype=float), where=total_hist != 0)
bin_centers = 0.5 * (PT_BINS[:-1] + PT_BINS[1:])

# === Print Summary ===
print("=== Summary (Reco-Based) ===")
print(f"Reco electrons matched to gen from Z: {total_from_Z}")
print(f"Passed Electron_passBDT:             {passed_bdt}")
overall_eff = passed_bdt / total_from_Z if total_from_Z > 0 else 0
print(f"Overall BDT Efficiency:              {overall_eff:.4f}")
print()

# === Plot ===
plt.figure(figsize=(8, 5))
plt.step(bin_centers, eff, where="mid", color="green", label="BDT Efficiency")
plt.xlabel("Reco Electron $p_T$ [GeV]")
plt.ylabel("Efficiency (Electron_passBDT)")
plt.title("Electron BDT Efficiency vs Reco $p_T$")
plt.grid(True)
plt.ylim(0, 1.05)
plt.legend()
plt.tight_layout()

plt.savefig(PLOT_FILENAME)
print(f"Plot saved as: {PLOT_FILENAME}")


