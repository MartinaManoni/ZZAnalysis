import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

file_name = "25c8f5ff-9de0-4a0c-9e2f-757332ad392f_Skim.root"
tree_name = "Events"

with uproot.open(file_name) as f:
    tree = f[tree_name]

    # Read branches as awkward arrays
    SF    = tree["Electron_dataMC"].array()
    dataMC    = tree["Electron_dataMCUnc"].array()
    RECO_stat = tree["Electron_RECO_statUnc"].array()
    RECO_syst = tree["Electron_RECO_systUnc"].array()
    ID_stat   = tree["Electron_ID_statUnc"].array()
    ID_syst   = tree["Electron_ID_systUnc"].array()

    # Compute total from components (keeping your original formula)
    total_from_components = np.sqrt(RECO_stat + RECO_syst + ID_stat + ID_syst)
    #error = total_from_components / SF

    # Print side by side per electron per event
    for i in range(len(dataMC)):
        print(f"Event {i}:")
        for j, (sf, dmc, tot, RECO_stat_, RECO_syst_, ID_stat_, ID_syst_) in enumerate(
            zip(SF[i], dataMC[i], total_from_components[i], RECO_stat[i], RECO_syst[i], ID_stat[i], ID_syst[i])
        ):
            print(f"  Electron {j:2d} | SF = {sf:.5f}| SFUnc/SF_original = {dmc:.5f} |"
                  f" SFUnc/SF = {tot:.5f} | RECO_stat={RECO_stat_:.5f} | RECO_syst={RECO_syst_:.5f} |"
                  f" ID_stat={ID_stat_:.5f} | ID_syst={ID_syst_:.5f}")

    # Mask electrons with SF != 1
    mask = SF != 1
    dataMC_masked = dataMC[mask]
    error_masked = total_from_components[mask]

    # Compute relative difference
    rel_diff_masked = abs(dataMC_masked - error_masked) / dataMC_masked

    # Flatten all electrons across all events to compute average
    rel_diff_flat = ak.flatten(rel_diff_masked)
    avg_rel_diff = np.mean(rel_diff_flat)

    print(f"Average relative difference (excluding SF=1): {avg_rel_diff:.6f}")

    # Plot
    plt.figure()
    plt.hist(rel_diff_flat, bins=50)
    plt.xlabel(r"$| \Delta | = |dataMCUnc - computed| / dataMCUnc$")
    plt.ylabel("Number of electrons")
    plt.title("Relative difference distribution (excluding SF = 1)")
    plt.savefig("relative_difference_distribution.png", dpi=300, bbox_inches="tight")

