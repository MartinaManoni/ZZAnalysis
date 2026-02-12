import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# Configuration
# ============================================================

datasets = [
    "WminusH125",
    "WplusH125",
    "ttH125",
    "VBFH125",
    "ggH125",
    "ZH125",
]

variables = [
    "ZZCand_Z1mass",
    "ZZCand_Z2mass",
    "ZZCand_mass",
    "ZCand_mass",
]

tree_name = "Events"
nbins = 50

file_official = (
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/"
    "HIG-25-015/RunIII_byZ1Z2/031125/2023preBPix_MC/"
    "{dataset}/ZZ4lAnalysis.root"
)

file_private = (
    "/eos/user/m/mmanoni/PrivateProd_test2023/"
    "PROD_samplesNano_2023preBPix_MC_0953f30c/"
    "{dataset}/ZZ4lAnalysis.root"
)

outdir = Path("plots_comparison")
outdir.mkdir(exist_ok=True)

# ============================================================
# Helper
# ============================================================

def read_branch(filename, branch):
    """Read and flatten a NanoAOD branch into a 1D NumPy array."""
    with uproot.open(filename) as f:
        arr = f[tree_name][branch].array()
        arr = ak.flatten(arr, axis=None)
        return ak.to_numpy(arr)

# ============================================================
# Main loop
# ============================================================

for dataset in datasets:
    for var in variables:

        print(f"→ {dataset} : {var}")

        # Read data
        data_off = read_branch(file_official.format(dataset=dataset), var)
        data_priv = read_branch(file_private.format(dataset=dataset), var)

        # Remove non-finite values
        data_off = data_off[np.isfinite(data_off)]
        data_priv = data_priv[np.isfinite(data_priv)]

        # Common binning
        xmin = min(data_off.min(), data_priv.min())
        xmax = max(data_off.max(), data_priv.max())
        bins = np.linspace(xmin, xmax, nbins + 1)

        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        bin_width = bins[1] - bins[0]

        # Histograms
        h_off, _ = np.histogram(data_off, bins=bins)
        h_priv, _ = np.histogram(data_priv, bins=bins)

        # Statistical uncertainties (Poisson)
        err_priv = np.sqrt(h_priv)

        # ====================================================
        # Plot: distributions + ratio
        # ====================================================

        fig, (ax_top, ax_bot) = plt.subplots(
            2, 1,
            figsize=(7, 8),
            sharex=True,
            gridspec_kw={"height_ratios": [3, 1]},
        )

        # --- Top panel ---
        ax_top.step(
            bins[:-1],
            h_off,
            where="post",
            linewidth=2,
            label="Official production",
        )

        ax_top.errorbar(
            bin_centers,
            h_priv,
            yerr=err_priv,
            xerr=bin_width / 2,
            fmt="o",
            capsize=2,
            label="Private production",
        )

        ax_top.set_ylabel("Events")
        ax_top.set_title(dataset)
        ax_top.legend()
        ax_top.grid(alpha=0.3)

        # --- Bottom panel (ratio) ---
        mask = h_off > 0

        ratio = np.zeros_like(h_priv, dtype=float)
        ratio_err = np.zeros_like(h_priv, dtype=float)

        ratio[mask] = h_priv[mask] / h_off[mask]
        ratio_err[mask] = err_priv[mask] / h_off[mask]

        ax_bot.errorbar(
            bin_centers[mask],
            ratio[mask],
            yerr=ratio_err[mask],
            xerr=bin_width / 2,
            fmt="o",
            capsize=2,
        )

        ax_bot.axhline(1.0, linestyle="--", color="black")
        ax_bot.set_ylabel("Private / Official")
        ax_bot.set_xlabel(var)
        ax_bot.set_ylim(0.0, 2.0)
        ax_bot.grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(outdir / f"{dataset}_{var}.png")
        plt.close()

print("All plots produced successfully.")


# ============================================================
# Extra plot: GenPart_iso (private production only, 0-1)
# ============================================================

gen_var = "GenPart_iso"
gen_nbins = 50
gen_min, gen_max = 0.0, 1.0  # limit range

for dataset in datasets:

    print(f"→ {dataset} : {gen_var} (private only, 0-1)")

    # Read GenPart_iso from private file
    data_gen = read_branch(
        file_private.format(dataset=dataset),
        gen_var
    )

    # Remove non-finite values
    data_gen = data_gen[np.isfinite(data_gen)]

    # Keep only values in [0, 1]
    data_gen = data_gen[(data_gen >= gen_min) & (data_gen <= gen_max)]

    if len(data_gen) == 0:
        print(f"  [WARNING] No entries in [0,1] for {dataset} {gen_var}")
        continue

    # Binning
    bins = np.linspace(gen_min, gen_max, gen_nbins + 1)

    # Histogram
    h_gen, _ = np.histogram(data_gen, bins=bins)

    # Plot
    plt.figure(figsize=(7, 6))
    plt.step(
        bins[:-1],
        h_gen,
        where="post",
        linewidth=2,
        label="Private production",
    )

    plt.yscale("log")

    plt.xlabel(gen_var)
    plt.ylabel("Events")
    plt.title(f"{dataset} – {gen_var} (private, 0-1)")
    plt.legend()
    plt.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(outdir / f"{dataset}_{gen_var}_private_0to1.png")
    plt.close()
