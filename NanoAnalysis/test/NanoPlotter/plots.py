import uproot
import matplotlib.pyplot as plt
import awkward as ak

# Input ROOT files and tree name
file1 = "/eos/user/m/mmanoni/HZZ_prod_170625/MC/2023preBPix/ggH125/ZZ4lAnalysis.root"
file2 = "/afs/cern.ch/user/m/mmanoni/NanoProd/tmp/output/nano_0.root"
tree_name = "Events"

# Open files
t1 = uproot.open(file1)[tree_name]
t2 = uproot.open(file2)[tree_name]

# --- File 1: FidDressedLeps_RelIso ---
fid_dressed_leps = t1["FidDressedLeps_RelIso"].array(library="ak")
fid_dressed_leps = ak.flatten(fid_dressed_leps)   # già leptoni, quindi niente filtro pdgId necessario

# --- File 2: GenPart_iso con selezione leptoni ---
pdg2 = t2["GenPart_pdgId"].array(library="ak")
iso2 = t2["GenPart_iso"].array(library="ak")

# mask: solo muoni (±13) ed elettroni (±11)
mask_leptons = (abs(pdg2) == 11) | (abs(pdg2) == 13) #abs(pdg2) == 11) | (abs(pdg2) == 13
genpart_iso_leptons = ak.flatten(iso2[mask_leptons])

# --- Plot ---
plt.figure(figsize=(8,6))
#plt.hist(fid_dressed_leps, bins=1000, range=(50,1000), alpha=0.5, label="FidDressedLeps_RelIso (file1)")
plt.hist(genpart_iso_leptons, bins=100, range=(-2,100), alpha=0.5, label="GenPart_iso leptons (file2)")

plt.xlabel("Relative isolation")
plt.ylabel("Entries")
plt.yscale("log")  
plt.title("GenPart_iso") #omparison: FidDressedLeps_RelIso vs GenPart_iso (leptons only
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("genIso_leptons.png")
plt.close()

