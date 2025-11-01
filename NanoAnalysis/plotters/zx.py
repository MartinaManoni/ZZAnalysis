import numpy as np
import pandas as pd
import uproot
from math import sqrt, log
import ROOT

# Configuration
years = ["2022", "2022EE", "2023preBPix", "2023postBPix"]
eos_path_template = '/eos/user/l/lurda/CMS/HZZ/XS_analysis/250303/FAKERATES/{}/'
branches_ZX = ['ZZMass', 'Z1Flav', 'Z2Flav', 'LepLepId', 'LepEta', 'LepPt', 'Z1Mass', 'Z2Mass']


def FindFinalState(z1_flav, z2_flav):
    if z1_flav == -121:
        if z2_flav == 121:
            return 0  # 4e
        if z2_flav == 169:
            return 2  # 2e2mu
    if z1_flav == -169:
        if z2_flav == 121:
            return 3  # 2mu2e
        if z2_flav == 169:
            return 1  # 4mu


def GetFakeRate(lep_Pt, lep_eta, lep_ID, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE):
    my_lep_Pt = min(lep_Pt, 79.)
    my_lep_ID = abs(lep_ID)

    if 5 < my_lep_Pt <= 7:
        bin = 0
    elif 7 < my_lep_Pt <= 10:
        bin = 1
    elif 10 < my_lep_Pt <= 20:
        bin = 2
    elif 20 < my_lep_Pt <= 30:
        bin = 3
    elif 30 < my_lep_Pt <= 40:
        bin = 4
    elif 40 < my_lep_Pt <= 50:
        bin = 5
    elif 50 < my_lep_Pt <= 80:
        bin = 6
    else:
        bin = 6  # fallback

    if my_lep_ID == 11:
        bin = bin - 1  # Electron FR has no bin for [5, 7]
        if abs(lep_eta) < 1.479:
            return g_FR_e_EB.GetY()[bin]
        else:
            return g_FR_e_EE.GetY()[bin]

    if my_lep_ID == 13:
        if abs(lep_eta) < 1.2:
            return g_FR_mu_EB.GetY()[bin]
        else:
            return g_FR_mu_EE.GetY()[bin]


def openFR(year):
    print("OpenFR")
    eos_path = eos_path_template.format(year)
    fnameFR = eos_path + f'FakeRates_SS_{year}.root'

    # Open file with ROOT to access TGraphErrors
    input_file_FR = ROOT.TFile(fnameFR)
    g_FR_mu_EB = input_file_FR.Get("FR_SS_muon_EB")
    g_FR_mu_EE = input_file_FR.Get("FR_SS_muon_EE")
    g_FR_e_EB = input_file_FR.Get("FR_SS_electron_EB")
    g_FR_e_EE = input_file_FR.Get("FR_SS_electron_EE")

    return g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE


def findFSZX(df):
    print("findFSZX")
    df['FinState'] = [FindFinalState(x, y) for x, y in zip(df['Z1Flav'], df['Z2Flav'])]
    return df


def comb(year):
    # cb_SS
    if year == "2022":
        return np.array([
            1.239, 
            1.093, 
            1.057, 
            1.254,
            ])  # 2022preEE
    if year == "2022EE":
        return np.array([
            1.067, # 4e
            1.015, # 4mu
            1.049, # 2e2mu
            0.905, # 2mu2e
            ])  # 2022EE
    if year == "2023preBPix":
        return np.array([
            1.116, 
            1.036, 
            0.989, 
            1.141,
            ])  # 2023preBPix
    else:
        return np.array([
            0.795, # 4e
            1.025, # 4mu
            1.074, # 2e2mu
            1.078, # 2mu2e
            ])  # 2023postBPix


def ratio(year):
    # fs_ROS_SS
    if year == "2022":
        return np.array([
            1.030, 
            1.165, 
            1.057,
            1.254,
            ])  # 2022preEE
    if year == "2022EE":
        return np.array([
            0.990,   # 4e
            0.997,  # 4mu
            1.039,   # 2e2mu
            1.016,  # 2mu2e
            ])  # 2022EE
    if year == "2023preBPix":
        return np.array([
            0.992, 
            1.024, 
            1.102, 
            1.024,
            ])  # 2023preBPix
    else:
        return np.array([
            0.795,   # 4e
            1.040,  # 4mu
            1.074,   # 2e2mu
            1.078,  # 2mu2e
            ])  # 2023postBPix


def ZXYield(df, year, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE):
    print("ZXYield")
    cb_SS = comb(year)
    fs_ROS_SS = ratio(year)

    #print("df", df)

    #vec = df.to_numpy()
    #Yield = np.zeros(len(vec), float)

    #print(len(vec))
    #print("vec", vec)

    Yield = np.zeros(len(df), float)

    for i in range(len(df)):
        row = df.iloc[i]
        finSt = row['FinState']
        lepPt = row['LepPt']   # This is a list/array
        lepEta = row['LepEta'] # list/array
        lepID = row['LepLepId']  # list/array of lepton IDs

        #print("lepPt", lepPt)
        #print("lepPt2", lepPt[2])
        #print("lepEta", lepEta)
        #print("lepID", lepID)

        fr1 = GetFakeRate(lepPt[2], lepEta[2], lepID[2], g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE)
        fr2 = GetFakeRate(lepPt[3], lepEta[3], lepID[3], g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE)

        Yield[i] = cb_SS[finSt] * fs_ROS_SS[finSt] * fr1 * fr2

    return Yield


def doZX(year, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE):
    dir_name = "CRZLLTree" #"CRZLLTree"       # TDirectoryFile
    tree_name = "candTree"       # Actual TTree inside the directory

    # Map year strings to file paths
    data_files = {
        "2022": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2022/Data_eraCD_preEE_SKIMMED_FR_ok.root",
        "2022EE": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2022/Data_eraEFG_postEE_FR_ok.root",
        "2023preBPix": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2023/Data_eraC_preBPix_FR_ok.root",
        "2023postBPix": "/eos/user/m/mmanoni/HZZ_prod_170625/Data/2023/Data_eraD_postBPix_FR_ok.root",
    }

    if year not in data_files:
        raise ValueError(f"Year '{year}' is not recognized. Available options: {list(data_files.keys())}")

    data = data_files[year]

    print(f"Using data file: {data}")

    print("Reading data:", data)

    with uproot.open(data) as file:
        # Navigate to directory first
        if dir_name not in file:
            raise KeyError(f"Directory '{dir_name}' not found in file.")

        subdir = file[dir_name]

        # Now access the TTree inside the subdirectory
        if tree_name not in subdir:
            raise KeyError(f"TTree '{tree_name}' not found inside '{dir_name}'.")

        ttreeZX = subdir[tree_name]

        # Load data as DataFrame
        dfZX = ttreeZX.arrays(branches_ZX, library="pd")
        print("dfZX", dfZX)

    dfZX = dfZX[dfZX.Z2Flav > 0]  # Keep only same-sign events
    #print("dfZX", dfZX)
    dfZX = findFSZX(dfZX)
    dfZX['weight1'] = ZXYield(dfZX, year, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE)

    return dfZX

def save_to_root(df, output_path, tree_name="candTree"):
    # uproot requires arrays, convert DataFrame columns to dict of arrays
    arrays = {col: df[col].values for col in df.columns}
    
    with uproot.recreate(output_path) as root_file:
        root_file[tree_name] = arrays
    print(f"Saved ROOT file to {output_path}")

# ------------------------- Main -------------------------
def zx():
    d_ZX = {}
    for year in years:
        g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE = openFR(year)
        d_ZX[year] = doZX(year, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE)

        # Define output path
        output_file = f"ZX_results_{year}.root"
        # Save the DataFrame with weights to ROOT file
        save_to_root(d_ZX[year], output_file)

        print(year, 'done')
        print("d_ZX", d_ZX[year])
    return d_ZX


# If script is run directly
if __name__ == "__main__":
    zx()
