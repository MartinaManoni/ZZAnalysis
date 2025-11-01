import numpy as np
import pandas as pd
import uproot3 as uproot
from math import sqrt, log
import ROOT
# from config import *

years = [2016, 2017, 2018]
# years = [2017]
eos_path = '/eos/user/a/atarabin/FRfiles_UL/'
branches_ZX = ['ZZMass', 'Z1Flav', 'Z2Flav', 'LepLepId', 'LepEta', 'LepPt', 'Z1Mass', 'Z2Mass', 'ZZPt', 'ZZEta', 
               'helcosthetaZ1','helcosthetaZ2', 'helphi', 'costhetastar', 'phistarZ1', 'ZZPhi',
               'pTj1', 'pTj2', 'absdetajj', 'TCjmax', 'TBjmax',
               'pTHj', 'pTHjj', 'mHj', 'mHjj', 'detajj', 'dphijj', 'mjj', 'njets_pt30_eta4p7', 'njets_pt30_eta2p5', 'ZZy',
               'D0m', 'Dcp', 'D0hp', 'Dint', 'DL1', 'DL1int', 'DL1Zg', 'DL1Zgint', 'Dbkg', 'Dbkg_kin']


def FindFinalState(z1_flav, z2_flav):
    if(z1_flav == -121):
        if(z2_flav == +121): return 0 # 4e
        if(z2_flav == +169): return 2 # 2e2mu
    if(z1_flav == -169):
        if(z2_flav == +121): return 3 # 2mu2e
        if(z2_flav == +169): return 1 # 4mu
        

def GetFakeRate(lep_Pt, lep_eta, lep_ID, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE):
    
    if(lep_Pt >= 80.): 
        my_lep_Pt = 79.
    else:
        my_lep_Pt = lep_Pt
        
    my_lep_ID = abs(lep_ID)
    
    if((my_lep_Pt > 5) & (my_lep_Pt <= 7)): bin = 0
    if((my_lep_Pt >  7) & (my_lep_Pt <= 10)): bin = 1
    if((my_lep_Pt > 10) & (my_lep_Pt <= 20)): bin = 2
    if((my_lep_Pt > 20) & (my_lep_Pt <= 30)): bin = 3
    if((my_lep_Pt > 30) & (my_lep_Pt <= 40)): bin = 4
    if((my_lep_Pt > 40) & (my_lep_Pt <= 50)): bin = 5
    if((my_lep_Pt > 50) & (my_lep_Pt <= 80)): bin = 6
        
    if(abs(my_lep_ID) == 11): bin = bin-1 # There is no [5, 7] bin in the electron fake rate
    
    if(my_lep_ID == 11):
        if(abs(lep_eta) < 1.479): return g_FR_e_EB.GetY()[bin]
        else: return g_FR_e_EE.GetY()[bin]
    
    if(my_lep_ID == 13):
        if(abs(lep_eta) < 1.2): return g_FR_mu_EB.GetY()[bin]
        else: return g_FR_mu_EE.GetY()[bin]
        

# Open Fake Rates files
def openFR(year):
    fnameFR = eos_path + 'FakeRates_SS_%i.root' %year
    file = uproot.open(fnameFR)
    
    # Retrieve FR from TGraphErrors
    input_file_FR = ROOT.TFile(fnameFR)

    g_FR_mu_EB = input_file_FR.Get("FR_SS_muon_EB")
    g_FR_mu_EE = input_file_FR.Get("FR_SS_muon_EE")
    g_FR_e_EB  = input_file_FR.Get("FR_SS_electron_EB")
    g_FR_e_EE  = input_file_FR.Get("FR_SS_electron_EE")
    
    return g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE
    

# Find final state
def findFSZX(df):
    df['FinState'] = [FindFinalState(x,y) for x,y in zip(df['Z1Flav'], df['Z2Flav'])]
    return df


# Define combination coefficients
def comb(year):
    if year == 2016:
        cb_SS = np.array([
            1.175,   # 4e
            0.975,   # 4mu
            1.052,    # 2e2mu
            1.148,    # 2mu2e
        ])
    elif year == 2017:
        cb_SS = np.array([
            1.094,   # 4e
            0.948,   # 4mu
            0.930,    # 2e2mu
            1.139,    # 2mu2e
        ])   
    else:
        cb_SS = np.array([
            1.157,   # 4e
            0.974,   # 4mu
            0.930,    # 2e2mu
            1.143,    # 2mu2e
        ])
        
    return cb_SS


# Define ration OppositeSign/SameSign
def ratio(year):
    if year == 2016:
        fs_ROS_SS = np.array([
            1.0039,   # 4e
            0.999103,  # 4mu
            1.0332,   # 2e2mu
            1.00216,  # 2mu2e
            ])
    elif year == 2017:
        fs_ROS_SS = np.array([
            0.990314,   # 4e
            1.02903,  # 4mu
            1.0262,   # 2e2mu
            1.00154,  # 2mu2e
            ])      
    else: 
        fs_ROS_SS = np.array([
            1.00322,   # 4e
            1.0187,  # 4mu
            1.04216,   # 2e2mu
            0.996253,  # 2mu2e
            ])
        
    return fs_ROS_SS


# Calculate yield for Z+X (data in CRZLL control region are scaled in signal region through yields)
def ZXYield(df, year, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE):
    cb_SS = comb(year)
    fs_ROS_SS = ratio(year)
    
    vec = df.to_numpy()
    Yield = np.zeros(len(vec), float) 
    for i in range(len(vec)):
        finSt  = vec[i][len(branches_ZX)] #Final state information is in the last column which is added afterthe last column of branches_ZX
        lepPt  = vec[i][5]
        lepEta = vec[i][4]
        lepID  = vec[i][3]
        Yield[i] = cb_SS[finSt] * fs_ROS_SS[finSt] * GetFakeRate(lepPt[2], lepEta[2], lepID[2], g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE) * GetFakeRate(lepPt[3], lepEta[3], lepID[3], g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE)

    return Yield


def doZX(year, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE):
    keyZX = 'CRZLL'
    
    data = '/eos/user/a/atarabin/Data/reducedTree_AllData_%i.root' %year
    ttreeZX = uproot.open(data)[keyZX]    

    dfZX = ttreeZX.pandas.df(branches_ZX, flatten = False)
    dfZX = dfZX[dfZX.Z2Flav > 0] #Keep just same-sign events
    dfZX = findFSZX(dfZX)
     
    dfZX['weight1'] = ZXYield(dfZX, year, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE)
    return dfZX

#------------------------- Main -------------------------
def zx(): 
    d_ZX = {}
    for year in years:
        g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE = openFR(year)
        d_ZX[year] = doZX(year, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE)
#         d_ZX[year] = add_lead_lep(d_ZX[year])
#         d_ZX[year] = add_rapidity(d_ZX[year])
        print(year, 'done')
    return d_ZX