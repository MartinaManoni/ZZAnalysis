#!/usr/bin/env python3
###
# Example for running the analysis locally, after customizing variables.
# Run with: 
# python runLocal.py
###
from __future__ import print_function
from ZZAnalysis.NanoAnalysis.tools import setConf, getConf, insertAfter

# Check that the checkout recipe has been properly updated 
from ZZAnalysis.AnalysisStep.validateCheckout import validateCheckout 
if not validateCheckout() :
    exit(1)

#SampleToRun = "MCsync_2018Rereco" # for mini vs nano sync
#SampleToRun = "MCsync_2017UL" # for mini vs nano sync
#SampleToRun = "Data2022"
SampleToRun = "MC2023"
#SampleToRun = "MELA_Test"
#SampleToRun = "ggh125_2018UL"
#SampleToRun = "forNanoDoc" # To prepare variable lists with inspectNanoFile.py


### Customize processing variables.
#setConf("runMELA", False)
#setConf("bestCandByMELA", False)
#setConf("APPLYMUCORR", False)
#setConf("APPLYELECORR", False)

## Force filling K factors and weights (default: all off)
#setConf("APPLY_K_NNLOQCD_ZZGG", 1) # 0:None; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
#setConf("APPLY_K_NNLOQCD_ZZQQB", True)
#setConf("APPLY_K_NNLOEW_ZZQQB", True)
#setConf("APPLY_QCD_GGF_UNCERT", True)

setConf("PROCESS_CR", True)
setConf("PROCESS_ZL", True)
setConf("DEBUG", False)
setConf("SYNCMODE", True) # Force muon resolution correction with fixed +1 sigma smearing
#setConf("ADD_ALLEVENTS", True) # Add extra tree of gen info for all events
#setConf("FILTER_EVENTS", 'Z') # Store all events which contain a good Z candidate
#setConf("FILTER_EVENTS", '3L_20_10') # for trigger studies
#setConf("FILTER_EVENTS", 'NoFilter') # don't skip events with no candidates
#setConf("TRIGPASSTHROUGH", True) #don't skip events failing triggers
#setConf("APPLYJETCORR", False)
#setConf("CANDSTOSTORE",'AllWithRelaxedMuId')

json = None #replace this if needed

################################################################################
if SampleToRun == "Data2022" :
    # 2022 data sample from /MuonEG/Run2022D-22Sep2023-v1/NANOAOD
    setConf("IsMC", False)
    setConf("LEPTON_SETUP", 2022)
    setConf("PD", "any")
    setConf("SAMPLENAME", "test")
    setConf("TRIGPASSTHROUGH", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/27453cd2-36d5-4b34-9cf6-1303480b5dcf.root",
        ])

elif SampleToRun == "Data2023" :
    # 2022 data sample from /MuonEG/Run2023D-22Sep2023_v1-v1/NANOAOD
    setConf("IsMC", False)
    setConf("LEPTON_SETUP", 2023)
    setConf("PD", "any")
    setConf("SAMPLENAME", "test")
    setConf("TRIGPASSTHROUGH", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/data/Run2023D/MuonEG/NANOAOD/22Sep2023_v1-v1/40000/180fca36-4680-4295-8364-bcc292910808.root",
        ])

elif SampleToRun == "Data2024" :
    # 2024 data sample from /EGamma1/Run2024I-MINIv6NANOv15_v2-v1/NANOAOD
    setConf("IsMC", False)
    setConf("LEPTON_SETUP", 2023)
    setConf("PD", "any")
    setConf("DATA_TAG", "post_BPix")
    setConf("NANOVERSION", 15)
    setConf("APPLYELECORR", True)
    setConf("APPLYMUCORR", True)
    setConf("APPLYJETCORR", False)
    setConf("SAMPLENAME", "test")
    setConf("TRIGPASSTHROUGH", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/data/Run2024I/EGamma1/NANOAOD/MINIv6NANOv15_v2-v1/2540000/01acd988-764d-4615-88f3-f385408d25f1.root",
        ])


################################################################################
elif SampleToRun == "ggh125_2018UL" : ### 2018 UL test sample
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("LEPTON_SETUP", 2018)
    setConf("NANOVERSION", 9)    
    setConf("DATA_TAG", "UL")
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/mc/RunIISummer20UL18NanoAODv2/WplusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/270000/3B6A5CB5-2B7C-924D-85B4-FC3B0C1F4909.root",
        ])

################################################################################
elif SampleToRun == "MCsync_2017UL" :
    # Custom-reprocessed Rereco nanoAOD file with updated FSR and electron MVA,
    # no packing for genparticle p3; 26000 events
    # corresponding to:/store/mc/RunIISummer20UL17MiniAODv2/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/130000/3E4E8D55-3993-2B43-AF3B-7AB45BBE0BDA.root
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("LEPTON_SETUP", 2017)
    setConf("NANOVERSION", 10) # variable defined as per nanoAOD v10 (notably electron_mvaHZZIso)
    setConf("DATA_TAG", "UL")
    setConf("store","")
    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_2017UL_fixedFSR.root"])
#    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_2017UL_fixedFSR_nopacking.root"]) # with no packing of muon eta, phi, mass


################################################################################
elif SampleToRun == "MCsync_2018Rereco" :
     # Custom-reprocessed Rereco nanoAOD file with updated FSR,
     # corresponding to:/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root
    setConf("APPLYMUCORR", True)
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("NANOVERSION", 9)
    setConf("store","")
    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_fixedFSR.root"])


################################################################################
elif SampleToRun == "MC2022" :
    # 2022 MC sample
    setConf("SAMPLENAME", "ggH125")
    setConf("DATA_TAG", "post_EE")
    setConf("XSEC", 52.23*0.0002745)
    setConf("LEPTON_SETUP", 2022)
    setConf("IsMC", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("APPLY_QCD_GGF_UNCERT", True) # for ggH
    setConf("fileNames",[
       "/store/mc/Run3Summer22NanoAODv12/GluGluHtoZZto4L_M-124p5_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/28b181a2-3ef5-4ffa-8e70-9eff44bcdc04.root",
        ])

 ################################################################################
#root file from /GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg-jhugen-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v3/NANOAODSIM
elif SampleToRun == "MC2023" :
    # 2022 MC sample
    setConf("SAMPLENAME", "ggH125")
    setConf("DATA_TAG", "pre_BPix")
    setConf("XSEC", 52.23*0.0002745)
    setConf("LEPTON_SETUP", 2023)
    setConf("IsMC", True)
    setConf("APPLYELECORR", True)
    setConf("APPLYMUCORR", True)
    setConf("APPLYJETCORR", False)
    setConf("FILTER_EVENTS", 'NoFilter')
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("APPLY_QCD_GGF_UNCERT", True) # for ggH
    setConf("fileNames",[
    "/store/mc/Run3Summer23NanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg-jhugen-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v3/50000/4cb201c2-3bfd-459f-8b1c-fb54ee1f8d3f.root",
    ])


################################################################################
elif SampleToRun == "forNanoDoc" :
    # Create a file with a complete set of variables to feed to inspectNanoFile to generate variable documentation
    setConf("SAMPLENAME", "ggH125")
    setConf("DATA_TAG", "post_EE")
    setConf("XSEC", 52.23*0.0002745)
    setConf("LEPTON_SETUP", 2022)
    setConf("IsMC", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("runMELA", True)
    setConf("APPLYMUCORR", True)
    setConf("APPLYELECORR", True)
    setConf("APPLYJETCORR", True)
    # setConf("APPLY_K_NNLOQCD_ZZGG", 1) # requires mcHistoryTools before weightFiller when AllEvents=true, which is not needed in practical cases
    # setConf("APPLY_K_NNLOQCD_ZZQQB", True) # ditto
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("APPLY_QCD_GGF_UNCERT", True)
    setConf("PROCESS_CR", True)
    setConf("PROCESS_ZL", True)
    setConf("ADD_ALLEVENTS", True)
    setConf("fileNames",["/store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2540000/25c8f5ff-9de0-4a0c-9e2f-757332ad392f.root"])


###################################################################################
elif SampleToRun == "MELA_Test" : 
    setConf("SAMPLENAME", "ggH125")
    setConf("LEPTON_SETUP", 2022)  
    setConf("XSEC", 290.58626*0.0002745)
    setConf("IsMC", True)
    setConf("ADD_ALLEVENTS", True)
    setConf("NANOVERSION", 15)
    setConf("store", "")
    setConf("fileNames", ["/eos/user/n/nipinto/old_CMSSW_13_3_3/src/ggH_test.root"]) # private reprocessing to add LHE mothers/daughters as in v15
    


#####################################################################
### This import should be done AFTER all customizations (setConf calls)
from ZZAnalysis.NanoAnalysis.nanoZZ4lAnalysis import *
######################################################################

### Tweak postprocessor parameters as necessary
p.prefetch=True # Prefetch remote files
p.longTermCache=True # keep prefetched files (useful for rerunning tests several times)
if len(p.inputFiles) == 1 :
    p.haddFileName = None # Skip final hadd
#p.maxEntries = 10000

### Select specific events to debug
#p.cut = "run==316239  && luminosityBlock==226 && event==284613817"

### Print out detailed candidate information for debug purposes
#from ZZAnalysis.NanoAnalysis.dumpEvents import dumpEvents
#p.cut = None # Remove preselction
#insertAfter(p.modules,"lepFiller",dumpEvents(level=-1),getConf("NANOVERSION", 11)) 

### Dump MC and LHE history for selected events
#from ZZAnalysis.NanoAnalysis.mcHistoryDump import mcHistoryDump
#p.modules.append(mcHistoryDump( printGen=True, printLHE=True))

#p.branchsel=None #Read all branches
#p.outputbranchsel=None #Output all branches

#replace JSON
p.json = json

### Run the postprocessor
p.run()
