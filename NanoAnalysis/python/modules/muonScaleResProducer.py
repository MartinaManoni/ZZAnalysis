import os

# Set up the NATModules muonScaleRes module
def getMuonScaleRes(era, tag, is_mc, overwritePt=True) :
#FIXME: module to be moved in NATModules
#    from PhysicsTools.NATModules.modules.muonScaleRes import muonScaleRes 
    from ZZAnalysis.NanoAnalysis.modules.muonScaleRes import muonScaleRes    

    if era not in [2022, 2023]:  # Add support for 2023
        raise ValueError(f"getMuonScaleRes: Era {era} is not supported")

    if era == 2022:
        if "pre_EE" in tag :
            fname = "2022_schemaV2.json.gz"
        else :
            fname = "2022EE_schemaV2.json.gz"
    elif era == 2023:
        if "pre_BPix" in tag:
            fname = "2023_schemaV2.json.gz"
        else:
            fname = "2023BPix_schemaV2.json.gz"

 
    json = "%s/src/ZZAnalysis/NanoAnalysis/data/MuonScale/%s" % (os.environ['CMSSW_BASE'], fname)

    print("***muonScaleRes: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "json:", json)
    return muonScaleRes(json, is_mc, overwritePt)

