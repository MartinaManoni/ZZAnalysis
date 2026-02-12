'''Add data/MC weights to leptons.
'''

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ROOT import LeptonSFHelper


class lepDataMCWeight(Module):
    def __init__(self, year, data_tag, muonIdByMVA = False):
        '''Add data/MC weights to leptons.'''
        
        print("***lepDataMCWeight: year:", year, "data_tag:", data_tag, "muonIdByMVA:", muonIdByMVA, flush=True)
        self.year = year
        self.lepSFHelper = LeptonSFHelper(year, data_tag+("MUON_ID_BYMVA" if muonIdByMVA else "")) # Squeeze muonIdByMVA in data_tag to avoid changing the C++ interface

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # Original SF branches
        self.out.branch("Muon_dataMC", "F", lenVar="nMuon", title="data/MC correction", limitedPrecision=12)
        self.out.branch("Muon_dataMCUnc", "F", lenVar="nMuon", title="data/MC correction relative uncertainty", limitedPrecision=12)
        self.out.branch("Electron_dataMC", "F", lenVar="nElectron", title="data/MC correction", limitedPrecision=12)
        self.out.branch("Electron_dataMCUnc", "F", lenVar="nElectron", title="data/MC correction relative uncertainty", limitedPrecision=12)

        # Decorrelated uncertainty branches (electrons only)
        self.out.branch("Electron_RECO_statUnc", "F", lenVar="nElectron", title="electron RECO statistical uncertainty", limitedPrecision=12)
        self.out.branch("Electron_RECO_systUnc", "F", lenVar="nElectron", title="electron RECO systematic uncertainty", limitedPrecision=12)
        self.out.branch("Electron_ID_statUnc", "F", lenVar="nElectron", title="electron ID statistical uncertainty", limitedPrecision=12)
        self.out.branch("Electron_ID_systUnc", "F", lenVar="nElectron", title="electron ID systematic uncertainty", limitedPrecision=12)

        self.out.branch("Electron_dataMC_RECO", "F", lenVar="nElectron", title="electron ID statistical uncertainty", limitedPrecision=12)
        self.out.branch("Electron_dataMC_ID", "F", lenVar="nElectron", title="electron ID systematic uncertainty", limitedPrecision=12)


    def analyze(self, event):
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")

        # --- Original SFs
        e_SFs = [1.]*event.nElectron
        e_SFsUnc = [1.]*event.nElectron
        e_SFsRECO= [1.]*event.nElectron
        e_SFsID= [1.]*event.nElectron
        for ie, ele in enumerate(electrons):
            e_SFs[ie], e_SFsUnc[ie],e_SFsRECO[ie],e_SFsID[ie] = self.getLepSF(ele)

        m_SFs = [1.]*event.nMuon
        m_SFsUnc = [1.]*event.nMuon
        m_SFsRECO = [1.]*event.nMuon
        m_SFsID = [1.]*event.nMuon
        for im, mu in enumerate(muons):
            m_SFs[im], m_SFsUnc[im], m_SFsRECO[im],m_SFsID[im]  = self.getLepSF(mu)

        self.out.fillBranch("Electron_dataMC", e_SFs)
        self.out.fillBranch("Electron_dataMC_RECO", e_SFsRECO)  
        self.out.fillBranch("Electron_dataMC_ID", e_SFsID)    
        self.out.fillBranch("Electron_dataMCUnc", e_SFsUnc)
        self.out.fillBranch("Muon_dataMC", m_SFs)
        self.out.fillBranch("Muon_dataMCUnc", m_SFsUnc)

        # --- Decorrelated electron uncertainties
        e_RECO_stat = [0.]*event.nElectron
        e_RECO_syst = [0.]*event.nElectron
        e_ID_stat   = [0.]*event.nElectron
        e_ID_syst   = [0.]*event.nElectron

        for ie, ele in enumerate(electrons):
            reco_stat, reco_syst, id_stat, id_syst = self.getLepSF_decorr(ele)
            e_RECO_stat[ie] = reco_stat
            e_RECO_syst[ie] = reco_syst
            e_ID_stat[ie]   = id_stat
            e_ID_syst[ie]   = id_syst

        self.out.fillBranch("Electron_RECO_statUnc", e_RECO_stat)
        self.out.fillBranch("Electron_RECO_systUnc", e_RECO_syst)
        self.out.fillBranch("Electron_ID_statUnc", e_ID_stat)
        self.out.fillBranch("Electron_ID_systUnc", e_ID_syst)

        return True

    
    def getLepSF(self, lep):
        '''Return lepton efficiency scale factor'''

        myLepID = abs(lep.pdgId)
        mySCeta = lep.eta
        isCrack = False # FIXME: isGap() is not available in nanoAODs, and cannot be recomputed easily based on eta, phi. We thus use the non-gap SFs for all electrons.
        isHoleBPix = False  # default

        if myLepID==11 :
            mySCeta = lep.eta + lep.deltaEtaSC # Use the SC eta and not the electron eta

        # Deal with very rare cases when SCeta is out of 2.5 bounds
        mySCeta = min(mySCeta,2.49)
        mySCeta = max(mySCeta,-2.49)

        SF, SFerror, SFReco, SFID = self.lepSFHelper.getSF(myLepID, lep.pt, lep.eta, mySCeta, lep.phi, isCrack)

        # Add a protection for leptons outside standard acceptance (pt<5/7 for mu/ele or |eta|>2.4 for mu) which get SF=0 and SFError = nan, since they may be still
        # be used for dedicated studies
        if SF==0 :
            SF, SFerror = 1., 0.5
        if SFReco == 0:
            SFReco = 1.
        if SFID == 0:
            SFID = 1.

        return SF, SFerror, SFReco, SFID
      

    def getLepSF_decorr(self, lep):
        """Return lepton SF uncertainties split into RECO/ID stat/syst (electrons only)."""

        myLepID = abs(lep.pdgId)
        mySCeta = lep.eta
        isCrack = False  # FIXME: cannot recompute from nanoAODs
        isHoleBPix = False

        # Only electrons have decorrelated uncertainties
        if myLepID != 11:
            return 0., 0., 0., 0.

        # Use SC eta for electrons
        mySCeta = lep.eta + lep.deltaEtaSC

        # Protect SCeta from out-of-bounds
        mySCeta = min(max(mySCeta, -2.49), 2.49)

        # Call the C++ function
        reco_stat, reco_syst, id_stat, id_syst = self.lepSFHelper.getSF_decorrUnc(
            myLepID,
            lep.pt,
            lep.eta,
            mySCeta,
            lep.phi,
            isCrack
        )

        # Get the corresponding SF
        SF, SFerror, SFReco, SFID= self.getLepSF(lep)

        # Convert absolute uncertainties to relative uncertainties
        if SF ==0:
            reco_stat, reco_syst, id_stat, id_syst = 0.5, 0.5, 0.5, 0.5

        return reco_stat, reco_syst, id_stat, id_syst

        # Protection for very rare edge cases (pt out-of-range etc.)
        #if reco_stat == 0 and reco_syst == 0 and id_stat == 0 and id_syst == 0:
            #reco_stat, reco_syst, id_stat, id_syst = 0., 0., 0., 0.

        #return reco_stat, reco_syst, id_stat, id_syst
