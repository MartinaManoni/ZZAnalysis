void check_leading_jet() {
    // Open the ROOT file
    TFile *file = TFile::Open("/eos/user/m/mmanoni/HZZ_samples_2022/MC_jetVeto/PROD_samplesNano_2022_MC_901ffb16/DYJetsToLL_forPOG/ZZ4lAnalysis.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Access the Events tree
    TTree *tree = (TTree*)file->Get("Events");
    if (!tree) {
        std::cerr << "Error accessing tree 'Events'!" << std::endl;
        file->Close();
        return;
    }

    // Set up branches for Jet_pt and Jet_eta as arrays
    const int maxJets = 100; // Assume a reasonable maximum number of jets per event
    float Jet_pt[maxJets];
    float Jet_eta[maxJets];
    int nJet;

    tree->SetBranchAddress("Jet_pt", Jet_pt);
    tree->SetBranchAddress("Jet_eta", Jet_eta);
    tree->SetBranchAddress("nJet", &nJet);  // Assuming 'nJet' exists and represents the number of jets in each event

    Long64_t nEntries = tree->GetEntries();
    bool allEventsMatch = true;
    int mismatchCount = 0; // Counter for mismatches

    // Loop over each event in the tree
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // Skip events with no jets
        if (nJet == 0) continue;

        // Find the index of the jet with the highest pt
        size_t leadingJetIndex = 0;
        float maxPt = Jet_pt[0];
        for (int j = 1; j < nJet; j++) {
            if (Jet_pt[j] > maxPt) {
                maxPt = Jet_pt[j];
                leadingJetIndex = j;
            }
        }

        // Check if Jet_eta[0] corresponds to the jet with highest Jet_pt
        if (leadingJetIndex != 0) {
            std::cout << "Mismatch in event " << i << ":\n";
            std::cout << "  - Jet_pt[0] = " << Jet_pt[0] << " (corresponds to Jet_eta[0] = " << Jet_eta[0] << ")\n";
            std::cout << "  - Actual leading Jet_pt = " << maxPt 
                      << " (found at index " << leadingJetIndex 
                      << " with Jet_eta[" << leadingJetIndex << "] = " << Jet_eta[leadingJetIndex] << ")\n";
            std::cout << "  - First 4 jets in the event:\n";

            // Print first 4 jets
            for (int j = 0; j < std::min(4, nJet); j++) {
                std::cout << "    Jet " << j << ": Jet_pt = " << Jet_pt[j] << ", Jet_eta = " << Jet_eta[j] << "\n";
            }
            allEventsMatch = false;
            mismatchCount++;
        }
    }

    if (allEventsMatch) {
        std::cout << "In all events, Jet_eta[0] corresponds to the jet with the highest Jet_pt.\n";
    } else {
        std::cout << "Total mismatches found: " << mismatchCount << "\n";
    }

    // Close the file
    file->Close();
}
