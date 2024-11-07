#include <iostream>
#include <vector>
#include <algorithm>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

void plot_leading_jet_eta_filtered() {
    // Open the ROOT file
    TFile *file = TFile::Open("/eos/user/m/mmanoni/HZZ_samples_2022/Data_jetVeto/PROD_samplesNano_2022_Data_901ffb16/Data_eraCD_preEE.root");
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

    // Set up branches for Jet properties
    const int maxJets = 100; // Define a maximum number of jets
    Float_t Jet_pt[maxJets];
    Float_t Jet_eta[maxJets];
    UChar_t Jet_jetId[maxJets]; // Change to UChar_t for Jet IDs
    int nJet;
    
    // Change the type of ZLCand_lepIdx to Short_t
    Short_t ZLCand_lepIdx;

    tree->SetBranchAddress("Jet_pt", Jet_pt);
    tree->SetBranchAddress("Jet_eta", Jet_eta);
    tree->SetBranchAddress("Jet_jetId", Jet_jetId); // Use UChar_t
    tree->SetBranchAddress("nJet", &nJet);  // Assuming 'nJet' exists and represents the number of jets in each event
    tree->SetBranchAddress("ZLCand_lepIdx", &ZLCand_lepIdx); // Use Short_t for ZLCand_lepIdx

    Long64_t nEntries = tree->GetEntries();

    // Create a histogram to store Jet_eta[0] of leading jets after filtering
    TH1F *hJetEta = new TH1F("hJetEta", "Leading Jet Eta (ID=6); Jet #eta; Entries", 100, -5, 5);

    // Loop over each event in the tree
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // Maintain only events with ZLCand_lepIdx >= 0
        if (ZLCand_lepIdx >= 0) {  // Keep the event if condition is met
            // Create vectors to hold the filtered jets
            std::vector<float> filteredJetPt;
            std::vector<float> filteredJetEta;

            // Filter jets based on Jet_jetId == 6
            for (int j = 0; j < nJet; j++) {
                if (Jet_jetId[j] == 6) {
                    filteredJetPt.push_back(Jet_pt[j]);
                    filteredJetEta.push_back(Jet_eta[j]);
                }
            }

            // Check if we have filtered jets
            if (filteredJetPt.empty()) continue; // Skip if no jets pass the filter

            // Find the leading jet (highest pT) among the filtered jets
            auto maxPtIter = std::max_element(filteredJetPt.begin(), filteredJetPt.end());
            int leadingJetIndex = std::distance(filteredJetPt.begin(), maxPtIter);

            // Fill histogram with Jet_eta of the leading filtered jet
            hJetEta->Fill(filteredJetEta[leadingJetIndex]);
        }
    }

    // Draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Leading Jet Eta (ID=6)", 800, 600);
    hJetEta->Draw();

    // Save the plot as a .png file (optional)
    canvas->SaveAs("leading_jet_eta_filtered_DATI.png");

    // Clean up
    delete hJetEta;
    file->Close();
}
