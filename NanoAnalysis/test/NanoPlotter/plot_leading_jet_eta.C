#include <iostream>
#include <algorithm>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

void plot_leading_jet_eta() {
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

    // Set up branches for Jet_pt and Jet_eta as fixed-size arrays
    const int maxJets = 100; // Define a maximum number of jets
    Float_t Jet_pt[maxJets];
    Float_t Jet_eta[maxJets];
    int nJet;

    tree->SetBranchAddress("Jet_pt", Jet_pt);
    tree->SetBranchAddress("Jet_eta", Jet_eta);
    tree->SetBranchAddress("nJet", &nJet);  // Assuming 'nJet' exists and represents the number of jets in each event

    Long64_t nEntries = tree->GetEntries();

    // Create a histogram to store Jet_eta[0] of leading jets
    TH1F *hJetEta = new TH1F("hJetEta", "Leading Jet Eta; Jet #eta; Entries", 100, -5, 5);

    // Loop over each event in the tree
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // Skip events with no jets
        if (nJet == 0) continue;

        // Find the index of the leading jet (highest Jet_pt)
        int leadingJetIndex = std::distance(Jet_pt, std::max_element(Jet_pt, Jet_pt + nJet));

        // Fill histogram with Jet_eta of the leading jet
        if (leadingJetIndex < nJet) { // Ensure the index is valid
            hJetEta->Fill(Jet_eta[leadingJetIndex]);
        }
    }

    // Draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Leading Jet Eta", 800, 600);
    hJetEta->Draw();

    // Save the plot as a .png file (optional)
    canvas->SaveAs("leading_jet_eta.png");

    // Clean up
    delete hJetEta;
    file->Close();
}
