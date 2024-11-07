#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <vector>
#include <iostream>

void FilterAndOrderJets() {
    // Hardcoded input file path
    std::string inputFile = "/eos/user/m/mmanoni/HZZ_samples_2022/MC_jetVeto/PROD_samplesNano_2022_MC_901ffb16/DYJetsToLL_forPOG/ZZ4lAnalysis.root";
    
    // Hardcoded output file path
    std::string outputFileName = "/eos/user/m/mmanoni/HZZ_samples_2022/MC_jetVeto/PROD_samplesNano_2022_MC_901ffb16/DYJetsToLL_forPOG/filtered_ZZ4lAnalysis.root";

    // Initialize RDataFrame for reading the input ROOT file
    ROOT::RDataFrame df("Events", inputFile);

    // Filter jets with JetId == 6 and debug the result
    auto filtered_df = df.Define("FilteredJetPt", "Jet_pt[Jet_jetId == 6]")
                         .Define("FilteredJetEta", "Jet_eta[Jet_jetId == 6]")
                         .Define("FilteredJetSmearUpPt", "Jet_smearUp_pt[Jet_jetId == 6]")
                         .Define("FilteredJetSmearDnPt", "Jet_smearDn_pt[Jet_jetId == 6]");

    // Debug: Check the number of filtered jets
    auto nFilteredJets = filtered_df.Filter("FilteredJetPt.size() > 0").Count();
    std::cout << "Number of events with at least one filtered jet: " << *nFilteredJets << std::endl;

    // Additional Debug: Print the first few entries of filtered jets
    filtered_df.Foreach([](const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta) {
        if (pt.size() > 0) {
            std::cout << "Filtered Jet Pt (first event): " << pt[0] << ", Eta: " << eta[0] << std::endl;
        }
    }, {"FilteredJetPt", "FilteredJetEta"});

    if (*nFilteredJets > 0) {
        // Sort the jets in descending order of pt and apply sorting indices to all quantities
        auto sorted_df = filtered_df.Define("SortedIndices", "Reverse(Argsort(FilteredJetPt))")
                                    .Define("OrderedJetPt", "FilteredJetPt[SortedIndices]")
                                    .Define("OrderedJetEta", "FilteredJetEta[SortedIndices]")
                                    .Define("OrderedJetSmearUpPt", "FilteredJetSmearUpPt[SortedIndices]")
                                    .Define("OrderedJetSmearDnPt", "FilteredJetSmearDnPt[SortedIndices]");

        // More Debugging: Print sorted jets before saving
        sorted_df.Foreach([](const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta) {
            if (pt.size() > 0) {
                std::cout << "Ordered Jet Pt (first event): " << pt[0] << ", Eta: " << eta[0] << std::endl;
            }
        }, {"OrderedJetPt", "OrderedJetEta"});

        // Define columns to save in the output file
        std::vector<std::string> columnsToSave = {"OrderedJetPt", "OrderedJetEta", "OrderedJetSmearUpPt", "OrderedJetSmearDnPt"};

        // Snapshot the filtered and ordered jets to a new ROOT file
        sorted_df.Snapshot("Events", outputFileName, columnsToSave);
        std::cout << "Filtered and ordered jets saved to " << outputFileName << std::endl;
    } else {
        std::cout << "No jets passed the filter criteria." << std::endl;
    }
}
