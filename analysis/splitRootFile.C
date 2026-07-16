void splitRootFile(const char* inputFileName, const char* treeName = "nTuple") {
    // Open the input ROOT file
    TFile* inputFile = TFile::Open(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open input file." << std::endl;
        return;
    }

    // Retrieve the tree from the input file
    TTree* inputTree = (TTree*)inputFile->Get(treeName);
    if (!inputTree) {
        std::cerr << "Error: Could not find tree " << treeName << " in file." << std::endl;
        return;
    }

    // Get the total number of entries in the tree
    Long64_t nEntries = inputTree->GetEntries();
    Long64_t entriesPerFile = nEntries / 8;

    // Loop to create 8 output files
    for (int i = 0; i < 8; ++i) {
        // Create a new output file
        TString outputFileName = TString::Format("output_%d.root", i);
        TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
        if (!outputFile || outputFile->IsZombie()) {
            std::cerr << "Error: Could not create output file " << outputFileName << std::endl;
            return;
        }

        // Clone the tree structure to the new file
        TTree* outputTree = inputTree->CloneTree(0); // Clone the tree structure, no entries yet

        // Copy the appropriate entries to the new tree
        Long64_t startEntry = i * entriesPerFile;
        Long64_t endEntry = (i == 7) ? nEntries : (i + 1) * entriesPerFile; // Last file may take the remaining entries

        for (Long64_t j = startEntry; j < endEntry; ++j) {
            inputTree->GetEntry(j); // Get the entry from the input tree
            outputTree->Fill();     // Fill the entry into the output tree
        }

        // Write the tree to the output file
        outputFile->cd();
        outputTree->Write();

        // Close the output file
        outputFile->Close();
        delete outputFile;
    }

    // Close the input file
    inputFile->Close();
    delete inputFile;

    std::cout << "Splitting completed successfully!" << std::endl;
}