#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <output_file_name>" << std::endl;
        return 1;
    }

    // Get the output file name from the command line arguments
    std::string filename = argv[1];

    // Open the ROOT file
    TFile *file = TFile::Open("output_t0.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file output_t0.root" << std::endl;
        return 1;
    }

    // Get the ntuple from the file
    TTree *tree = nullptr;
    file->GetObject("nTuple", tree);
    if (!tree) {
        std::cerr << "Error: ntuple 'nTuple' not found in file" << std::endl;
        file->Close();
        return 1;
    }

    // Define variables to hold the data
    Int_t Evn, ParticleID, ParticleTag, ParentID, VolumeNumber, pdgID_hits, currentTrackID;
    Double_t x_hits, y_hits, z_hits, EnergyDeposit, VolumeTraslX, VolumeTraslY, VolumeTraslZ;
    Double_t px_particle, py_particle, pz_particle, tracklen_hits;
    Char_t Nucleus[20];
    Char_t ParticleName[20];
    Char_t CreationProcess[20];

    // Set branch addresses
    tree->SetBranchAddress("EventNumber", &Evn);
    tree->SetBranchAddress("ParticleName", &ParticleName);
    tree->SetBranchAddress("ParticleID", &ParticleID);
    tree->SetBranchAddress("ParticleTag", &ParticleTag);
    tree->SetBranchAddress("ParentID", &ParentID);
    tree->SetBranchAddress("x_hits", &x_hits);
    tree->SetBranchAddress("y_hits", &y_hits);
    tree->SetBranchAddress("z_hits", &z_hits);
    tree->SetBranchAddress("EnergyDeposit", &EnergyDeposit);
    tree->SetBranchAddress("VolumeNumber", &VolumeNumber);
    tree->SetBranchAddress("VolumeTraslX", &VolumeTraslX);
    tree->SetBranchAddress("VolumeTraslY", &VolumeTraslY);
    tree->SetBranchAddress("VolumeTraslZ", &VolumeTraslZ);
    tree->SetBranchAddress("Nucleus", &Nucleus);
    tree->SetBranchAddress("ProcessType", &CreationProcess);
    tree->SetBranchAddress("px_particle", &px_particle);
    tree->SetBranchAddress("py_particle", &py_particle);
    tree->SetBranchAddress("pz_particle", &pz_particle);
    tree->SetBranchAddress("pdgID_hits", &pdgID_hits);
    tree->SetBranchAddress("tracklen_hits", &tracklen_hits);
    tree->SetBranchAddress("currentTrackID", &currentTrackID);

    // Loop over the entries in the ntuple and print them
    Long64_t nEntries = tree->GetEntries();

    // Initialize variables for output tree
    Int_t Out_event;
    Int_t nhits_out;
    Double_t ETotal;
    Double_t ETotal_NR;
    vector<Int_t> pdgID;
    vector<Double_t> tracklen;
    //Double_t px_part;
    //Double_t py_part;
    //Double_t pz_part;
    vector<Double_t> px_part;
    vector<Double_t> py_part;
    vector<Double_t> pz_part;
    vector<Double_t> EdepHits_out;
    vector<Double_t> x_hits_out;
    vector<Double_t> y_hits_out;
    vector<Double_t> z_hits_out;
    //others
    std::string nucl;


    // Create a new ROOT file to store output data
    TFile* f_out = new TFile(Form("%s.root", filename.c_str()), "RECREATE");
    TTree* outTree = new TTree("nTuple", "nTuple");

    // Define branches for the output tree
    outTree->Branch("eventnumber", &Out_event);
    outTree->Branch("numhits", &nhits_out);
    outTree->Branch("energyDep", &ETotal);
    outTree->Branch("energyDep_NR", &ETotal_NR);
    outTree->Branch("pdgID_hits", &pdgID);
    outTree->Branch("tracklen_hits", &tracklen);
    outTree->Branch("px_particle", &px_part);
    outTree->Branch("py_particle", &py_part);
    outTree->Branch("pz_particle", &pz_part);
    outTree->Branch("energyDep_hits", &EdepHits_out);
    outTree->Branch("x_hits", &x_hits_out);
    outTree->Branch("y_hits", &y_hits_out);
    outTree->Branch("z_hits", &z_hits_out);

    // Read the first entry to initialize variables.
    tree->GetEntry(0);
    Out_event=Evn;
    nucl=Nucleus;
    ETotal=0;
    nhits_out=0;

    // Loop over the entries in the ntuple and fill the output tree
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        // Print progress for every 10000 entries processed.
        if(i%10000==0)std::cout << i << "/" << tree->GetEntries() << std::endl;

        // Check if the current hit belongs to the same event and nucleus as previously processed.
        if(Out_event == Evn && nucl == Nucleus){
        // Accumulate total energy deposited.
        ETotal += EnergyDeposit*1000;//digi wants it in keV
        nhits_out++;
        // Store hit data.
        px_part.push_back(px_particle);
        py_part.push_back(py_particle);
        pz_part.push_back(pz_particle);
        x_hits_out.push_back(x_hits);
        y_hits_out.push_back(y_hits);
        z_hits_out.push_back(z_hits);
        EdepHits_out.push_back(EnergyDeposit*1000);
        pdgID.push_back(pdgID_hits);
        tracklen.push_back(tracklen_hits);
        }
        else{
        // If a new event or nucleus is encountered, save the data from the previous event.
        Out_event=Evn;
        //nhits_out
        //ETotal
        ETotal_NR=0;
        //px_part=px_particle;
        //py_part=py_particle;
        //pz_part=pz_particle;
        //FILL
        outTree->Fill();
        //RESET
        Out_event=Evn;
        nucl=Nucleus;
        ETotal=0;
        nhits_out=1;
        x_hits_out.clear(); y_hits_out.clear(); z_hits_out.clear(); EdepHits_out.clear();
        px_part.clear(); py_part.clear(); pz_part.clear();
        // Start accumulating new event data.
        ETotal+=EnergyDeposit*1000;//digi wants it in keV
        x_hits_out.push_back(x_hits);
        y_hits_out.push_back(y_hits);
        z_hits_out.push_back(z_hits);
        px_part.push_back(px_particle);
        py_part.push_back(py_particle);
        pz_part.push_back(pz_particle);
        EdepHits_out.push_back(EnergyDeposit*1000);
        }





    }

    // Write the output tree to the new file and close it
    outTree->Write();
    f_out->Close();

    // Close the input file
    file->Close();

    return 0;
}