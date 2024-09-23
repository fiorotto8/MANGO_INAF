#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <vector>

using namespace std;

double InnerSourceContThick = 5;
double GasRadius = 36.9; // Radius of the gas within the cylinder in mm
double GasDistanceFromCollim = 10;
double CollimatorDepth = 2;
double CollimatorDistance = 0;
double GasThickness = 50;
double containment_off=5;//mm
/*Positive z-direction is at 0 radians.
Positive x-direction would be at π/2 radians.
Negative z-direction  would be at π radians. //so this one
Negative x-direction would be at 3π/2 radians.*/
double acceptedAngleStart=M_PI-(30*M_PI/180);
double acceptedAngleEnd=M_PI+(30*M_PI/180);
// Calculate the z-coordinate of the cylinder's center
double cyl_center_z = InnerSourceContThick / 2 + CollimatorDepth + CollimatorDistance + GasRadius + GasDistanceFromCollim;
double zedge_min= InnerSourceContThick / 2 + CollimatorDepth + CollimatorDistance+ GasDistanceFromCollim;
double zedge_max=InnerSourceContThick / 2 + CollimatorDepth + CollimatorDistance + 2*GasRadius + GasDistanceFromCollim;



bool areAllPointsInsideCylinder(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z) {
    // Iterate over all points
    for (size_t i = 0; i < x.size(); ++i) {
        // Calculate the distance from the center of the cylinder (assumed to be at (0, cylCenterZ))
        double dx = x[i]; // x[i] - 0 is simply x[i]
        double dz = z[i] - cyl_center_z;
        double r = std::sqrt(dx * dx + dz * dz); // radial distance in the x-z plane

         // Calculate the angle in radians from the positive z-axis
        double angle = atan2(dx, dz); // atan2 returns the angle in radians between [-π, π]
        // Normalize the angle to be within [0, 2π]
        if (angle < 0) angle += 2 * M_PI;

        // Check if the point is within the accepted angle range
        bool isWithinAcceptedAngle = (angle >= acceptedAngleStart && angle <= acceptedAngleEnd);

        // Check if the point is inside the circle or within the accepted angle
        if (!(r <= (GasRadius - containment_off) || isWithinAcceptedAngle)) {
            // Optionally, for debugging, print which point is outside the acceptable conditions
            // std::cout << "Point at index " << i << " is outside the acceptable area" << std::endl;
            return false;  // If any point is outside the acceptable conditions, return false immediately
        }
    }
    return true;  // If all points are inside the acceptable conditions, return true
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <output_file_name> <fill_option: 1=check, 0=fill_all>" << std::endl;
        return 1;
    }

    // Get the output file name from the command line arguments
    std::string filename = argv[1];

    // Get the fill option from the command line arguments (1=check, 0=fill everything)
    bool check_points = std::stoi(argv[2]) == 1;  // If the argument is 1, check points

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
        // Check points if the flag is set, otherwise fill unconditionally
        if (!check_points || areAllPointsInsideCylinder(x_hits_out, y_hits_out, z_hits_out)) {
            outTree->Fill();
        }
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