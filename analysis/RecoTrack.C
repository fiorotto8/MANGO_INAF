#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Constants defining the cylinder
double InnerSourceContThick = 5;
double GasRadius = 36.9; // Radius of the gas within the cylinder in mm
double GasDistanceFromCollim = 10;
double CollimatorDepth = 2;
double CollimatorDistance = 0;
double GasThickness = 50;
double containment_off=5;//mm
double cyl_center_y = 5;
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

// Function to determine if a point (x, y, z) is fully contained within a cylindrical volume.
// The cylinder is defined with specific dimensions and positional parameters.
bool is_fully_contained(double x, double y, double z) {
    // Thickness of an inner source container component i.e. the thickness of the source base
    double InnerSourceContThick = 5;
    // Radius of the gas within the cylinder.
    double GasRadius = 36.9;//in mm, from DetectorConstruciton instead of previous 40mm;
    // Distance from the collimator to the start of the gas region i.e. radialRingThickness?
    double GasDistanceFromCollim = 10;
    // Depth of the collimator.
    double CollimatorDepth = 2; //instead of previos 0.6
    // Distance from the collimator to another reference point, not used as it's zero.
    double CollimatorDistance = 0;
    // Thickness of the gas region in z
    double GasThickness = 50;

    // Calculating the z-coordinate of the cylinder's center. This includes half the thickness of the inner source,
    // half the depth of the collimator, and additional distances and radius to place the center correctly.
    //double cyl_center_z = InnerSourceContThick / 2 + CollimatorDepth / 2 + CollimatorDistance + GasRadius + GasDistanceFromCollim;
    double cyl_center_z = InnerSourceContThick / 2 + CollimatorDepth + CollimatorDistance + GasRadius + GasDistanceFromCollim;
    //std::cout<<cyl_center_z<<std::endl;

    // Calculate the distance from the point (x, z) to the center of the cylinder in the xz-plane.
    // This uses the Pythagorean theorem to find the distance from the point to the cylinder's central axis.
    double distance_from_center = std::sqrt(std::pow(x, 2) + std::pow(cyl_center_z - z, 2));
    // Calculate the absolute value of y to determine the vertical distance from the horizontal mid-plane of the cylinder.
    double vertical_distance = abs(y);
    // Check if the point is within the vertical limits (half the gas thickness minus a safety margin of 1 unit)
    // and radial limits (gas radius minus a safety margin of 1 unit) of the cylinder.
    if (vertical_distance < (GasThickness / 2 - 1) && distance_from_center < (GasRadius - 1)) {
        return true;  // The point is fully contained within the cylinder.
    } else {
        return false; // The point is outside the cylinder.
    }
}

bool areAllPointsInsideCylinder(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z) {
    // Iterate over all points
    for (size_t i = 0; i < x.size(); ++i) {
        if (std::abs(y[i] - cyl_center_y) > (GasThickness / 2 - containment_off)) {
            return false;
        }

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

// Function to create a 3D plot from vectors of x, y, and z data.
void Draw3DPlot(const std::vector<double>& xData, const std::vector<double>& yData, const std::vector<double>& zData, int i) {
    // Check that the vectors have the same size
    if (xData.size() != yData.size() || yData.size() != zData.size()) {
        std::cerr << "Error: xData, yData, and zData must be of the same size." << std::endl;
        return;
    }

    // Create a new TGraph2D object
    TGraph2D *graph = new TGraph2D();

    // Fill the graph with data
    for (size_t i = 0; i < xData.size(); ++i) {
        graph->SetPoint(i, xData[i], yData[i], zData[i]);
    }

    // Set graph title and axis labels
    graph->SetTitle("3D Scatter Plot;X axis;Y axis;Z axis");

    // Create a new canvas to draw the graph
    TCanvas *c1 = new TCanvas("c1", "3D Scatter Plot", 1000, 1000);
    c1->cd();

    // Draw the graph with points
    graph->Draw("P0");  // "P0" option to draw the graph with points only.

    // Update the canvas to display the graph
    c1->Update();

    // Optionally save the canvas to a file
    c1->SaveAs(Form("3DScatterPlot_%d.png", i));
    c1->Write();
}

void DrawGraphAndLineZY(const std::vector<double>& zData, const std::vector<double>& yData, int i) {
    if (zData.size() != yData.size()) {
        std::cerr << "Error: zData and yData must be of the same size." << std::endl;
        return;
    }

    // Create a TGraph from the z and y data
    TGraph *graph = new TGraph(zData.size());
    for (size_t i = 0; i < zData.size(); ++i) {
        graph->SetPoint(i, zData[i], yData[i]);
    }

    // Setup graph appearance
    graph->SetTitle("Graph of Z vs Y;Z axis;Y axis");
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue);

    // Create a new canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas for Graph and Line", 800, 600);
    c1->cd();

    // Set fixed ranges for the axes before drawing
    graph->Draw("AP");  // Draw graph to access its histogram
    graph->GetHistogram()->GetXaxis()->SetLimits(0, 100); // Set X-axis range
    graph->GetHistogram()->SetMinimum(-100);   // Set minimum Y-axis
    graph->GetHistogram()->SetMaximum(100); // Set maximum Y-axis

    // Redraw the graph with updated axes
    graph->Draw("AP");

    double xmin,xmax, yup,ydown;
    xmin=zedge_min;
    xmax=zedge_max;
    yup=GasThickness;
    ydown=-1*GasThickness;
    // Draw a vertical line at x = 51.9
    TLine *line = new TLine(xmax, ydown, xmax, yup);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw();
    // Draw a vertical line at x = 14.5
    TLine *line3 = new TLine(xmin, ydown, xmin, yup);
    line3->SetLineColor(kRed);
    line3->SetLineWidth(2);
    line3->Draw();
    // Draw a horizontal line at y = 50
    TLine *line1 = new TLine(xmin, yup, xmax, yup);
    line1->SetLineColor(kRed);
    line1->SetLineWidth(2);
    line1->Draw();
    // Draw a horizontal line at y = -50
    TLine *line2 = new TLine(xmin, ydown, xmax, ydown);
    line2->SetLineColor(kRed);  // You can change the color if needed
    line2->SetLineWidth(2);
    line2->Draw();

    // Update the canvas to show everything
    c1->Update();
    // Optionally save the canvas to a file
    c1->SaveAs(Form("ZYplot_%d.png", i));
    c1->Write();
}

void DrawGraphAndLineZX(const std::vector<double>& zData, const std::vector<double>& xData, int i) {
    if (zData.size() != xData.size()) {
        std::cerr << "Error: zData and xData must be of the same size." << std::endl;
        return;
    }

    // Create a TGraph from the z and y data
    TGraph *graph = new TGraph(zData.size());
    for (size_t i = 0; i < zData.size(); ++i) {
        graph->SetPoint(i, zData[i], xData[i]);
    }

    // Setup graph appearance
    graph->SetTitle("Graph of Z vs X;Z axis;X axis");
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue);

    // Create a new canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas for Graph and Line", 800, 600);
    c1->cd();

    // Set fixed ranges for the axes before drawing
    graph->Draw("AP");  // Draw graph to access its histogram
    graph->GetHistogram()->GetXaxis()->SetLimits(0, 100); // Set X-axis range
    graph->GetHistogram()->SetMinimum(-100);   // Set minimum Y-axis
    graph->GetHistogram()->SetMaximum(100); // Set maximum Y-axis

    // Redraw the graph with updated axes
    graph->Draw("AP");


    double xmin,xmax, yup,ydown;
    xmin=zedge_min;
    xmax=zedge_max;
    yup=GasThickness;
    ydown=-1*GasThickness;
    // Draw a vertical line at x = 51.9
    TLine *line = new TLine(xmax, ydown, xmax, yup);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw();
    // Draw a vertical line at x = 14.5
    TLine *line3 = new TLine(xmin, ydown, xmin, yup);
    line3->SetLineColor(kRed);
    line3->SetLineWidth(2);
    line3->Draw();
    // Draw a horizontal line at y = 50
    TLine *line1 = new TLine(xmin, yup, xmax, yup);
    line1->SetLineColor(kRed);
    line1->SetLineWidth(2);
    line1->Draw();
    // Draw a horizontal line at y = -50
    TLine *line2 = new TLine(xmin, ydown, xmax, ydown);
    line2->SetLineColor(kRed);  // You can change the color if needed
    line2->SetLineWidth(2);
    line2->Draw();

    // Update the canvas to show everything
    c1->Update();
    // Optionally save the canvas to a file
    c1->SaveAs(Form("ZXplot_%d.png", i));
    c1->Write();
}

double ComputeTrackLength(const std::vector<Double_t>& xs,
                          const std::vector<Double_t>& ys,
                          const std::vector<Double_t>& zs) {
    double L = 0.0;
    if (xs.size() < 2) return 0.0;

    for (size_t i = 1; i < xs.size(); ++i) {
        double dx = xs[i] - xs[i-1];
        double dy = ys[i] - ys[i-1];
        double dz = zs[i] - zs[i-1];
        L += std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    return L;
}


// Function to map ParticleTag to ParticleName_out
    std::string GetParticleName(int ParticleTag) {
        if (ParticleTag == 0) {
            return "e-";
        } else if (ParticleTag == 1) {
            return "e+";
        } else if (ParticleTag == 2) {
            return "gamma";
        } else if (ParticleTag == 3) {
            return "alpha";
        } else {
            return "unknown";
        }
    }

// Reduce Geant4 hit steps into one entry per physical track.
// EnergyDeposit is stored by Geant4 in MeV.
void RecoTrack(std::string filename,
               std::string particle = "e-",
               bool primaryOnly = true) {
    const double W_factor = 38E-6; // MeV per ion pair

    TFile* inputFile = TFile::Open(filename.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Cannot open input ROOT file: " << filename << std::endl;
        return;
    }

    TTree* tree = nullptr;
    inputFile->GetObject("nTuple", tree);
    if (!tree) {
        std::cerr << "Tree 'nTuple' is missing in " << filename << std::endl;
        inputFile->Close();
        return;
    }

    const char* requiredBranches[] = {
        "EventNumber", "ParticleName", "x_hits", "y_hits", "z_hits",
        "EnergyDeposit", "Nucleus", "ProcessType",
        "tracklen_hits", "currentTrackID"
    };
    for (const char* branchName : requiredBranches) {
        if (!tree->GetBranch(branchName)) {
            std::cerr << "Required branch is missing: " << branchName << std::endl;
            inputFile->Close();
            return;
        }
    }

    Int_t eventNumber = 0;
    Int_t trackId = 0;
    Double_t xHit = 0.;
    Double_t yHit = 0.;
    Double_t zHit = 0.;
    Double_t energyDeposit = 0.;
    Double_t stepLength = 0.;
    Char_t particleName[64] = {};
    Char_t nucleus[64] = {};
    Char_t processType[64] = {};

    tree->SetBranchStatus("*", 0);
    for (const char* branchName : requiredBranches) {
        tree->SetBranchStatus(branchName, 1);
    }

    tree->SetBranchAddress("EventNumber", &eventNumber);
    tree->SetBranchAddress("ParticleName", particleName);
    tree->SetBranchAddress("x_hits", &xHit);
    tree->SetBranchAddress("y_hits", &yHit);
    tree->SetBranchAddress("z_hits", &zHit);
    tree->SetBranchAddress("EnergyDeposit", &energyDeposit);
    tree->SetBranchAddress("Nucleus", nucleus);
    tree->SetBranchAddress("ProcessType", processType);
    tree->SetBranchAddress("tracklen_hits", &stepLength);
    tree->SetBranchAddress("currentTrackID", &trackId);

    TString inputBase = gSystem->BaseName(filename.c_str());
    TString outputName = "elab_" + inputBase;
    TFile outputFile(outputName, "RECREATE");
    TTree outputTree("elabHits", "elabHits");

    Int_t outputEvent = -1;
    Int_t outputTrack = -1;
    Int_t numberOfHits = 0;
    Double_t totalEnergyDeposit = 0.;
    Double_t trackLength = 0.;
    Double_t primaries = 0.;
    Bool_t fullyContained = true;
    std::string outputParticle;
    std::string outputNucleus;
    std::vector<Double_t> xHits;
    std::vector<Double_t> yHits;
    std::vector<Double_t> zHits;
    std::vector<Double_t> energyDeposits;
    std::vector<Double_t> clusterPrimaries;

    outputTree.Branch("EventNumber", &outputEvent);
    outputTree.Branch("TrackID", &outputTrack);
    outputTree.Branch("TotalEDep", &totalEnergyDeposit);
    outputTree.Branch("nhits", &numberOfHits);
    outputTree.Branch("x_hits", &xHits);
    outputTree.Branch("y_hits", &yHits);
    outputTree.Branch("z_hits", &zHits);
    outputTree.Branch("Edep_hits", &energyDeposits);
    outputTree.Branch("FullyContained", &fullyContained);
    outputTree.Branch("NameParticle", &outputParticle);
    outputTree.Branch("Nucleus", &outputNucleus);
    outputTree.Branch("TrackLength", &trackLength);
    outputTree.Branch("ClusterPrimaries", &clusterPrimaries);
    outputTree.Branch("Primaries", &primaries);

    bool hasTrack = false;

    auto flushTrack = [&]() {
        if (!hasTrack) {
            return;
        }
        primaries = totalEnergyDeposit / W_factor;
        fullyContained = areAllPointsInsideCylinder(xHits, yHits, zHits);
        outputTree.Fill();
    };

    const Long64_t numberOfEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < numberOfEntries; ++entry) {
        if (entry % 10000 == 0) {
            std::cout << entry << "/" << numberOfEntries << std::endl;
        }

        tree->GetEntry(entry);

        if (particleName != particle) {
            continue;
        }
        if (primaryOnly && std::string(processType) != "RadioactiveDecay") {
            continue;
        }

        const bool sameTrack =
            hasTrack && outputEvent == eventNumber && outputTrack == trackId;

        if (!sameTrack) {
            flushTrack();

            outputEvent = eventNumber;
            outputTrack = trackId;
            outputParticle = particleName;
            outputNucleus = nucleus;
            numberOfHits = 0;
            totalEnergyDeposit = 0.;
            trackLength = 0.;
            primaries = 0.;
            fullyContained = true;
            xHits.clear();
            yHits.clear();
            zHits.clear();
            energyDeposits.clear();
            clusterPrimaries.clear();
            hasTrack = true;
        }

        totalEnergyDeposit += energyDeposit;
        trackLength += stepLength;
        xHits.push_back(xHit);
        yHits.push_back(yHit);
        zHits.push_back(zHit);
        energyDeposits.push_back(energyDeposit);
        clusterPrimaries.push_back(energyDeposit / W_factor);
        ++numberOfHits;
    }

    flushTrack();
    outputFile.cd();
    outputTree.Write();
    outputFile.Close();
    inputFile->Close();
}
