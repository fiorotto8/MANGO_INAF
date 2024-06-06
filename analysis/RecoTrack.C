// Constants defining the cylinder
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

// Function to process track data from a ROOT file and write selected information to a new ROOT file.
void RecoTrack(std::string filename){
    // Open the input ROOT file.
    TFile* f = TFile::Open(filename.c_str());
    // Retrieve a TTree named "Hits" from the file.
    TTree* tree = (TTree*)f->Get("Hits");

    // Define variables to hold data from the tree.
    Int_t Evn,ParticleID,ParticleTag,ParentID,VolumeNumber;
    Double_t x_hits,y_hits,z_hits,VolumeTraslX,VolumeTraslY,VolumeTraslZ,EnergyDeposit;
    Char_t Nucleus[20];
    Char_t ParticleName[20];
    Char_t CreationProcess[20];

    // Set the addresses of local variables where tree data will be stored during the reading process.
    tree->SetBranchAddress("EventNumber",&Evn);
    tree->SetBranchAddress("ParticleID",&ParticleID);
    tree->SetBranchAddress("ParticleTag",&ParticleTag);
    tree->SetBranchAddress("ParticleName",&ParticleName);
    tree->SetBranchAddress("ParentID",&ParentID);
    tree->SetBranchAddress("VolumeNumber",&VolumeNumber);
    tree->SetBranchAddress("x_hits",&x_hits);
    tree->SetBranchAddress("y_hits",&y_hits);
    tree->SetBranchAddress("z_hits",&z_hits);
    tree->SetBranchAddress("EnergyDeposit",&EnergyDeposit);
    tree->SetBranchAddress("VolumeNumber",&VolumeNumber);
    tree->SetBranchAddress("VolumeTraslX",&VolumeTraslX);
    tree->SetBranchAddress("VolumeTraslY",&VolumeTraslY);
    tree->SetBranchAddress("VolumeTraslZ",&VolumeTraslZ);
    tree->SetBranchAddress("Nucleus",&Nucleus);
    tree->SetBranchAddress("ProcessType",&CreationProcess);

    // Initialize variables for output data.
    Int_t Out_event;
    Double_t ETotal;
    Int_t nhits_out;
    std::vector<Double_t> x_hits_out, y_hits_out, z_hits_out, EdepHits_out;
    std::string nucl;
    bool fullyCont;

    // Create a new ROOT file to store output data.
    TFile* f_out = new TFile(Form("elab_%s",filename.c_str()),"recreate");
    TTree* outTree = new TTree("elabHits","elabHits");

    // Define branches for the output tree.
    outTree->Branch("EventNumber",&Out_event);
    outTree->Branch("TotalEDep",&ETotal);
    outTree->Branch("nhits",&nhits_out);
    outTree->Branch("x_hits",&x_hits_out);
    outTree->Branch("y_hits",&y_hits_out);
    outTree->Branch("z_hits",&z_hits_out);
    outTree->Branch("Edep_hits",&EdepHits_out);
    outTree->Branch("FullyContained",&fullyCont);

    // Read the first entry to initialize variables.
    tree->GetEntry(0);
    Out_event=Evn;
    nucl=Nucleus;
    ETotal=0;
    nhits_out=0;
    fullyCont=true;

    // Process each entry in the tree.
    for(int i=0; i<tree->GetEntries(); i++){
    //for(int i=0; i<200; i++){

    // Print progress for every 10000 entries processed.
    if(i%10000==0)std::cout << i << "/" << tree->GetEntries() << std::endl;

    // Read the current entry.
    tree->GetEntry(i);

    // Check if the current hit belongs to the same event and nucleus as previously processed.
    if(Out_event == Evn && nucl == Nucleus){
        // Accumulate total energy deposited.
        ETotal += EnergyDeposit;
        // Store hit data.
        x_hits_out.push_back(x_hits);
        y_hits_out.push_back(y_hits);
        z_hits_out.push_back(z_hits);
        EdepHits_out.push_back(EnergyDeposit);
        nhits_out++;
        // Check containment if the number of hits exceeds 50.
        //if(fullyCont) fullyCont = is_fully_contained(x_hits,y_hits,z_hits);
        /* if(nhits_out>10){
        if(fullyCont) fullyCont = is_fully_contained(x_hits,y_hits,z_hits);
        } */
    }
    else{
        // If a new event or nucleus is encountered, save the data from the previous event.
        fullyCont=areAllPointsInsideCylinder(x_hits_out, y_hits_out, z_hits_out);
        //if (fullyCont) outTree->Fill();
        outTree->Fill();

        // Reset variables for the new event.
        Out_event=Evn;
        nucl=Nucleus;
        ETotal=0;
        nhits_out=1;
        // Draw the plot
        //Draw3DPlot(x_hits_out, y_hits_out, z_hits_out,i);
        //DrawGraphAndLineZY(z_hits_out,y_hits_out,i);
        //DrawGraphAndLineZX(z_hits_out,x_hits_out,i);
        x_hits_out.clear(); y_hits_out.clear(); z_hits_out.clear(); EdepHits_out.clear();

        // Start accumulating new event data.
        ETotal+=EnergyDeposit;
        x_hits_out.push_back(x_hits);
        y_hits_out.push_back(y_hits);
        z_hits_out.push_back(z_hits);
        EdepHits_out.push_back(EnergyDeposit);
        fullyCont=true;
    }//chioudo if

    }//chiudo for entrate

    // Write the remaining data to the file.
    f_out->cd();
    outTree->Write();
    f_out->Save();
    f_out->Close();

}
