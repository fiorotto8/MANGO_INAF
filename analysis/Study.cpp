#include <TTree.h>
#include <TFile.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMath.h>
#include <TH1D.h>
#include <memory>
#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>
#include <cmath>

int main(int argc, char** argv){
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <elaborated.root>" << std::endl;
    return EXIT_FAILURE;
  }

  // Open a ROOT file specified by the command line argument
  auto f = std::unique_ptr<TFile>( TFile::Open(Form("%s",argv[1]),"r") );
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open input file: " << argv[1] << std::endl;
    return EXIT_FAILURE;
  }

  // Retrieve a TTree object named "elabHits" from the file
  TTree* mytree = nullptr;
  f->GetObject("elabHits", mytree);
  if (!mytree) {
    std::cerr << "Tree 'elabHits' is missing in " << argv[1] << std::endl;
    return EXIT_FAILURE;
  }

  // Setup TTreeReader to read the TTree efficiently
  TTreeReader reader(mytree);
  TTreeReaderValue<Int_t> EventNumber(reader, "EventNumber");
  TTreeReaderValue<Int_t> nhits(reader, "nhits");
  TTreeReaderValue<Double_t> TotalEDep(reader, "TotalEDep");
  TTreeReaderValue<bool> FullyContained(reader, "FullyContained");
  TTreeReaderValue<std::vector<double>> x_hits(reader, "x_hits");
  TTreeReaderValue<std::vector<double>> y_hits(reader, "y_hits");
  TTreeReaderValue<std::vector<double>> z_hits(reader, "z_hits");
  TTreeReaderValue<std::vector<double>> edep_hits(reader, "Edep_hits");

  // Initialize a graph and a linear fit function
  std::unique_ptr<TGraph> graph(new TGraph());
  std::unique_ptr<TF1> fit( new TF1("fit","pol1") );

  std::vector<Double_t> X,Z;

  //1.2cm is the fitted distance
  // Histogram to store angle calculations
  std::unique_ptr<TH1D> histoThr( new TH1D("histoIntrRes(rad)","histoIntrRes(rad)",150,-TMath::Pi(),TMath::Pi()) );
  std::unique_ptr<TH1D> histoThd( new TH1D("histoIntrRes(deg)","histoIntrRes(deg)",150,-180,180 ));

  std::vector<TGraph*> TGvec;
  double dist=0;

  Double_t angDeg,angRad,EDep;
  bool FContained;
  // Prepare file for output
  auto fout = std::unique_ptr<TFile>( new TFile("outHisto.root","recreate") );
  TTree* outTree = new TTree("angles","angles");
  fout.get()->mkdir("exampleTracks");
  fout.get()->cd("exampleTracks");

  // Define branches for the output tree.
  outTree->Branch("AngleDegree",&angDeg);
  outTree->Branch("AngleRadians",&angRad);
  outTree->Branch("EDep",&EDep);
  outTree->Branch("FContained",&FContained);

  // Iterate over all entries in the TTree
  for (auto entry : reader){
    EDep = *TotalEDep;  // Assign the TotalEDep from the input TTree to EDep
    FContained=*FullyContained;
    const std::size_t numberOfPoints = std::min((*x_hits).size(), (*z_hits).size());
    if (numberOfPoints < 2) {
      continue;
    }

    // Loop through hits and calculate distances
    for(std::size_t i = 0; i + 1 < numberOfPoints; ++i){
      if(dist<5){
        X.push_back((*x_hits)[i]);
        Z.push_back((*z_hits)[i]);
        //std::cout << (*x) << " "<< (*z) << std::endl;
      }
      else {
        break;
      }//chiudo else

      // Calculate the Euclidean distance between consecutive points
      const double dx = (*x_hits)[i] - (*x_hits)[i + 1];
      const double dz = (*z_hits)[i] - (*z_hits)[i + 1];
      dist += std::sqrt(dx*dx + dz*dz);
      //std::cout << dist << std::endl;

    }//chiudo for on coordinates

    dist=0; // Reset distance for next entry

    graph->Set(0);

    if (X.size() < 2) {
      X.clear();
      Z.clear();
      continue;
    }

    for(std::size_t i = 0; i < X.size(); ++i){
      graph.get()->SetPoint(i,X[i],Z[i]);
    }
    // Fit the graph with a linear function
    graph.get()->Fit(fit.get(),"Q");

    // Save the graph every 300 entries
    if(entry%300 == 0) {
      graph.get()->SetName(Form("Plot%lli",entry));
      graph.get()->SetMarkerStyle(8);
      graph.get()->SetTitle(Form("Example Plot%lli;X axis;Z axis", entry));
      graph.get()->Write();
    }

    double radangl=atan(fit.get()->GetParameter(1));

    if (radangl>=0) radangl = radangl - TMath::Pi()/2;
    else radangl = radangl + TMath::Pi()/2;

    double degangl= radangl * 180.0 / TMath::Pi();
    angDeg=degangl;
    angRad=radangl;
    //if(degangl>0) degangl=degangl-180;
    //else degangl=degangl+180;
    // Fill histogram with the arctangent of the slope from the fit
    if (FContained){
      histoThr->Fill(radangl);
      histoThd->Fill(degangl);
      outTree->Fill();
    }
    //std:: cout << fit.get()->GetParameter(1) << "  " << atan(fit.get()->GetParameter(1)) << "\n";
    // Clear X and Z vectors for the next iteration
    X.clear();
    Z.clear();
  }//chiudo for entries
  fout.get()->cd();
  // Write histogram to file and close the file
  histoThr.get()->Write();
  histoThd.get()->Write();
  outTree->Write();
  fout.get()->Save();
  fout.get()->Close();

  return EXIT_SUCCESS;
}
