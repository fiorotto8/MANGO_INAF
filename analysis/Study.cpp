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
  // Open a ROOT file specified by the command line argument
  auto f = std::unique_ptr<TFile>( TFile::Open(Form("%s",argv[1]),"r") );

  // Retrieve a TTree object named "elabHits" from the file
  auto mytree = std::unique_ptr<TTree>( (TTree*)f.get()->Get("elabHits") );

  // Setup TTreeReader to read the TTree efficiently
  TTreeReader reader(mytree.get());
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
    // Loop through hits and calculate distances
    for(auto x = (*x_hits).begin(), z = (*z_hits).begin(); x!= (*x_hits).end()-1 || z!=(*z_hits).end()-1; ++x, ++z   ){
      if(dist<5){
        X.push_back((*x));
        Z.push_back((*z));
        //std::cout << (*x) << " "<< (*z) << std::endl;
      }
      else {
        break;
      }//chiudo else

      // Calculate the Euclidean distance between consecutive points
      dist+=sqrt( ((*x)-(*x+1)) * ((*x)-(*x+1)) + ((*z)-(*z+1)) * ((*z)-(*z+1)) );
      //std::cout << dist << std::endl;

    }//chiudo for on coordinates

    dist=0; // Reset distance for next entry

    for(int i=0;i<graph->GetN();i++){
      graph.get()->RemovePoint(i);
    }

    for(int i=0;i<X.size();i++){
      graph.get()->SetPoint(i,X[i],Z[i]);
    }
    // Fit the graph with a linear function
    graph.get()->Fit(fit.get(),"Q");

    // Save the graph every 300 entries
    if(entry%300 == 0) {
      graph.get()->SetName(Form("Plot%lli",entry));
      graph.get()->SetMarkerStyle(8);
      graph.get()->SetTitle("Example Plot%lli;X axis;Z axis");
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
