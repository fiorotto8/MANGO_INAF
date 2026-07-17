//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(PrimaryGeneratorAction* kin)
:G4UserRunAction(),
 fPrimary(kin),
 fRun(0) 
{
  fRunMessenger = new G4GenericMessenger(this, "/output/","Output file");
  fRunMessenger->DeclareProperty("OutFile", fOutFileName, "Output file name");

  fOutFileName="output";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run();
  return fRun;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{ 
  // keep run condition
  if (fPrimary) { 
    G4ParticleDefinition* particle 
      = fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
    fRun->SetPrimary(particle, energy);
  }    
      
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //if ( analysisManager->IsActive() ) {
  analysisManager->OpenFile("outfiles/"+fOutFileName+".root");
    //}

  analysisManager->CreateNtuple("Hits","Hits");
  analysisManager->CreateNtupleIColumn("EventNumber");
  analysisManager->CreateNtupleSColumn("ParticleName");
  analysisManager->CreateNtupleIColumn("ParticleID");
  analysisManager->CreateNtupleIColumn("ParticleTag");
  analysisManager->CreateNtupleIColumn("ParentID");
  analysisManager->CreateNtupleDColumn("x_hits");
  analysisManager->CreateNtupleDColumn("y_hits");
  analysisManager->CreateNtupleDColumn("z_hits");
  analysisManager->CreateNtupleDColumn("EnergyDeposit");
  //analysisManager->CreateNtuple("TotalEnergyDeposit","TotalEnergyDeposit");
  analysisManager->CreateNtupleIColumn("VolumeNumber");
  analysisManager->CreateNtupleDColumn("VolumeTraslX");
  analysisManager->CreateNtupleDColumn("VolumeTraslY");
  analysisManager->CreateNtupleDColumn("VolumeTraslZ");
  analysisManager->CreateNtupleSColumn("Nucleus");
  analysisManager->CreateNtupleSColumn("ProcessType");
  
  analysisManager->FinishNtuple(0);

  // One compact row per beta electron created directly by radioactive decay.
  // Directions use the manuscript convention: x is the drift/GEM normal and
  // y-z is the GEM/camera plane.  Energies are stored in keV and angles in rad.
  analysisManager->CreateNtuple("Acceptance", "Sr90 beta acceptance");
  analysisManager->CreateNtupleIColumn("EventID");
  analysisManager->CreateNtupleIColumn("TrackID");
  analysisManager->CreateNtupleIColumn("ParentID");
  analysisManager->CreateNtupleDColumn("GeneratedEnergy_keV");
  analysisManager->CreateNtupleDColumn("GeneratedUx");
  analysisManager->CreateNtupleDColumn("GeneratedUy");
  analysisManager->CreateNtupleDColumn("GeneratedUz");
  analysisManager->CreateNtupleDColumn("GeneratedTheta_rad");
  analysisManager->CreateNtupleDColumn("GeneratedPhi_rad");
  analysisManager->CreateNtupleDColumn("GeneratedSinTheta");
  analysisManager->CreateNtupleDColumn("GeneratedAbsUx");
  analysisManager->CreateNtupleIColumn("GasEntered");
  analysisManager->CreateNtupleDColumn("EntryEnergy_keV");
  analysisManager->CreateNtupleDColumn("Ux");
  analysisManager->CreateNtupleDColumn("Uy");
  analysisManager->CreateNtupleDColumn("Uz");
  analysisManager->CreateNtupleDColumn("Theta_rad");
  analysisManager->CreateNtupleDColumn("Phi_rad");
  analysisManager->CreateNtupleDColumn("SinTheta");
  analysisManager->CreateNtupleDColumn("AbsUx");
  analysisManager->CreateNtupleDColumn("DepositedEnergy_keV");
  analysisManager->CreateNtupleIColumn("HitCount");
  analysisManager->CreateNtupleIColumn("FullyContained");
  analysisManager->CreateNtupleIColumn("Selected");
  analysisManager->CreateNtupleIColumn("RecoAvailable");
  analysisManager->CreateNtupleDColumn("RecoUx");
  analysisManager->CreateNtupleDColumn("RecoUy");
  analysisManager->CreateNtupleDColumn("RecoUz");
  analysisManager->CreateNtupleDColumn("RecoTheta_rad");
  analysisManager->CreateNtupleDColumn("RecoPhi_rad");
  analysisManager->CreateNtupleDColumn("RecoSinTheta");
  analysisManager->CreateNtupleDColumn("RecoAbsUx");
  analysisManager->FinishNtuple(1);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  if (isMaster) fRun->EndOfRun();
            
 //save histograms
 //
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 //if ( analysisManager->IsActive() ) {
 analysisManager->Write();
 analysisManager->CloseFile();
  //} 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
