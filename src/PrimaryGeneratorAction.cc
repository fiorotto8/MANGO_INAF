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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "PrimaryGeneratorAction.hh"
#include "G4LogicalVolume.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4VSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* Detector):
  G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fDetector(Detector)
{
  
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("geantino");  
    
  fParticleGun->SetParticleEnergy(0*eV);
  fParticleGun->SetParticlePosition(GetPointOnSource());
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,0.)); //forse (1,0,0)          
  fParticleGun->SetParticleDefinition(particle);

  fPrimaryMessenger = new G4GenericMessenger(this, "/isotope/","Radioactive isotope");
  fPrimaryMessenger->DeclareProperty("AtomicNumber", fZIsotope, "Select atomic number");
  fPrimaryMessenger->DeclareProperty("MassNumber", fAIsotope, "Select mass number");
    
  fZIsotope = 38, fAIsotope = 90;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {  
    
    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.*keV;
    
    G4ParticleDefinition* ion
      = G4IonTable::GetIonTable()->GetIon(fZIsotope,fAIsotope,excitEnergy);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
  }    
    
  //create vertex
  //   
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4ThreeVector PrimaryGeneratorAction::GetPointOnSource(){

  G4double width = (fDetector->GetSourceWidth());
    
  G4VPhysicalVolume* vol = fDetector->GetVolumeStored()->GetVolume("Source"); // get the corresponding random physical volume

  G4VSolid* solid = vol->GetLogicalVolume()->GetSolid(); // get the corresponding solid
  
  G4ThreeVector PointOnSurface = solid->GetPointOnSurface(); // get a random point on the surface

  G4ThreeVector Normal = solid->SurfaceNormal(PointOnSurface); // get the normal to the surface in that point
  
  G4ThreeVector TranslationVolume = vol->GetObjectTranslation(); // get the translation vector of that physical volume

  G4ThreeVector NormDepth = (G4UniformRand()*width)*Normal;
  
  G4ThreeVector Point = TranslationVolume + PointOnSurface - NormDepth; //random point in the random volume as the translation vector + a point on the surface + a random depth 

  /*
  G4cout<< "ElementNumber_____ " << nEl << "\n";
  G4cout<< "Element " << Elements[nEl] << "\n";
  
  G4cout <<"PointOnSurface " << PointOnSurface << "\n";
  G4cout <<"TranslationVolume " << TranslationVolume << G4endl;
  */
  
  return Point;
  
}
