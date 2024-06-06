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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GenericMessenger.hh"

#include "SensitiveDetector.hh"
#include "globals.hh"

class G4UniformElectricField;
class G4EqMagElectricField;
class G4MagIntegratorStepper;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
  DetectorConstruction();
  ~DetectorConstruction();
  // Other member function declarations
  void SetupElectricField();

  virtual     
  G4VPhysicalVolume* Construct();

  void SetSourceWidth(G4double a) {fSourceWidth=a;}
  G4double GetSourceWidth() {return fSourceWidth;}
  G4PhysicalVolumeStore* GetVolumeStored() {return fPhysVolStore;}
  SensitiveDetector* GetSensitiveDetector(){return fSensitiveDetector;}
  
  private:

  G4double fWorldSize_x;
  G4double fWorldSize_y;
  G4double fWorldSize_z;

  G4double fSourceWidth;

  G4double fModeratorDepth;
  G4double fCollimatorDepth;
  G4double fCollimatorHoleRadius;
  G4double fCollimatorDistance;
  
  G4PhysicalVolumeStore* fPhysVolStore;
  G4LogicalVolume* fLogicalGasVolume;
  SensitiveDetector* fSensitiveDetector;
  virtual void ConstructSDandField();

  G4GenericMessenger* fDetectorMessenger;
  // Member variables related to electric field setup
  G4UniformElectricField* fElectricField;
  G4EqMagElectricField* fEquation;
  G4MagIntegratorStepper* fStepper;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

