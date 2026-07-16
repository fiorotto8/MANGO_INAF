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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVReplica.hh"
#include "G4VPVParameterisation.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Tubs.hh"
#include "SensitiveDetector.hh"
#include "G4Torus.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fSourceWidth(0.),
    fSourceMaterialZ(95),
    fModeratorDepth(0.01*mm),
    fCollimatorDepth(2.0*mm),
    fCollimatorHoleRadius(1.0*mm),
    fCollimatorDistance(0.0*mm),
    fPhysVolStore(nullptr),
    fLogicalGasVolume(nullptr),
    fLogicalSource(nullptr),
    fSensitiveDetector(nullptr),
    fDetectorMessenger(nullptr),
    fElectricField(nullptr),
    fEquation(nullptr),
    fStepper(nullptr)
{
  fWorldSize_x = 0.3*m;
  fWorldSize_y = 0.3*m;
  fWorldSize_z = 0.3*m;  

  fDetectorMessenger = new G4GenericMessenger(this, "/detector/","Element of the detector");
  fDetectorMessenger->DeclareProperty("ModeratorThickness", fModeratorDepth, "Select moderator thickness");
  fDetectorMessenger->DeclareProperty("CollimatorThickness", fCollimatorDepth, "Select moderator thickness");
  fDetectorMessenger->DeclareProperty("CollimatorHoleRadius", fCollimatorHoleRadius, "Select moderator thickness");
  fDetectorMessenger->DeclareProperty("CollimatorDistance", fCollimatorDistance, "Select moderator thickness");
  fDetectorMessenger->DeclareMethod("SourceMaterialZ",
                                    &DetectorConstruction::SetSourceMaterialZ,
                                    "Set the elemental source material by atomic number");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

void DetectorConstruction::SetSourceMaterialZ(G4int atomicNumber)
{
  if (atomicNumber < 1 || atomicNumber > 98) {
    G4cerr << "Unsupported source material atomic number: "
           << atomicNumber << G4endl;
    return;
  }

  fSourceMaterialZ = atomicNumber;
  if (fLogicalSource) {
    G4Material* material =
      G4NistManager::Instance()->FindOrBuildSimpleMaterial(fSourceMaterialZ);
    fLogicalSource->SetMaterial(material);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  //
  // define a material
  //

  G4NistManager* nist = G4NistManager::Instance();
  
  G4Material* Air =
    nist->FindOrBuildMaterial("G4_AIR"); 

  G4Material* Alluminium =
    nist->FindOrBuildMaterial("G4_Al");

  G4Material* SourceMaterial =
    nist->FindOrBuildSimpleMaterial(fSourceMaterialZ);

  G4Material* Tungsten =
    nist->FindOrBuildMaterial("G4_W");

  // Define Silver using NIST material database
  G4Material* Silver = nist->FindOrBuildMaterial("G4_Ag");

  G4Material* PMMA = nist->FindOrBuildMaterial("G4_PLEXIGLASS");

  
  //
  //defining detector gas mixture
  //

  G4double aHe = 4.002602*g/mole;
  G4Element* elHe = new G4Element("Helium","He", 2, aHe);

  G4double aC = 12.0107*g/mole;
  G4Element* elC = new G4Element("Carbon", "C", 6, aC);
  
  G4double aF=18.998*g/mole;
  G4Element* elF = new G4Element("Fluorine"  ,"F" , 9., aF);

  G4double aS = 32.065*g/mole;
  G4Element* elS = new G4Element("Sulfur","S",16,aS);

  // User knobs
G4double gasPressure    = 0.65 * atmosphere;  // <-- you change only this
G4double gasTemperature = 300*kelvin;

// Fractions (sum to 1.0 in volume/partial pressure sense)
//G4double He_frac  = 0.6;
//G4double CF4_frac = 0.4;
//G4double Ar_frac  = 0.0;
//G4double SF6_frac = 0.0;
//Fraction for NID
G4double He_frac  = 0.590551;
G4double CF4_frac = 0.393701;
G4double SF6_frac = 0.015748;


// Reference densities at 1 atm, 300 K for *pure* gas
// (your numbers were these *times* the fraction already; we separate it)
G4double rhoHe_ref   = 162.488  * g/m3;   // He @1 atm
G4double rhoCF4_ref  = 3574.736 * g/m3;   // CF4 @1 atm
G4double rhoSF6_ref  = 6010.368 * g/m3;   // SF6 @1 atm

// Scale by total pressure and by fraction
// ρ_component = ρ_ref * (gasPressure / 1 atm) * fraction
G4double densityHe   = rhoHe_ref  * (gasPressure/atmosphere) * He_frac;
G4double densityCF4  = rhoCF4_ref * (gasPressure/atmosphere) * CF4_frac;
G4double densitySF6  = rhoSF6_ref * (gasPressure/atmosphere) * SF6_frac;

// Partial pressures (still useful to store)
G4double pressureHe   = gasPressure * He_frac;
G4double pressureCF4  = gasPressure * CF4_frac;
G4double pressureSF6  = gasPressure * SF6_frac;

// Build each component material
auto* He_gas = new G4Material("He_gas",  densityHe, 1,
                              kStateGas, gasTemperature, pressureHe);
He_gas->AddElement(elHe,1);

auto* CF4_gas = new G4Material("CF4_gas", densityCF4, 2,
                               kStateGas, gasTemperature, pressureCF4);
CF4_gas->AddElement(elC,1);
CF4_gas->AddElement(elF,4);

auto* SF6_gas = new G4Material("SF6_gas", densitySF6, 2,
                               kStateGas, gasTemperature, pressureSF6);
SF6_gas->AddElement(elS,1);
SF6_gas->AddElement(elF,6);

// Now the mixture
G4double densityMix   = densityHe + densityCF4 + densitySF6;
G4double pressureMix  = pressureHe + pressureCF4 + pressureSF6;

auto* CYGNO_gas = new G4Material("CYGNO_gas", densityMix, 3,
                                 kStateGas, gasTemperature, pressureMix);

// Add components by MASS FRACTION = componentMass / totalMass
CYGNO_gas->AddMaterial(He_gas,  densityHe  / densityMix);
CYGNO_gas->AddMaterial(CF4_gas, densityCF4 / densityMix);
CYGNO_gas->AddMaterial(SF6_gas, densitySF6 / densityMix);
  //
  //defining kapton
  //

  G4double a = 1.01*g/mole;
  G4Element* elH  = new G4Element("Hydrogen", "H", 1., a);
  
  a = 16.00*g/mole;
  G4Element* elO  = new G4Element("Oxygen"  , "O", 8., a);

  a = 14.01*g/mole;
  G4Element* elN  = new G4Element("Nitrogen", "N", 7., a);

  
  
  G4double densityKapt = 1.413*g/cm3;
  G4Material* Kapton = new G4Material("Kapton", densityKapt, 4);
  Kapton->AddElement(elO,5);
  Kapton->AddElement(elC,22);
  Kapton->AddElement(elN,2);
  Kapton->AddElement(elH,10);


  //     
  // World
  //
  
  G4Box*  
  solidWorld = new G4Box("World",                          //its name
                   fWorldSize_x/2,fWorldSize_y/2,fWorldSize_z/2);//its size
                   
  G4LogicalVolume*                         
  logicWorld = new G4LogicalVolume(solidWorld,             //its solid
                                   CYGNO_gas,                    //its material
                                   "World");               //its name
  G4VPhysicalVolume*                                   
  physiWorld = new G4PVPlacement(0,                      //no rotation
                                 G4ThreeVector(0,0,0),        //at (0,0,0)
                                 logicWorld,             //its logical volume
                                 "World",                //its name
                                 0,                      //its mother  volume
                                 false,                  //no boolean operation
                                 0);                     //copy number


  G4Colour BaseColor(0.5, 0.5, 0.5,0.3);
  G4VisAttributes* baseVisAttributes = new G4VisAttributes(BaseColor);
  baseVisAttributes->SetForceSolid(true);
  
  //
  //defining the source base
  //
  
  G4double outerRad = 25./2*mm;
  G4double InnerSourceContThick = 5*mm;
    
  G4Tubs *OuterTube = new G4Tubs("outerSupport",0,outerRad,InnerSourceContThick/2,0,360);

  G4double innerRad = 5*mm;
  G4double depth = 0.01*mm;
  
  G4Tubs *InnerTube = new G4Tubs("outerSupport",0,innerRad,depth/2,0,360);

  G4VSolid* solidSourceContainer = new G4SubtractionSolid("SourceContainer", OuterTube,InnerTube, 0, G4ThreeVector(0.,0.,InnerSourceContThick/2));

  G4LogicalVolume*                         
    logicSourceContainer = new G4LogicalVolume(solidSourceContainer,             //its solid
					       Alluminium,                    //its material
					       "SourceContainer");               //its name

  logicSourceContainer->SetVisAttributes(baseVisAttributes);
  
  G4VPhysicalVolume*                                   
    physSourceContainer = new G4PVPlacement(0,                      //no rotation
				   G4ThreeVector(),        //at (0,0,0)
				   logicSourceContainer,             //its logical volume
				   "SourceContainer",                //its name
				   logicWorld,                      //its mother  volume
				   false,                  //no boolean operation
				   0);                     //copy number
  
  
  //
  //defining source region
  //

  G4Colour SourceColor = G4Colour::Yellow();
  G4VisAttributes* SourceVisAttributes = new G4VisAttributes(SourceColor);
  SourceVisAttributes->SetForceSolid(true);
  
  //G4double fSourceWidth= depth/2;
  SetSourceWidth(depth/2);
  
  G4Tubs *SolidSource = new G4Tubs("Source",0,innerRad,fSourceWidth/2,0,360);  
  
  fLogicalSource = new G4LogicalVolume(SolidSource,
                                      SourceMaterial,
                                      "Source");

  fLogicalSource->SetVisAttributes(SourceVisAttributes);
  
  G4VPhysicalVolume*                                   
    physiSource = new G4PVPlacement(0,                      //no rotation
				    G4ThreeVector(0,0,InnerSourceContThick/2-fSourceWidth/2),        //at (0,0,0)
				   fLogicalSource,             //its logical volume
				   "Source",                //its name
				   logicWorld,                      //its mother  volume
				   false,                  //no boolean operation
				   0);                     //copy number

  
  //
  //defining moderator
  //
/* 
  //  G4Colour ModeratorColor(0.5, 0.5, 0.5, 0.04);
  //  G4VisAttributes* moderatorVisAttributes = new G4VisAttributes(ModeratorColor);
  //  moderatorVisAttributes->SetForceSolid(true);

  
  G4Tubs *SolidModerator = new G4Tubs("Moderator",0,outerRad,fModeratorDepth/2,0,360);
  
  G4LogicalVolume*                         
    logicModerator = new G4LogicalVolume(SolidModerator,             //its solid
					 Alluminium,                    //its material
					 "Moderator");               //its name
 */
  //logicModerator->SetVisAttributes(moderatorVisAttributes);

  /*
  G4VPhysicalVolume*                                   
    physModerator = new G4PVPlacement(0,                      //no rotation
				      G4ThreeVector(0,0,InnerSourceContThick/2+fModeratorDepth/2),        //at (0,0,0)
				      logicModerator,             //its logical volume
				      "Moderator",                //its name
				      logicWorld,                      //its mother  volume
				      false,                  //no boolean operation
				      0);                     //copy number



  */


/*
  //
  //defining collimator
  //

  G4Colour CollimatorColor(0.45,0.25,0.0,0.5);
  G4VisAttributes* collimatorVisAttributes = new G4VisAttributes(CollimatorColor);
  collimatorVisAttributes->SetForceSolid(true);

  G4double CollimatorRadius = 16./2*mm; 
      
  G4Tubs *OuterCollimator = new G4Tubs("OutCollimator",0,CollimatorRadius,fCollimatorDepth/2,0,360);
  G4Tubs *InnerCollimator = new G4Tubs("InnCollimator",0,fCollimatorHoleRadius,fCollimatorDepth,0,360);

  G4VSolid* SolidCollimator = new G4SubtractionSolid("Collimator", OuterCollimator,InnerCollimator, 0, G4ThreeVector(0.,0.,0.));  
  
  G4LogicalVolume*
    logicCollimator = new G4LogicalVolume(SolidCollimator,             //its solid
					  Tungsten,                    //its material
					 "Collimator");               //its name

  logicCollimator->SetVisAttributes(collimatorVisAttributes);
  
  G4VPhysicalVolume*                                   
    physCollimator = new G4PVPlacement(0,                      //no rotation
				      G4ThreeVector(0,0,InnerSourceContThick/2+fCollimatorDepth/2+fCollimatorDistance),        //at (0,0,0)
				      logicCollimator,             //its logical volume
				      "Collimator",                //its name
				      logicWorld,                      //its mother  volume
				      false,                  //no boolean operation
				       0);                     //copy number

*/
  
  
/*WORKING!!
  //
  //Creating the RING
  //

  G4double radialRingThickness=1.0*cm;
  G4double RingThicknessAlongDrift=0.5*cm;  
  G4double GasRadius = 36.9*mm;
  G4double GasThickness = 50*mm;

  G4RotationMatrix* rotX = new G4RotationMatrix();
  rotX->rotateX(90*degree);

  G4Colour darkGreyColor(0.3, 0.3, 0.3, 1.0);  // Opaque dark grey
  G4VisAttributes* RingVisAttributes = new G4VisAttributes(darkGreyColor);
  RingVisAttributes->SetForceSolid(true);

  G4Tubs* solidRing = new G4Tubs("solidRing",GasRadius,GasRadius+radialRingThickness,RingThicknessAlongDrift/2,0,360);

  G4LogicalVolume*
    logicRing = new G4LogicalVolume(solidRing,             //its solid
					  PMMA,                    //its material
					  "Ring");               //its name
  G4int Nrings = 5;
  G4double ringSpacing = (GasThickness-5*RingThicknessAlongDrift)/4;

  for(int i =-2;i<3;i++){

    G4VPhysicalVolume* phisicalRings = new G4PVPlacement(rotX,
                  G4ThreeVector(0,-i*ringSpacing-i*RingThicknessAlongDrift+RingThicknessAlongDrift,InnerSourceContThick/2+fCollimatorDepth+fCollimatorDistance+GasRadius+radialRingThickness),
                  logicRing,
                  "Ring",
                  logicWorld,
                  false,
                  0);

  }
  logicRing->SetVisAttributes(RingVisAttributes);
 */


//
//Creating the RING with a Trench and filling the trench with Silver
//

// Define dimensions
G4double radialRingThickness = 1.0 * cm;
G4double RingThicknessAlongDrift = 0.5 * cm;  
G4double GasRadius = 36.9 * mm;
G4double GasThickness = 50 * mm;
G4double trenchHeight = 1 * mm;
G4double trenchMinRadius = GasRadius + 1 * mm;

// Rotation matrix for the ring
G4RotationMatrix* rotX = new G4RotationMatrix();
rotX->rotateX(90 * degree);

// Define visualization attributes for the ring
G4Colour darkGreyColor(0.3, 0.3, 0.3, 1.0);  // Opaque dark grey
G4VisAttributes* RingVisAttributes = new G4VisAttributes(darkGreyColor);
RingVisAttributes->SetForceSolid(true);
// Define visualization attributes for the ring
G4Colour lightGreyColor(0.7, 0.7, 0.7, 1.0);  // Opaque light grey
G4VisAttributes* FieldRingAttributes = new G4VisAttributes(lightGreyColor);
FieldRingAttributes->SetForceSolid(true);

// Create solid geometry for the ring
G4Tubs* solidRing = new G4Tubs("solidRing", GasRadius, GasRadius + radialRingThickness, RingThicknessAlongDrift / 2, 0, 360);

// Create solid geometry for the trench
G4Tubs* solidTrench = new G4Tubs("solidTrench", GasRadius - 10 * mm, trenchMinRadius, trenchHeight / 2, 0, 360);

// Combine the ring and the trench
G4VSolid* solidWithTrench = new G4SubtractionSolid("solidWithTrench", solidRing, solidTrench, 0, G4ThreeVector(0, 0, 0));

// Create logical volume for the ring with trench
G4LogicalVolume* logicRingWithTrench = new G4LogicalVolume(solidWithTrench, PMMA, "RingWithTrench");

// Set visualization attributes for the ring with trench
logicRingWithTrench->SetVisAttributes(RingVisAttributes);

// Create solid geometry for Silver filling the trench
G4Tubs* solidSilverFill = new G4Tubs("solidSilverFill", GasRadius , trenchMinRadius, trenchHeight / 2, 0, 360);

// Create logical volume for Silver filling
G4LogicalVolume* logicSilverFill = new G4LogicalVolume(solidSilverFill, Silver, "SilverFill");
// Set visualization attributes for the ring with trench
logicSilverFill->SetVisAttributes(FieldRingAttributes);

// Loop to place multiple instances of the ring with trench in the simulation world
G4int Nrings = 5;
G4double ringSpacing = (GasThickness - 5 * RingThicknessAlongDrift) / 4;

for (int i = -2; i < 3; i++) {
    G4ThreeVector ringPosition(
      0,
      -i * ringSpacing - i * RingThicknessAlongDrift + RingThicknessAlongDrift,
      InnerSourceContThick / 2 + fCollimatorDepth + fCollimatorDistance
        + GasRadius + radialRingThickness);

    new G4PVPlacement(rotX,
                      ringPosition,
                      logicRingWithTrench,
                      "RingWithTrench",
                      logicWorld,
                      false,
                      i + 2);

    new G4PVPlacement(rotX,
                      ringPosition,
                      logicSilverFill,
                      "SilverFill",
                      logicWorld,
                      false,
                      i + 2);
}


  //
  //CF4 sensitive volume
  //
    
    
  G4Colour GasColor(0.0,0.0,1.0,0.5);
  G4VisAttributes* GasVisAttributes = new G4VisAttributes(GasColor);
  //GasVisAttributes->SetForceSolid(true);
  
  G4Tubs* solidGasVolume = new G4Tubs("GasVolume",0,GasRadius,GasThickness/2,0,360);

  fLogicalGasVolume = new G4LogicalVolume(solidGasVolume,
					  CYGNO_gas,
					  "GasVolume"
					  );
  
  fLogicalGasVolume->SetVisAttributes(GasVisAttributes);
  
  
  G4VPhysicalVolume* GasVolume = new G4PVPlacement(rotX,
                  G4ThreeVector(0,RingThicknessAlongDrift,InnerSourceContThick/2+fCollimatorDepth+fCollimatorDistance+GasRadius+radialRingThickness),
                  fLogicalGasVolume,
                  "GasVolume",
                  logicWorld,
                  false,
                  0);
  
    

  
  //
  //Creating the physicalvolumestore
  //
  
  fPhysVolStore = G4PhysicalVolumeStore::GetInstance();
  
  
  //
  //always return the physical World
  //  
  
  return physiWorld;
}


void DetectorConstruction::ConstructSDandField()
{

  fSensitiveDetector = new SensitiveDetector("SensitiveDetector");  
  G4SDManager::GetSDMpointer()->AddNewDetector(fSensitiveDetector);
  fLogicalGasVolume->SetSensitiveDetector(fSensitiveDetector);
  
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
