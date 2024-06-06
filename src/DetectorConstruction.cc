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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():G4VUserDetectorConstruction()
{
  fWorldSize_x = 0.3*m;
  fWorldSize_y = 0.3*m;
  fWorldSize_z = 0.3*m;  

  fDetectorMessenger = new G4GenericMessenger(this, "/detector/","Element of the detector");
  fDetectorMessenger->DeclareProperty("ModeratorThickness", fModeratorDepth, "Select moderator thickness");
  fDetectorMessenger->DeclareProperty("CollimatorThickness", fCollimatorDepth, "Select moderator thickness");
  fDetectorMessenger->DeclareProperty("CollimatorHoleRadius", fCollimatorHoleRadius, "Select moderator thickness");
  fDetectorMessenger->DeclareProperty("CollimatorDistance", fCollimatorDistance, "Select moderator thickness");
  
  fModeratorDepth = 0.01*mm;//Moderator is commented
  fCollimatorDepth = 2.0*mm;
  fCollimatorHoleRadius = 2. / 2 *mm;
  fCollimatorDistance = 0.0*mm;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

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

  G4Material* Strontium =
    nist->FindOrBuildMaterial("G4_Sr");

  G4Material* Tungsten =
    nist->FindOrBuildMaterial("G4_W");

  // Define stainless steel using NIST material database
  G4Material* StainlessSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");

  std::vector<G4int> natoms;
  std::vector<G4String> elements;

  elements.push_back("O");
  natoms.push_back(2);

  G4double PMMADensity = 1.190*g/cm3;

  G4Material* PMMA = nist->ConstructNewMaterial("PMMA", elements, natoms, PMMADensity);

  
  //
  //defining detector gas mixture
  //

  G4double aHe = 4.002602*g/mole;
  G4Element* elHe = new G4Element("Helium","He", 2, aHe);

  //density 1394
  G4double aAr = 39.948*g/mole;
  G4Element* elAr = new G4Element("Argon","Ar", 18, aAr);

  G4double aC = 12.0107*g/mole;
  G4Element* elC = new G4Element("Carbon", "C", 6, aC);
  
  G4double aF=18.998*g/mole;
  G4Element* elF = new G4Element("Flourine"  ,"F" , 9., aF);

  G4double He_frac = 0.6;
  G4double CF4_frac = 0.4;
  G4double Ar_frac = 0.;

  //Ar gas
  G4double densityAr = 1394*Ar_frac*g/m3;
  G4double pressureAr = 1*Ar_frac*atmosphere;
  G4double temperatureAr = 300*kelvin;
  G4Material* Ar_gas = new G4Material("Ar_gas", densityAr, 1, kStateGas, temperatureAr, pressureAr);
  Ar_gas->AddElement(elAr, 1);

  //He gas
  G4double densityHe = 162.488*He_frac*g/m3;
  G4double pressureHe = 1*He_frac*atmosphere;
  G4double temperatureHe = 300*kelvin;
  G4Material* He_gas = new G4Material("He_gas", densityHe, 1, kStateGas, temperatureHe, pressureHe);
  He_gas->AddElement(elHe, 1);

  //CF4_gas
  G4double densityCF4 = 3574.736*CF4_frac*g/m3;
  G4double pressureCF4 = 1*CF4_frac*atmosphere;
  G4double temperatureCF4 = 300*kelvin;
  G4Material* CF4_gas = new G4Material("CF4_gas", densityCF4, 2, kStateGas, temperatureCF4, pressureCF4);
  CF4_gas->AddElement(elC, 1);
  CF4_gas->AddElement(elF, 4);

  //CYGNO_gas
  G4double densityMix = He_gas->GetDensity()+CF4_gas->GetDensity()+Ar_gas->GetDensity();
  G4double pressureMix = He_gas->GetPressure()+CF4_gas->GetPressure()+Ar_gas->GetPressure();
  G4double temperatureMix = 300*kelvin;
  G4Material* CYGNO_gas = new G4Material("CYGNO_gas", densityMix, 3, kStateGas, temperatureMix, pressureMix);
  CYGNO_gas->AddMaterial(He_gas, He_gas->GetDensity()/densityMix*100*perCent);
  CYGNO_gas->AddMaterial(CF4_gas, CF4_gas->GetDensity()/densityMix*100*perCent);
  CYGNO_gas->AddMaterial(Ar_gas, Ar_gas->GetDensity()/densityMix*100*perCent);

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
  
  G4LogicalVolume*                         
    logicSource = new G4LogicalVolume(SolidSource,             //its solid
				      Strontium,                    //its material
				      "Source");               //its name

  logicSource->SetVisAttributes(SourceVisAttributes);
  
  G4VPhysicalVolume*                                   
    physiSource = new G4PVPlacement(0,                      //no rotation
				    G4ThreeVector(0,0,InnerSourceContThick/2-fSourceWidth/2),        //at (0,0,0)
				   logicSource,             //its logical volume
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
//Creating the RING with a Trench and filling the trench with Stainless Steel
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

// Create solid geometry for stainless steel filling the trench
G4Tubs* solidSteelFill = new G4Tubs("solidSteelFill", GasRadius , trenchMinRadius, trenchHeight / 2, 0, 360);

// Create logical volume for stainless steel filling
G4LogicalVolume* logicSteelFill = new G4LogicalVolume(solidSteelFill, StainlessSteel, "SteelFill");
// Set visualization attributes for the ring with trench
logicSteelFill->SetVisAttributes(FieldRingAttributes);

// Place the stainless steel fill in the trench
new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicSteelFill, "SteelFill", logicRingWithTrench, false, 0);

// Loop to place multiple instances of the ring with trench in the simulation world
G4int Nrings = 5;
G4double ringSpacing = (GasThickness - 5 * RingThicknessAlongDrift) / 4;

for (int i = -2; i < 3; i++) {
    new G4PVPlacement(rotX,
                      G4ThreeVector(0, -i * ringSpacing - i * RingThicknessAlongDrift + RingThicknessAlongDrift, InnerSourceContThick / 2 + fCollimatorDepth + fCollimatorDistance + GasRadius + radialRingThickness),
                      logicRingWithTrench,
                      "RingWithTrench",
                      logicWorld,
                      false,
                      0);
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
  fLogicalGasVolume->SetSensitiveDetector(fSensitiveDetector);
  
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
