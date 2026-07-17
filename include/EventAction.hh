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
/// \file EventAction.hh
/// \brief Definition of the EventAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4GenericMessenger.hh" 

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <map>
#include <vector>

class G4Step;
class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
  EventAction(PrimaryGeneratorAction* );
  ~EventAction();
  
public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);
  
  void AddDecayChain(G4String val) {fDecayChain += val;};
  void AddEvisible(G4double val)   {fEvisTot    += val;};
  void SetElement(G4String el) {fElement = el;}
  void SetLastDecay(G4String lastEl) {fLastDecay=lastEl;}
  G4String GetLastDecsay() {return fLastDecay;}

  // Record beta electrons created directly by radioactive decay, then attach
  // their steps in the sensitive gas.  One acceptance-ntuple row is written
  // per generated beta electron at EndOfEventAction.
  void RecordIon(const G4Track*);
  void RecordBetaGenerated(const G4Track*);
  void RecordSensitiveStep(const G4Step*);
  
private:
  struct BetaRecord {
    G4int trackID = -1;
    G4int parentID = -1;
    G4double generatedEnergy = 0.;
    G4ThreeVector generatedDirection;
    G4bool gasEntered = false;
    G4double entryEnergy = 0.;
    G4ThreeVector entryDirection;
    G4double depositedEnergy = 0.;
    std::vector<G4ThreeVector> hits;
  };

  PrimaryGeneratorAction* fPrimary;
  G4String        fDecayChain;                   
  G4double        fEvisTot;
  G4String        fElement;
  G4String        fLastDecay;
  std::map<G4int, G4String> fIonNames;
  std::map<G4int, BetaRecord> fBetaRecords;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
