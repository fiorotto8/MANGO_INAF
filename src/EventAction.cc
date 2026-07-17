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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "HistoManager.hh"
#include "Run.hh"
#include "PrimaryGeneratorAction.hh" 

#include "G4Event.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>

namespace {

// The Geant4 geometry uses global y as the drift/GEM normal.  The acceptance
// analysis convention requested for the manuscript is x normal and y-z in the
// GEM/camera plane, hence (x,y,z)_analysis = (y,x,z)_Geant4.
G4ThreeVector ToAcceptanceCoordinates(const G4ThreeVector& vector)
{
  return G4ThreeVector(vector.y(), vector.x(), vector.z());
}

G4double ClampUnit(G4double value)
{
  return std::max(-1., std::min(1., value));
}

G4double Theta(const G4ThreeVector& direction)
{
  return std::acos(ClampUnit(direction.x()));
}

G4double Phi(const G4ThreeVector& direction)
{
  return std::atan2(direction.y(), direction.z());
}

// This is the areAllPointsInsideCylinder definition from analysis/RecoTrack.C,
// applied to the hit sequence of one primary beta electron.  Coordinates here
// remain in the original Geant4 frame because this is a geometry cut.
G4bool IsFullyContained(const std::vector<G4ThreeVector>& hits)
{
  constexpr G4double gasRadius = 36.9*mm;
  constexpr G4double containmentOffset = 5.*mm;
  constexpr G4double cylinderCenterZ = (2.5 + 2.0 + 36.9 + 10.0)*mm;
  constexpr G4double acceptedAngleStart = CLHEP::pi - 30.*CLHEP::pi/180.;
  constexpr G4double acceptedAngleEnd = CLHEP::pi + 30.*CLHEP::pi/180.;

  for (const auto& hit : hits) {
    const G4double dx = hit.x();
    const G4double dz = hit.z() - cylinderCenterZ;
    const G4double radius = std::hypot(dx, dz);
    G4double angle = std::atan2(dx, dz);
    if (angle < 0.) angle += CLHEP::twopi;
    const G4bool inAcceptedAngle =
      angle >= acceptedAngleStart && angle <= acceptedAngleEnd;
    if (!(radius <= gasRadius - containmentOffset || inAcceptedAngle)) {
      return false;
    }
  }
  return !hits.empty();
}

// Orthogonal straight-line fit to the first four hit positions projected onto
// the requested y-z camera plane.  The eigenvector is oriented from hit 1 to
// hit 4 so that atan2(uy,uz) retains the track direction rather than only its
// unoriented axis.
G4bool ReconstructFirstFour(const std::vector<G4ThreeVector>& g4Hits,
                            G4ThreeVector& direction)
{
  if (g4Hits.size() < 4) return false;

  G4double meanY = 0.;
  G4double meanZ = 0.;
  for (std::size_t i = 0; i < 4; ++i) {
    const auto hit = ToAcceptanceCoordinates(g4Hits[i]);
    meanY += hit.y();
    meanZ += hit.z();
  }
  meanY /= 4.;
  meanZ /= 4.;

  G4double syy = 0.;
  G4double szz = 0.;
  G4double syz = 0.;
  for (std::size_t i = 0; i < 4; ++i) {
    const auto hit = ToAcceptanceCoordinates(g4Hits[i]);
    const G4double dy = hit.y() - meanY;
    const G4double dz = hit.z() - meanZ;
    syy += dy*dy;
    szz += dz*dz;
    syz += dy*dz;
  }
  if (syy + szz <= std::numeric_limits<G4double>::epsilon()) return false;

  const G4double phi = 0.5*std::atan2(2.*syz, szz - syy);
  G4double uy = std::sin(phi);
  G4double uz = std::cos(phi);

  const auto first = ToAcceptanceCoordinates(g4Hits[0]);
  const auto fourth = ToAcceptanceCoordinates(g4Hits[3]);
  if (uy*(fourth.y() - first.y()) + uz*(fourth.z() - first.z()) < 0.) {
    uy = -uy;
    uz = -uz;
  }
  direction = G4ThreeVector(0., uy, uz);
  return true;
}

}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(PrimaryGeneratorAction* primary)
:G4UserEventAction(),
 fDecayChain(),
 fEvisTot(0.),
 fPrimary(primary)
{
  // Set default print level 
  G4RunManager::GetRunManager()->SetPrintProgress(10000);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
 fDecayChain = " ";
 fEvisTot = 0.;
 fIonNames.clear();
 fBetaRecords.clear();
 
 fPrimary->GetParticleGun()->SetParticlePosition(fPrimary->GetPointOnSource());

 //fPrimary->GetParticleGun()->GeneratePrimaryVertex(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID(); 
 G4int printProgress = G4RunManager::GetRunManager()->GetPrintProgress();
 //printing survey
 //

 if(evtNb%500==0){
   G4cout << evtNb << G4endl;
 }
 
 if (evtNb%printProgress == 0) 
   G4cout << "    End of event. Decay chain:" << fDecayChain 
          << G4endl << G4endl;

 //total visible energy
 /*G4AnalysisManager::Instance()->FillH1(9, fEvisTot);
 Run* run 
  = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 run->EvisEvent(fEvisTot);
  */

 G4AnalysisManager* analysis = G4AnalysisManager::Instance();
 const G4double nan = std::numeric_limits<G4double>::quiet_NaN();
 for (const auto& item : fBetaRecords) {
   const BetaRecord& record = item.second;
   const G4ThreeVector generated =
     ToAcceptanceCoordinates(record.generatedDirection).unit();
   const G4ThreeVector entry = record.gasEntered
     ? ToAcceptanceCoordinates(record.entryDirection).unit()
     : G4ThreeVector();

   G4ThreeVector reconstructed;
   const G4bool recoAvailable =
     record.gasEntered && ReconstructFirstFour(record.hits, reconstructed);
   const G4bool fullyContained =
     record.gasEntered && IsFullyContained(record.hits);
   // The retained simulation analysis applies containment and then requires a
   // reconstructed direction.  No additional simulation-side selection exists.
   const G4bool selected = fullyContained && recoAvailable;

   const G4double generatedTheta = Theta(generated);
   const G4double generatedPhi = Phi(generated);
   const G4double entryTheta = record.gasEntered ? Theta(entry) : nan;
   const G4double entryPhi = record.gasEntered ? Phi(entry) : nan;
   const G4double recoTheta = recoAvailable ? Theta(reconstructed) : nan;
   const G4double recoPhi = recoAvailable ? Phi(reconstructed) : nan;

   analysis->FillNtupleIColumn(1, 0, evtNb);
   analysis->FillNtupleIColumn(1, 1, record.trackID);
   analysis->FillNtupleIColumn(1, 2, record.parentID);
   analysis->FillNtupleDColumn(1, 3, record.generatedEnergy/keV);
   analysis->FillNtupleDColumn(1, 4, generated.x());
   analysis->FillNtupleDColumn(1, 5, generated.y());
   analysis->FillNtupleDColumn(1, 6, generated.z());
   analysis->FillNtupleDColumn(1, 7, generatedTheta);
   analysis->FillNtupleDColumn(1, 8, generatedPhi);
   analysis->FillNtupleDColumn(1, 9, std::sin(generatedTheta));
   analysis->FillNtupleDColumn(1, 10, std::abs(generated.x()));
   analysis->FillNtupleIColumn(1, 11, record.gasEntered);
   analysis->FillNtupleDColumn(1, 12,
     record.gasEntered ? record.entryEnergy/keV : nan);
   analysis->FillNtupleDColumn(1, 13, record.gasEntered ? entry.x() : nan);
   analysis->FillNtupleDColumn(1, 14, record.gasEntered ? entry.y() : nan);
   analysis->FillNtupleDColumn(1, 15, record.gasEntered ? entry.z() : nan);
   analysis->FillNtupleDColumn(1, 16, entryTheta);
   analysis->FillNtupleDColumn(1, 17, entryPhi);
   analysis->FillNtupleDColumn(1, 18,
     record.gasEntered ? std::sin(entryTheta) : nan);
   analysis->FillNtupleDColumn(1, 19,
     record.gasEntered ? std::abs(entry.x()) : nan);
   analysis->FillNtupleDColumn(1, 20, record.depositedEnergy/keV);
   analysis->FillNtupleIColumn(1, 21, record.hits.size());
   analysis->FillNtupleIColumn(1, 22, fullyContained);
   analysis->FillNtupleIColumn(1, 23, selected);
   analysis->FillNtupleIColumn(1, 24, recoAvailable);
   analysis->FillNtupleDColumn(1, 25,
     recoAvailable ? reconstructed.x() : nan);
   analysis->FillNtupleDColumn(1, 26,
     recoAvailable ? reconstructed.y() : nan);
   analysis->FillNtupleDColumn(1, 27,
     recoAvailable ? reconstructed.z() : nan);
   analysis->FillNtupleDColumn(1, 28, recoTheta);
   analysis->FillNtupleDColumn(1, 29, recoPhi);
   analysis->FillNtupleDColumn(1, 30,
     recoAvailable ? std::sin(recoTheta) : nan);
   analysis->FillNtupleDColumn(1, 31,
     recoAvailable ? std::abs(reconstructed.x()) : nan);
   analysis->AddNtupleRow(1);
 }
}

void EventAction::RecordIon(const G4Track* track)
{
  fIonNames[track->GetTrackID()] =
    track->GetParticleDefinition()->GetParticleName();
}

void EventAction::RecordBetaGenerated(const G4Track* track)
{
  const auto parent = fIonNames.find(track->GetParentID());
  if (parent == fIonNames.end() ||
      (parent->second != "Sr90" && parent->second != "Y90")) {
    return;
  }

  BetaRecord record;
  record.trackID = track->GetTrackID();
  record.parentID = track->GetParentID();
  record.generatedEnergy = track->GetKineticEnergy();
  record.generatedDirection = track->GetMomentumDirection();
  fBetaRecords[record.trackID] = record;
}

void EventAction::RecordSensitiveStep(const G4Step* step)
{
  const G4Track* track = step->GetTrack();
  auto found = fBetaRecords.find(track->GetTrackID());
  if (found == fBetaRecords.end()) return;

  BetaRecord& record = found->second;
  const G4StepPoint* pre = step->GetPreStepPoint();
  if (!record.gasEntered) {
    record.gasEntered = true;
    record.entryEnergy = pre->GetKineticEnergy();
    record.entryDirection = pre->GetMomentumDirection();
  }
  record.depositedEnergy += step->GetTotalEnergyDeposit();
  record.hits.push_back(pre->GetPosition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
