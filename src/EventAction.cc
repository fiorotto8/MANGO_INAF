/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
#include "EventAction.hh"
#include "HistoManager.hh"
#include "Run.hh"
#include "PrimaryGeneratorAction.hh" 
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"
#include <iomanip>


EventAction::EventAction(PrimaryGeneratorAction* primary)
:G4UserEventAction(),
fDecayChain(),
fEvisTot(0.),
fPrimary(primary)
{
  // Set default print level 
  G4RunManager::GetRunManager()->SetPrintProgress(10000);
}

EventAction::~EventAction()
{}


void EventAction::BeginOfEventAction(const G4Event*)
{
  fDecayChain = " ";
  fEvisTot = 0.;
  fIonNames.clear();

  fPrimary->GetParticleGun()->SetParticlePosition(fPrimary->GetPointOnSource());

 //fPrimary->GetParticleGun()->GeneratePrimaryVertex(); 
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID(); 
  G4int printProgress = G4RunManager::GetRunManager()->GetPrintProgress();
  //printing survey
  //

  if(evtNb%1000==0){
    G4cout << evtNb << G4endl;
  }

  if (evtNb%printProgress == 0) 
    G4cout << "    End of event. Decay chain:" << fDecayChain 
            << G4endl << G4endl;

  //total visible energy
  Run* run 
    = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->EvisEvent(fEvisTot);
}

void EventAction::RegisterIon(G4int trackId, const G4String& name)
{
  fIonNames[trackId] = name;
}

G4String EventAction::GetIonName(G4int trackId) const
{
  auto found = fIonNames.find(trackId);
  return found == fIonNames.end() ? G4String() : found->second;
}
