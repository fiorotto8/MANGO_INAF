// TrackingAction.cc - Implementation of the TrackingAction class

#include "TrackingAction.hh"
#include "HistoManager.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "TrackingMessenger.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4IonTable.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

// Constructor for the TrackingAction class
// Initializes event action, detector, and time windows
TrackingAction::TrackingAction(EventAction* EA, DetectorConstruction* det)
  : G4UserTrackingAction(),
    fEvent(EA),        // Event action
    fDetector(det),    // Detector construction
    fTrackingMessenger(new TrackingMessenger(this)),
    fFullChain(false)
{
  fTimeWindow1 = fTimeWindow2 = 0.;  // Initialize time windows to 0
}

// Destructor for the TrackingAction class
TrackingAction::~TrackingAction()
{
  delete fTrackingMessenger;
}

// Sets the time window for the tracking action
void TrackingAction::SetTimeWindow(G4double t1, G4double dt)
{
  fTimeWindow1 = t1;         // Set start of the time window
  fTimeWindow2 = t1 + dt;    // Set end of the time window
}

// Function called before the tracking of a particle begins
void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  // Get current run object
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  
  // Get particle properties
  G4ParticleDefinition* particle = track->GetDefinition();
  G4String name = particle->GetParticleName();
  fCharge = particle->GetPDGCharge(); // Particle charge
  fMass = particle->GetPDGMass();     // Particle mass

  // Get track properties
  G4double Ekin = track->GetKineticEnergy(); // Kinetic energy
  G4int ID = track->GetTrackID();            // Track ID

  // Check particle lifetime
  G4double meanLife = particle->GetPDGLifeTime();

  // Count particles and record energy spectrum
  run->ParticleCount(name, Ekin, meanLife);

  G4int ih = 0;  // Histogram index
  if (particle == G4Electron::Electron() || particle == G4Positron::Positron()) {
    ih = 1;
  } else if (particle == G4NeutrinoE::NeutrinoE() || particle == G4AntiNeutrinoE::AntiNeutrinoE()) {
    ih = 2;
  } else if (particle == G4Gamma::Gamma()) {
    ih = 3;
  } else if (particle == G4Alpha::Alpha()) {
    ih = 4;
  } else if (fCharge > 2.) {
    ih = 5;
  }

  // Handle ion particles (charge > 2)
  if (fCharge > 2.) {
    fEvent->RegisterIon(ID, name);

    if (ID == 1) {
      fEvent->AddDecayChain(name);  // Add to decay chain
    } else {
      fEvent->AddDecayChain(" ---> " + name);  // Add subsequent particles
    }
    
    // Set the last decay name in event and detector
    fEvent->SetLastDecay(name);

    // Manage the track status depending on full chain setting
    G4Track* tr = const_cast<G4Track*>(track);
    if (fFullChain) {
      tr->SetKineticEnergy(0.); // Stop the particle
      tr->SetTrackStatus(fStopButAlive);
    } else if (ID > 1) {
      tr->SetTrackStatus(fStopAndKill); // Kill secondary particles
    }

    fTime_birth = track->GetGlobalTime();  // Record birth time of particle
  }
}

// Function called after the tracking of a particle ends
void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  // Only process ions (charge >= 3)
  if (fCharge < 3.) return;

  // Get current run object
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();

  // Get global time and track ID
  G4double time = track->GetGlobalTime();
  G4int ID = track->GetTrackID();
  
  // Record the primary ion's lifetime if it's the primary particle
  if (ID == 1) run->PrimaryTiming(time);

  fTime_end = time;  // Record the end time

  // Process secondary particles for energy and momentum balance
  const std::vector<const G4Track*>* secondaries = track->GetStep()->GetSecondaryInCurrentStep();
  size_t nbtrk = secondaries->size();
  if (nbtrk > 0) {
    G4double EkinTot = 0., EkinVis = 0.;
    G4ThreeVector Pbalance = -track->GetMomentum();
    
    // Loop through secondaries to calculate energy and momentum balance
    for (size_t itr = 0; itr < nbtrk; ++itr) {
      const G4Track* trk = (*secondaries)[itr];
      G4ParticleDefinition* particle = trk->GetDefinition();
      G4double Ekin = trk->GetKineticEnergy();
      
      EkinTot += Ekin;  // Total kinetic energy
      
      G4bool visible = !(particle == G4NeutrinoE::NeutrinoE() || particle == G4AntiNeutrinoE::AntiNeutrinoE());
      if (visible) EkinVis += Ekin;  // Visible kinetic energy

      // Exclude gamma from momentum balance
      if (particle != G4Gamma::Gamma()) Pbalance += trk->GetMomentum();
    }

    G4double Pbal = Pbalance.mag();  // Magnitude of momentum balance
    run->Balance(EkinTot, Pbal);     // Record energy and momentum balance
    fEvent->AddEvisible(EkinVis);    // Record visible energy
  }

  // If no secondaries, end of chain; record total lifetime
  if (nbtrk == 0) {
    run->EventTiming(time);  // Record total lifetime of the particle
    fTime_end = DBL_MAX;     // No end time
  }

  // Check if the particle is within the defined time windows
  run->SetTimeWindow(fTimeWindow1, fTimeWindow2);
  G4String name = track->GetDefinition()->GetParticleName();
  G4bool life1 = false, life2 = false, decay = false;
  
  // Check conditions for time windows and decay
  if ((fTime_birth <= fTimeWindow1) && (fTime_end > fTimeWindow1)) life1 = true;
  if ((fTime_birth <= fTimeWindow2) && (fTime_end > fTimeWindow2)) life2 = true;
  if ((fTime_end > fTimeWindow1) && (fTime_end < fTimeWindow2)) decay = true;
  
  if (life1 || life2 || decay) {
    run->CountInTimeWindow(name, life1, life2, decay);  // Count particle within time window
  }
}
