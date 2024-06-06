#include "PhysicsList.hh"
#include "G4EmStandardPhysics.hh" 
#include "G4Radioactivation.hh"
#include "G4VMultipleScattering.hh" 
#include "G4ParticleTypes.hh"
#include "G4IonConstructor.hh"
#include "G4NuclideTable.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"  
#include "G4NuclearLevelData.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4ProcessManager.hh"
#include "G4PhysListFactory.hh"

PhysicsList::PhysicsList(): G4VModularPhysicsList()
{

  RegisterPhysics(new G4EmStandardPhysics_option4());
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  
  // mandatory for G4NuclideTable
  //
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(1*ps*std::log(2.));
  G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);

  //read new PhotonEvaporation data set 
  //
  
  G4DeexPrecoParameters* deex = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  deex->SetCorrelatedGamma(false);
  deex->SetStoreAllLevels(true);
  deex->SetIsomerProduction(true);  
  deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()/std::log(2.));


  
  
}

PhysicsList::~PhysicsList()
{}

void PhysicsList::ConstructParticle()
{
  // pseudo-particles
  //G4Geantino::GeantinoDefinition();
  
  // gamma
  G4Gamma::GammaDefinition();

  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  
  // baryons
  G4Proton::ProtonDefinition();
  G4Neutron::NeutronDefinition();  

  // ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();  
}

void PhysicsList::ConstructProcess()
{
  AddTransportation();

  //G4VModularPhysicsList::ConstructProcess();
  
  G4Radioactivation* radioactiveDecay = new G4Radioactivation();
  
  G4bool ARMflag = false;
  radioactiveDecay->SetARM(ARMflag);        //Atomic Rearangement

  // need to initialize atomic deexcitation
  //
  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* deex = man->AtomDeexcitation();

  if (!deex) {
     ///G4EmParameters::Instance()->SetFluo(true);
    G4EmParameters::Instance()->SetAugerCascade(ARMflag);
    G4EmParameters::Instance()->SetDeexcitationIgnoreCut(ARMflag);    
    deex = new G4UAtomicDeexcitation();
    deex->InitialiseAtomicDeexcitation();
    man->SetAtomDeexcitation(deex);
  }

  // register radioactiveDecay
  //
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());
  
  //printout
  //
  G4cout << "\n  Set atomic relaxation mode " << ARMflag << G4endl;
}

