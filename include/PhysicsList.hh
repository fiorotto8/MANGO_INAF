#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4VUserPhysicsList.hh"  
#include "G4RadioactiveDecayPhysics.hh"
#include "G4DecayPhysics.hh"

class PhysicsList: public G4VModularPhysicsList
{
  
public:
  PhysicsList();
  ~PhysicsList(); 
  
protected:
  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();
};


#endif
