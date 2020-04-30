#ifndef SourceCo60_h
#define SourceCo60_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

#include <map>

class G4ParticleGun;
class G4Event;
class MG4DetectorConstruction;

class SourceCo60 : public G4VUserPrimaryGeneratorAction
{
public:
  SourceCo60();    
  virtual ~SourceCo60();
  
  // static access method
  static SourceCo60* GetInstance();
  
  G4double GetEnergy();

  // method from the base class
  virtual void GeneratePrimaries(G4Event*);         
  
  // method to access particle gun
  //  const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
  
private:
  G4ThreeVector GetPosition();
  G4ThreeVector GetDirection(G4ThreeVector& pos);

  static SourceCo60* fgInstance;
  
  //  G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
  G4ParticleDefinition* particle;
  std::map<double, double> m_lines;
  G4double energy;
};

#endif


