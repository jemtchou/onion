#include "SourceCo60.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

SourceCo60* SourceCo60::fgInstance = 0;

SourceCo60* SourceCo60::GetInstance()
{
  return fgInstance;
}      

SourceCo60::SourceCo60()
: G4VUserPrimaryGeneratorAction()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  particle = particleTable->FindParticle("gamma");

  fgInstance = this;
}

SourceCo60::~SourceCo60()
{
  fgInstance = 0;
}


G4ThreeVector SourceCo60::GetPosition()
{
   G4double r2 = G4UniformRand();
   G4double r = sqrt(r2)*cm;
   G4double phi = G4UniformRand()*2*CLHEP::pi;
   G4double cost = G4UniformRand()*2-1;
   G4double sint = sqrt(1-cost*cost);
  
   G4ThreeVector position(r*sint*cos(phi), r*sint*sin(phi), r*cost); //x,y,z
   return position; 
}


G4ThreeVector SourceCo60::GetDirection(G4ThreeVector& pos)
{
   // phi = [0,2pi], theta=[0,pi]
  G4double cost = G4UniformRand()*0.348-0.174; // 90+- 10 deg

  G4double sint = sqrt(1-cost*cost);

  G4double phi = G4UniformRand()*0.349-0.175+1.571; // 90 +- 10 deg
  G4ThreeVector dir(sint*cos(phi), sint*sin(phi), cost);

  return dir;
}


void SourceCo60::GeneratePrimaries(G4Event* evt)
{
  //Co-60 : Q=2824 keV
  // Co^60_27 -> Ni^60_28 + (e-) + gamma[1173] + gamma[1332] (100% beta-decay)
  //                        direct transition p ~ 2e-6, deltaJ = 4
  // Ni^60_28 - stable
  double lines[] = {1173.2*keV, 1332.5*keV};
  double intensity[] = {0.999, 1.0};

  int i;
  for (i=0; i<2; i++)
  {
    G4double rand = G4UniformRand();
    if(rand>intensity[i]) continue; // skip lines with intensity <100%

    G4double ptime=0.0;
    G4ThreeVector position = GetPosition();

    G4PrimaryVertex* vertex = new G4PrimaryVertex(position, ptime); 

    G4double mass =  particle->GetPDGMass();

    G4ThreeVector direction = GetDirection(position);

    G4PrimaryParticle* p =
        new G4PrimaryParticle(particle);
    p->SetKineticEnergy( lines[i] );
    p->SetMass( mass );
    p->SetMomentumDirection( direction );
    p->SetCharge( 0. );
    p->SetPolarization(0., 0., 0.);
    vertex->SetPrimary( p );

    evt->AddPrimaryVertex( vertex );
  }
}

