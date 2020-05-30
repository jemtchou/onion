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

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Threading.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction)
{
   G4cout << "THREAD "  << G4Threading::G4GetThreadId() << G4endl;
   G4int id = G4Threading::G4GetThreadId();
   if(id<0) id = 0;
 
   char name[20];
   sprintf(name,"fout_%d.csv",id);
   fout.open(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
   fout.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
 
  G4String name = volume->GetName();

  auto postvolume
    = step->GetPostStepPoint()->GetTouchableHandle()
      ->GetVolume();

  if(postvolume)
   {
      G4String postname = postvolume->GetLogicalVolume()->GetName();
     /* if(name(0,11)!="Detector_1_" && postname(0,11)=="Detector_1_")
      {
        G4int pdg = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
        G4cout << "Enter detector 1 " << pdg << G4endl; 

        if(pdg == 22) // for gamma only
        {
 	   auto energy = step->GetTrack()->GetKineticEnergy();
	   auto cosx = step->GetTrack()->GetMomentumDirection().x();
           auto cosy = step->GetTrack()->GetMomentumDirection().y();
           auto cosz = step->GetTrack()->GetMomentumDirection().z();
           cosx = acos(cosx)*180./CLHEP::pi - 90.;
           cosy = acos(cosy)*180./CLHEP::pi;
           cosz = acos(cosz)*180./CLHEP::pi - 90.;
//           G4cout << energy << "," << cosx << "," << cosy << "," << cosz << G4endl;

     //      fout << energy << ", " << cosx << ", " << cosy << ", " << cosz << std::endl; 
        }
      }*/
     if(name(0,5)=="World" && postname(0,4)=="ICRU")
      {
        G4int pdg = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
        if(pdg == 22) // for gamma only
         { 
            auto energy = step->GetTrack()->GetKineticEnergy();
	    fout << energy << std::endl;
        }
      }
   }

  if(name(0,6)!="Sensor")
     return;

  G4int nCopy = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();

 // G4cout << volume->GetName() << " " << nCopy << G4endl;
      
  // collect energy deposited in this step
  if(nCopy > 0)
  {
   G4double edepStep = step->GetTotalEnergyDeposit();
   if(edepStep > 0.)   
       fEventAction->AddEdep(nCopy, edepStep);  
  }
}

