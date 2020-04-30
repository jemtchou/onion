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
#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4AccumulableManager.hh"

#include <iostream>

RunAction::RunAction() : G4UserRunAction()
{ 
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(1000);     

  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  for(auto i = 1; i<=20; i++)
  {
    fEdep[i] = new G4Accumulable<G4double>(0);
    fEdep2[i] = new G4Accumulable<G4double>(0);
     
    accumulableManager->RegisterAccumulable(fEdep[i]);
    accumulableManager->RegisterAccumulable(fEdep2[i]); 
  }
}

RunAction::~RunAction()
{
  for(auto i = 1; i<=20; i++)
  {
    if (fEdep[i]) delete fEdep[i];
    if (fEdep2[i]) delete fEdep2[i];
  }
}

void RunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

// reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();
}

void RunAction::EndOfRunAction(const G4Run* run)
{
 G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  double edep[21];
  double edep2[21];
  double rms[21];

  for(auto i = 1; i<=20; i++)
  { 
    edep[i]=0;
    edep2[i]=0;
    for(auto k = 1; k<=i; k++) 
    {
//    G4cout << "---> " << i << " " << edep[i] << " " << k << " " << fEdep[k]->GetValue() << G4endl;
      edep[i]  = edep[i] + fEdep[k]->GetValue();
      edep2[i] = edep2[i] + fEdep2[k]->GetValue();
    }
   
    rms[i] = edep2[i] - edep[i]*edep[i]/nofEvents;
    if (rms[i] > 0.) rms[i] = std::sqrt(rms[i]); else rms[i] = 0.;
    if(IsMaster())
       G4cout << "Detector " << i << " " << edep[i] << " " << rms[i] << G4endl;
    else
       G4cout << "  Local: Detector " << i << " " << edep[i] << " " << rms[i] << G4endl;
  }
}

void RunAction::AddEdep(G4int n, G4double edep)
{
  //  std::cout << "XXX " << n << " " <<  fEdep[n] << std::flush << std::endl;
    (*fEdep[n])+=(edep);
    (*fEdep2[n])+=(edep*edep);
}



