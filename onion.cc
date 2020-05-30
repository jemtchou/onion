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
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "nwlPhysicsList.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4EmLivermorePhysics.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  
  int runNumber = 0;
  
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }
  else
  {
    if(argc>1)
    {
       if(argc == 3)  runNumber = atol(argv[2]);
    }
  }

  //CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED  
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(7);
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new DetectorConstruction());

  G4VModularPhysicsList* physicsList = new nwlPhysicsList;
  physicsList->SetVerboseLevel(0);
  physicsList->ReplacePhysics(new G4EmLivermorePhysics);
  runManager->SetUserInitialization(physicsList);
    
  // Set user action classes
  runManager->SetUserInitialization(new ActionInitialization());

  G4cout << "Run Number " << runNumber << G4endl;

  CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
  G4Random::setTheEngine( &defaultEngine ); 
  G4int seed = G4int(runNumber*1000000 + 123456); 
  G4Random::setTheSeed( seed ); 

  runManager->Initialize();
  
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;
}
