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

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(), 
 fStepLimit(NULL),
 fCheckOverlaps(true)
{

}

DetectorConstruction::~DetectorConstruction()
{
  delete fStepLimit;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

//  G4Material* Tungsten = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

  G4Material* Water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

  G4Element* W = nist->FindOrBuildElement("W");
  G4Element* Ni = nist->FindOrBuildElement("Ni");
  G4Element* Cu = nist->FindOrBuildElement("Cu");
  G4Material* VND = new G4Material("Wol", 19.25*g/cm3, 3);
  VND->AddElement(W,95*perCent);
  VND->AddElement(Ni, 3*perCent);
  VND->AddElement(Cu, 2*perCent);

  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World
  G4Box* worldS
    = new G4Box("world", 2*m, 2*m, 2*m);                                   
  
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,   //its solid
                 Air,      //its material
                 "World"); //its name
  
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 worldLV,         // its logical volume
                 "World",         // its name
                 0,               // its mother  volume
                 false,           // no boolean operations
                 0,               // copy number
                 fCheckOverlaps); // checking overlaps 

  // Shielding
  // 2.25 - 10
  // 4.5  - 100
  // 6.75  - 1000
  // 9.00  - 10000
  // 11.35 - 100000
  G4Sphere* shieldS 
    = new G4Sphere("shield", 1*cm, 12.35*cm, 0, 2*CLHEP::pi, 0, CLHEP::pi);

  G4LogicalVolume* shieldLV
    = new G4LogicalVolume(
                 shieldS,   //its solid
                 VND,      //its material
                 "Shield"); //its name
/*
  G4VPhysicalVolume* shieldPV
    = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 shieldLV,         // its logical volume
                 "Shield",         // its name
                 worldLV,               // its mother  volume
                 false,           // no boolean operations
                 0,               // copy number
                 fCheckOverlaps); // checking overlaps 
*/
  // Detector
  //     
  // Phantom
  //  
  G4Element* O = nist->FindOrBuildElement("O");
  G4Element* C = nist->FindOrBuildElement("C");
  G4Element* H = nist->FindOrBuildElement("H");
  G4Element* N = nist->FindOrBuildElement("N");

  G4Material* icru = new G4Material("ICRU", 1*g/cm3, 4);
  icru->AddElement(O, 0.762);
  icru->AddElement(C, 0.111);
  icru->AddElement(H, 0.101);
  icru->AddElement(N, 0.026);

  G4ThreeVector posICRU = G4ThreeVector(0*cm, 100*cm, 0*cm);
        
  // ICRU sphere       
  G4Sphere* sICRU =    
    new G4Sphere("ICRUsphere", 0., 15.*cm, 0, 360*deg, 0, 180*deg);
                      
  G4LogicalVolume* lICRU =                         
    new G4LogicalVolume(sICRU,         //its solid
                        icru,          //its material
                        "ICRUSphere");           //its name

    new G4PVPlacement(0,                       //no rotation
                    posICRU,                    //at position
                    lICRU,             //its logical volume
                    "ICRUSphere",                //its name
                    worldLV,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);          //overlaps checking

  //     
  // Sensor
  //
  G4ThreeVector posThin = G4ThreeVector(1*cm, -14*cm, 0);
  G4ThreeVector posThick = G4ThreeVector(-1*cm, -14*cm, 0);

  G4Box* sSensor1 =    
    new G4Box("SensorThin", 0.5*cm, 0.5*0.5*cm, 0.5*cm); 
                
  G4LogicalVolume* lSensor1 =                         
    new G4LogicalVolume(sSensor1, icru, "SensorThin");

  new G4PVPlacement(0,                       //no rotation
                    posThin,                    //at position
                    lSensor1,             //its logical volume
                    "SensorThin",                //its name
                    lICRU,                //its mother  volume
                    false,                   //no boolean operation
                    1,                       //copy number
                    fCheckOverlaps);          //overlaps checking

  G4Box* sSensor2 =
    new G4Box("SensorThin", 0.5*cm, 0.5*cm, 0.5*cm);     

  G4LogicalVolume* lSensor2 =   
    new G4LogicalVolume(sSensor2, icru, "SensorThick");

  new G4PVPlacement(0,                       //no rotation
                    posThick,                    //at position
                    lSensor2,             //its logical volume
                    "SensorThick",                //its name
                    lICRU,                //its mother  volume
                    false,                   //no boolean operation
                    2,                       //copy number
                    fCheckOverlaps);          //overlaps checking

  G4cout << "Detector " << lSensor1->GetName()  
         <<  " mass " << lSensor1->GetMass()/kg << " kg " << G4endl;

  G4cout << "Detector " << lSensor2->GetName()             
         <<  " mass " << lSensor2->GetMass()/kg << " kg " << G4endl;



/* 
  G4Sphere* detS 
    = new G4Sphere("Container", 100*cm, 120*cm, 0, 2*CLHEP::pi, 0, CLHEP::pi);

  G4LogicalVolume* detLV = new G4LogicalVolume(detS, Water, "Container");

  G4VPhysicalVolume* detPV = new G4PVPlacement(
                 0, G4ThreeVector(), detLV, "Container",
                 worldLV, false, 0, fCheckOverlaps);  

  // Segments
  auto i = 0;

  auto r_inner = 100*cm;
  auto r_outer = 101*cm;

  for(i=1; i<=20; i++)
  {
     auto angle_start = 90.*deg - 5.73/2*deg;

     char name[15];
     sprintf(name, "Detector_%d_", i);

     G4Sphere* detS
         = new G4Sphere(name, r_inner, r_outer, angle_start, 5.73*deg, angle_start, 5.73*deg);

      G4LogicalVolume* tempLV = new G4LogicalVolume(detS, Water, name);

      G4VPhysicalVolume* detPV = new G4PVPlacement(
                 0, G4ThreeVector(), tempLV, name,
                 detLV, false, i, fCheckOverlaps); 
    
      G4cout << "Detector " << i << " inner " << r_inner << " outer " << r_outer 
	     <<  " mass " << tempLV->GetMass()/kg << " kg " << G4endl;

    r_inner += 1*cm;
    r_outer += 1*cm;
  }
*/

  return worldPV;
}

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors
/*
  G4String trackerChamberSDname = "B2/TrackerChamberSD";
  B2TrackerSD* aTrackerSD = new B2TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name 
  // of "Chamber_LV".
  SetSensitiveDetector("Chamber_LV", aTrackerSD, true);
*/
}
