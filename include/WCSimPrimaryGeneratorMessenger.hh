#ifndef WCSimPrimaryGeneratorMessenger_h
#define WCSimPrimaryGeneratorMessenger_h 1

class WCSimPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithADouble;

#include "G4UImessenger.hh"
#include "globals.hh"

class WCSimPrimaryGeneratorMessenger: public G4UImessenger
{
 public:
  WCSimPrimaryGeneratorMessenger(WCSimPrimaryGeneratorAction* mpga);
  ~WCSimPrimaryGeneratorMessenger();
  
 public:
  void     SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);
  
 private:
  WCSimPrimaryGeneratorAction* myAction;
  
 private: //commands
  G4UIdirectory*      mydetDirectory;
  G4UIcmdWithAString* genCmd;
  G4UIcmdWithAString* fileNameCmd;
  G4UIdirectory*      ArbDistDirectory;
  G4UIcmdWithAString* ArbDistFileNameCmd;
  G4UIcmdWithADoubleAndUnit* smearTimeCmd;
  G4UIcmdWithADoubleAndUnit* smearWlCmd;
  G4UIcmdWith3VectorAndUnit* translateCmd;
  G4UIcmdWith3Vector* rotationAxisCmd;
  G4UIcmdWithADoubleAndUnit* rotationAngleCmd;
  G4UIcmdWithADouble* numPhotonsCmd;
  G4UIcmdWithADouble* numPhotonsSigmaCmd;
  
};

#endif


