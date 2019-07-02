#ifndef WCSimPrimaryGeneratorAction_h
#define WCSimPrimaryGeneratorAction_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "jhfNtuple.h"
#include <vector>
#include <random>

#include "TTree.h"
#include "TFile.h"

#include <fstream>

#include "WCSimRootOptions.hh"

class WCSimDetectorConstruction;
class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class WCSimPrimaryGeneratorMessenger;
class G4Generator;
class ArbDistFile;

class WCSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  WCSimPrimaryGeneratorAction(WCSimDetectorConstruction*);
  ~WCSimPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);

  // Gun, laser & gps setting calls these functions to fill jhfNtuple and Root tree
  void SetVtx(G4ThreeVector i)     { vtxs[0] = i; nvtxs = 1; };
  void SetBeamEnergy(G4double i, G4int n = 0)   { beamenergies[n] = i;};
  void SetBeamDir(G4ThreeVector i, G4int n = 0) { beamdirs[n] = i;};
  void SetBeamPDG(G4int i, G4int n = 0)         { beampdgs[n] = i;};
  void SetNvtxs(G4int i)     { nvtxs = i; };
  void SetVtxs(G4int i, G4ThreeVector v)     { vtxs[i] = v; };

  // These go with jhfNtuple
  G4int GetVecRecNumber(){return vecRecNumber;}
  G4int GetMode() {return mode;};
  G4int GetNvtxs() {return nvtxs;};
  G4int GetVtxVol(G4int n = 0) {return vtxsvol[n];};
  G4ThreeVector GetVtx(G4int n = 0) {return vtxs[n];}
  G4int GetNpar() {return npar;};
  G4int GetBeamPDG(G4int n = 0) {return beampdgs[n];};
  G4double GetBeamEnergy(G4int n = 0) {return beamenergies[n];};
  G4ThreeVector GetBeamDir(G4int n = 0) {return beamdirs[n];};
  G4int GetTargetPDG(G4int n = 0) {return targetpdgs[n];};
  G4double GetTargetEnergy(G4int n = 0) {return targetenergies[n];};
  G4ThreeVector GetTargetDir(G4int n = 0) {return targetdirs[n];};

  // older ...
  G4double GetNuEnergy() {return nuEnergy;};
  G4double GetEnergy() {return energy;};
  G4double GetXPos() {return xPos;};
  G4double GetYPos() {return yPos;};
  G4double GetZPos() {return zPos;};
  G4double GetXDir() {return xDir;};
  G4double GetYDir() {return yDir;};
  G4double GetZDir() {return zDir;};

  G4String GetGeneratorTypeString();
  
  void SaveOptionsToOutput(WCSimRootOptions * wcopt);

private:
  WCSimDetectorConstruction*      myDetector;
  G4ParticleGun*                  particleGun;
  G4GeneralParticleSource*        MyGPS;  //T. Akiri: GPS to run Laser
  WCSimPrimaryGeneratorMessenger* messenger;

  // Variables set by the messenger
  G4bool   useMulineEvt;
  G4bool   useGunEvt;
  G4bool   useLaserEvt;  //T. Akiri: Laser flag
  G4bool   useGPSEvt;
  G4bool   useArbDistEvt;
  std::fstream inputFile;
  G4String vectorFileName;
  G4bool   GenerateVertexInRock;

  // These go with jhfNtuple
  G4int mode;
  G4int nvtxs;
  G4int vtxsvol[MAX_N_PRIMARIES];
  G4ThreeVector vtxs[MAX_N_PRIMARIES];
  G4int npar;
  G4int beampdgs[MAX_N_PRIMARIES], targetpdgs[MAX_N_PRIMARIES];
  G4ThreeVector beamdirs[MAX_N_PRIMARIES], targetdirs[MAX_N_PRIMARIES];
  G4double beamenergies[MAX_N_PRIMARIES], targetenergies[MAX_N_PRIMARIES];
  G4int vecRecNumber;

  G4double nuEnergy;
  G4double energy;
  G4double xPos, yPos, zPos;
  G4double xDir, yDir, zDir;

  G4int    _counterRock; 
  G4int    _counterCublic;

  G4String arbDistFileName;
  G4double arbDistTimeSmearSigma;
  G4double arbDistWavelengthSmearSigma;
  G4ThreeVector arbDistTranslate;
  G4ThreeVector arbDistRotationAxis;
  G4double arbDistRotationAngle;
  G4double arbDistNumPhotons;
  G4double arbDistNumPhotonsSigma;
  
public:

  inline void SetMulineEvtGenerator(G4bool choice) { useMulineEvt = choice; }
  inline G4bool IsUsingMulineEvtGenerator() { return useMulineEvt; }

  inline void SetGunEvtGenerator(G4bool choice) { useGunEvt = choice; }
  inline G4bool IsUsingGunEvtGenerator()  { return useGunEvt; }

  //T. Akiri: Addition of function for the laser flag
  inline void SetLaserEvtGenerator(G4bool choice) { useLaserEvt = choice; }
  inline G4bool IsUsingLaserEvtGenerator()  { return useLaserEvt; }

  inline void SetGPSEvtGenerator(G4bool choice) { useGPSEvt = choice; }
  inline G4bool IsUsingGPSEvtGenerator()  { return useGPSEvt; }

  inline void SetArbDistEvtGenerator(G4bool choice) { useArbDistEvt = choice; }
  inline G4bool IsUsingArbDistEvtGenerator() { return useArbDistEvt; }
  
  inline void OpenVectorFile(G4String fileName) 
  {
    if ( inputFile.is_open() ) 
      inputFile.close();

    vectorFileName = fileName;
    inputFile.open(vectorFileName, std::fstream::in);

    if ( !inputFile.is_open() ) {
      G4cout << "Vector file " << vectorFileName << " not found" << G4endl;
      exit(-1);
    }
  }
  inline G4bool IsGeneratingVertexInRock() { return GenerateVertexInRock; }
  inline void SetGenerateVertexInRock(G4bool choice) { GenerateVertexInRock = choice; }

  inline void SetArbDistFileName(G4String choice) { arbDistFileName = choice; }
  inline G4String GetArbDistFileName() { return arbDistFileName; }

  inline void SetArbDistTimeSmearSigma(G4double choice) { arbDistTimeSmearSigma = choice; }
  inline G4double GetArbDistTimeSmearSigma() { return arbDistTimeSmearSigma; }

  inline void SetArbDistWavelengthSmearSigma(G4double choice) { arbDistWavelengthSmearSigma = choice; }
  inline G4double GetArbDistWavelengthSmearSigma() { return arbDistWavelengthSmearSigma; }

  inline void SetArbDistTranslateXYZ(G4ThreeVector choice) { arbDistTranslate = choice; }
  inline G4ThreeVector GetArbDistTranslateXYZ() { return arbDistTranslate; }

  inline void SetArbDistRotationAxis(G4ThreeVector choice) { arbDistRotationAxis = choice; }
  inline G4ThreeVector GetArbDistRotationAxis() { return arbDistRotationAxis; }

  inline void SetArbDistRotationAngle(G4double choice) { arbDistRotationAngle = choice; }
  inline G4double GetArbDistRotationAngle() { return arbDistRotationAngle; }

  inline void SetArbDistNumPhotons(G4double choice) { arbDistNumPhotons = choice; }
  inline G4double GetArbDistNumPhotons() { return arbDistNumPhotons; }

  inline void SetArbDistNumPhotonsSigma(G4double choice) { arbDistNumPhotonsSigma = choice; }
  inline G4double GetArbDistNumPhotonsSigma() { return arbDistNumPhotonsSigma; }
  
private:
  ArbDistFile *ad_file;
  G4bool arbDistFileRead;

  
  
};

class ArbDistFile {
public:
  
  struct Data {
    unsigned int index;
    int status, pdg;
    double xpos, ypos, zpos, xdir, ydir, zdir;
    double time, energy, ncdf, polang;
  };

  ArbDistFile(G4String fname);
  G4bool IsValid() { return isValid; }
  G4bool SampleDistribution(double random, Data & result);
  
private:
  G4bool isValid;
  std::vector<Data> dvec;
  size_t dvec_size;
  
};

#endif


