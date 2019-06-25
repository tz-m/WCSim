#include "WCSimPrimaryGeneratorMessenger.hh"
#include "WCSimPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4ios.hh"

WCSimPrimaryGeneratorMessenger::WCSimPrimaryGeneratorMessenger(WCSimPrimaryGeneratorAction* pointerToAction)
:myAction(pointerToAction)
{
  mydetDirectory = new G4UIdirectory("/mygen/");
  mydetDirectory->SetGuidance("WCSim detector control commands.");

  genCmd = new G4UIcmdWithAString("/mygen/generator",this);
  genCmd->SetGuidance("Select primary generator.");
  //T. Akiri: Addition of laser
  genCmd->SetGuidance(" Available generators : muline, gun, laser, gps, arbdist");
  genCmd->SetParameterName("generator",true);
  genCmd->SetDefaultValue("muline");
  //T. Akiri: Addition of laser
  genCmd->SetCandidates("muline gun laser gps arbdist");

  fileNameCmd = new G4UIcmdWithAString("/mygen/vecfile",this);
  fileNameCmd->SetGuidance("Select the file of vectors.");
  fileNameCmd->SetGuidance(" Enter the file name of the vector file");
  fileNameCmd->SetParameterName("fileName",true);
  fileNameCmd->SetDefaultValue("inputvectorfile");

  ArbDistDirectory = new G4UIdirectory("/mygen/arbdist/");
  ArbDistDirectory->SetGuidance("Parameters for controlling the arbitrary distribution generator.");
  
  ArbDistFileNameCmd = new G4UIcmdWithAString("/mygen/arbdist/fileName",this);
  ArbDistFileNameCmd->SetGuidance("Select the ROOT file containing a TTree called 'source'. The branches and entries of the tree are described in the relevant section of WCSimPrimaryGeneratorAction.cc.");
  ArbDistFileNameCmd->SetParameterName("fileName",true);
  ArbDistFileNameCmd->SetDefaultValue("file.root");

  smearTimeCmd = new G4UIcmdWithADoubleAndUnit("/mygen/arbdist/smearTime",this);
  smearTimeCmd->SetGuidance("(Optional) Smear vertex time by a Gaussian with given sigma.");
  smearTimeCmd->SetParameterName("smearTime",true);
  smearTimeCmd->SetDefaultValue(0.0);
  smearTimeCmd->SetDefaultUnit("ns");
  smearTimeCmd->SetUnitCategory("Time");
  smearTimeCmd->SetUnitCandidates("s ns");

  smearWlCmd = new G4UIcmdWithADoubleAndUnit("/mygen/arbdist/smearWavelength",this);
  smearWlCmd->SetGuidance("(Optional) If simulating optical photons, smear the wavelength by a Gaussian with given sigma.");
  smearWlCmd->SetParameterName("smearWavelength",true);
  smearWlCmd->SetDefaultValue(0.0);
  smearWlCmd->SetDefaultUnit("nm");
  smearWlCmd->SetUnitCategory("Length");
  smearWlCmd->SetUnitCandidates("m nm Ang");

  translateCmd = new G4UIcmdWith3VectorAndUnit("/mygen/arbdist/translateVtx",this);
  translateCmd->SetGuidance("(Optional) Translate vertex in file by a distance x,y,z.");
  translateCmd->SetParameterName("translateX","translateY","translateZ",true);
  translateCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
  translateCmd->SetUnitCategory("Length");
  translateCmd->SetDefaultUnit("cm");

  rotationAxisCmd = new G4UIcmdWith3Vector("/mygen/arbdist/rotationAxis",this);
  rotationAxisCmd->SetGuidance("(Optional) Axis around which to rotate direction vector specified in file.");
  rotationAxisCmd->SetParameterName("l","m","n",true);
  rotationAxisCmd->SetDefaultValue(G4ThreeVector(0.,0.,1.));

  rotationAngleCmd = new G4UIcmdWithADoubleAndUnit("/mygen/arbdist/rotationAngle",this);
  rotationAngleCmd->SetGuidance("(Optional) Angle around which to rotate direction vector specified in file.");
  rotationAngleCmd->SetGuidance("  Angle is defined as counterclockwise when viewing the origin from (l,m,n).");
  rotationAngleCmd->SetParameterName("theta",true);
  rotationAngleCmd->SetDefaultValue(0.0);
  rotationAngleCmd->SetDefaultUnit("deg");
  rotationAngleCmd->SetUnitCategory("Angle");
  rotationAngleCmd->SetUnitCandidates("rad deg");

  numPhotonsCmd = new G4UIcmdWithADouble("/mygen/arbdist/numPhotonsMean",this);
  numPhotonsCmd->SetGuidance("Mean number of opticalphotons to simulate per event.");
  numPhotonsCmd->SetParameterName("Num",true);
  numPhotonsCmd->SetDefaultValue(100.0);

  numPhotonsSigmaCmd = new G4UIcmdWithADouble("/mygen/arbdist/numPhotonsSigma",this);
  numPhotonsSigmaCmd->SetGuidance("Gaussian sigma of number of optical photons to simulate per event.");
  numPhotonsSigmaCmd->SetParameterName("Sigma",true);
  numPhotonsSigmaCmd->SetDefaultValue(10.0);
  
}

WCSimPrimaryGeneratorMessenger::~WCSimPrimaryGeneratorMessenger()
{
  //delete genCmd;
  //delete mydetDirectory;
}

void WCSimPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==genCmd )
  {
    if (newValue == "muline")
    {
      myAction->SetMulineEvtGenerator(true);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetArbDistEvtGenerator(false);
    }
    else if ( newValue == "gun")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(true);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetArbDistEvtGenerator(false);
    }
    else if ( newValue == "laser")   //T. Akiri: Addition of laser
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(true);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetArbDistEvtGenerator(false);
    }
    else if ( newValue == "gps")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(true);
      myAction->SetArbDistEvtGenerator(false);
    }
    else if ( newValue == "arbdist")
      {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetArbDistEvtGenerator(true);
      }
  }

  if( command == fileNameCmd )
  {
    myAction->OpenVectorFile(newValue);
    G4cout << "Input vector file set to " << newValue << G4endl;
  }

  if ( command==ArbDistFileNameCmd )
    { myAction->SetArbDistFileName(newValue); }
  if ( command==smearTimeCmd )
    { myAction->SetArbDistTimeSmearSigma(smearTimeCmd->GetNewDoubleValue(newValue)); }
  if ( command==smearWlCmd )
    { myAction->SetArbDistWavelengthSmearSigma(smearWlCmd->GetNewDoubleValue(newValue)); }
  if ( command==translateCmd )
    { myAction->SetArbDistTranslateXYZ(translateCmd->GetNew3VectorValue(newValue)); }
  if ( command==rotationAxisCmd )
    { myAction->SetArbDistRotationAxis(rotationAxisCmd->GetNew3VectorValue(newValue)); }
  if ( command==rotationAngleCmd )
    { myAction->SetArbDistRotationAngle(rotationAngleCmd->GetNewDoubleValue(newValue)); }
  if ( command==numPhotonsCmd )
    { myAction->SetArbDistNumPhotons(numPhotonsCmd->GetNewDoubleValue(newValue)); }
  if ( command==numPhotonsSigmaCmd )
    { myAction->SetArbDistNumPhotonsSigma(numPhotonsSigmaCmd->GetNewDoubleValue(newValue)); }

}

G4String WCSimPrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String cv;
  
  if( command==genCmd )
  {
    if(myAction->IsUsingMulineEvtGenerator())
      { cv = "muline"; }
    else if(myAction->IsUsingGunEvtGenerator())
      { cv = "gun"; }
    else if(myAction->IsUsingLaserEvtGenerator())
      { cv = "laser"; }   //T. Akiri: Addition of laser
    else if(myAction->IsUsingGPSEvtGenerator())
      { cv = "gps"; }
    else if(myAction->IsUsingArbDistEvtGenerator())
      { cv = "arbdist"; }
  }

  if (myAction->IsUsingArbDistEvtGenerator())
    {
      if ( command==smearTimeCmd )
	{ cv = smearTimeCmd->ConvertToString(myAction->GetArbDistTimeSmearSigma(),"ns"); }
      else if ( command==smearWlCmd )
	{ cv = smearWlCmd->ConvertToString(myAction->GetArbDistWavelengthSmearSigma(),"nm"); }
      else if ( command==translateCmd )
	{ cv = translateCmd->ConvertToString(myAction->GetArbDistTranslateXYZ(),"cm"); }
      else if ( command==rotationAxisCmd )
	{ cv = rotationAxisCmd->ConvertToString(myAction->GetArbDistRotationAxis()); }
      else if ( command==rotationAngleCmd )
	{ cv = rotationAngleCmd->ConvertToString(myAction->GetArbDistRotationAngle(),"deg"); }
      else if ( command==numPhotonsCmd )
	{ cv = numPhotonsCmd->ConvertToString(myAction->GetArbDistNumPhotons()); }
      else if ( command==numPhotonsSigmaCmd )
	{ cv = numPhotonsSigmaCmd->ConvertToString(myAction->GetArbDistNumPhotonsSigma()); }
    }
  
  return cv;
}

