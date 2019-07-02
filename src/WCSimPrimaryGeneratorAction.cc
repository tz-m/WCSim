#include "WCSimPrimaryGeneratorAction.hh"
#include "WCSimDetectorConstruction.hh"
#include "WCSimPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <fstream>
#include <vector>
#include <string>

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using std::vector;
using std::string;
using std::fstream;

vector<string> tokenize( string separators, string input );

inline vector<string> readInLine(fstream& inFile, int lineSize, char* inBuf)
{
  // Read in line break it up into tokens
  inFile.getline(inBuf,lineSize);
  return tokenize(" $\r", inBuf);
}

inline float atof( const string& str ) {return std::atof( str.c_str() );}
inline int   atoi( const string& str ) {return std::atoi( str.c_str() );}

WCSimPrimaryGeneratorAction::WCSimPrimaryGeneratorAction(
					  WCSimDetectorConstruction* myDC)
  :myDetector(myDC), vectorFileName("")
{
  //T. Akiri: Initialize GPS to allow for the laser use 
  MyGPS = new G4GeneralParticleSource();

  // Initialize to zero
  mode = 0;
  nvtxs = 0;
  for( Int_t u=0; u<50; u++){
    vtxsvol[u] = 0;
    vtxs[u] = G4ThreeVector(0.,0.,0.);
  }
  nuEnergy = 0.;
  _counterRock=0; // counter for generated in Rock
  _counterCublic=0; // counter generated
  
  //---Set defaults. Do once at beginning of session.
  
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.0));

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //G4IonTable* ionTable = G4IonTable::GetIonTable();
  G4String particleName;
  particleGun->
    SetParticleDefinition(particleTable->FindParticle(particleName="mu+"));

  particleGun->
    SetParticlePosition(G4ThreeVector(0.*m,0.*m,0.*m));
    
  messenger = new WCSimPrimaryGeneratorMessenger(this);
  useMulineEvt = true;
  useGunEvt    = false;
  useLaserEvt  = false;
  useGPSEvt    = false;
  useArbDistEvt = false;

  arbDistFileName = "";
  arbDistTimeSmearSigma = 0.0;
  arbDistWavelengthSmearSigma = 0.0;
  arbDistTranslate = G4ThreeVector(0.,0.,0.);
  arbDistRotationAxis = G4ThreeVector(0.,0.,0.);
  arbDistRotationAngle = 0.0;
  arbDistFileRead = false;
  arbDistNumPhotons = 100.0;
  arbDistNumPhotonsSigma = 10.0;
}

WCSimPrimaryGeneratorAction::~WCSimPrimaryGeneratorAction()
{
  if (IsGeneratingVertexInRock()){
    G4cout << "Fraction of Rock volume is : " << G4endl;
      G4cout << " Random number generated in Rock / in Cublic = " 
             << _counterRock << "/" << _counterCublic 
             << " = " << _counterRock/(G4double)_counterCublic << G4endl;
  }
  inputFile.close();
  delete particleGun;
  delete MyGPS;   //T. Akiri: Delete the GPS variable
  delete messenger;
}

void WCSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  // We will need a particle table
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ionTable = G4IonTable::GetIonTable();

  // Temporary kludge to turn on/off vector text format 

  G4bool useNuanceTextFormat = true;


  // Do for every event

  if (useMulineEvt)
  { 

    if ( !inputFile.is_open() )
    {
      G4cout << "Set a vector file using the command /mygen/vecfile name"
	     << G4endl;
      exit(-1);
    }

    //
    // Documentation describing the nuance text format can be found here: 
    // http://neutrino.phy.duke.edu/nuance-format/
    //
    // The format must be strictly adhered to for it to be processed correctly.
    // The lines and their meanings from begin through info are fixed, and then
    // a variable number of tracks may follow.
    //
    if (useNuanceTextFormat)
      {
	const int lineSize=100;
	char      inBuf[lineSize];
	vector<string> token(1);
	
	token = readInLine(inputFile, lineSize, inBuf);
	  
        if (token.size() == 0) 
	  {
	    G4cout << "end of nuance vector file!" << G4endl;
	  }
	else if (token[0] != "begin")
	  {
	    G4cout << "unexpected line begins with " << token[0] << G4endl;
	  }
	else   // normal parsing begins here
	  {
	    // Read the nuance line (ignore value now)

	    token = readInLine(inputFile, lineSize, inBuf);
	    mode = atoi(token[1]);

	    // Read the Vertex line
	    token = readInLine(inputFile, lineSize, inBuf);
	    vtxs[0] = G4ThreeVector(atof(token[1])*cm,
				    atof(token[2])*cm,
				    atof(token[3])*cm);
	    
            // true : Generate vertex in Rock , false : Generate vertex in WC tank
            SetGenerateVertexInRock(false);

	    // Next we read the incoming neutrino and target
	    
	    // First, the neutrino line

	    token=readInLine(inputFile, lineSize, inBuf);
	    beampdgs[0] = atoi(token[1]);
	    beamenergies[0] = atof(token[2])*MeV;
	    beamdirs[0] = G4ThreeVector(atof(token[3]),
					atof(token[4]),
					atof(token[5]));

	    // Now read the target line

	    token=readInLine(inputFile, lineSize, inBuf);
	    targetpdgs[0] = atoi(token[1]);
	    targetenergies[0] = atof(token[2])*MeV;
	    targetdirs[0] = G4ThreeVector(atof(token[3]),
					  atof(token[4]),
					  atof(token[5]));

	    // Read the info line, basically a dummy
	    token=readInLine(inputFile, lineSize, inBuf);
	    G4cout << "Vector File Record Number " << token[2] << G4endl;
            vecRecNumber = atoi(token[2]);
	    
	    // Now read the outgoing particles
	    // These we will simulate.


	    while ( token=readInLine(inputFile, lineSize, inBuf),
		    token[0] == "track" )
	      {
		// We are only interested in the particles
		// that leave the nucleus, tagged by "0"


		if ( token[6] == "0")
		  {
		    G4int pdgid = atoi(token[1]);
		    G4double en = atof(token[2])*MeV;
		    G4ThreeVector dir = G4ThreeVector(atof(token[3]),
						      atof(token[4]),
						      atof(token[5]));
		    std::cout<<"PDGcode "<<pdgid<<"\n";
		    //must handle the case of an ion speratly from other particles
		    //check PDG code if we have an ion.
		    //PDG code format for ions ±10LZZZAAAI
		    char strPDG[11];
		    char strA[10]={0};
		    char strZ[10]={0};
		    

		    long int A=0,Z=0;
		    //		    A=strotl(strPDG,&str);
		    if(abs(pdgid) >= 1000000000)
		      {
			//ion
			sprintf(strPDG,"%i",abs(pdgid));
			strncpy(strZ, &strPDG[3], 3);
			strncpy(strA, &strPDG[6], 3);
			strA[3]='\0';
			strZ[3]='\0';
			A=atoi(strA);
			Z=atoi(strZ);
			G4ParticleDefinition* ion;
			ion =  ionTable->GetIon(Z, A, 0.);
			particleGun->SetParticleDefinition(ion);
			particleGun->SetParticleCharge(0);
		      }
		    else {
		      //not ion
		      particleGun->
			SetParticleDefinition(particleTable->
		      FindParticle(pdgid));
		    }
		    G4double mass = 
		      particleGun->GetParticleDefinition()->GetPDGMass();

		    G4double ekin = en - mass;

		    particleGun->SetParticleEnergy(ekin);
		    //G4cout << "Particle: " << pdgid << " KE: " << ekin << G4endl;
		    particleGun->SetParticlePosition(vtxs[0]);
		    particleGun->SetParticleMomentumDirection(dir);
		    particleGun->GeneratePrimaryVertex(anEvent);
		    SetVtx(vtxs[0]);
		    SetBeamEnergy(ekin);
		    SetBeamDir(dir);
		  }
	      }
	  }
      }
    else 
      {    // old muline format  
	inputFile >> nuEnergy >> energy >> xPos >> yPos >> zPos 
		  >> xDir >> yDir >> zDir;
	
	G4double random_z = ((myDetector->GetWaterTubePosition())
			     - .5*(myDetector->GetWaterTubeLength()) 
			     + 1.*m + 15.0*m*G4UniformRand())/m;
	zPos = random_z;
	G4ThreeVector vtx = G4ThreeVector(xPos, yPos, random_z);
	G4ThreeVector dir = G4ThreeVector(xDir,yDir,zDir);

	particleGun->SetParticleEnergy(energy*MeV);
	particleGun->SetParticlePosition(vtx);
	particleGun->SetParticleMomentumDirection(dir);
	particleGun->GeneratePrimaryVertex(anEvent);
      }
  }

  else if (useGunEvt)
  {      // manual gun operation
    particleGun->GeneratePrimaryVertex(anEvent);

    //To prevent occasional seg fault from an un assigned targetpdg 
    targetpdgs[0] = 2212; //ie. proton

    G4ThreeVector P  =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
    G4ThreeVector vtx=anEvent->GetPrimaryVertex()->GetPosition();
    G4double mass       =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
    G4int pdg        =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();

    char strPDG[11];
    char strA[10]={0};
    char strZ[10]={0};
    
    
    long int A=0,Z=0;
    //		    A=strotl(strPDG,&str);
    if(abs(pdg) >= 1000000000)
      {
	//ion
	sprintf(strPDG,"%i",abs(pdg));
	strncpy(strZ, &strPDG[3], 3);
	strncpy(strA, &strPDG[6], 3);
	strA[3]='\0';
	strZ[3]='\0';
	A=atoi(strA);
	Z=atoi(strZ);

	G4ParticleDefinition* ion   = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
	ion->SetPDGStable(false);
	ion->SetPDGLifeTime(0.);
	
	G4ParticleDefinition* ion2   = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
	std::cout<<"ion2 "<<ion2->GetPDGLifeTime()<<"\n";
      }
    
    
    G4ThreeVector dir  = P.unit();
    G4double E         = std::sqrt((P.dot(P))+(mass*mass));

//     particleGun->SetParticleEnergy(E);
//     particleGun->SetParticlePosition(vtx);
//     particleGun->SetParticleMomentumDirection(dir);

    SetVtx(vtx);
    SetBeamEnergy(E);
    SetBeamDir(dir);
    SetBeamPDG(pdg);
  }
  else if (useLaserEvt)
    {
      targetpdgs[0] = 2212; //ie. proton 
      //T. Akiri: Create the GPS LASER event
      MyGPS->GeneratePrimaryVertex(anEvent);
      
      G4ThreeVector P   =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
      G4ThreeVector vtx =anEvent->GetPrimaryVertex()->GetPosition();
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      //     G4ThreeVector dir  = P.unit();
      G4double E         = std::sqrt((P.dot(P)));
      
      //SetVtx(vtx);
      SetBeamEnergy(E);
      //SetBeamDir(dir);
      SetBeamPDG(pdg);
    }
  else if (useGPSEvt)
    {
      MyGPS->GeneratePrimaryVertex(anEvent);
      
      G4ThreeVector P   =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
      G4ThreeVector vtx =anEvent->GetPrimaryVertex()->GetPosition();
      G4double mass        =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      G4ThreeVector dir  = P.unit();
      G4double E         = std::sqrt((P.dot(P))+(mass*mass));
      
      SetVtx(vtx);
      SetBeamEnergy(E);
      SetBeamDir(dir);
      SetBeamPDG(pdg);
    }
  else if (useArbDistEvt)
    {
      // Use arbitrary (normalised cumulative) probability distribution to define particle
      // release parameters. The input file will be a ROOT TTree (called "source") with branches for:
      //    -- "index/i"  ,  for bookkeeping purposes
      //    -- "status/I" ,  nuance-style status code
      //    -- "pdg/I"    ,  PDG code of particle to release
      //    -- "xpos/D"   ,  x position of release vertex (cm)
      //    -- "ypos/D"   ,  y position of release vertex (cm)
      //    -- "zpos/D"   ,  z position of release vertex (cm)
      //    -- "xdir/D"   ,  x direction cosine (unitless)
      //    -- "ydir/D"   ,  y direction cosine (unitless)
      //    -- "zdir/D"   ,  z direction cosine (unitless)
      //    -- "time/D"   ,  release time (ns)
      //    -- "energy/D" ,  particle total energy (MeV)
      //    -- "ncdf/D"   ,  normalised cumulative probability (unitless)
      //    -- "polang/D" ,  polarisation angle (rad), see Geant4 Physics Reference Manual in section for Optical Photons
      // Each entry in the tree corresponds to one "bin" of the cumulative probability distribution.
      // The process for event generation will go like this:
      //    1. Generate random number, r, on [0,1],
      //    2. Step through normalised cumulative probability distribution bins, i,  searching
      //       for the bin which corresponds to ncdf_i <= r < ncdf_(i+1),
      //    3. Generate particle based on properties in the tree entry which corresponds to
      //       the bin which was found in step 2.
      //
      // Get options from PrimaryGeneratorMessenger:
      //   -- choose ArbDist generator
      //   -- specify filename containing "source" tree with above branches
      //   -- gaussian smear time by sigma(ns) (optional)
      //   -- gaussian smear wavelength by sigma(nm) (optional, if using optical photons)
      //   -- vertex translation x,y,z(cm)
      //   -- direction rotation
      //           --> Rotate the vector (xdir,ydir,zdir) by theta(deg) around the vector (l,m,n).
      //           --> Direction is counterclockwise when looking at the origin from (l,m,n).
      //           --> See https://en.wikipedia.org/wiki/Transformation_matrix#Rotation_2
      //   -- Mean number of photons to simulate per event
      //   -- Sigma number of photons to simulate per event

      std::cout << "Read parameters:" << std::endl;
      std::cout << "  arbDistFileName=" << arbDistFileName << std::endl;
      std::cout << "  arbDistTimeSmearSigma=" << arbDistTimeSmearSigma/ns << " ns,   arbDistWavelengthSmearSigma=" << arbDistWavelengthSmearSigma/nm << " nm" << std::endl;
      std::cout << "  arbDistTranslate=(" << arbDistTranslate.x()/cm << "," << arbDistTranslate.y()/cm << "," << arbDistTranslate.z()/cm << ") cm" << std::endl;
      std::cout << "  arbDistRotationAxis=(" << arbDistRotationAxis.x() << "," << arbDistRotationAxis.y() << "," << arbDistRotationAxis.z() << ")  arbDistRotationAngle=" << arbDistRotationAngle/deg << " deg" << std::endl;
      std::cout << "  arbDistNumPhotonsMean=" << arbDistNumPhotons << "  arbDistNumPhotonsSigma=" << arbDistNumPhotonsSigma << std::endl;
      
      if (!arbDistFileRead)
	{
	  ad_file = new ArbDistFile(arbDistFileName);
	  arbDistFileRead = true;
	}
      if (ad_file->IsValid())
	{
	  G4double numPhots = G4RandGauss::shoot(arbDistNumPhotons,arbDistNumPhotonsSigma);
	  G4int numPhots_Int = std::max(int(0),int(std::round(numPhots)));

	  std::cout << "Firing " << numPhots_Int << " photons." << std::endl;
	  
	  for (G4int i_phot = 0; i_phot < numPhots_Int; ++i_phot)
	    {
	      G4double rand_num = G4UniformRand();
	      ArbDistFile::Data ad_data;
	      G4bool ret = ad_file->SampleDistribution(rand_num,ad_data);
	      if (ret)
		{
		  // do generation here
		  //std::cout << "Particle to create using random number " << rand_num << ":" << std::endl;
		  //std::cout << "  index=" << ad_data.index << "  status=" << ad_data.status << "  pdg=" << ad_data.pdg << std::endl;
		  //std::cout << "  x,y,z=(" << ad_data.xpos << "," << ad_data.ypos << "," << ad_data.zpos << ") cm,   xd,yd,zd=(" << ad_data.xdir << "," << ad_data.ydir << "," << ad_data.zdir << ")" << std::endl;
		  //std::cout << "  time=" << ad_data.time << " ns,  energy=" << ad_data.energy << " MeV,  ncdf=" << ad_data.ncdf << std::endl;
		  
		  G4double t_smear = G4RandGauss::shoot(0,arbDistTimeSmearSigma/ns);
		  G4double l_smear = G4RandGauss::shoot(0,arbDistWavelengthSmearSigma/nm);
		  
		  const G4double hc = 1.239842075191944e-3;// nm*MeV
		  G4double en_smear = hc/l_smear;
		  
		  G4double shoot_time = ad_data.time+t_smear;// ns
		  G4double shoot_energy = hc/((hc/ad_data.energy)+l_smear);// MeV
		  
		  //std::cout << "Smeared particle: " << std::endl;
		  //std::cout << "  time=" << ad_data.time << "+" << t_smear << "=" << ad_data.time+t_smear << " ns,  energy=" << ad_data.energy << "+" << en_smear << "=" << hc/((hc/ad_data.energy)+l_smear) << " MeV,  wavelength=" << hc/ad_data.energy << "+" << l_smear << "=" << (hc/ad_data.energy)+l_smear << " nm" << std::endl;
		  
		  G4ThreeVector dirvec(ad_data.xdir,ad_data.ydir,ad_data.zdir);
		  G4ThreeVector origdirvec = dirvec;
		  G4ThreeVector rotaxis(arbDistRotationAxis.x(),arbDistRotationAxis.y(),arbDistRotationAxis.z());
		  G4double rotangle = arbDistRotationAngle/rad;
		  G4ThreeVector newdir = dirvec.rotate(rotangle,rotaxis);
		  
		  //std::cout << "Rotate axes: " << std::endl;
		  //std::cout << "  old direction=(" << origdirvec.x() << "," << origdirvec.y() << "," << origdirvec.z() << ")  new direction=(" << newdir.x() << "," << newdir.y() << "," << newdir.z() << ")" << std::endl;
		  
		  G4ThreeVector oldvertex(ad_data.xpos,ad_data.ypos,ad_data.zpos);
		  G4ThreeVector newvertex = oldvertex + arbDistTranslate;
		  
		  //std::cout << "Translate vertex: " << std::endl;
		  //std::cout << "  old vertex=(" << oldvertex.x() << "," << oldvertex.y() << "," << oldvertex.z() << ")  new vertex=(" << newvertex.x() << "," << newvertex.y() << "," << newvertex.z() << ")" << std::endl;
		  
		  particleGun->SetParticleTime(shoot_time);
		  if (ad_data.pdg == 30)
		    particleGun->SetParticleDefinition(particleTable->FindParticle("opticalphoton"));
		  else
		    particleGun->SetParticleDefinition(particleTable->FindParticle(ad_data.pdg));
		  G4double mass = particleGun->GetParticleDefinition()->GetPDGMass();
		  G4double ekin = shoot_energy - mass;
		  
		  particleGun->SetParticleEnergy(ekin);
		  particleGun->SetParticleMomentumDirection(newdir);
		  particleGun->SetParticlePosition(newvertex);
		  
		  if (ad_data.pdg == 30)
		    {
		      // I don't pretend to know exactly how this snippet works. It's part of the "standard" code for creating linear polarization based on an angle.
		      G4ThreeVector normal (1., 0., 0.);
		      G4ThreeVector kphoton = particleGun->GetParticleMomentumDirection();
		      G4ThreeVector product = normal.cross(kphoton);
		      G4double modul2       = product*product;
		      
		      G4ThreeVector e_perpend (0., 0., 1.);
		      if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
		      G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
		      
		      G4ThreeVector polar = std::cos(ad_data.polang)*e_paralle + std::sin(ad_data.polang)*e_perpend;
		      particleGun->SetParticlePolarization(polar);
		      //std::cout << "Polarization angle: " << ad_data.polang << " rad (" << ad_data.polang*180.0/M_PI << " deg)" << std::endl;
		      //std::cout << "    Polarization vector: " << polar << std::endl;
		    }
     
		  particleGun->GeneratePrimaryVertex(anEvent);	  
		}
	    }
	} 
    }
}

ArbDistFile::ArbDistFile(G4String fname)
  : isValid(false)
{
  TFile * file = new TFile(fname.c_str(),"READ");
  isValid = file->IsOpen();
  if (isValid)
    {
      std::cout << "Open File Fresh" << std::endl;
      TTree * sourcetree = (TTree*)file->Get("source");
      unsigned int r_index;
      int r_status, r_pdg;
      double r_xpos, r_ypos, r_zpos;
      double r_xdir, r_ydir, r_zdir;
      double r_time, r_energy, r_ncdf, r_polang;
      sourcetree->SetBranchAddress("index",&r_index);
      sourcetree->SetBranchAddress("status",&r_status);
      sourcetree->SetBranchAddress("pdg",&r_pdg);
      sourcetree->SetBranchAddress("xpos",&r_xpos);
      sourcetree->SetBranchAddress("ypos",&r_ypos);
      sourcetree->SetBranchAddress("zpos",&r_zpos);
      sourcetree->SetBranchAddress("xdir",&r_xdir);
      sourcetree->SetBranchAddress("ydir",&r_ydir);
      sourcetree->SetBranchAddress("zdir",&r_zdir);      
      sourcetree->SetBranchAddress("time",&r_time);
      sourcetree->SetBranchAddress("energy",&r_energy);
      sourcetree->SetBranchAddress("ncdf",&r_ncdf);
      sourcetree->SetBranchAddress("polang",&r_polang);
      for (Long64_t i=0; i<sourcetree->GetEntries();i++)
	{
	  sourcetree->GetEntry(i);
	  Data d;
	  d.index = r_index; d.status = r_status; d.pdg = r_pdg;
	  d.xpos = r_xpos; d.ypos = r_ypos; d.zpos = r_zpos;
	  d.xdir = r_xdir; d.ydir = r_ydir; d.zdir = r_zdir;
	  d.time = r_time; d.energy = r_energy; d.ncdf = r_ncdf;
	  d.polang = r_polang;
	  dvec.push_back(d);
	}
      dvec_size = dvec.size();
      std::sort(dvec.begin(),dvec.end(),[](const Data &d1, const Data &d2) {return d1.ncdf<d2.ncdf;});
      file->Close();
    }
}

G4bool ArbDistFile::SampleDistribution(double random, Data & result)
{
  if (!isValid) return false;

  // Binary search algorithm.
  size_t min = 0, max = dvec_size-1;
  if (dvec[min].ncdf <= random && dvec[max].ncdf > random)
    {
      while (max-min != 1)
	{
	  size_t mid = (min+max)/2;
	  if (dvec[mid].ncdf < random)
	    min = mid;
	  else
	    max = mid;
	}
      result = dvec[min];
      return true;
    }
  
  return false;
}

void WCSimPrimaryGeneratorAction::SaveOptionsToOutput(WCSimRootOptions * wcopt)
{
  if(useMulineEvt)
    wcopt->SetVectorFileName(vectorFileName);
  else
    wcopt->SetVectorFileName("");
  wcopt->SetGeneratorType(GetGeneratorTypeString());
}

G4String WCSimPrimaryGeneratorAction::GetGeneratorTypeString()
{
  if(useMulineEvt)
    return "muline";
  else if(useGunEvt)
    return "gun";
  else if(useGPSEvt)
    return "gps";
  else if(useLaserEvt)
    return "laser";
  else if(useArbDistEvt)
    return "arbdist";
  return "";
}

// Returns a vector with the tokens
vector<string> tokenize( string separators, string input ) 
{
  std::size_t startToken = 0, endToken; // Pointers to the token pos
  vector<string> tokens;  // Vector to keep the tokens
  
  if( separators.size() > 0 && input.size() > 0 ) 
    {
    
      while( startToken < input.size() )
	{
	  // Find the start of token
	  startToken = input.find_first_not_of( separators, startToken );
      
	  // If found...
	  if( startToken != input.npos ) 
	    {
	      // Find end of token
	      endToken = input.find_first_of( separators, startToken );
	      if( endToken == input.npos )
		// If there was no end of token, assign it to the end of string
		endToken = input.size();
        
	      // Extract token
	      tokens.push_back( input.substr( startToken, endToken - startToken ) );
        
	      // Update startToken
	      startToken = endToken;
	    }
	}
    }
  
  return tokens;
}

