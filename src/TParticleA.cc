////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Jun  8 14:35:28 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TParticleA.h"
#include "TOSCARSSR.h"

#include <algorithm>
#include <cmath>



TParticleA::TParticleA ()
{
  // Default constructor
}




TParticleA::TParticleA (std::string const& Type)
{
  // Constructor.  This requires a valid type name
  this->SetParticleType(Type);
}




TParticleA::TParticleA (std::string const& Type, TVector3D const& X0, TVector3D const& B0, double const T0)
{
  // Constructor.  Requires a valid type name, also sets initial conditions for this particle
  // UPDATE: check this is needed and not superceeded by ParticleBeam

  this->SetParticleType(Type);
  this->SetX0(X0);
  this->SetB0(B0);
  this->SetT0(T0);

  // Set the "gamma" variable
  SetGamma();
}




TParticleA::~TParticleA ()
{
  // Destroy me
}




void TParticleA::SetParticleType (std::string const& Type)
{
  // If it's a supported / built-in type I should know the charge and mass

  // Save name you give it
  fType = Type;

  // Get a lowercase version
  std::string type = Type;
  std::transform(type.begin(), type.end(), type.begin(), ::tolower);


  // built-in particle types
  // Leptons first
  // UPDATE: values
  if (type == "electron" || type == "anti-positron") {
    this->SetQM( -TOSCARSSR::Qe(), TOSCARSSR::Me() );
  } else if (type == "positron" || type == "anti-electron") {
    this->SetQM(  TOSCARSSR::Qe(), TOSCARSSR::Me() );

  } else if (type == "muon") {
    this->SetQM( -TOSCARSSR::Qe(), TOSCARSSR::GeVTokg(0.1056583715) );
  } else if (type == "anti-muon") {
    this->SetQM(  TOSCARSSR::Qe(), TOSCARSSR::GeVTokg(0.1056583715) );


  } else if (type == "proton") {
    this->SetQM(  TOSCARSSR::Qe(), TOSCARSSR::GeVTokg(0.938272046) );
  } else if (type == "anti-proton") {
    this->SetQM( -TOSCARSSR::Qe(), TOSCARSSR::GeVTokg(0.938272046) );

  } else if (type == "pi+") {
    this->SetQM(  TOSCARSSR::Qe(), TOSCARSSR::GeVTokg(0.13957018) );
  } else if (type == "pi-") {
    this->SetQM( -TOSCARSSR::Qe(), TOSCARSSR::GeVTokg(0.13957018) );

  } else if (type == "custom") {
    // Don't do anything.  Q and M are set elsewhere for custom particles
  } else {
    std::cerr << "TParticleA::SetParticleType type not found: " << Type << std::endl;
    throw;
  }

  return;
}




void TParticleA::SetParticleTypeCustom (std::string const& Type, double const Charge, double const Mass)
{
  // Set a custom particle type.
  // Type - whatever name you want to call it
  // Charge - in SI units [C]
  // Mass   - in SI units [kg]

  // Save name you give it
  fType = Type;


  // Set charge and mass
  this->SetQM(Charge, Mass);

  return;
}




void TParticleA::SetParticleTypeFromPDGID (int const ID)
{
  // UPDATE: If I find an easy way to implement a list for this
  throw;
  return;
}




void TParticleA::SetQ (double const Q)
{
  // Set the charge of this particle
  fQ = Q;

  // Need to adjust variable for fast readback
  SetQoverMGamma();

  return;
}




void TParticleA::SetM (double const M)
{
  // Set mass variable for this particle
  fM = M;

  // Need to adjust variable for fast readback
  SetQoverMGamma();

  return;
}




void TParticleA::SetQM (double const Q, double const M)
{
  // Set charge and mass for this particle
  fQ = Q;
  fM = M;

  // Need to adjust variable for fast readback
  SetQoverMGamma();

  return;
}




void TParticleA::SetX0 (TVector3D const& X0)
{
  // Set initial position in space
  fX0 = X0;

  return;
}




void TParticleA::SetB0 (TVector3D const& B0)
{
  // Set initial Beta (V/c) for this particle
  fB0 = B0;

  SetGamma();

  return;
}




void TParticleA::SetT0 (double const T0)
{
  // Set initial time for this particle
  fT0 = T0;

  return;
}




double TParticleA::GetQ () const
{
  // Return Charge
  return fQ;
}




double TParticleA::GetM () const
{
  // Return mass
  return fM;
}




double TParticleA::GetGamma () const
{
  // Return gamma
  // UPDATE: This should be "inline"
  return fGamma;
}




double TParticleA::GetQoverMGamma () const
{
  // Return computed Q / M / gamma
  // UPDATE: this should be "inline"
  return fQoverMGamma;
}




void TParticleA::SetCurrent (double const Current)
{
  // Set the current.  This is mostly utilized by TParticleBeam
  // UPDATE: Check for negative current?
  if (Current != 0) {
    fCurrent = Current;
  } else {
    fCurrent = this->GetQ();
  }


  return;
}




double TParticleA::GetCurrent () const
{
  // Set the current.  This is mostly utilized by TParticleBeam
  // and SRS functions that need current.  This was just a convenient
  // place to put it
  return fCurrent;
}




TVector3D const& TParticleA::GetX0 () const
{
  // Get initial position
  return fX0;
}




TVector3D const& TParticleA::GetB0 () const
{
  // Get initial Beta (V/c)
  return fB0;
}




double TParticleA::GetE0 () const
{
  // Get initial energy for this particle
  return TOSCARSSR::kgToGeV(fM) * fGamma;
}




double TParticleA::GetT0 () const
{
  // Get initial time
  return fT0;
}




std::string const& TParticleA::GetType () const
{
  // Get reference to type string
  return fType;
}




void TParticleA::SetInitialParticleConditions (TVector3D const& X0, TVector3D const& B0, double const T0)
{
  // Set the initial conditions for this particle

  fX0 = X0;
  fB0 = B0;
  fT0 = T0;

  // Need to adjust variable for fast readback
  SetGamma();

  return;
}




TParticleTrajectoryPoints& TParticleA::GetTrajectory ()
{
  // Get reference to the trajectory member
  return fTrajectory;
}




void TParticleA::SetGamma ()
{
  // Set gamma variable for fast readback
  double const Beta2 = fB0.Mag2();
  if (Beta2 == 1) {
    return;
  }

  fGamma = 1.0 / sqrt(1.0 - Beta2);

  // Need to adjust variable for fast readback
  SetQoverMGamma();

  return;
}




void TParticleA::SetQoverMGamma ()
{
  // Set internal variable for Q / M / gamma for fast readback
  if (GetM() == 0 || GetGamma() == 0) {
    return;
  }

  fQoverMGamma = GetQ() / GetM() / GetGamma();

  return;
}
