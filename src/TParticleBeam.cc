////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Jun  8 17:18:39 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TParticleBeam.h"

#include "TOSCARSSR.h"
#include "TRandomA.h"

#include <cmath>
#include <algorithm>


// External global random generator
extern TRandomA* gRandomA;



TParticleBeam::TParticleBeam ()
{
  // Default constructor
}




TParticleBeam::TParticleBeam (std::string const& PredefinedBeamType, std::string const& Name, double const Weight)
{
  // Constructor given a predefined beam name
  this->SetPredefinedBeam(PredefinedBeamType);
  this->SetName(Name);
  this->SetWeight(Weight);

}




TParticleBeam::TParticleBeam (std::string const& ParticleType, std::string const& Name, double const Energy, double const Current, double const Weight)
{
  // Constructor given a particle type.

  this->SetParticleType(ParticleType);
  this->SetName(Name);
  fE0 = Energy;

  this->SetCurrent(Current);
  this->SetWeight(Weight);

  SetBetaEmittance(TVector3D(0, 0, 0), TVector2D(0, 0), TVector2D(0, 0), TVector3D(0, 0, 0), 0);
}




TParticleBeam::TParticleBeam (std::string const& ParticleType, std::string const& Name, TVector3D const& X0, TVector3D const& D0, double const Energy, double const Current, double const Weight)
{
  // Constructor given a particle type.
  // Sets the initial time to 0

  this->SetParticleType(ParticleType);
  this->SetName(Name);

  fX0 = X0;
  fU0 = D0.UnitVector();
  fE0 = Energy;
  fT0 = 0;

  this->SetCurrent(Current);
  this->SetWeight(Weight);

  SetBetaEmittance(TVector3D(0, 0, 0), TVector2D(0, 0), TVector2D(0, 0), TVector3D(0, 0, 0), 0);
}




TParticleBeam::TParticleBeam (std::string const& ParticleType, std::string const& Name, TVector3D const& X0, TVector3D const& D0, double const Energy, double const T0, double const Current, double const Charge, double const Mass, double const Weight)
{
  // Constructor given a particle type.

  if (ParticleType == "custom") {
    this->SetParticleTypeCustom(ParticleType, Charge, Mass);
  } else {
    this->SetParticleType(ParticleType);
  }

  this->SetName(Name);
  fX0 = X0;
  fU0 = D0.UnitVector();
  fE0 = Energy;
  fT0 = T0;

  this->SetCurrent(Current);
  this->SetWeight(Weight);

  SetBetaEmittance(TVector3D(0, 0, 0), TVector2D(0, 0), TVector2D(0, 0), TVector3D(0, 0, 0), 0);
}




TParticleBeam::~TParticleBeam ()
{
  // Destruction!
}




void TParticleBeam::SetPredefinedBeam (std::string const& Beam)
{
  // Set the beam to a predefined beam

  std::string BeamU = Beam;
  std::transform(BeamU.begin(), BeamU.end(), BeamU.begin(), ::toupper);

  if (BeamU == "NSLSII") {
    this->SetParticleType("electron");
    this->SetCurrent(0.500);
    this->SetE0(3);
    this->SetU0(TVector3D(0, 0, 1));
    this->SetT0(0);
    this->SetX0(0);
    this->SetBetaEmittance(TVector3D(1, 0, 0), TVector2D(0, 0), TVector2D(0.55e-9, 0.008e-9), TVector3D(0, 0, 0), 3 * 0.001);

  } else if (BeamU == "NSLSII-LONGSTRAIGHT") {
    this->SetParticleType("electron");
    this->SetCurrent(0.500);
    this->SetE0(3);
    this->SetU0(TVector3D(0, 0, 1));
    this->SetT0(0);
    this->SetX0(0);
    this->SetBetaEmittance(TVector3D(1, 0, 0), TVector2D(18, 3.1), TVector2D(0.55e-9, 0.008e-9), TVector3D(0, 0, 0), 3 * 0.001);

  } else if (BeamU == "NSLSII-SHORTSTRAIGHT") {
    this->SetParticleType("electron");
    this->SetCurrent(0.500);
    this->SetE0(3);
    this->SetU0(TVector3D(0, 0, 1));
    this->SetT0(0);
    this->SetX0(0);
    this->SetBetaEmittance(TVector3D(1, 0, 0), TVector2D(1.5, 0.8), TVector2D(0.55e-9, 0.008e-9), TVector3D(0, 0, 0), 3 * 0.001);
  }

  return;
}




void TParticleBeam::SetInitialConditions (double const X, double const Y, double const Z, double const Dx, double const Dy, double const Dz, double const E0, double const T)
{
  // This function takes the initial position, initial direction (normalized or not), initial time, and Energy of a particle

  this->SetInitialConditions( TVector3D(X, Y, Z), TVector3D(Dx, Dy, Dz), E0, T );
  return;
}




void TParticleBeam::SetInitialConditions (TVector3D const& X, TVector3D const& D, double const E0, double const T0)
{
  // Set the initial conditions variables for this particle
  // X - Initial position
  // D - 3D vector of the initial direction of velocity (Magnitude is arbitrary)
  // E0 - The initial energy
  // T0 - The initial time (in [m])

  // Rescale the direction to a unit vector

  this->fX0 = X;
  this->fU0 = D.UnitVector();
  this->fE0 = E0;
  this->fT0 = T0;

  return;
}





void TParticleBeam::SetBetaEmittance (TVector3D const& HorizontalDirection, TVector2D const& Beta, TVector2D const& Emittance, TVector3D const& SigmaAt, double const SigmaEnergyGeV)
{
  // Set beam parameters.  This function likely to be updated with better twiss functions
  // UPDATE: Twiss?

  fBeta = Beta;
  fEmittance = Emittance;


  fSigmaU[0]  = sqrt(Emittance[0] * Beta[0]);
  fSigmaU[1]  = sqrt(Emittance[1] * Beta[1]);
  fSigmaUP[0] = sqrt(Emittance[0] / Beta[0]);
  fSigmaUP[1] = sqrt(Emittance[1] / Beta[1]);

  fHorizontalDirection = HorizontalDirection.UnitVector();
  fSigmaAt             = SigmaAt;
  fSigmaEnergyGeV      = SigmaEnergyGeV;

  // The vertical direction has to be orthogonal to the two other directions
  fVerticalDirection = fU0.Cross(fHorizontalDirection).UnitVector();

  if (fabs(fHorizontalDirection.Dot(this->GetU0())) > 0.00000001) {
    // We're supposed to be doing precision science here!
    throw;
  }


  return;
}




void TParticleBeam::SetSigma (TVector3D const& HorizontalDirection, TVector2D const& SigmaU, TVector2D const& SigmaUP, TVector3D const& SigmaAt, double const SigmaEnergyGeV)
{
  // Set beam parameters.  This function likely to be updated with better twiss functions
  // UPDATE: Twiss?

  fHorizontalDirection = HorizontalDirection.UnitVector();
  fSigmaU              = SigmaU;
  fSigmaUP             = SigmaUP;
  fSigmaAt             = SigmaAt;
  fSigmaEnergyGeV      = SigmaEnergyGeV;

  // The vertical direction has to be orthogonal to the two other directions
  fVerticalDirection = fU0.Cross(fHorizontalDirection).UnitVector();

  if (fabs(fHorizontalDirection.Dot(this->GetU0())) > 0.00000001) {
    // We're supposed to be doing precision science here!
    throw;
  }


  return;
}




TVector3D const& TParticleBeam::GetX0 () const
{
  // Return const reference to the initial position
  return fX0;
}




TVector3D const& TParticleBeam::GetU0 () const
{
  // Get unit vector in the direction of initial velocity
  return fU0;
}




double TParticleBeam::GetE0 () const
{
  // Return initial energy
  return fE0;
}




double TParticleBeam::GetT0 () const
{
  // Return initial time
  return fT0;
}




std::string const& TParticleBeam::GetName () const
{
  // Return the internal name of the particle beam
  return fName;
}




double TParticleBeam::GetWeight () const
{
  // Return the internal weight for the particle beam
  return fWeight;
}




TVector2D TParticleBeam::GetBeta () const
{
  // Return the beta function values
  return fBeta;
}




TVector2D TParticleBeam::GetEmittance () const
{
  // Return the emittance values
  return fEmittance;
}




TVector3D TParticleBeam::GetHorizontalDirection () const
{
  // Return the horizontal direction
  return fHorizontalDirection;
}




TVector3D TParticleBeam::GetVerticalDirection () const
{
  // Return the vertical direction
  return fVerticalDirection;
}




TVector3D TParticleBeam::GetSigmaAt () const
{
  // Return the lattice reference
  return fSigmaAt;
}




void TParticleBeam::SetX0 (TVector3D const& X)
{
  // Set the initial position of particle beam
  fX0 = X;
  return;
}




void TParticleBeam::SetU0 (TVector3D const& U)
{
  // Set the initial direction of particle beam (use a unit vector)
  fU0 = U.UnitVector();
  return;
}




void TParticleBeam::SetE0 (double const En)
{
  // Set energy of the beam [GeV]
  fE0 = En;
  return;
}




void TParticleBeam::SetT0 (double const Time)
{
  // Set initial time of beam (int [m])
  fT0 = Time;
  return;
}





void TParticleBeam::SetSigmaEnergyGeV (double const Sigma)
{
  // Set initial time of beam (int [m])
  fSigmaEnergyGeV = Sigma;
  return;
}




void TParticleBeam::SetName (std::string const& Name)
{
  // Set the name of this beam
  fName = Name;
  return;
}




void TParticleBeam::SetWeight (double const Weight)
{
  // Set the weight for this beam
  if (Weight <= 0) {
    std::cerr << "Weight cannot be <= 0" << std::endl;
    throw;
  }

  fWeight = Weight;
  return;
}




TParticleA TParticleBeam::GetNewParticle (std::string const& IdealOrRandom)
{
  // If you want the ideal particle from the beam definition, otherwise return a random

  // I'll take either case
  std::string idor = IdealOrRandom;
  std::transform(idor.begin(), idor.end(), idor.begin(), ::tolower);

  // The ideal trajectory
  if (idor == "ideal") {
    // Calculate Beta from original beam E0
    double const Gamma = fE0 / TOSCARSSR::kgToGeV(this->GetM());
    double const Beta = sqrt(1.0 - 1.0 / (Gamma * Gamma));

    // Copy this particle and set ideal conditions
    TParticleA NewParticle = (TParticleA) *this;
    NewParticle.SetInitialParticleConditions(fX0, Beta * fU0, fT0);
    return NewParticle;
  }

  // If it's not above, return the default
  return GetNewParticle();
}





TParticleA TParticleBeam::GetNewParticle ()
{
  // Intended to get you a new random particle based on the beam parameters
  // given.

  // UPDATE: not just below, but all
  this->GetTrajectory().Clear();

  // UPDATE: Needs rand for twiss, or other beam configurations...
  // UPDATE: Could also take a python function

  double    ENew = fE0 + fSigmaEnergyGeV * gRandomA->Normal(); // correlated with BNew, not sure how to handle this yet

  double const Gamma = ENew / TOSCARSSR::kgToGeV(this->GetM());
  double const Beta = sqrt(1.0 - 1.0 / (Gamma * Gamma));

  // Distance from t0 location to lattice midpoint
  double const DistanceToMidpoint = (fSigmaAt - fX0).Dot(this->GetU0());
  // UPDATE: ME
  TVector3D XNew = this->GetX0();
  XNew += fHorizontalDirection * fSigmaU[0] * (1 + DistanceToMidpoint) * gRandomA->Normal();
  XNew += fVerticalDirection   * fSigmaU[1] * (1 + DistanceToMidpoint) * gRandomA->Normal();

  TVector3D BetaNew = this->GetU0() * Beta;

  // UPDATE: Rotate about the horizontal and vertical beam axes (arbitrary)
  BetaNew.RotateSelfY(fSigmaUP[0] * gRandomA->Normal());
  BetaNew.RotateSelfX(-fSigmaUP[1] * gRandomA->Normal());

  double    TNew = fT0;


  TParticleA NewParticle = (TParticleA) *this;
  NewParticle.SetInitialParticleConditions(XNew, BetaNew, TNew);


  return NewParticle;
}
