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

  fBeamDistribution = kBeamDistribution_None;
}




TParticleBeam::TParticleBeam (std::string const& PredefinedBeamType, std::string const& Name, double const Weight)
{
  fBeamDistribution = kBeamDistribution_None;

  // Constructor given a predefined beam name
  this->SetPredefinedBeam(PredefinedBeamType);
  this->SetName(Name);
  this->SetWeight(Weight);

}




TParticleBeam::TParticleBeam (std::string const& ParticleType, std::string const& Name, double const Energy, double const Current, double const Weight)
{
  // Constructor given a particle type.

  fBeamDistribution = kBeamDistribution_None;

  this->SetParticleType(ParticleType);
  this->SetName(Name);
  fE0 = Energy < TOSCARSSR::kgToGeV(this->GetM()) ? this->GetM() : Energy;

  double const Gamma = this->GetE0() / TOSCARSSR::kgToGeV(this->GetM());
  double const Beta = sqrt(1.0 - 1.0 / (Gamma * Gamma));

  this->SetU0(TVector3D(0, 0, 1));
  this->SetB0(this->GetU0() * Beta);

  this->SetCurrent(Current);
  this->SetWeight(Weight);
}




TParticleBeam::TParticleBeam (std::string const& ParticleType, std::string const& Name, TVector3D const& X0, TVector3D const& D0, double const Energy, double const Current, double const Weight)
{
  // Constructor given a particle type.
  // Sets the initial time to 0

  fBeamDistribution = kBeamDistribution_None;

  this->SetParticleType(ParticleType);
  this->SetName(Name);

  fX0 = X0;
  fU0 = D0.Mag2() > 0 ? D0.UnitVector() : TVector3D(0, 0, 0);
  fE0 = Energy < TOSCARSSR::kgToGeV(this->GetM()) ? this->GetM() : Energy;
  fT0 = 0;

  double const Gamma = this->GetE0() / TOSCARSSR::kgToGeV(this->GetM());
  double const Beta = sqrt(1.0 - 1.0 / (Gamma * Gamma));

  this->SetU0(TVector3D(0, 0, 1));
  this->SetB0(this->GetU0() * Beta);

  this->SetCurrent(Current);
  this->SetWeight(Weight);
}




TParticleBeam::TParticleBeam (std::string const& ParticleType, std::string const& Name, TVector3D const& X0, TVector3D const& D0, double const Energy, double const T0, double const Current, double const Charge, double const Mass, double const Weight)
{
  // Constructor given a particle type.

  fBeamDistribution = kBeamDistribution_None;

  if (ParticleType == "custom") {
    this->SetParticleTypeCustom(ParticleType, Charge, Mass);
  } else {
    this->SetParticleType(ParticleType);
  }

  this->SetName(Name);
  fX0 = X0;
  fU0 = D0.Mag2() > 0 ? D0.UnitVector() : TVector3D(0, 0, 0);
  fE0 = Energy < TOSCARSSR::kgToGeV(this->GetM()) ? this->GetM() : Energy;
  fT0 = T0;

  double const Gamma = this->GetE0() / TOSCARSSR::kgToGeV(this->GetM());
  double const Beta = sqrt(1.0 - 1.0 / (Gamma * Gamma));

  this->SetU0(TVector3D(0, 0, 1));
  this->SetB0(this->GetU0() * Beta);

  this->SetCurrent(Current);
  this->SetWeight(Weight);
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

  if (BeamU == "NSLSII" || BeamU == "NSLS2" || BeamU == "NSLS-II") {
    this->SetParticleType("electron");
    this->SetCurrent(0.500);
    this->SetE0(3);

    double const Gamma = this->GetE0() / TOSCARSSR::kgToGeV(this->GetM());
    double const Beta = sqrt(1.0 - 1.0 / (Gamma * Gamma));

    this->SetU0(TVector3D(0, 0, 1));
    this->SetB0(this->GetU0() * Beta);
    this->SetT0(0);
    this->SetX0(0);
    this->SetEmittance(TVector2D(0.55e-9, 0.008e-9));
    this->SetTwissLatticeReference(TVector3D(0, 0, 0));
    this->SetTwissBetaAlpha(TVector2D(1.5, 0.8), TVector2D(0, 0));
    this->SetSigmaEnergyGeV(0);
    this->SetVerticalDirection(TVector3D(0, 1, 0));
    this->SetBeamDistribution(kBeamDistribution_Filament);

  } else if (BeamU == "NSLSII-LONGSTRAIGHT" || BeamU == "NSLS2-LONGSTRAIGHT" || BeamU == "NSLS-II-LONGSTRAIGHT") {
    this->SetParticleType("electron");
    this->SetCurrent(0.500);
    this->SetE0(3);

    double const Gamma = this->GetE0() / TOSCARSSR::kgToGeV(this->GetM());
    double const Beta = sqrt(1.0 - 1.0 / (Gamma * Gamma));

    this->SetU0(TVector3D(0, 0, 1));
    this->SetB0(this->GetU0() * Beta);
    this->SetT0(0);
    this->SetX0(TVector3D(0, 0, 0));
    this->SetEmittance(TVector2D(0.55e-9, 0.008e-9));
    this->SetTwissLatticeReference(TVector3D(0, 0, 0));
    this->SetTwissBetaAlpha(TVector2D(18, 3.1), TVector2D(0, 0));
    this->SetSigmaEnergyGeV(3. * 0.00089);
    this->SetVerticalDirection(TVector3D(0, 1, 0));
    this->SetBeamDistribution(kBeamDistribution_Gaussian);

  } else if (BeamU == "NSLSII-SHORTSTRAIGHT" || BeamU == "NSLS2-SHORTSTRAIGHT" || BeamU == "NSLS-II-SHORTSTRAIGHT") {
    this->SetParticleType("electron");
    this->SetCurrent(0.500);
    this->SetE0(3);

    double const Gamma = this->GetE0() / TOSCARSSR::kgToGeV(this->GetM());
    double const Beta = sqrt(1.0 - 1.0 / (Gamma * Gamma));

    this->SetU0(TVector3D(0, 0, 1));
    this->SetB0(this->GetU0() * Beta);
    this->SetT0(0);
    this->SetX0(0);
    this->SetEmittance(TVector2D(0.55e-9, 0.008e-9));
    this->SetTwissLatticeReference(TVector3D(0, 0, 0));
    this->SetTwissBetaAlpha(TVector2D(1.5, 0.8), TVector2D(0, 0));
    this->SetSigmaEnergyGeV(3. * 0.00089);
    this->SetVerticalDirection(TVector3D(0, 1, 0));
    this->SetBeamDistribution(kBeamDistribution_Gaussian);
  } else {
    throw std::invalid_argument("no beam by that name found");
  }

  return;
}




void TParticleBeam::SetInitialConditions (double const X, double const Y, double const Z, double const Dx, double const Dy, double const Dz, double const E0, double const T)
{
  // This function takes the initial position, initial direction (normalized or not), initial time, and Energy of a particle

  this->SetInitialConditions( TVector3D(X, Y, Z), TVector3D(Dx, Dy, Dz), E0, T );
  return;
}




void TParticleBeam::SetTwissBetaAlpha (TVector2D const& Beta,
                                       TVector2D const& Alpha,
                                       TVector3D const& Reference,
                                       bool const HasReference)
{
  // Beta and Gamma values at the reference point
  // index as follows: 0-Horizontal, 1-vertical

  // Check that beta is not zero
  if (Beta[0] <= 0 || Beta[1] <= 0) {
    throw std::out_of_range("Beta cannot be <= 0");
  }

  fTwissBeta = Beta;
  fTwissGamma = TVector2D( (1. + Alpha[0] * Alpha[0]) / Beta[0], (1. + Alpha[1] * Alpha[1]) / Beta[1] );
  fTwissAlpha = Alpha;

  if (HasReference) {
    fTwissLatticeReference = Reference;
  }

  SetTwissParametersAtX0();

  return;
}




void TParticleBeam::SetTwissBetaGamma (TVector2D const& Beta,
                                       TVector2D const& Gamma,
                                       TVector3D const& Reference,
                                       bool const HasReference)
{
  // Beta and Gamma values at the reference point
  // index as follows: 0-Horizontal, 1-vertical

  fTwissBeta = Beta;
  fTwissAlpha = TVector2D( sqrt(Gamma[0] * Beta[0] - 1), sqrt(Gamma[1] * Beta[1] - 1) );
  fTwissGamma = Gamma;

  if (HasReference) {
    fTwissLatticeReference = Reference;
  }

  SetTwissParametersAtX0();

  return;
}





void TParticleBeam::SetTwissAlphaGamma (TVector2D const& Alpha,
                                        TVector2D const& Gamma,
                                        TVector3D const& Reference,
                                        bool const HasReference)
{
  // Beta and Gamma values at the reference point
  // index as follows: 0-Horizontal, 1-vertical

  // Check that gamma is not 0
  if (Gamma[0] == 0 || Gamma[1] == 0) {
    throw std::out_of_range("Gamma cannot be <= 0");
  }
  fTwissBeta = TVector2D( (1. + Alpha[0] * Alpha[0]) / Gamma[0], (1. + Alpha[1] * Alpha[1]) / Gamma[1] );
  fTwissAlpha = Alpha;
  fTwissGamma = Gamma;

  if (HasReference) {
    fTwissLatticeReference = Reference;
  }

  SetTwissParametersAtX0();

  return;
}




void TParticleBeam::SetEmittance (TVector2D const& Emittance)
{
  // Set the horizontal and vertical emittance

  fEmittance = Emittance;
  return;
}




void TParticleBeam::SetTwissParametersAtX0 ()
{
  // Transform twiss parameters according to drift space.

  // Calculate distance to transform
  double const D = (fX0 - fTwissLatticeReference).Mag();
  double const L = (fX0 - fTwissLatticeReference).Dot(fU0) >= 0 ? D : -1. * D;

  fTwissBetaX0 = fTwissBeta - 2 * L * fTwissAlpha + L * L * fTwissGamma;
  fTwissAlphaX0 = fTwissAlpha - L * fTwissGamma;
  fTwissGammaX0 = fTwissGamma;

  return;
}




void TParticleBeam::SetTwissParameters (TVector2D const& Beta,
                                        TVector2D const& Alpha,
                                        TVector2D const& Gamma,
                                        TVector3D const& Reference,
                                        bool const HasReference)
{
  fTwissBeta = Beta;
  fTwissAlpha = Alpha;
  fTwissGamma = Gamma;

  if (HasReference) {
    fTwissLatticeReference = Reference;
  }

  SetTwissParametersAtX0();

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
  this->fU0 = D.Mag2() > 0.001 ? D.UnitVector() : TVector3D(0, 0, 1);
  this->fE0 = E0 < TOSCARSSR::kgToGeV(this->GetM()) ? this->GetM() : E0;
  this->fT0 = T0;

  return;
}





void TParticleBeam::SetHorizontalDirection (TVector3D const& HorizontalDirection)
{
  // Set the horizontal direction.  When you do this it resets the vertical direction
  fHorizontalDirection = HorizontalDirection.UnitVector();
  fVerticalDirection = - HorizontalDirection.Cross(fU0).UnitVector();
  return;
}




void TParticleBeam::SetTwissLatticeReference (TVector3D const& L)
{
  fTwissLatticeReference = L;

  // If you change the reference need to update twiss parameters at X0
  SetTwissParametersAtX0();

  return;
}




TVector3D TParticleBeam::GetTwissLatticeReference () const
{
  return fTwissLatticeReference;
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




TVector2D TParticleBeam::GetTwissBeta () const
{
  // Return the beta function values
  return fTwissBeta;
}




TVector2D TParticleBeam::GetTwissAlpha () const
{
  // Return the alpha function values
  return fTwissAlpha;
}




TVector2D TParticleBeam::GetTwissGamma () const
{
  // Return the gamma function values
  return fTwissGamma;
}




TVector2D TParticleBeam::GetTwissBetaX0 () const
{
  // Return the beta function values at X0
  return fTwissBetaX0;
}




TVector2D TParticleBeam::GetTwissAlphaX0 () const
{
  // Return the alpha function values at X0
  return fTwissAlphaX0;
}




TVector2D TParticleBeam::GetTwissGammaX0 () const
{
  // Return the gamma function values at X0
  return fTwissGammaX0;
}




TVector2D TParticleBeam::GetEmittance () const
{
  // Return the emittance values
  return fEmittance;
}




void TParticleBeam::SetEta (TVector2D const& Eta)
{
  // Set eta value
  fEta = Eta;
  return;
}




TVector2D TParticleBeam::GetEta () const
{
  // Return the eta values
  return fEta;
}




TVector3D TParticleBeam::GetHorizontalDirection () const
{
  // Return the horizontal direction
  return fHorizontalDirection;
}




void TParticleBeam::SetVerticalDirection (TVector3D const& VerticalDirection)
{
  fVerticalDirection = VerticalDirection.UnitVector();
  fHorizontalDirection = fVerticalDirection.Cross(fU0).UnitVector();

  return;
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

  // Update TwissX0 for this change
  SetTwissParametersAtX0();

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
  fE0 = En < TOSCARSSR::kgToGeV(this->GetM()) ? this->GetM() : En;
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
  // Set beam energy spread in [GeV]
  fSigmaEnergyGeV = Sigma;
  return;
}




double TParticleBeam::GetSigmaEnergyGeV () const
{
  return fSigmaEnergyGeV;
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
    throw std::runtime_error("Particle weight cannot be negative");
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
    double const Gamma = fE0 / TOSCARSSR::kgToGeV(this->GetM()) < 1 ? 1 : fE0 / TOSCARSSR::kgToGeV(this->GetM());
    double const Beta = Gamma != 1 ? sqrt(1.0 - 1.0 / (Gamma * Gamma)) : 0;

    // Copy this particle and set ideal conditions
    TParticleA NewParticle((TParticleA) *this);
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
  // UPDATE: Could also take a python function

  // Clear current trajectory
  this->GetTrajectory().Clear();

  // If this is a filament beam return the ideal case
  if (this->GetBeamDistribution() == kBeamDistribution_Filament) {
    return this->GetNewParticle("ideal");
  }

  // Constant PI
  double const PI = TOSCARSSR::Pi();

  // The new particle to return
  TParticleA NewParticle((TParticleA) *this);

  // For a gaussian beam
  if (this->GetBeamDistribution() == kBeamDistribution_Gaussian) {

    // New energy
    double ENew = fE0 + fSigmaEnergyGeV * gRandomA->Normal();
    if (ENew < TOSCARSSR::kgToGeV(this->GetM())) {
      std::cerr << "WARNING in TParticleBeam::GetNewParticle(): ENew < mc^2.  Setting to mc^2" << std::endl;
      std::cerr << "  ENew fSigmaEnergyGeV: " << ENew << "  " << fSigmaEnergyGeV << std::endl;
      ENew = TOSCARSSR::kgToGeV(this->GetM());
    }

    // Gamma and Beta magnitude from energy distribution
    double const GammaE = ENew / TOSCARSSR::kgToGeV(this->GetM()) < 1 ? 1 : ENew / TOSCARSSR::kgToGeV(this->GetM());
    double const Beta = GammaE != 1 ? sqrt(1.0 - 1.0 / (GammaE * GammaE)) : 0;

    // Converging or diverging
    int const hConverge = fTwissAlphaX0[0] > 0 ? -1 : 1;
    int const vConverge = fTwissAlphaX0[1] > 0 ? -1 : 1;

    // For sampling
    double const hSigX  = sqrt(fEmittance[0] *  fTwissBetaX0[0]);
    double const hSigXP = sqrt(fEmittance[0] * fTwissGammaX0[0]);
    double const hRho   = sqrt(1 - pow(fEmittance[0] / (hSigX * hSigXP), 2));

    double const vSigX  = sqrt(fEmittance[1] *  fTwissBetaX0[1]);
    double const vSigXP = sqrt(fEmittance[1] * fTwissGammaX0[1]);
    double const vRho   = sqrt(1 - pow(fEmittance[1] / (vSigX * vSigXP), 2));


    // Horizontal x and xp
    double const hu  = gRandomA->Uniform();
    double const hv  = gRandomA->Uniform();
    double const hX  = hSigX * sqrt(-2 * log(hu)) * (sqrt(1 - hRho*hRho) * cos(2 * PI * hv) + hRho * sin(2 * PI * hv));
    double const hXP = hConverge * hSigXP * sqrt(-2 * log(hu)) * sin(2 * PI * hv);

    // Vertical x and xp
    double const vu  = gRandomA->Uniform();
    double const vv  = gRandomA->Uniform();
    double const vX  = vSigX * sqrt(-2 * log(vu)) * (sqrt(1 - vRho*vRho) * cos(2 * PI * vv) + vRho * sin(2 * PI * vv));
    double const vXP = vConverge * vSigXP * sqrt(-2 * log(vu)) * sin(2 * PI * vv);

    // New X0 location for this particle
    TVector3D XNew = this->GetX0();
    XNew += fHorizontalDirection * hX;
    XNew += fVerticalDirection   * vX;

    // Beta vector from randomization and energy spread
    TVector3D BetaNew = fHorizontalDirection.Cross(fVerticalDirection) * Beta;
    BetaNew.RotateSelf(hXP, fVerticalDirection);
    BetaNew.RotateSelf(vXP, -fHorizontalDirection);

    // Possibility to shift the time
    double const TNew = fT0;

    NewParticle.SetInitialParticleConditions(XNew, BetaNew, TNew);

  } else if (this->GetBeamDistribution() == kBeamDistribution_KV) {
    // We have a KV distribution - uniform within the ellipse

    // New energy
    double ENew = fE0 + fSigmaEnergyGeV * gRandomA->Normal();
    if (ENew < TOSCARSSR::kgToGeV(this->GetM())) {
      std::cerr << "WARNING in TParticleBeam::GetNewParticle(): ENew < mc^2.  Setting to mc^2" << std::endl;
      std::cerr << "  ENew fSigmaEnergyGeV: " << ENew << "  " << fSigmaEnergyGeV << std::endl;
      ENew = TOSCARSSR::kgToGeV(this->GetM());
    }

    // Gamma and Beta magnitude from energy distribution
    double const GammaE = ENew / TOSCARSSR::kgToGeV(this->GetM()) < 1 ? 1 : ENew / TOSCARSSR::kgToGeV(this->GetM());
    double const Beta = GammaE != 1 ? sqrt(1.0 - 1.0 / (GammaE * GammaE)) : 0;

    // Converging or diverging
    int const hConverge = fTwissAlphaX0[0] > 0 ? -1 : 1;
    int const vConverge = fTwissAlphaX0[1] > 0 ? -1 : 1;

    double const hTheta = gRandomA->Uniform() * 2 * PI;
    double const vTheta = gRandomA->Uniform() * 2 * PI;
    
    double const hAX = sqrt(fEmittance[0] * fTwissBetaX0[0]) * cos(hTheta);
    double const vAX = sqrt(fEmittance[1] * fTwissBetaX0[1]) * cos(vTheta);

    double const hAXP = -sqrt(fEmittance[0] / fTwissBetaX0[0]) * (fTwissAlphaX0[0] * cos(hTheta) + sin(hTheta));
    double const vAXP = -sqrt(fEmittance[1] / fTwissBetaX0[1]) * (fTwissAlphaX0[1] * cos(vTheta) + sin(vTheta));
    
    double const hR = sqrt(hAX*hAX + hAXP*hAXP);
    double const vR = sqrt(vAX*vAX + vAXP*vAXP);

    double hRealTheta = atan(hAXP/hAX);
    double vRealTheta = atan(vAXP/vAX);

    if (hAX < 0 && hAXP > 0) {
      hRealTheta += PI;
    } else if (hAX < 0 && hAXP < 0) {
      hRealTheta -= PI;
    }
    if (vAX < 0 && vAXP > 0) {
      vRealTheta += PI;
    } else if (vAX < 0 && vAXP < 0) {
      vRealTheta -= PI;
    }

    double const hNewR = sqrt(gRandomA->Uniform()) * hR;
    double const vNewR = sqrt(gRandomA->Uniform()) * vR;

    double const hX = hNewR * cos(hRealTheta);
    double const vX = vNewR * cos(vRealTheta);

    double const hXP = hNewR * sin(hRealTheta) * hConverge;
    double const vXP = vNewR * sin(vRealTheta) * vConverge;

    // New X0 location for this particle
    TVector3D XNew = this->GetX0();
    XNew += fHorizontalDirection * hX;
    XNew += fVerticalDirection   * vX;

    // Beta vector from randomization and energy spread
    TVector3D BetaNew = fHorizontalDirection.Cross(fVerticalDirection) * Beta;
    BetaNew.RotateSelf(hXP, fVerticalDirection);
    BetaNew.RotateSelf(vXP, -fHorizontalDirection);

    // Possibility to shift the time
    double const TNew = fT0;

    NewParticle.SetInitialParticleConditions(XNew, BetaNew, TNew);

  }





  return NewParticle;
}



void TParticleBeam::SetBeamDistribution (TParticleBeam_BeamDistribution const D)
{
  // Set the beam distribution type
  fBeamDistribution = D;
  return;
}



void TParticleBeam::SetBeamDistribution (std::string const& Name)
{
  // Set the beam distribution type
  fBeamDistribution = this->GetBeamDistribution(Name);
  return;
}



TParticleBeam::TParticleBeam_BeamDistribution TParticleBeam::GetBeamDistribution (std::string const& Name) const
{
  // Get the TParticleBeam_BeamDistribution for a named beam
  
  std::string L = Name;
  std::transform(L.begin(), L.end(), L.begin(), ::tolower);

  if (L == "") {
    return kBeamDistribution_None;
  } else if (L == "filament") {
    return kBeamDistribution_Filament;
  } else if (L == "gaussian") {
    return kBeamDistribution_Gaussian;
  } else if (L == "kv") {
    return kBeamDistribution_KV;
  } else {
    throw std::invalid_argument(("Unknown beam distribution name: " + Name).c_str());
  }

  return kBeamDistribution_None;
}







TParticleBeam::TParticleBeam_BeamDistribution const TParticleBeam::GetBeamDistribution () const
{
  return fBeamDistribution;
}




std::string TParticleBeam::GetBeamDistributionName () const
{
  switch (fBeamDistribution) {
    case kBeamDistribution_None:
      return std::string("none");
    case kBeamDistribution_Filament:
      return std::string("filament");
    case kBeamDistribution_Gaussian:
      return std::string("gaussian");
    case kBeamDistribution_KV:
      return std::string("kv");
  }

  return std::string("");
}


