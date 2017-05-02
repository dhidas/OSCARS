#ifndef GUARD_TParticleBeam_h
#define GUARD_TParticleBeam_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jun  7 18:01:19 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TParticleA.h"

#include <map>
#include <string>
#include <vector>

#include "TVector3D.h"
#include "TVector2D.h"

class TParticleBeam : public TParticleA
{
  public:
    TParticleBeam ();
    TParticleBeam (std::string const& PredefinedBeamType,
                   std::string const& Name,
                   double const Weight = 1);

    TParticleBeam (std::string const& ParticleType,
                   std::string const& Name,
                   double const Energy,
                   double const Current,
                   double const Weight = 1);

    TParticleBeam (std::string const& ParticleType,
                   std::string const& Name,
                   TVector3D const& X0,
                   TVector3D const& D0,
                   double const Energy,
                   double const Current,
                   double const Weight = 1);

    TParticleBeam (std::string const& ParticleType,
                   std::string const& Name,
                   TVector3D const& X0,
                   TVector3D const& D0,
                   double const Energy,
                   double const T0,
                   double const Current,
                   double const Charge = 0,
                   double const Mass = 0,
                   double const Weight = 1);

    ~TParticleBeam ();

    void SetPredefinedBeam (std::string const& Beam);

    void SetInitialConditions (double const X,
                               double const Y,
                               double const Z,
                               double const Dx,
                               double const Dy,
                               double const Dz,
                               double const E0,
                               double const T);

    void SetInitialConditions (TVector3D const& X,
                               TVector3D const& D,
                               double    const E0,
                               double    const T0);

    void SetBetaEmittance (TVector3D const& HorizontalDirection,
                           TVector2D const& Beta,
                           TVector2D const& Emittance,
                           TVector3D const& SigmaAt,
                           double const SigmaEnergyGeV);

    void SetSigma (TVector3D const& HorizontalDirection,
                   TVector2D const& SigmaU,
                   TVector2D const& SigmaUP,
                   TVector3D const& SigmaAt,
                   double const SigmaEnergyGeV);


    TVector3D const& GetX0 () const;
    TVector3D const& GetU0 () const;
    double           GetE0 () const;
    double           GetT0 () const;

    std::string const& GetName () const;
    double             GetWeight () const;

    TVector2D GetBeta () const;
    TVector2D GetEmittance () const;

    TVector3D GetVerticalDirection () const;
    TVector3D GetHorizontalDirection () const;
    TVector3D GetSigmaAt () const;

    void SetX0 (TVector3D const&);
    void SetU0 (TVector3D const&);
    void SetE0 (double const);
    void SetT0 (double const);
    void SetSigmaEnergyGeV (double const);

    void SetName (std::string const&);
    void SetWeight (double const);


    TParticleA GetNewParticle ();
    TParticleA GetNewParticle (std::string const&);

  private:
    std::string fName;
    double      fWeight;

    TVector3D fX0;  // Coordinates of initial conditions [m]
    TVector3D fU0;  // Direction at fX0 (stored as a unit vector) [unitless]
    double    fE0;  // Energy [GeV]
    double    fT0;  // Time at initial conditions [s]

    // Horizontal and vertical beta function and emittance
    TVector2D fBeta;
    TVector2D fEmittance;

    // Source size, calculated from Beta and Emittance
    TVector2D fSigmaU;
    TVector2D fSigmaUP;

    // Energy spread (in units of [GeV])
    double    fSigmaEnergyGeV;

    // Reference point for Beta and Emittance (and SigmaU and SigmaUP)
    TVector3D fSigmaAt;

    TVector3D fHorizontalDirection;
    TVector3D fVerticalDirection;
};




inline std::ostream& operator << (std::ostream& os, TParticleBeam const& o)
{
  os << "Name:       " << o.GetName() << "\n"
     << "Weight:     " << o.GetWeight() << "\n"
     << "X0:         " << o.GetX0() << "\n"
     << "U0:         " << o.GetU0() << "\n"
     << "T0:         " << o.GetT0() << "\n"
     << "E0:         " << o.GetE0() << "\n"
     << "Current     " << o.GetCurrent() << "\n"
     << "Beta        " << o.GetBeta() << "\n"
     << "Emittance   " << o.GetEmittance() << "\n"
     << "V-direction " << o.GetVerticalDirection() << "\n"
     << "H-direction " << o.GetHorizontalDirection() << "\n"
     << "Lattice Ref " << o.GetSigmaAt() << "\n";

  return os;
}

















#endif
