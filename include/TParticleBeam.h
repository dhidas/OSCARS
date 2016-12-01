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
    TParticleBeam (std::string const&);
    TParticleBeam (std::string const&, double const, double const);
    TParticleBeam (std::string const&, TVector3D const&, TVector3D const&, double const, double const);
    TParticleBeam (std::string const&, TVector3D const&, TVector3D const&, double const, double const, double const, double const Charge = 0, double const Mass = 0);
    ~TParticleBeam ();


    void SetInitialConditions (double const, double const, double const, double const, double const, double const, double const, double const);
    void SetInitialConditions (TVector3D const&, TVector3D const&, double const, double const);

    void SetSigma(TVector3D const&, TVector2D const&, TVector2D const&, TVector3D const&, double const);

    TVector3D const& GetX0 () const;
    TVector3D const& GetU0 () const;
    double           GetE0 () const;
    double           GetT0 () const;


    TParticleA GetNewParticle ();
    TParticleA GetNewParticle (std::string const&);

  private:

    TVector3D fX0;  // Coordinates of initial conditions [m]
    TVector3D fU0;  // Direction at fX0 (stored as a unit vector) [unitless]
    double    fE0;  // Energy [GeV]
    double    fT0;  // Time at initial conditions [s]

    TVector2D fSigmaU;
    TVector2D fSigmaUP;
    TVector3D fSigmaAt;
    double    fSigmaE;

    TVector3D fHorizontalDirection;
    TVector3D fVerticalDirection;
};




inline std::ostream& operator << (std::ostream& os, TParticleBeam const& o)
{
  os << "X0: " << o.GetX0() << "\n"
     << "U0: " << o.GetU0() << "\n"
     << "T0: " << o.GetT0() << "\n"
     << "E0: " << o.GetE0() << "\n"
     << "Current " << o.GetCurrent() << "\n";

  return os;
}

















#endif
