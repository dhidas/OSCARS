#ifndef GUARD_TParticleA_h
#define GUARD_TParticleA_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Jun  8 14:35:28 EDT 2016
//
// Class describing a particle, or at least its relevent properties
// for calculations in this package
//
////////////////////////////////////////////////////////////////////

#include <string>

#include "TVector3D.h"
#include "TParticleTrajectoryPoints.h"

class TParticleA
{
  public:
    TParticleA ();
    TParticleA (std::string const&);
    TParticleA (std::string const&, TVector3D const&, TVector3D const&, double const);
    ~TParticleA ();

    void SetParticleType (std::string const&);
    void SetParticleTypeCustom (std::string const&, double const, double const);
    void SetParticleTypeFromPDGID (int const);

    void SetQ  (double const);
    void SetM  (double const);
    void SetQM (double const, double const);

    void SetX0 (TVector3D const&);
    void SetB0 (TVector3D const&);
    void SetT0 (double const);

    double GetQ () const;
    double GetM () const;
    double GetGamma () const;
    double GetQoverMGamma () const;

    void   SetCurrent (double const);
    double GetCurrent () const;

    TVector3D const& GetX0 () const;
    TVector3D const& GetB0 () const;
    double           GetE0 () const;
    double           GetT0 () const;

    std::string const& GetType () const;

    void SetInitialParticleConditions (TVector3D const&, TVector3D const&, double const);

    TParticleTrajectoryPoints& GetTrajectory ();



  private:
    void SetGamma ();
    void SetQoverMGamma ();

    std::string fType;    // Type name
    double fQ;            // Charge
    double fM;            // Mass
    double fGamma;        // Relativistic gamma factor (derived quantity)
    double fQoverMGamma;  // Factor computed often in differential equations of motion

    TVector3D fX0;  // Coordinates of initial conditions
    TVector3D fB0;  // Initial Beta (velocity / c)
    double    fT0;  // Time at initial conditions

    TParticleTrajectoryPoints fTrajectory;

    // This is a funny one so I'll explain it here.
    // This is here because TParticleBeam inherits this class
    // but it is more convenient for SRS to have this as a TParticleA
    // class member.
    double fCurrent;

};








inline std::ostream& operator << (std::ostream& os, TParticleA const& o)
{
  os << "Q:            " << o.GetQ()  << "\n"
     << "M:            " << o.GetM()  << "\n"
     << "Gamma:        " << std::scientific << o.GetGamma()  << "\n"
     << "QoverMGamma:  " << o.GetQoverMGamma()  << "\n"
     << "X0:           " << o.GetX0() << "\n"
     << "U0:           " << o.GetB0() << "\n"
     << "T0:           " << o.GetT0() << "\n";

  return os;
}



















#endif
