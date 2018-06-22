#ifndef GUARD_TParticleTrajectoryPoint_h
#define GUARD_TParticleTrajectoryPoint_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Sat Aug 26 12:29:14 PDT 2017
//
// Class for one point in a particle trajectory
//
////////////////////////////////////////////////////////////////////


#include "TVector3D.h"

class TParticleTrajectoryPoint
{
  public:
    TParticleTrajectoryPoint ();

    TParticleTrajectoryPoint (double const);

    TParticleTrajectoryPoint (TVector3D const& X,
          TVector3D const& B,
          TVector3D const& AoverC);

    ~TParticleTrajectoryPoint ();

    TVector3D const& GetX () const;
    TVector3D const& GetB () const;
    void             SetB (TVector3D const&);
    TVector3D const& GetAoverC () const;
    void             SetAoverC (TVector3D const&);

    TParticleTrajectoryPoint operator  + (TParticleTrajectoryPoint const&) const;
    TParticleTrajectoryPoint operator  - (TParticleTrajectoryPoint const&) const;
    TParticleTrajectoryPoint operator  - () const;
    TParticleTrajectoryPoint operator  / (double const) const;
    TParticleTrajectoryPoint operator  * (double const) const;

    // This one is special for spline calculation.
    TParticleTrajectoryPoint  operator  * (TParticleTrajectoryPoint const&) const;
    TParticleTrajectoryPoint  operator  / (TParticleTrajectoryPoint const&) const;
    TParticleTrajectoryPoint  operator  - (double const&) const;


  private:
    TVector3D fX;
    TVector3D fB;
    TVector3D fAoverC;
};


inline std::ostream& operator << (std::ostream& os, TParticleTrajectoryPoint const& o)
{
  // For easy printing
  os << "(" << o.GetX() << ", " << o.GetB() << ", " << o.GetAoverC() << ")";
  return os;
}




inline TParticleTrajectoryPoint operator * (double const V, TParticleTrajectoryPoint const& R)
{
  // Multiply vector by some scalar
  return TParticleTrajectoryPoint(R.GetX() * V, R.GetB() * V, R.GetAoverC() * V);
}




inline TParticleTrajectoryPoint operator / (double const V, TParticleTrajectoryPoint const& R)
{
  // Multiply vector by some scalar
  return TParticleTrajectoryPoint(R.GetX() / V, R.GetB() / V, R.GetAoverC() / V);
}









#endif
