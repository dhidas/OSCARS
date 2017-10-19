#ifndef GUARD_TParticleTrajectoryInterpolatedPoints_h
#define GUARD_TParticleTrajectoryInterpolatedPoints_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Sep 13 09:51:20 EDT 2017
//
// For getting "points" for a certain "level" from the interpolating
// structure.  This mimics TParticleTrajectoryPoints, the difference
// being that it creates points as needed.  This is slower and
// should only be used in cases where memory is an issue and
// where cpu is not an issue.  The returns are NOT of reference
// type, saving memory, but costing cpu.
//
////////////////////////////////////////////////////////////////////

#include "TParticleTrajectoryInterpolated.h"

class TParticleTrajectoryInterpolatedPoints
{
  public:
    TParticleTrajectoryInterpolatedPoints ();

    TParticleTrajectoryInterpolatedPoints (TParticleTrajectoryInterpolated* TPTI,
                                           int const Level);

    ~TParticleTrajectoryInterpolatedPoints ();

    void Set (TParticleTrajectoryInterpolated* TPTI,
              int const Level);

    TParticleTrajectoryPoint GetTrajectoryPoint (int const i) const;

    double GetT (int const i) const;
    double GetDeltaT () const;
    int    GetNPoints () const;

  private:
    double fDeltaT;
    int    fNPoints;
    double fTStart;
    TParticleTrajectoryInterpolated* fTPTI;
};










#endif
