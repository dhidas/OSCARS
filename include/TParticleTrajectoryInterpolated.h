#ifndef GUARD_TParticleTrajectoryInterpolated_h
#define GUARD_TParticleTrajectoryInterpolated_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Aug 31 08:40:13 EDT 2017
//
// TParticleTrajectoryInterpolated
//
//   This class is for interpolating a particle trajectory.  All
//   coefficients are stored internally
//
////////////////////////////////////////////////////////////////////

#include <vector>

#include "TOMATH.h"
#include "TParticleTrajectoryPoint.h"
#include "TParticleTrajectoryPoints.h"

class TParticleTrajectoryInterpolated
{
  public:
    TParticleTrajectoryInterpolated ();

    TParticleTrajectoryInterpolated (TParticleTrajectoryPoints const& TPTP);

    TParticleTrajectoryInterpolated (std::vector<double> const& T,
                                     std::vector<TParticleTrajectoryPoint> const& P);

    ~TParticleTrajectoryInterpolated ();

    void Set (TParticleTrajectoryPoints const& TPTP);

    void Set (std::vector<double> const& T,
              std::vector<TParticleTrajectoryPoint> const& P);

    TParticleTrajectoryPoint GetTrajectoryPoint (double const T) const;

    void Clear ();

    void FillTParticleTrajectoryPointsLevel (TParticleTrajectoryPoints& TPTP,
                                             int const Level);

    void FillTParticleTrajectoryPoints (TParticleTrajectoryPoints& TPTP,
                                        double const TStart,
                                        double const TStop,
                                        int const NPoints);

    int    GetNPointsInclusiveToLevel (int const Level) const;
    double GetDeltaTInclusiveToLevel  (int const Level) const;
    double GetDeltaTThisLevel         (int const Level) const;

  private:
    void LevelCheck (int const Level) const;

    TOMATH::TSpline1D3<TParticleTrajectoryPoint> fP;

    // First and last point in time for input trajectory
    double fTStart;
    double fTStop;


};





#endif
