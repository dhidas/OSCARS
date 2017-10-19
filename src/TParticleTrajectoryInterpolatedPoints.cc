////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Sep 13 09:51:20 EDT 2017
//
////////////////////////////////////////////////////////////////////

#include "TParticleTrajectoryInterpolatedPoints.h"



TParticleTrajectoryInterpolatedPoints::TParticleTrajectoryInterpolatedPoints ()
{
}




TParticleTrajectoryInterpolatedPoints::TParticleTrajectoryInterpolatedPoints (TParticleTrajectoryInterpolated* TPTI,
                                                                             int const Level)
{
  this->Set(TPTI, Level);
}




TParticleTrajectoryInterpolatedPoints::~TParticleTrajectoryInterpolatedPoints ()
{
}




void TParticleTrajectoryInterpolatedPoints::Set (TParticleTrajectoryInterpolated* TPTI, int const Level)
{
  fTPTI = TPTI;

  fDeltaT = fTPTI->GetDeltaTThisLevel(Level);
  fNPoints = fTPTI->GetNPointsThisLevel(Level);
  fTStart = fTPTI->GetTStartThisLevel(Level);

}



TParticleTrajectoryPoint TParticleTrajectoryInterpolatedPoints::GetTrajectoryPoint (int const i) const
{
  double const Time = fTStart + fDeltaT * (double) i;

  return fTPTI->GetTrajectoryPoint(Time);
}



double TParticleTrajectoryInterpolatedPoints::GetT (int const i) const
{
  return fTStart + fDeltaT * (double) i;
}



double TParticleTrajectoryInterpolatedPoints::GetDeltaT () const
{
  return fDeltaT;
}



int TParticleTrajectoryInterpolatedPoints::GetNPoints () const
{
  return fNPoints;
}
