////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Aug 31 08:40:13 EDT 2017
//
////////////////////////////////////////////////////////////////////

#include "TParticleTrajectoryInterpolated.h"



TParticleTrajectoryInterpolated::TParticleTrajectoryInterpolated ()
{
  // Default constructor
}




TParticleTrajectoryInterpolated::TParticleTrajectoryInterpolated (std::vector<double> const& T,
                                                                  std::vector<TParticleTrajectoryPoint> const& P)
{
  this->Set(T, P);
}




TParticleTrajectoryInterpolated::~TParticleTrajectoryInterpolated ()
{
  // Destruction
}




void TParticleTrajectoryInterpolated::Set (std::vector<double> const& T,
                                           std::vector<TParticleTrajectoryPoint> const& P)
{
  fP.Set(T, P);
  return;
}



void TParticleTrajectoryInterpolated::Clear ()
{
  fP.Clear();

  return;
}




TParticleTrajectoryPoint TParticleTrajectoryInterpolated::GetTrajectoryPoint (double const T) const
{
  // Get the trajectory values at time T
  return fP.GetValue(T);
}




void TParticleTrajectoryInterpolated::FillTParticleTrajectoryPoints (TParticleTrajectoryPoints& TPTP,
                                                                     double const TStart,
                                                                     double const TStop,
                                                                     int    const NPoints)
{
  // Fill a TParticleTrajectoryPoints object with trajectory points from this interpolation

  if (TStop <= TStart) {
    throw;
  }

  if (NPoints < 2) {
    throw;
  }

  // Time step
  double const DeltaT = (TStop - TStart) / ( ((double) NPoints) - 1 );

  // Set the deltaT of the particle trajectory points
  TPTP.SetDeltaT(DeltaT);

  for (int i = 0; i < NPoints; ++i) {
    double const T = TStart + DeltaT * (double) i;
    TPTP.AddPoint( this->GetTrajectoryPoint(T) );
  }

  return;
}

