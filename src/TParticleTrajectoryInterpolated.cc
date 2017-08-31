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
  return;
}



void TParticleTrajectoryInterpolated::Clear ()
{
  //fT.Clear();
  //fP.Clear();

  return;
}
