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

  // Set start and stop to zero
  fTStart = 0;
  fTStop = 0;
}




TParticleTrajectoryInterpolated::TParticleTrajectoryInterpolated (TParticleTrajectoryPoints const& TPTP)
{
  // Constructor

  // Number of points in inpt
  size_t const NPoints = TPTP.GetNPoints();
  if (NPoints < 2) {
    throw std::length_error("TParticleTrajectoryInterpolated::TParticleTrajectoryInterpolated not enough points in input TPTP trajectory");
  }

  this->Set(TPTP.GetTimePoints(), TPTP.GetTrajectory());
}




TParticleTrajectoryInterpolated::TParticleTrajectoryInterpolated (std::vector<double> const& T,
                                                                  std::vector<TParticleTrajectoryPoint> const& P)
{
  this->Set(T, P);
}




TParticleTrajectoryInterpolated::~TParticleTrajectoryInterpolated ()
{
  this->Clear();
  // Destruction
}




void TParticleTrajectoryInterpolated::Set (TParticleTrajectoryPoints const& TPTP)
{
  // Check that there are at least two points
  if (TPTP.GetNPoints() < 2) {
    throw std::length_error("TParticleTrajectoryInterpolated::Set NPoints is too small");
  }

  // Set the interpolating structure
  fP.Set(TPTP.GetTimePoints(), TPTP.GetTrajectory());

  // Get the start and stop times (first and last).  These times will be "inclusive"
  fTStart = TPTP.GetTStart();
  fTStop  = TPTP.GetTStop();

  return;
}




void TParticleTrajectoryInterpolated::Set (std::vector<double> const& T,
                                           std::vector<TParticleTrajectoryPoint> const& P)
{
  // Check that there are at least two points and that the number of points is
  // the same as the number of time points
  if (T.size() < 2 || T.size() != P.size()) {
    throw std::length_error("TParticleTrajectoryInterpolated::Set NPoints is too small or T and P do not match");
  }

  // Set the interpolating structure
  fP.Set(T, P);

  // Get the start and stop times (first and last).  These times will be "inclusive"
  fTStart = T.front();
  fTStop = T.back();

  return;
}



void TParticleTrajectoryInterpolated::Clear ()
{
  // Clear interpolating structure
  fP.Clear();

  // Set start and stop back to zero
  fTStart = 0;
  fTStop = 0;

  return;
}




TParticleTrajectoryPoint TParticleTrajectoryInterpolated::GetTrajectoryPoint (double const T) const
{
  // Get the trajectory values at time T
  return fP.GetValue(T);
}




void TParticleTrajectoryInterpolated::FillTParticleTrajectoryPointsLevel (TParticleTrajectoryPoints& TPTP,
                                                                          int const Level) const
{
  // Fill a TParticleTrajectoryPoints object with trajectory points from this interpolation

  if (fTStop <= fTStart) {
    throw std::logic_error("TParticleTrajectoryInterpolated::FillTParticleTrajectoryPointsLevel throwing because fTStop <= fTStart.  Internal logic error.  Please report this.");
  }

  // Level checking
  this->LevelCheck(Level);


  // Number of points in this level
  int const NPoints = pow(2, Level);

  // Spacing of points in THIS level
  double const ThisTSpacing = (fTStop - fTStart) / pow(2., Level);

  // DeltaT including all points up to and including this level
  //double const DeltaTAllLevels = (fTStop - fTStart) / pow(2., Level + 1);



  // Set the deltaT of the particle trajectory points
  TPTP.SetDeltaT(ThisTSpacing);

  // First point of this trajectory is at:
  double const ThisTStart = this->GetTStartThisLevel(Level);

  for (int i = 0; i < NPoints; ++i) {
    double const T = ThisTStart + ThisTSpacing * (double) i;
    TPTP.AddPoint( this->GetTrajectoryPoint(T), T);
  }

  return;
}





void TParticleTrajectoryInterpolated::FillTParticleTrajectoryPoints (TParticleTrajectoryPoints& TPTP,
                                                                     int    const NPoints) const
{
  // Fill a TParticleTrajectoryPoints object with trajectory points from this interpolation

  this->FillTParticleTrajectoryPoints(TPTP, NPoints, fTStart, fTStop);

  return;
}




void TParticleTrajectoryInterpolated::FillTParticleTrajectoryPoints (TParticleTrajectoryPoints& TPTP,
                                                                     int    const NPoints,
                                                                     double const TStart,
                                                                     double const TStop) const
{
  // Fill a TParticleTrajectoryPoints object with trajectory points from this interpolation

  if (TStop <= TStart) {
    throw std::logic_error("TParticleTrajectoryInterpolated::FillTParticleTrajectoryPoints finding TStop <= TStart.  Please report this");
  }

  if (NPoints < 2) {
    throw std::logic_error("TParticleTrajectoryInterpolated::FillTParticleTrajectoryPoints finding NPoints < 2.  Please report this");
  }

  // Time step
  double const DeltaT = (TStop - TStart) / ( ((double) NPoints) - 1 );

  // Set the deltaT of the particle trajectory points
  TPTP.SetDeltaT(DeltaT);

  for (int i = 0; i < NPoints; ++i) {
    double const T = TStart + DeltaT * (double) i;
    TPTP.AddPoint( this->GetTrajectoryPoint(T), T );
  }

  return;
}




int TParticleTrajectoryInterpolated::GetNPointsThisLevel (int const Level)
{
  return pow(2, Level);
}




int TParticleTrajectoryInterpolated::GetNPointsInclusiveToLevel (int const Level)
{
  return pow(2, Level + 1) - 1;
}




double TParticleTrajectoryInterpolated::GetDeltaTInclusiveToLevel (int const Level) const
{
  // DeltaT for all points from level 0 up to this level

  // Level checking
  this->LevelCheck(Level);
    
  return (fTStop - fTStart) / pow(2., Level + 1);

}




double TParticleTrajectoryInterpolated::GetDeltaTThisLevel (int const Level) const
{
  // DeltaT for all points from level 0 up to this level

  // Level checking
  this->LevelCheck(Level);
    
  return (fTStop - fTStart) / pow(2., Level); 

}




size_t TParticleTrajectoryInterpolated::GetNPoints () const
{
  return fP.GetNPoints();
}




double TParticleTrajectoryInterpolated::GetTStartThisLevel (int const Level) const
{
  // TSTart for this level

  return fTStart + (fTStop - fTStart) / pow(2., Level + 1); 

}




double TParticleTrajectoryInterpolated::GetTStart () const
{
  // TStart for interpolation

  return fTStart; 

}




double TParticleTrajectoryInterpolated::GetTStop () const
{
  // TStop for interpolation

  return fTStop; 

}




TOMATH::TSpline1D3<TParticleTrajectoryPoint> const& TParticleTrajectoryInterpolated::GetSpline () const
{
  return fP;
}




void TParticleTrajectoryInterpolated::LevelCheck (int const Level) const
{
  // Common place to check the level and throw exception if incorrect

  if (Level < 0) {
    throw std::logic_error("TParticleTrajectoryInterpolated::LevelCheck seeing Level < 0.  Please report this");
  }
  return;
}
