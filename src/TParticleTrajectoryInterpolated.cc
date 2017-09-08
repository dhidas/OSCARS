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

  std::cout << "in constructor" << std::endl;
  // Number of points in inpt
  size_t const NPoints = TPTP.GetNPoints();
  if (NPoints < 2) {
    std::cerr << "throwing npoints too low" << std::endl;
    throw;
  }

  double const DeltaT = TPTP.GetDeltaT();
  fTStart = 0;
  fTStop = fTStart + DeltaT * NPoints;

  std::vector<double> T;
  T.reserve(NPoints);
  for (size_t i = 0; i != NPoints; ++i) {
    T.push_back(DeltaT * (double) i);
    //std::cout << "DeltaT * (double) i: " << DeltaT * (double) i << "  " << TPTP.GetX(i) << std::endl;
  }

  this->Set(T, TPTP.GetTrajectory());
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




void TParticleTrajectoryInterpolated::Set (std::vector<double> const& T,
                                           std::vector<TParticleTrajectoryPoint> const& P)
{
  // Check that there are at least two points and that the number of points is
  // the same as the number of time points
  if (T.size() < 2 || T.size() != P.size()) {
    std::cerr << "T is too small" << std::endl;
    throw;
  }

  std::cout << "trying to Set T, P" << std::endl;
  // Set the interpolating structure
  fP.Set(T, P);
  std::cout << "done Set T, P" << std::endl;

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
                                                                          int const Level)
{
  // Fill a TParticleTrajectoryPoints object with trajectory points from this interpolation

  if (fTStop <= fTStart) {
    throw;
  }

  // Level checking
  this->LevelCheck(Level);


  // Number of points in this level
  int const NPoints = pow(2, Level);

  // Offset for first point
  double const Offset = (fTStop - fTStart) / pow(2., Level + 1);

  // Spacing of points in THIS level
  double const ThisTSpacing = (fTStop - fTStart) / pow(2., Level);

  // DeltaT including all points up to and including this level
  double const DeltaTAllLevels = Offset;



  // Set the deltaT of the particle trajectory points
  TPTP.SetDeltaT(ThisTSpacing);

  // First point of this trajectory is at:
  double const ThisTStart = fTStart + Offset;

  for (int i = 0; i < NPoints; ++i) {
    double const T = ThisTStart + ThisTSpacing * (double) i;
    TPTP.AddPoint( this->GetTrajectoryPoint(T) );
    std::cout << "Added point: T, P: " << T << "  " << this->GetTrajectoryPoint(T) << std::endl;
  }

  return;
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




int TParticleTrajectoryInterpolated::GetNPointsInclusiveToLevel (int const Level) const
{
  // Level checking
  this->LevelCheck(Level);

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




void TParticleTrajectoryInterpolated::LevelCheck (int const Level) const
{
  // Common place to check the level and throw exception if incorrect

  if (Level < 0) {
    throw;
  }
  return;
}
