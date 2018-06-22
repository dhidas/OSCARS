////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Sat Aug 26 12:29:14 PDT 2017
//
////////////////////////////////////////////////////////////////////

#include "TParticleTrajectoryPoint.h"



TParticleTrajectoryPoint::TParticleTrajectoryPoint ()
{
  // Default constructor
}




TParticleTrajectoryPoint::TParticleTrajectoryPoint (double const V)
{
  fX = TVector3D(V);
  fB = TVector3D(V);
  fAoverC = TVector3D(V);
}




TParticleTrajectoryPoint::TParticleTrajectoryPoint (TVector3D const& X,
            TVector3D const& B,
            TVector3D const& AoverC)
{
  fX = X;
  fB = B;
  fAoverC = AoverC;
}




TParticleTrajectoryPoint::~TParticleTrajectoryPoint ()
{
  // Default destructor
}



TVector3D const& TParticleTrajectoryPoint::GetX () const
{
  return fX;
}




TVector3D const& TParticleTrajectoryPoint::GetB () const
{
  return fB;
}




void TParticleTrajectoryPoint::SetB (TVector3D const& B)
{
  fB = B;
  return;
}




TVector3D const& TParticleTrajectoryPoint::GetAoverC () const
{
  return fAoverC;
}




void TParticleTrajectoryPoint::SetAoverC (TVector3D const& AoverC)
{
  fAoverC = AoverC;
  return;
}




TParticleTrajectoryPoint TParticleTrajectoryPoint::operator + (TParticleTrajectoryPoint const& T) const
{
  // Vector addition, add components and return a vector
  return TParticleTrajectoryPoint(fX + T.GetX(), fB + T.GetB(), fAoverC + T.GetAoverC());
}




TParticleTrajectoryPoint TParticleTrajectoryPoint::operator - (TParticleTrajectoryPoint const& T) const
{
  // Vector addition, add components and return a vector
  return TParticleTrajectoryPoint(fX - T.GetX(), fB - T.GetB(), fAoverC - T.GetAoverC());
}




TParticleTrajectoryPoint TParticleTrajectoryPoint::operator - () const
{
  // Vector addition, add components and return a vector
  return TParticleTrajectoryPoint(-fX, -fB, -fAoverC);
}




TParticleTrajectoryPoint TParticleTrajectoryPoint::operator / (double const V) const
{
  // Divide vector by some scalar
  return TParticleTrajectoryPoint(fX / V, fB / V, fAoverC / V);
}




TParticleTrajectoryPoint TParticleTrajectoryPoint::operator * (double const V) const
{
  // Divide vector by some scalar
  return TParticleTrajectoryPoint(fX * V, fB * V, fAoverC * V);
}



TParticleTrajectoryPoint TParticleTrajectoryPoint::operator * (TParticleTrajectoryPoint const& V) const
{
  // Vector addition, add components and return a vector
  return TParticleTrajectoryPoint(fX * V.GetX(), fB * V.GetB(), fAoverC * V.GetAoverC());
}


TParticleTrajectoryPoint TParticleTrajectoryPoint::operator - (double const& V) const
{
  // Vector addition, add components and return a vector
  return TParticleTrajectoryPoint(fX - V, fB - V, fAoverC - V);
}



TParticleTrajectoryPoint TParticleTrajectoryPoint::operator / (TParticleTrajectoryPoint const& V) const
{
  // Vector addition, add components and return a vector
  return TParticleTrajectoryPoint(fX / V.GetX(), fB / V.GetB(), fAoverC / V.GetAoverC());
}




