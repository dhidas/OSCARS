////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Jun 29 08:52:22 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBField3D_Gaussian.h"
#include "TSRS.h"

#include <cmath>

TBField3D_Gaussian::TBField3D_Gaussian ()
{
  // Constructor
}



TBField3D_Gaussian::TBField3D_Gaussian (TVector3D const& BField, TVector3D const& Center, TVector3D const& Sigma, TVector3D const& Rotations)
{
  // Constructor you should use.. just a suggestion...

  fBField = BField;
  fBField.RotateSelfXYZ(Rotations);

  fCenter = Center;
  fSigma = Sigma;
  fRotated = Rotations;

  fIgnoreAxisX = false;
  fIgnoreAxisY = false;
  fIgnoreAxisZ = false;

  if (fSigma.GetX() <= 0) {
    fIgnoreAxisX = true;
  }
  if (fSigma.GetY() <= 0) {
    fIgnoreAxisY = true;
  }
  if (fSigma.GetZ() <= 0) {
    fIgnoreAxisZ = true;
  }
}



double TBField3D_Gaussian::GetBx (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z)).GetX();
}




double TBField3D_Gaussian::GetBy (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z)).GetY();
}




double TBField3D_Gaussian::GetBz (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z)).GetZ();
}




TVector3D TBField3D_Gaussian::GetB (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z));
}




TVector3D TBField3D_Gaussian::GetB (TVector3D const& X) const
{
  // Get the magnetic field at a point in space.

  // If you rotate the object the field is rotated in fBField and the coordinate rotation is done here

  // Translate back into box frame
  TVector3D XInBoxCoordinates = X;
  XInBoxCoordinates.RotateSelfXYZ(fRotated);

  // Position in the box frame with respect to the center
  TVector3D const RX = XInBoxCoordinates - fCenter;

  double Fraction = 1;
  if (fSigma.GetX() > 0) {
    Fraction *= exp(-pow((RX.GetX() - fCenter.GetX()) / fSigma.GetX(), 2) / 2.);
  }
  if (fSigma.GetY() > 0) {
    Fraction *= exp(-pow((RX.GetY() - fCenter.GetY()) / fSigma.GetY(), 2) / 2.);
  }
  if (fSigma.GetZ() > 0) {
    Fraction *= exp(-pow((RX.GetZ() - fCenter.GetZ()) / fSigma.GetZ(), 2) / 2.);
  }

  return Fraction * fBField;
}
