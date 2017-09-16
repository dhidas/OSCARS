////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Jun 29 08:52:22 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TField3D_Gaussian.h"
#include "TOSCARSSR.h"

#include <cmath>

TField3D_Gaussian::TField3D_Gaussian (std::string const Name)
{
  // Constructor

  // Set the name and default scale factors
  this->SetName(Name);
  this->SetScaleFactorMinimumMaximum();
}



TField3D_Gaussian::TField3D_Gaussian (TVector3D const& PeakField,
                                      TVector3D const& Center,
                                      TVector3D const& Sigma,
                                      TVector3D const& Rotations,
                                      std::string const& Name)
{
  // Constructor you should use.. just a suggestion...

  // Set the name and default scale factors
  this->SetName(Name);
  this->SetScaleFactorMinimumMaximum();

  fPeakField = PeakField;
  fPeakField.RotateSelfXYZ(Rotations);

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



TField3D_Gaussian::~TField3D_Gaussian ()
{
  // Anti-Constructor
}



TVector3D TField3D_Gaussian::GetF (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z));
}




TVector3D TField3D_Gaussian::GetF (TVector3D const& X) const
{
  // Get the magnetic field at a point in space.

  // If you rotate the object the field is rotated in fPeakField and the coordinate rotation is done here

  // Translate back into box frame
  TVector3D XInBoxCoordinates = X;
  XInBoxCoordinates.RotateSelfXYZ(fRotated);

  // Position in the box frame with respect to the center
  TVector3D const RX = XInBoxCoordinates - fCenter;

  double Fraction = 1;
  if (fSigma.GetX() > 0) {
    Fraction *= exp(-pow((RX.GetX()) / fSigma.GetX(), 2) / 2.);
  }
  if (fSigma.GetY() > 0) {
    Fraction *= exp(-pow((RX.GetY()) / fSigma.GetY(), 2) / 2.);
  }
  if (fSigma.GetZ() > 0) {
    Fraction *= exp(-pow((RX.GetZ()) / fSigma.GetZ(), 2) / 2.);
  }

  return Fraction * fPeakField;
}





TVector3D const& TField3D_Gaussian::GetPeakField () const
{
  return fPeakField;
}





TVector3D const& TField3D_Gaussian::GetCenter () const
{
  return fCenter;
}





TVector3D const& TField3D_Gaussian::GetSigma () const
{
  return fSigma;
}





TVector3D const& TField3D_Gaussian::GetRotations () const
{
  return fRotated;
}




void TField3D_Gaussian::Print (std::ostream& os) const
{
  os << *this << std::endl;
  return;
}
