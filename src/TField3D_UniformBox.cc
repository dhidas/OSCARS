////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jun 30 08:09:53 EDT 2016
//
// UPDATE: Comments
//
////////////////////////////////////////////////////////////////////

#include "TField3D_UniformBox.h"

#include <cmath>

TField3D_UniformBox::TField3D_UniformBox (double const Fx, double const Fy, double const Fz)
{
  fField = TVector3D(Fx, Fy, Fz);
  fWidth  = TVector3D(0, 0, 0);
  fCenter = TVector3D(0, 0, 0);
  fRotated = TVector3D(0, 0, 0);

  fIgnoreAxisX = true;
  fIgnoreAxisY = true;
  fIgnoreAxisZ = true;
}




TField3D_UniformBox::TField3D_UniformBox (TVector3D const& Field, TVector3D const& Width, TVector3D const& Center, TVector3D const& Rotations)
{
  fField = Field;
  fField.RotateSelfXYZ(Rotations);

  fWidth  = Width;
  fCenter = Center;
  fRotated = Rotations;

  fIgnoreAxisX = false;
  fIgnoreAxisY = false;
  fIgnoreAxisZ = false;

  if (fWidth.GetX() <= 0) {
    fIgnoreAxisX = true;
  }
  if (fWidth.GetY() <= 0) {
    fIgnoreAxisY = true;
  }
  if (fWidth.GetZ() <= 0) {
    fIgnoreAxisZ = true;
  }
}


TField3D_UniformBox::~TField3D_UniformBox ()
{
  // Destruction!
}


double TField3D_UniformBox::GetFx (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z)).GetX();
}




double TField3D_UniformBox::GetFy (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z)).GetY();
}




double TField3D_UniformBox::GetFz (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z)).GetZ();
}




TVector3D TField3D_UniformBox::GetF (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z));
}




TVector3D TField3D_UniformBox::GetF (TVector3D const& X) const
{
  // Get the magnetic field at a point in space.

  // If you rotate the object the field is rotated in fField and the coordinate rotation is done here

  // Translate back into box frame
  TVector3D XInBoxCoordinates = X;
  XInBoxCoordinates.RotateSelfXYZ(fRotated);

  // Position in the box frame with respect to the center
  TVector3D const RX = XInBoxCoordinates - fCenter;

  if ((!fIgnoreAxisX && fabs(RX.GetX()) > fabs(fWidth.GetX() / 2.)) || (!fIgnoreAxisY && fabs(RX.GetY()) > fabs(fWidth.GetY() / 2.)) || (!fIgnoreAxisZ && fabs(RX.GetZ()) > fabs(fWidth.GetZ() / 2.))) {
    return TVector3D(0, 0, 0);
  }

  return fField;
}






