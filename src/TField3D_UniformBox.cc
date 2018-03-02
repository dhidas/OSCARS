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
#include "TOSCARSSR.h"

#include <cmath>

TField3D_UniformBox::TField3D_UniformBox (double      const  Fx,
                                          double      const  Fy,
                                          double      const  Fz,
                                          double      const  Frequency,
                                          double      const  FrequencyPhase,
                                          double      const  TimeOffset,
                                          std::string const& Name)
{
  // Set the name and default scale factors
  this->SetName(Name);
  this->SetScaleFactorMinimumMaximum();

  fField = TVector3D(Fx, Fy, Fz);
  fWidth  = TVector3D(0, 0, 0);
  fCenter = TVector3D(0, 0, 0);
  fRotated = TVector3D(0, 0, 0);

  fFrequency = Frequency;
  fFrequencyPhase = FrequencyPhase;
  fTimeOffset = TimeOffset;

  fIgnoreAxisX = true;
  fIgnoreAxisY = true;
  fIgnoreAxisZ = true;
}



TField3D_UniformBox::TField3D_UniformBox (TVector3D   const& Field,
                                          TVector3D   const& Width,
                                          TVector3D   const& Center,
                                          TVector3D   const& Rotations,
                                          double      const  Frequency,
                                          double      const  FrequencyPhase,
                                          double      const  TimeOffset,
                                          std::string const& Name)
{
  // Set the name and default scale factors
  this->SetName(Name);
  this->SetScaleFactorMinimumMaximum();

  fField = Field;
  fField.RotateSelfXYZ(Rotations);

  fWidth  = Width;
  fCenter = Center;
  fRotated = Rotations;

  fFrequency = Frequency;
  fFrequencyPhase = FrequencyPhase;
  fTimeOffset = TimeOffset;

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


TVector3D TField3D_UniformBox::GetF (double const X, double const Y, double const Z, double const T) const
{
  return this->GetF(TVector3D(X, Y, Z), T);
}




TVector3D TField3D_UniformBox::GetF (TVector3D const& X, double const T) const
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

  if (fFrequency == 0) {
    // It has no time dependence
    return fField;
  }
  return fField * cos(TOSCARSSR::TwoPi() * fFrequency * (T + fTimeOffset) + fFrequencyPhase);
}




TVector3D TField3D_UniformBox::GetField () const
{
  // Return the field
  return fField;
}




TVector3D TField3D_UniformBox::GetWidth () const
{
  // Return the width
  return fWidth;
}




TVector3D TField3D_UniformBox::GetRotated () const
{
  // Return the rotations
  return fRotated;
}




TVector3D TField3D_UniformBox::GetCenter () const
{
  // Return the Center postion
  return fCenter;
}




double TField3D_UniformBox::GetFrequency () const
{
  // Return the frequency
  return fFrequency;
}




double TField3D_UniformBox::GetFrequencyPhase () const
{
  // Return the frequency
  return fFrequencyPhase;
}




double TField3D_UniformBox::GetTimeOffset () const
{
  // Return the time offset for the frequency
  return fTimeOffset;
}




void TField3D_UniformBox::Print (std::ostream& os) const
{
  os << *this << std::endl;
  return;
}

