////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed May  4 08:43:33 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TVector3DC.h"

#include <cmath>



TVector3DC::TVector3DC ()
{
  // Default constructor
}




TVector3DC::TVector3DC (TVector3D const& V)
{
  // Default constructor
  fX = V.GetX();
  fY = V.GetY();
  fZ = V.GetZ();
}




TVector3DC::TVector3DC (std::complex<double> const& X, std::complex<double> const& Y, std::complex<double> const& Z)
{
  // Probably most used and useful constructor
  fX = X;
  fY = Y;
  fZ = Z;
}




TVector3DC::~TVector3DC ()
{
  // Destroy me
}




std::complex<double> TVector3DC::GetX () const
{
  // Return the X-component
  return fX;
}




std::complex<double> TVector3DC::GetY () const
{
  // Return the Y-component
  return fY;
}




std::complex<double> TVector3DC::GetZ () const
{
  // Return the Z-component
  return fZ;
}




void TVector3DC::SetX (std::complex<double> const& X)
{
  // Set the X component
  fX = X;
  return;
}




void TVector3DC::SetY (std::complex<double> const& Y)
{
  // Set the Y component
  fY = Y;
  return;
}




void TVector3DC::SetZ (std::complex<double> const& Z)
{
  // Set the Z component
  fZ = Z;
  return;
}




void TVector3DC::SetXYZ (std::complex<double> const& X, std::complex<double> const& Y, std::complex<double> const& Z)
{
  // Set the X, Y, and Z components
  fX = X;
  fY = Y;
  fZ = Z;

  return;
}




std::complex<double> TVector3DC::Dot (TVector3DC const& V) const
{
  // Get the dot product of this dot V
  return fX * V.GetX() + fY * V.GetY() + fZ * V.GetZ();
}





TVector3DC TVector3DC::Cross (TVector3DC const& V) const
{
  // Get the cross product of this cross V using the right hand convention
  return TVector3DC(fY * V.GetZ() - V.GetY() * fZ, fZ * V.GetX() - V.GetZ() * fX, fX * V.GetY() - V.GetX() * fY);
}




TVector3DC TVector3DC::UnitVector () const
{
  // Get the cross product of this cross V using the right hand convention
  return *this / this->Mag();
}




TVector3DC TVector3DC::CC () const
{
  // Get the cross product of this cross V using the right hand convention
  return TVector3DC(std::conj(fX), std::conj(fY), std::conj(fZ));
}




double TVector3DC::Mag2 () const
{
  // Magnitude squared of 3DC vector
  return this->Dot(this->CC()).real();
}




double TVector3DC::Mag () const
{
  // Magnitude squared of 3DC vector
  return sqrt(this->Mag2());
}



std::complex<double> TVector3DC::MagC2 () const
{
  // Magnitude squared of 3DC vector
  return this->Dot(this->CC());
}




std::complex<double> TVector3DC::MagC () const
{
  // Magnitude squared of 3DC vector
  return sqrt(this->Mag2());
}

















TVector3DC TVector3DC::operator + (TVector3DC const& V) const
{
  // Vector addition, add components and return a vector
  return TVector3DC(fX + V.GetX(), fY + V.GetY(), fZ + V.GetZ());
}




TVector3DC TVector3DC::operator - (TVector3DC const& V) const
{
  // Vector subtraction, subtract components and return a vector
  return TVector3DC(fX - V.GetX(), fY - V.GetY(), fZ - V.GetZ());
}




TVector3DC TVector3DC::operator / (double const& V) const
{
  // Divide vector by some scalar
  return TVector3DC(fX / V, fY / V, fZ / V);
}




TVector3DC& TVector3DC::operator = (TVector3DC const& V)
{
  // Assignment operator
  fX = V.GetX();
  fY = V.GetY();
  fZ = V.GetZ();
  return *this;
}




TVector3DC TVector3DC::operator - ()
{
  // Negative vector
  return TVector3DC(-fX, -fY, -fZ);
}




TVector3DC& TVector3DC::operator += (TVector3DC const& V)
{
  // Add a vector to this vector by components
  fX += V.GetX();
  fY += V.GetY();
  fZ += V.GetZ();
  return *this;
}




TVector3DC& TVector3DC::operator -= (TVector3DC const& V)
{
  // Subtract a vector from this vector by components
  fX -= V.GetX();
  fY -= V.GetY();
  fZ -= V.GetZ();
  return *this;
}




TVector3DC& TVector3DC::operator *= (double const& V)
{
  // Multiply this vector by a scalar
  fX *= V;
  fY *= V;
  fZ *= V;
  return *this;
}




TVector3DC& TVector3DC::operator /= (double const& V)
{
  // Divide this vector by a scalar
  fX /= V;
  fY /= V;
  fZ /= V;
  return *this;
}




TVector3DC& TVector3DC::operator *= (std::complex<double> const& V)
{
  // Multiply this vector by a scalar
  fX *= V;
  fY *= V;
  fZ *= V;
  return *this;
}




TVector3DC& TVector3DC::operator /= (std::complex<double> const& V)
{
  // Divide this vector by a scalar
  fX /= V;
  fY /= V;
  fZ /= V;
  return *this;
}




bool TVector3DC::operator == (TVector3DC const& V) const
{
  // Is this vector equal to V by components
  return fX == V.GetX() && fY == V.GetY() && fZ == V.GetZ();
}




bool TVector3DC::operator != (TVector3DC const& V) const
{
  // Is any component of this vector not equal to the equivalent component of V
  return fX != V.GetX() || fY != V.GetY() || fZ != V.GetZ();
}









