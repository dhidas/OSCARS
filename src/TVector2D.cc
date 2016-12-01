////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jul 28 11:45:52 EDT 2016
//
// A 2D Vector class.
//
////////////////////////////////////////////////////////////////////



#include "TVector2D.h"

#include <cmath>



TVector2D::TVector2D ()
{
  // Default constructor
}




TVector2D::TVector2D (double const X, double const Y)
{
  // Probably most used and useful constructor
  fX = X;
  fY = Y;
}




TVector2D::~TVector2D ()
{
  // Destroy me
}




void TVector2D::SetX (double const X)
{
  // Set the X component
  fX = X;
  return;
}




void TVector2D::SetY (double const Y)
{
  // Set the Y component
  fY = Y;
  return;
}




void TVector2D::SetXY (double const X, double const Y)
{
  // Set the X, Y components
  fX = X;
  fY = Y;

  return;
}




double TVector2D::Mag() const
{
  // Get the magnitude
  return sqrt(Mag2());
}




double TVector2D::Mag2() const
{
  // Get the magnitude squared
  return fX * fX + fY * fY;
}




double TVector2D::Dot(TVector2D const& V) const
{
  // Get the dot product of this dot V
  return fX * V.GetX() + fY * V.GetY();
}




double TVector2D::Perp2(TVector2D const& p)  const {
  double const tot = p.Mag2();
  double const ss  = Dot(p);
  double per = Mag2();
  if (tot > 0.0) per -= ss*ss/tot;
  if (per < 0)   per = 0;
  return per;
}





TVector2D TVector2D::UnitVector () const
{
  // Get a unit vector in the direction of this
  return TVector2D(*this / Mag());
}



void TVector2D::RotateSelf (double const Angle) {
  // Rotate vector
  double const s = sin(Angle);
  double const c = cos(Angle);
  double const xx = fX;

  fX = fX * c - fY * s,
  fY = xx * s + fY * c;

  return;
}






TVector2D TVector2D::operator + (TVector2D const& V) const
{
  // Vector addition, add components and return a vector
  return TVector2D(fX + V.GetX(), fY + V.GetY());
}




TVector2D TVector2D::operator - (TVector2D const& V) const
{
  // Vector subtraction, subtract components and return a vector
  return TVector2D(fX - V.GetX(), fY - V.GetY());
}




TVector2D TVector2D::operator / (double const V) const
{
  // Divide vector by some scalar
  return TVector2D(fX / V, fY / V);
}




TVector2D& TVector2D::operator = (TVector2D const& V)
{
  // Assignment operator
  fX = V.GetX();
  fY = V.GetY();
  return *this;
}




TVector2D TVector2D::operator - ()
{
  // Negative vector
  return TVector2D(-fX, -fY);
}




TVector2D& TVector2D::operator += (TVector2D const& V)
{
  // Add a vector to this vector by components
  fX += V.GetX();
  fY += V.GetY();
  return *this;
}




TVector2D& TVector2D::operator -= (TVector2D const& V)
{
  // Subtract a vector from this vector by components
  fX -= V.GetX();
  fY -= V.GetY();
  return *this;
}




TVector2D& TVector2D::operator *= (double const V)
{
  // Multiply this vector by a scalar
  fX *= V;
  fY *= V;
  return *this;
}




TVector2D& TVector2D::operator /= (double const V)
{
  // Divide this vector by a scalar
  fX /= V;
  fY /= V;
  return *this;
}




bool TVector2D::operator == (TVector2D const& V) const
{
  // Is this vector equal to V by components
  return fX == V.GetX() && fY == V.GetY();
}




bool TVector2D::operator != (TVector2D const& V) const
{
  // Is any component of this vector not equal to the equivalent component of V
  return fX != V.GetX() || fY != V.GetY();
}




double TVector2D::operator [] (int const i) const
{
  // An operator to use an index like a vector
  // For getting a value

  switch (i) {
    case 0:
      return fX;
    case 1:
      return fY;
    default:
      std::cerr << "ERROR: TVector2D operator []" << std::endl;
      throw;
  }
  return 0.;
}





double& TVector2D::operator [] (int const i)
{
  // An operator to use an index like a vector
  // For setting a value

  switch (i) {
    case 0:
      return fX;
    case 1:
      return fY;
    default:
      std::cerr << "ERROR: TVector2D operator []" << std::endl;
      throw;
  }
  return fX;
}


