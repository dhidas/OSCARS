////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Mar 25 10:24:13 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"

#include <cmath>



TVector3D::TVector3D ()
{
  // Default constructor
}




TVector3D::TVector3D (double const X, double const Y, double const Z)
{
  // Probably most used and useful constructor
  fX = X;
  fY = Y;
  fZ = Z;
}




TVector3D::~TVector3D ()
{
  // Destroy me
}




void TVector3D::SetX (double const X)
{
  // Set the X component
  fX = X;
  return;
}




void TVector3D::SetY (double const Y)
{
  // Set the Y component
  fY = Y;
  return;
}




void TVector3D::SetZ (double const Z)
{
  // Set the Z component
  fZ = Z;
  return;
}




void TVector3D::SetXYZ (double const X, double const Y, double const Z)
{
  // Set the X, Y, and Z components
  fX = X;
  fY = Y;
  fZ = Z;

  return;
}




void TVector3D::SetXYZ (TVector3D const& V)
{
  // Set the X, Y, and Z components
  fX = V.GetX();
  fY = V.GetY();
  fZ = V.GetZ();

  return;
}




double TVector3D::Mag() const
{
  // Get the magnitude
  return sqrt(Mag2());
}




double TVector3D::Mag2() const
{
  // Get the magnitude squared
  return fX * fX + fY * fY + fZ * fZ;
}




double TVector3D::Dot(TVector3D const& V) const
{
  // Get the dot product of this dot V
  return fX * V.GetX() + fY * V.GetY() + fZ * V.GetZ();
}




double TVector3D::Perp2(TVector3D const& p)  const {
  double const tot = p.Mag2();
  double const ss  = Dot(p);
  double per = Mag2();
  if (tot > 0.0) per -= ss*ss/tot;
  if (per < 0)   per = 0;
  return per;
}




TVector3D TVector3D::Cross (TVector3D const& V) const
{
  // Get the cross product of this cross V using the right hand convention
  return TVector3D(fY * V.GetZ() - V.GetY() * fZ, fZ * V.GetX() - V.GetZ() * fX, fX * V.GetY() - V.GetX() * fY);
}




TVector3D TVector3D::UnitVector () const
{
  // Get a unit vector in the direction of this
  return TVector3D(*this / Mag());
}



void TVector3D::RotateSelfX(double const Angle) {
  // Rotate vector around X
  double const s = sin(Angle);
  double const c = cos(Angle);
  double const yy = fY;

  fY = c * yy - s * fZ;
  fZ = s * yy + c * fZ;

  return;
}




void TVector3D::RotateSelfY (double const Angle) {
  // Rotate vector around Y
  double const s = sin(Angle);
  double const c = cos(Angle);
  double const zz = fZ;

  fZ = c * zz - s * fX;
  fX = s * zz + c * fX;

  return;
}




void TVector3D::RotateSelfZ (double const Angle) {
  // Rotate vector around Z
  double const s = sin(Angle);
  double const c = cos(Angle);
  double const xx = fX;

  fX = c * xx - s * fY;
  fY = s * xx + c * fY;

  return;
}




void TVector3D::RotateSelfXYZ (TVector3D const& V)
{
  // Rotate a vector about X, then Y then Z.
  this->RotateSelfX(V.GetX());
  this->RotateSelfY(V.GetY());
  this->RotateSelfZ(V.GetZ());

  return;
}




void TVector3D::RotateSelf (double const A, TVector3D const& VIN) {
  // Rotate vector by an angle A around the vector V (which is VIN, but normalized)

  // Normalized vector
  TVector3D const V = VIN.UnitVector();

  // Rotation Matrix
  double M[3][3];

  M[0][0] = cos(A) + V[0] * V[0] * (1 - cos(A));
  M[0][1] = V[0] * V[1] * (1 - cos(A)) - V[2] * sin(A);
  M[0][2] = V[0] * V[2] * (1 - cos(A)) + V[1] * sin(A);

  M[1][0] = V[1] * V[0] * (1 - cos(A)) + V[2] * sin(A);
  M[1][1] = cos(A) + V[1] * V[1] * (1 - cos(A));
  M[1][2] = V[1] * V[2] * (1 - cos(A)) - V[0] * sin(A);

  M[2][0] = V[2] * V[0] * (1 - cos(A)) - V[1] * sin(A);
  M[2][1] = V[2] * V[1] * (1 - cos(A)) + V[0] * sin(A);
  M[2][2] = cos(A) + V[2] * V[2] * (1 - cos(A));

  // Get original vector and rotate it
  TVector3D O(fX, fY, fZ);

  this->SetXYZ(M[0][0] * O[0] + M[0][1] * O[1] + M[0][2] * O[2],
               M[1][0] * O[0] + M[1][1] * O[1] + M[1][2] * O[2],
               M[2][0] * O[0] + M[2][1] * O[1] + M[2][2] * O[2]);

  return;
}






TVector3D TVector3D::operator + (TVector3D const& V) const
{
  // Vector addition, add components and return a vector
  return TVector3D(fX + V.GetX(), fY + V.GetY(), fZ + V.GetZ());
}




TVector3D TVector3D::operator - (TVector3D const& V) const
{
  // Vector subtraction, subtract components and return a vector
  return TVector3D(fX - V.GetX(), fY - V.GetY(), fZ - V.GetZ());
}




TVector3D TVector3D::operator / (double const V) const
{
  // Divide vector by some scalar
  return TVector3D(fX / V, fY / V, fZ / V);
}




TVector3D& TVector3D::operator = (TVector3D const& V)
{
  // Assignment operator
  fX = V.GetX();
  fY = V.GetY();
  fZ = V.GetZ();
  return *this;
}




TVector3D TVector3D::operator - ()
{
  // Negative vector
  return TVector3D(-fX, -fY, -fZ);
}




TVector3D& TVector3D::operator += (TVector3D const& V)
{
  // Add a vector to this vector by components
  fX += V.GetX();
  fY += V.GetY();
  fZ += V.GetZ();
  return *this;
}




TVector3D& TVector3D::operator -= (TVector3D const& V)
{
  // Subtract a vector from this vector by components
  fX -= V.GetX();
  fY -= V.GetY();
  fZ -= V.GetZ();
  return *this;
}




TVector3D& TVector3D::operator *= (double const V)
{
  // Multiply this vector by a scalar
  fX *= V;
  fY *= V;
  fZ *= V;
  return *this;
}




TVector3D& TVector3D::operator /= (double const V)
{
  // Divide this vector by a scalar
  fX /= V;
  fY /= V;
  fZ /= V;
  return *this;
}




bool TVector3D::operator == (TVector3D const& V) const
{
  // Is this vector equal to V by components
  return fX == V.GetX() && fY == V.GetY() && fZ == V.GetZ();
}




bool TVector3D::operator != (TVector3D const& V) const
{
  // Is any component of this vector not equal to the equivalent component of V
  return fX != V.GetX() || fY != V.GetY() || fZ != V.GetZ();
}




double TVector3D::operator [] (int const i) const
{
  // An operator to use an index like a vector
  // For getting a value

  switch (i) {
    case 0:
      return fX;
    case 1:
      return fY;
    case 2:
      return fZ;
    default:
      std::cerr << "ERROR: TVector3D operator []" << std::endl;
      throw;
  }
  return 0.;
}





double& TVector3D::operator [] (int const i)
{
  // An operator to use an index like a vector
  // For setting a value

  switch (i) {
    case 0:
      return fX;
    case 1:
      return fY;
    case 2:
      return fZ;
    default:
      std::cerr << "ERROR: TVector3D operator []" << std::endl;
      throw;
  }
  return fX;
}


