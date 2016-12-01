////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jul  7 17:56:06 EDT 2016
//
// Surface of a rectangle
//
////////////////////////////////////////////////////////////////////

#include "TSurfacePoints_Rectangle.h"

#include <algorithm>

// UPDATE: Comments

TSurfacePoints_Rectangle::TSurfacePoints_Rectangle ()
{
  // Default Constructor
  // Do nothing

}




TSurfacePoints_Rectangle::TSurfacePoints_Rectangle (int const NX1, int const NX2, TVector3D const& X0, TVector3D const& X1, TVector3D const& X2, int const Normal)
{
  // Constructor taking three points to define a rectangle, or parallelagram as the case may be.

  // NX1    - number of points in X1
  // NX2    - number of points in X2
  // X0     - Starting point
  // X1     - Point that defines X1 axis (from X0)
  // X2     - Point that defines X2 axis (from X0)
  // Normal - If -1 reverse the direction of the calculated normal vector

  this->Init(NX1, NX2, X0, X1, X2, Normal);
}






TSurfacePoints_Rectangle::TSurfacePoints_Rectangle (std::string const& Plane, int const NX1, int const NX2, double const WidthX1, double const WidthX2, TVector3D const& Rotations, TVector3D const& Translation, int const Normal)
{
  // Constructor taking a string representing where to orient the plane to start, eg in the "XY" plane.
  // Widths are given such that the rectangle with be centered at 0, 0, 0 and have +/- width/2 in either direction.
  // Rotations are done before translations.  If the Normal int is -1 the direction of the normal is reversed.
  // Direction for the normal is taken from a right handed system, eg "XY" would give +Z direction, "YX" would
  // have a normal in the -Z direction

  // Plane       - The plane you want to start in (eg XY, YX, XZ, ZX, YZ, ZX)
  // NX1         - number of points in X1
  // NX2         - number of points in X2
  // WidthX1     - Width in the X1 direction
  // WidthX2     - Width in the X2 direction
  // Rotations   - Rotation angles about the X, Y, and Z axis.  Rotations are also performed in that order
  // Translation - Take this object and move it anywhere in sace you like, done after rotations
  // Normal      - If -1 reverse the direction of the calculated normal vector

  this->Init(Plane, NX1, NX2, WidthX1, WidthX2, Rotations, Translation, Normal);
}



TSurfacePoints_Rectangle::~TSurfacePoints_Rectangle ()
{
  // Destructor
}




void TSurfacePoints_Rectangle::Init (int const NX1, int const NX2, TVector3D const& X0, TVector3D const& X1, TVector3D const& X2, int const Normal)
{
  // Initializer taking three points in space to define a rectangle, the number of points along each axis,
  // and "Normal" in case you want to reverse the direction of the surface normal

  // NX1    - number of points in X1
  // NX2    - number of points in X2
  // X0     - Starting point
  // X1     - Point that defines X1 axis (from X0)
  // X2     - Point that defines X2 axis (from X0)
  // Normal - If -1 reverse the direction of the calculated normal vector

  // X0 is the first point and starting point for stepping
  fStartVector = X0;

  // The normal to this plane.  If it is -1 the normal vector is reversed
  fNormal = Normal;

  // Number of points to use in X1 and X2
  fNX1 = NX1;
  fNX2 = NX2;
  fNPoints = (size_t) (fNX1 * fNX2);

  // Size of X1 and X2
  double const WidthX1 = (X1 - X0).Mag();
  double const WidthX2 = (X2 - X0).Mag();

  // Stepsize for X1 and X2
  fX1StepSize = WidthX1 / (NX1 - 1);
  fX2StepSize = WidthX2 / (NX2 - 1);

  // The stepping vectors for X1 and X2
  fX1Vector = (X1 - X0) / (NX1 - 1);
  fX2Vector = (X2 - X0) / (NX2 - 1);

  // The normal vector is taken as the cross product
  fNormalVector = fX1Vector.Cross(fX2Vector).UnitVector();

  // If the input normal is -1 reverse the normal vector
  if (fNormal == -1) {
    fNormalVector *= -1;
  } else if (fNormal != 0 && fNormal != 1) {
    throw;
  }


  return;
}






void TSurfacePoints_Rectangle::Init (std::string const& Plane, int const NX1, int const NX2, double const WidthX1, double const WidthX2, TVector3D const& Rotations, TVector3D const& Translation, int const Normal)
{
  // Initialization taking a string representing where to orient the plane to start, eg in the "XY" plane.
  // Widths are given such that the rectangle with be centered at 0, 0, 0 and have +/- width/2 in either direction.
  // Rotations are done before translations.  If the Normal int is -1 the direction of the normal is reversed.
  // Direction for the normal is taken from a right handed system, eg "XY" would give +Z direction, "YX" would
  // have a normal in the -Z direction

  // Plane       - The plane you want to start in (eg XY, YX, XZ, ZX, YZ, ZX)
  // NX1         - number of points in X1
  // NX2         - number of points in X2
  // WidthX1     - Width in the X1 direction
  // WidthX2     - Width in the X2 direction
  // Rotations   - Rotation angles about the X, Y, and Z axis.  Rotations are also performed in that order
  // Translation - Take this object and move it anywhere in sace you like, done after rotations
  // Normal      - If -1 reverse the direction of the calculated normal vector


  // I will accept lower-case
  std::string P = Plane;
  std::transform(P.begin(), P.end(), P.begin(), ::toupper);

  // Number of points in each dimension and total
  fNX1 = NX1;
  fNX2 = NX2;
  fNPoints = (size_t) (fNX1 * fNX2);

  // Normal
  fNormal = Normal;

  // Stepsize in each dimension
  fX1StepSize = WidthX1 / (fNX1 - 1);
  fX2StepSize = WidthX2 / (fNX2 - 1);

  // Which plane did you pick to start!?
  if (P == "XY") {
    fStartVector.SetXYZ(-WidthX1 / 2., -WidthX2 / 2., 0);
    fX1Vector.SetXYZ(fX1StepSize, 0, 0);
    fX2Vector.SetXYZ(0, fX2StepSize, 0);
  } else if (P == "YX") {
    fStartVector.SetXYZ(-WidthX2 / 2., -WidthX1 / 2., 0);
    fX2Vector.SetXYZ(fX2StepSize, 0, 0);
    fX1Vector.SetXYZ(0, fX1StepSize, 0);
  } else if (P == "XZ") {
    fStartVector.SetXYZ(-WidthX1 / 2., 0, -WidthX2 / 2.);
    fX1Vector.SetXYZ(fX1StepSize, 0, 0);
    fX2Vector.SetXYZ(0, 0, fX2StepSize);
  } else if (P == "ZX") {
    fStartVector.SetXYZ(-WidthX2 / 2., 0, -WidthX1 / 2.);
    fX2Vector.SetXYZ(fX2StepSize, 0, 0);
    fX1Vector.SetXYZ(0, 0, fX1StepSize);
  } else if (P == "YZ") {
    fStartVector.SetXYZ(0, -WidthX1 / 2., -WidthX2 / 2.);
    fX1Vector.SetXYZ(0, fX1StepSize, 0);
    fX2Vector.SetXYZ(0, 0, fX2StepSize);
  } else if (P == "ZY") {
    fStartVector.SetXYZ(0, -WidthX2 / 2., -WidthX1 / 2.);
    fX2Vector.SetXYZ(0, fX2StepSize, 0);
    fX1Vector.SetXYZ(0, 0, fX1StepSize);
  } else {
    throw std::invalid_argument("not a valid surface string: XY YX XZ ZX YZ ZY");
  }

  // Rotate and translate the arms of our rectangle
  fStartVector.RotateSelfXYZ(Rotations);
  fX1Vector.RotateSelfXYZ(Rotations);
  fX2Vector.RotateSelfXYZ(Rotations);
  fStartVector += Translation;

  // Calculate a normal unit vector
  fNormalVector = fX1Vector.Cross(fX2Vector).UnitVector();

  // Reverse the direction of the normal if requested
  if (fNormal == -1) {
    fNormalVector *= -1;
  } else if (fNormal != 0 && fNormal != 1) {
    throw std::invalid_argument("normal must be -1, 0, or 1");
  }

  return;
}


TSurfacePoint const TSurfacePoints_Rectangle::GetPoint (size_t const i) const
{
  // Get the ith surface point
  return TSurfacePoint(this->GetXYZ(i), fNormalVector);
}




TVector3D TSurfacePoints_Rectangle::GetXYZ (size_t const i) const
{
  // Get the XYZ coordinate for the ith point

  int const ix1 = i / fNX2;
  int const ix2 = i % fNX2;

  return (fStartVector + ix1 * fX1Vector + ix2 * fX2Vector);
}




size_t TSurfacePoints_Rectangle::GetNPoints () const
{
  // Get the number of points
  return fNPoints;
}




double TSurfacePoints_Rectangle::GetX1 (size_t const i) const
{
  // Get the X1 coordinate in this frame
  int const ix1 = i / fNX2;
  return fX1StepSize * ix1 - fX1Vector.Mag() * (fNX1 - 1)/ 2.;
}




double TSurfacePoints_Rectangle::GetX2 (size_t const i) const
{
  // Get the X2 coordinate in this frame
  int const ix2 = i % fNX2;
  return fX2StepSize * ix2 - fX2Vector.Mag() * (fNX2 - 1)/ 2.;
}




double TSurfacePoints_Rectangle::GetElementArea () const
{
  // Get the elemental area
  // UPDATE: calculate area from vectors for skew
  return fX1StepSize * fX2StepSize;
}
