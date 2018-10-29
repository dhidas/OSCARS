#include "TTriangle3D.h"

#include <stdexcept>


TTriangle3D::TTriangle3D (TVector3D const& A,
                          TVector3D const& B,
                          TVector3D const& C,
                          TVector3D const& N)
{
  this->Set(A, B, C, N);
  this->SetValue(0);
  this->SetCompensation(0);
}



TTriangle3D::TTriangle3D (double const Ax, double const Ay, double const Az,
                          double const Bx, double const By, double const Bz,
                          double const Cx, double const Cy, double const Cz,
                          double const Nx, double const Ny, double const Nz)
{
  this->Set(TVector3D(Ax, Ay, Az), TVector3D(Bx, By, Bz), TVector3D(Cx, Cy, Cz), TVector3D(Nx ,Ny, Nz));
  this->SetValue(0);
  this->SetCompensation(0);
}



TTriangle3D::~TTriangle3D ()
{
}


void TTriangle3D::Set (TVector3D const& A,
                       TVector3D const& B,
                       TVector3D const& C,
                       TVector3D const& N)
{
  fA = A;
  fB = B;
  fC = C;
  fN = N;

  fCenter = this->CalculateCenter();

  return;
}

void TTriangle3D::Set (double const Ax, double const Ay, double const Az,
                       double const Bx, double const By, double const Bz,
                       double const Cx, double const Cy, double const Cz,
                       double const Nx, double const Ny, double const Nz)
{
  this->Set(TVector3D(Ax, Ay, Az), TVector3D(Bx, By, Bz), TVector3D(Cx, Cy, Cz), TVector3D(Nx, Ny, Nz));
  return;
}

void TTriangle3D::Translate (TVector3D const& T)
{
  fA += T;
  fB += T;
  fC += T;

  fCenter += T;

  return;
}


void TTriangle3D::RotateSelfXYZ (TVector3D const& R)
{
  // Rotate about x, y, then z axes
  fA.RotateSelfXYZ(R);
  fB.RotateSelfXYZ(R);
  fC.RotateSelfXYZ(R);
  fN.RotateSelfXYZ(R);

  fCenter.RotateSelfXYZ(R);

  return;
}




TVector3D TTriangle3D::CalculateCenter () const
{
  return (fA + fB + fC) / 3.;
}




TVector3D TTriangle3D::GetCenter () const
{
  return fCenter;
}




TVector3D TTriangle3D::GetNormal () const
{
  return fN;
}




void TTriangle3D::SetValue (double const V)
{
  fValue = V;
  this->SetCompensation(0);
  return;
}



double TTriangle3D::GetValue () const
{
  return fValue;
}




void TTriangle3D::SetCompensation (double const V)
{
  fCompensation = V;
  return;
}



double TTriangle3D::GetCompensation () const
{
  return fCompensation;
}



void TTriangle3D::AddToValue (double const V)
{
  // Compensated summation

  double Sum = fValue;
  double y = V - fCompensation;
  double t = Sum + y;
  fCompensation = (t - Sum) - y;
  fValue = V + Sum;

  return;
}




double TTriangle3D::RayIntersectionDistance (TVector3D const& RayOrigin, TVector3D const& RayDirectionIn) const
{
  // Returns the distance to intersection in units of ray direction (normalized)
  //
  // This is the method of Moller and Trumbore described in "Fast, Minimum Storage Ray, Triangle Intersection"

  TVector3D const RayDirection = RayDirectionIn.UnitVector();
  TVector3D edge1, edge2, tvec, pvec, qvec;
  double det, inv_det, t, u, v;

  // Find vectors for two edges sharing fA
  edge1 = fB - fA;
  edge2 = fC - fA;

  // Begin calculating determinant - also used to calculate U parameter
  pvec = RayDirection.Cross(edge2);

  // If determinant is near zero, ray lies in plane of triangle
  det = edge1.Dot(pvec);


  // EPSILON 0.000001
  if (det > -1e-12 && det < 1e-12) {
    return -1;
  }
  inv_det = 1.0 / det;

  // Calculate distance from fA to ray RayOriginin
  tvec = RayOrigin - fA;

  // Calculate U parameter and test bounds
  u = tvec.Dot(pvec) * inv_det;
  if (u < 0.0 || u > 1.0) {
    return -1;
  }

  // Prepare to test V parameter
  qvec = tvec.Cross(edge1);

  // Calculate V parameter and test bounds
  v = RayDirection.Dot(qvec) * inv_det;
  if (v < 0.0 || u + v > 1.0) {
    return -1;
  }

  // Calculate t, ray intersects triangle
  t = edge2.Dot(qvec) * inv_det;

  if (t < 0) {
    return -1;
  }

  return t;
}





TVector3D TTriangle3D::operator [] (int const i) const
{
  // An operator to use an index like a vector
  // For getting a point

  switch (i) {
    case 0:
      return fA;
    case 1:
      return fB;
    case 2:
      return fC;
    case 3:
      return fN;
    default:
      throw std::out_of_range("TTriangle3D::operator [] index out of range");
  }
}





TVector3D& TTriangle3D::operator [] (int const i)
{
  // An operator to use an index like a vector
  // For setting a point

  switch (i) {
    case 0:
      return fA;
    case 1:
      return fB;
    case 2:
      return fC;
    case 3:
      return fN;
    default:
      throw std::out_of_range("TTriangle3D::operator [] index out of range");
  }
}





TTriangle3D& TTriangle3D::operator += (TVector3D const& V)
{
  // Add a vector to this vector by components
  fA += V;
  fB += V;
  fC += V;

  fCenter += V;

  return *this;
}




