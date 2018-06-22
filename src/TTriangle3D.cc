#include "TTriangle3D.h"

#include <stdexcept>


TTriangle3D::TTriangle3D (TVector3D const& A,
                          TVector3D const& B,
                          TVector3D const& C,
                          TVector3D const& N)
{
  this->Set(A, B, C, N);
}



TTriangle3D::TTriangle3D (double const Ax, double const Ay, double const Az,
                          double const Bx, double const By, double const Bz,
                          double const Cx, double const Cy, double const Cz,
                          double const Nx, double const Ny, double const Nz)
{
  this->Set(TVector3D(Ax, Ay, Az), TVector3D(Bx, By, Bz), TVector3D(Cx, Cy, Cz), TVector3D(Nx ,Ny, Nz));
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
  return;
}


void TTriangle3D::RotateSelfXYZ (TVector3D const& R)
{
  // Rotate about x, y, then z axes
  fA.RotateSelfXYZ(R);
  fB.RotateSelfXYZ(R);
  fC.RotateSelfXYZ(R);
  fN.RotateSelfXYZ(R);

  return;
}




TVector3D TTriangle3D::GetCenter () const
{
  return (fA + fB + fC) / 3.;
}




TVector3D TTriangle3D::GetNormal () const
{
  return fN;
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
  return *this;
}




