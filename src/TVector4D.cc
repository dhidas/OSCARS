#include "TVector4D.h"


TVector4D::TVector4D ()
{
  // Default constructor
}




TVector4D::TVector4D (double const X, double const Y, double const Z, double const E)
{
  // Constructor given the individual arguments
  fP.SetXYZ(X, Y, Z);
  fE = E;
}




TVector4D::TVector4D (TVector3D const& X, double const E)
{
  // Constructor given the 3-vector and E/t
  fP = X;
  fE = E;
}




TVector4D::~TVector4D ()
{
  // I own nothing
}
